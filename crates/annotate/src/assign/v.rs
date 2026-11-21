use super::ties::tiebreak_cmp;
use super::{primary_by_tiebreak, Oriented, SegmentIndex};
use crate::align::{align_semiglobal, collapse_gap_pairs, place_gaps, score_semiglobal, Alignment, Op, Scoring};
use abstar_germline::{Allele, Locus};

/// The V call for one oriented read.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct VHit {
    /// Database index of the primary allele (record 0008); `ties[0]`.
    pub primary: usize,
    /// Database indices of the same-locus V alleles tied for the fewest mismatches over the
    /// read's V-derived bases (record 0029), primary first, then record 0008 order.
    pub ties: Vec<usize>,
    /// The read against the primary allele's ungapped sequence, cut to the V-derived bases.
    pub aln: Alignment,
    /// Read position -> IMGT gapped index (0-based) of the primary; None = not aligned to it
    /// (outside the alignment, or inserted).
    pub imgt: Vec<Option<usize>>,
    /// The tie set's mismatch count over the read's V-derived bases.
    pub mismatches: u32,
}

/// Up to this many k-mer-ranked V genes are scored by alignment (record 0029) ...
const TOP: usize = 3;
/// ... those whose k-mer count is at least this fraction of the leading gene's ...
const RATIO: f64 = 0.9;
/// ... and up to this many when none of them reaches `MIN_V_SCORE`.
const FALLBACK: usize = 12;
/// Diagonal bin width for the fallback ranking (`KmerIndex::diagonal_hits`).
const DIAG_BIN: usize = 8;
/// Below this many diagonal-banded k-mer hits for every V allele, the fallback is skipped (no V
/// hit on this strand).
const MIN_DIAG: u32 = 8;
/// Projection refinement rounds (see `assign_v`).
const REFINE: usize = 2;
/// J alleles scored for J evidence, and the J score that opens the search of every V gene of
/// the J's locus when the fallback finds no V (record 0029).
const J_TOP: usize = 3;
const MIN_J_SCORE: i32 = 30;
/// ... or the best V score so far that opens it on the best V's locus.
const NEAR_V_SCORE: i32 = 20;

/// Minimum V alignment score (`score_semiglobal`) for a V hit (record 0029).
pub const MIN_V_SCORE: i32 = 50;

/// The V's 5' end extends to the read's or the allele's start when that is at most this many
/// bases away on the alignment's first diagonal (record 0029) ...
const MAX_START_EXT: usize = 20;
/// ... and the extension is at most this long or matches the allele on at least a quarter of
/// its bases, as often as random sequence would.
const FREE_START_EXT: usize = 3;

/// The read against one V allele.
struct Projected {
    /// The overlap alignment's DP score (`score_semiglobal`'s), for ranking alleles.
    score: i32,
    /// The alignment cut to the read's V-derived bases.
    aln: Alignment,
    imgt: Vec<Option<usize>>,
    /// The V-derived bases as (IMGT gapped index, read base), in read order.
    pairs: Vec<(usize, u8)>,
    /// The alignment has a gap.
    gapped: bool,
}

/// `aln` cut to its maximum-scoring contiguous stretch of columns (a gap run is one column): the
/// read's V-derived bases. Columns score on a symmetric scale: a match `+s.matched`, a mismatch
/// `-s.matched` (±2 under record 0029's scoring, not `s.mismatch`), a gap run its affine cost under
/// `s`. Among equal scores the longest stretch wins, then the leftmost. With `keep_end` (J) only
/// the start is cut: the stretch ends at the last column, and among equal scores the shortest
/// wins. The score becomes that stretch's score on this scale (at least 0), so J's `MIN_SCORE`
/// floor of 30 is on it too.
pub(super) fn trim(aln: &Alignment, s: &Scoring, keep_end: bool) -> Alignment {
    let gap = |n: u32| s.gap_open as i32 + (n as i32 - 1) * s.gap_extend as i32;
    let mut cols: Vec<(Op, u32, i32)> = Vec::new();
    for &(op, n) in &aln.ops {
        match op {
            Op::Match => cols.extend((0..n).map(|_| (op, 1, s.matched as i32))),
            Op::Mismatch => cols.extend((0..n).map(|_| (op, 1, -(s.matched as i32)))),
            Op::Ins | Op::Del => cols.push((op, n, gap(n))),
        }
    }
    // Kadane over prefix sums: for each end, the start at the earliest minimal prefix.
    let (mut prefix, mut low, mut low_at) = (0, 0, 0);
    let (mut best, mut start, mut end) = (i32::MIN, 0, 0);
    for (c, col) in cols.iter().enumerate() {
        prefix += col.2;
        let (sum, len) = (prefix - low, c + 1 - low_at);
        if (!keep_end || c + 1 == cols.len()) && (sum > best || (sum == best && len > end - start)) {
            (best, start, end) = (sum, low_at, c + 1);
        }
        if prefix < low || (keep_end && prefix == low) {
            (low, low_at) = (prefix, c + 1);
        }
    }
    let advance = |cols: &[(Op, u32, i32)]| {
        cols.iter().fold((0, 0), |(q, g), &(op, n, _)| match op {
            Op::Match | Op::Mismatch => (q + n as usize, g + n as usize),
            Op::Ins => (q + n as usize, g),
            Op::Del => (q, g + n as usize),
        })
    };
    let (dq, dg) = advance(&cols[..start]);
    let (kq, kg) = advance(&cols[start..end]);
    let mut ops: Vec<(Op, u32)> = Vec::new();
    for &(op, n, _) in &cols[start..end] {
        match ops.last_mut() {
            Some((last, m)) if *last == op => *m += n,
            _ => ops.push((op, n)),
        }
    }
    let (q_start, g_start) = (aln.q_start + dq, aln.g_start + dg);
    Alignment { score: best.max(0), q_start, q_end: q_start + kq, g_start, g_end: g_start + kg, ops }
}

/// `aln` ending before its first gap that starts at germline position `from` or later. V has no
/// indels on the CDR3 side of Cys104 (the generator's rule), so such a gap marks where V ended.
fn before_gap_in(aln: &Alignment, from: usize, s: &Scoring) -> Alignment {
    let (mut qi, mut gi) = (aln.q_start, aln.g_start);
    for (i, &(op, n)) in aln.ops.iter().enumerate() {
        if matches!(op, Op::Ins | Op::Del) && gi >= from {
            let ops = aln.ops[..i].to_vec();
            let score = trim(&Alignment { ops: ops.clone(), q_end: qi, g_end: gi, ..aln.clone() }, s, false).score;
            return Alignment { score, q_start: aln.q_start, q_end: qi, g_start: aln.g_start, g_end: gi, ops };
        }
        match op {
            Op::Match | Op::Mismatch => (qi, gi) = (qi + n as usize, gi + n as usize),
            Op::Ins => qi += n as usize,
            Op::Del => gi += n as usize,
        }
    }
    aln.clone()
}

/// `aln` extended 5' along its first diagonal to the read's or the allele's start, whichever
/// comes first (record 0029): the trim cuts mutated first V bases, which reaching the read's or
/// the allele's start shows are V. The extension is taken when it is at most `MAX_START_EXT`
/// bases long and either at most `FREE_START_EXT` long or matching on a quarter of its bases.
/// The score adds the new columns as `trim` scores them.
fn extend_start(aln: Alignment, q: &[u8], g: &[u8], s: &Scoring) -> Alignment {
    let k = aln.q_start.min(aln.g_start);
    if k == 0 || k > MAX_START_EXT {
        return aln;
    }
    let (q_start, g_start) = (aln.q_start - k, aln.g_start - k);
    let hit = |x: usize| matches!(q[q_start + x], b'A' | b'C' | b'G' | b'T') && q[q_start + x] == g[g_start + x];
    let cols: Vec<Op> = (0..k).map(|x| if hit(x) { Op::Match } else { Op::Mismatch }).collect();
    if k > FREE_START_EXT && 4 * cols.iter().filter(|&&op| op == Op::Match).count() < k {
        return aln;
    }
    let score = aln.score + cols.iter().map(|&op| if op == Op::Match { s.matched as i32 } else { -(s.matched as i32) }).sum::<i32>();
    let mut ops: Vec<(Op, u32)> = Vec::new();
    for (op, n) in cols.into_iter().map(|op| (op, 1)).chain(aln.ops.iter().copied()) {
        match ops.last_mut() {
            Some((last, m)) if *last == op => *m += n,
            _ => ops.push((op, n)),
        }
    }
    Alignment { score, q_start, g_start, ops, ..aln }
}

/// Aligns `q` to V allele `k` of `v`, collapses spurious insertion/deletion pairs, cuts the
/// alignment to the read's V-derived bases (`trim`), places its gaps outside the CDR3 side of
/// Cys104 (leftmost among equals), and projects the read onto the allele's IMGT gapped
/// positions.
fn project(q: &[u8], alleles: &[Allele], v: &SegmentIndex, k: usize, s: &Scoring) -> Projected {
    let (a, g) = (&alleles[v.ids[k]], &v.seqs[k]);
    let raw = align_semiglobal(q, g, s);
    let protect = a.anchor.map(|c| c.saturating_sub(2)..g.len());
    let mut cut = trim(&collapse_gap_pairs(&raw, q, g, s), s, false);
    if let Some(zone) = &protect {
        cut = before_gap_in(&cut, zone.start, s);
    }
    let aln = extend_start(place_gaps(&cut, q, g, protect, s), q, g, s);
    // Ungapped allele index -> IMGT gapped index.
    let imgt_of: Vec<usize> = a.sequence_gapped.bytes().enumerate().filter(|&(_, b)| b != b'.').map(|(i, _)| i).collect();
    let mut imgt = vec![None; q.len()];
    let mut pairs = Vec::with_capacity(aln.q_end - aln.q_start);
    let (mut qi, mut gi) = (aln.q_start, aln.g_start);
    for &(op, n) in &aln.ops {
        let n = n as usize;
        match op {
            Op::Match | Op::Mismatch => {
                for _ in 0..n {
                    imgt[qi] = Some(imgt_of[gi]);
                    pairs.push((imgt_of[gi], q[qi]));
                    qi += 1;
                    gi += 1;
                }
            }
            Op::Ins => qi += n,
            Op::Del => gi += n,
        }
    }
    let gapped = aln.ops.iter().any(|&(op, _)| matches!(op, Op::Ins | Op::Del));
    Projected { score: raw.score, aln, imgt, pairs, gapped }
}

/// A read base differs from the candidate's (uppercase) base at its IMGT gapped index: a
/// non-ACGT read base, a gap or a position past the candidate's end is a mismatch.
fn differs((g, b): (usize, u8), gapped: &[u8]) -> bool {
    !matches!(b, b'A' | b'C' | b'G' | b'T') || gapped.get(g) != Some(&b)
}

/// Mismatch counts of every same-locus V allele (database index, count) over `p`'s V-derived
/// bases.
fn counts(p: &Projected, v: &SegmentIndex, locus: Locus) -> Vec<(usize, u32)> {
    (0..v.ids.len())
        .filter(|&i| v.loci[i] == locus)
        .map(|i| (v.ids[i], p.pairs.iter().filter(|&&x| differs(x, &v.gapped[i])).count() as u32))
        .collect()
}

/// Alleles with the fewest mismatches, and that count.
pub(super) fn fewest(c: &[(usize, u32)]) -> (Vec<usize>, u32) {
    let min = c.iter().map(|x| x.1).min().unwrap_or(0);
    (c.iter().filter(|x| x.1 == min).map(|x| x.0).collect(), min)
}

/// Assigns V to an oriented read (record 0029).
///
/// 1. **Candidates.** V genes are ranked by their best allele's k-mer count. The best allele of
///    each of the first `TOP` genes within `RATIO` of the leader is scored with
///    `score_semiglobal` (a lone leader is not scored). When the best scores below
///    `MIN_V_SCORE`, up to `FALLBACK` genes are; when that fails too and a J allele scores at
///    least `MIN_J_SCORE`, the best allele of every V gene of that J's locus is. Below
///    `MIN_V_SCORE` there is no V hit (`None`).
/// 2. **Projection.** The best allele is aligned, spurious insertion/deletion pairs are
///    collapsed, the alignment is cut to the read's V-derived bases, and the read is placed on
///    IMGT gapped positions through it.
/// 3. **Refinement.** While that alignment has a gap, or another allele has fewer mismatches
///    than the projecting allele, the fewest-mismatch alleles not yet scored (up to `TOP`) are
///    scored; a higher score projects instead (at most `REFINE` rounds).
/// 4. **Ties.** Every V allele of the projecting allele's locus is counted for mismatches over
///    the V-derived bases at their IMGT positions (inserted bases skipped; a base the allele
///    doesn't cover is a mismatch). The record 0008 pick among the fewest is then aligned and
///    cuts the V-derived bases itself, and the count is repeated once over its cut: the tie
///    set is the alleles with the fewest mismatches there, the primary their record 0008 pick.
/// 5. The primary's own alignment and projection are reported.
pub fn assign_v(read: &Oriented, alleles: &[Allele], v: &SegmentIndex, j: &SegmentIndex, scoring: &Scoring) -> Option<VHit> {
    let q = &read.seq;
    let mut hits = v.index.hits(q);
    hits.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
    // The best-ranked allele of each gene, in k-mer order.
    let mut genes: Vec<&str> = Vec::new();
    let gene_best: Vec<(usize, u32)> = hits
        .iter()
        .filter(|&&(k, _)| {
            let gene = alleles[v.ids[k]].gene.as_str();
            let new = !genes.contains(&gene);
            if new {
                genes.push(gene);
            }
            new
        })
        .copied()
        .collect();
    let lead = gene_best.first()?.1;
    let close = gene_best.iter().take(TOP).take_while(|g| g.1 as f64 >= RATIO * lead as f64).count();
    let mut scored: Vec<usize> = Vec::new();
    // Scores the not-yet-scored alleles among `ks`; the best (index, score).
    let score_best = |ks: &mut dyn Iterator<Item = usize>, scored: &mut Vec<usize>| {
        let mut best: Option<(usize, i32)> = None;
        for k in ks {
            if !scored.contains(&k) {
                scored.push(k);
                let sc = score_semiglobal(q, &v.seqs[k], scoring).0;
                if best.is_none_or(|(_, b)| sc > b) {
                    best = Some((k, sc));
                }
            }
        }
        best
    };
    let mut k = if close == 1 { gene_best[0].0 } else { score_best(&mut gene_best.iter().take(close).map(|g| g.0), &mut scored)?.0 };
    scored.push(k);
    let mut p = project(q, alleles, v, k, scoring);
    if p.score < MIN_V_SCORE {
        // Genes ranked by their best allele's diagonal-banded k-mer hits, which tolerate heavy
        // mutation better than plain counts; then the genes without any.
        let mut dh = v.index.diagonal_hits(q, DIAG_BIN);
        dh.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
        if dh.first().is_none_or(|h| h.1 < MIN_DIAG) {
            return None;
        }
        let mut dgenes: Vec<&str> = Vec::new();
        let ranked: Vec<usize> = dh
            .iter()
            .map(|h| h.0)
            .chain(0..v.ids.len())
            .filter(|&i| {
                let gene = alleles[v.ids[i]].gene.as_str();
                let new = !dgenes.contains(&gene);
                if new {
                    dgenes.push(gene);
                }
                new
            })
            .collect();
        let fallback = score_best(&mut ranked.iter().take(FALLBACK).copied(), &mut scored);
        // The best V so far, the projected allele included.
        let near = fallback.filter(|b| b.1 > p.score).unwrap_or((k, p.score));
        let mut found = Some(near).filter(|b| b.1 >= MIN_V_SCORE);
        if found.is_none() {
            // Every V gene of one locus: that of a J with real evidence, else that of the best V
            // if it scores above the null.
            let mut jh = j.index.hits(q);
            jh.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
            let jbest = jh.iter().take(J_TOP).map(|&(jk, _)| (jk, score_semiglobal(q, &j.seqs[jk], scoring).0)).max_by_key(|&(jk, sc)| (sc, std::cmp::Reverse(jk)));
            let locus = match jbest {
                Some((jk, sc)) if sc >= MIN_J_SCORE => Some(j.loci[jk]),
                _ => (near.1 >= NEAR_V_SCORE).then(|| v.loci[near.0]),
            };
            if let Some(locus) = locus {
                found = score_best(&mut ranked.iter().copied().filter(|&i| v.loci[i] == locus), &mut scored).filter(|b| b.1 >= MIN_V_SCORE);
            }
        }
        k = found?.0;
        p = project(q, alleles, v, k, scoring);
    }
    let locus = v.loci[k];
    let mut c = counts(&p, v, locus);
    for _ in 0..REFINE {
        let (best, _) = fewest(&c);
        if !p.gapped && best.contains(&v.ids[k]) {
            break;
        }
        let mut order: Vec<(u32, usize)> = c.iter().map(|&(id, n)| (n, id)).collect();
        order.sort_unstable();
        let mut fresh =
            order.iter().map(|&(_, id)| v.ids.binary_search(&id).expect("V allele")).filter(|k2| !scored.contains(k2)).take(TOP).collect::<Vec<_>>().into_iter();
        match score_best(&mut fresh, &mut scored) {
            Some((k2, sc)) if sc > p.score => {
                k = k2;
                p = project(q, alleles, v, k, scoring);
                c = counts(&p, v, locus);
            }
            _ => break,
        }
    }
    // The fewest-mismatch pick cuts the V-derived bases itself; the exact rule is applied once
    // over its cut.
    let pick = primary_by_tiebreak(&fewest(&c).0, alleles);
    if pick != v.ids[k] {
        k = v.ids.binary_search(&pick).expect("tie is a V allele");
        p = project(q, alleles, v, k, scoring);
        c = counts(&p, v, locus);
    }
    let (mut ties, mismatches) = fewest(&c);
    ties.sort_by(|&x, &y| tiebreak_cmp(&alleles[x], &alleles[y]).then(x.cmp(&y)));
    let primary = ties[0];
    if primary != v.ids[k] {
        p = project(q, alleles, v, v.ids.binary_search(&primary).expect("tie is a V allele"), scoring);
    }
    Some(VHit { primary, ties, aln: p.aln, imgt: p.imgt, mismatches })
}
