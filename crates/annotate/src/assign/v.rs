use super::ties::tiebreak_cmp;
use super::{Oriented, SegmentIndex};
use crate::align::{align_semiglobal, place_gaps, score_semiglobal, Alignment, Op, Scoring};
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
    /// The fewest mismatches over the read's V-derived bases (the primary's count).
    pub mismatches: u32,
}

/// Up to this many k-mer-ranked V genes are scored by alignment (record 0029) ...
const TOP: usize = 3;
/// ... those whose k-mer count is at least this fraction of the leading gene's ...
const RATIO: f64 = 0.9;
/// ... and up to this many when none of them reaches `MIN_V_SCORE`.
const FALLBACK: usize = 12;
/// Projection refinement rounds (see `assign_v`).
const REFINE: usize = 2;
/// Each end of the V-derived bases may lie this many bases either side of its estimate.
const END_SLACK: usize = 3;

/// Minimum V alignment score (`score_semiglobal`) for a V hit (record 0029).
pub const MIN_V_SCORE: i32 = 50;

/// The read against one V allele.
struct Projected {
    /// The overlap alignment's DP score (`score_semiglobal`'s), for ranking alleles.
    score: i32,
    /// The alignment cut to the read's V-derived bases.
    aln: Alignment,
    imgt: Vec<Option<usize>>,
    /// The V-derived bases as (IMGT gapped index, read base), in read order.
    pairs: Vec<(usize, u8)>,
    /// Up to `END_SLACK` read bases just before / after the V-derived bases, continuing the
    /// alignment's diagonal (nearest first), as (IMGT gapped index, read base).
    head: Vec<(usize, u8)>,
    tail: Vec<(usize, u8)>,
    /// The alignment has a gap.
    gapped: bool,
}

/// `aln` cut to its best-scoring stretch under the boundary score (match +1, mismatch -1, gaps
/// as `s`): the 3' end is the last column where the running score peaks, the 5' end the first
/// column where the score from there to that end peaks. The score becomes the kept ops' score
/// under `s`.
fn trim(aln: &Alignment, s: &Scoring) -> Alignment {
    let gap = |n: u32| s.gap_open as i32 + (n as i32 - 1) * s.gap_extend as i32;
    // Columns: (op, length); a gap run is one column.
    let mut cols: Vec<(Op, u32)> = Vec::new();
    for &(op, n) in &aln.ops {
        match op {
            Op::Match | Op::Mismatch => cols.extend((0..n).map(|_| (op, 1))),
            Op::Ins | Op::Del => cols.push((op, n)),
        }
    }
    let boundary = |&(op, n): &(Op, u32)| match op {
        Op::Match => 1,
        Op::Mismatch => -1,
        Op::Ins | Op::Del => gap(n),
    };
    let (mut run, mut peak, mut end) = (0, i32::MIN, 0);
    for (c, col) in cols.iter().enumerate() {
        run += boundary(col);
        if run >= peak {
            (peak, end) = (run, c + 1);
        }
    }
    let (mut run, mut peak, mut start) = (0, i32::MIN, end);
    for c in (0..end).rev() {
        run += boundary(&cols[c]);
        if run >= peak {
            (peak, start) = (run, c);
        }
    }
    let advance = |cols: &[(Op, u32)]| {
        cols.iter().fold((0, 0), |(q, g), &(op, n)| match op {
            Op::Match | Op::Mismatch => (q + n as usize, g + n as usize),
            Op::Ins => (q + n as usize, g),
            Op::Del => (q, g + n as usize),
        })
    };
    let (dq, dg) = advance(&cols[..start]);
    let (kq, kg) = advance(&cols[start..end]);
    let mut ops: Vec<(Op, u32)> = Vec::new();
    for &(op, n) in &cols[start..end] {
        match ops.last_mut() {
            Some((last, m)) if *last == op => *m += n,
            _ => ops.push((op, n)),
        }
    }
    let score = ops
        .iter()
        .map(|&(op, n)| match op {
            Op::Match => n as i32 * s.matched as i32,
            Op::Mismatch => n as i32 * s.mismatch as i32,
            Op::Ins | Op::Del => gap(n),
        })
        .sum();
    let (q_start, g_start) = (aln.q_start + dq, aln.g_start + dg);
    Alignment { score, q_start, q_end: q_start + kq, g_start, g_end: g_start + kg, ops }
}

/// Aligns `q` to V allele `k` of `v`, cuts the alignment to the read's V-derived bases
/// (`trim`), places its gaps outside the CDR3 side of Cys104 (leftmost among equals), and
/// projects the read onto the allele's IMGT gapped positions.
fn project(q: &[u8], alleles: &[Allele], v: &SegmentIndex, k: usize, s: &Scoring) -> Projected {
    let (a, g) = (&alleles[v.ids[k]], &v.seqs[k]);
    let raw = align_semiglobal(q, g, s);
    let protect = a.anchor.map(|c| c.saturating_sub(2)..g.len());
    let aln = place_gaps(&trim(&raw, s), q, g, protect, s);
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
    let empty = aln.ops.is_empty();
    let head = (1..=END_SLACK)
        .take_while(|&i| !empty && i <= aln.q_start && i <= aln.g_start)
        .map(|i| (imgt_of[aln.g_start - i], q[aln.q_start - i]))
        .collect();
    let tail = (0..END_SLACK)
        .take_while(|&i| !empty && aln.q_end + i < q.len() && aln.g_end + i < g.len())
        .map(|i| (imgt_of[aln.g_end + i], q[aln.q_end + i]))
        .collect();
    let gapped = aln.ops.iter().any(|&(op, _)| matches!(op, Op::Ins | Op::Del));
    Projected { score: raw.score, aln, imgt, pairs, head, tail, gapped }
}

/// A read base differs from the candidate's (uppercase) base at its IMGT gapped index: a
/// non-ACGT read base, a gap or a position past the candidate's end is a mismatch.
fn differs((g, b): (usize, u8), gapped: &[u8]) -> bool {
    !(matches!(b, b'A' | b'C' | b'G' | b'T') && gapped.get(g) == Some(&b))
}

/// Mismatch counts of one candidate over the V-derived bases and over each span with one end
/// moved by up to `END_SLACK` bases: [estimate, 5' end in by 1..=n, 5' end out by 1..=n,
/// 3' end in by 1..=n, 3' end out by 1..=n]. Which variants exist depends only on the
/// projection, so the vectors line up across candidates.
fn span_counts(p: &Projected, gapped: &[u8]) -> Vec<u32> {
    let core = p.pairs.iter().filter(|&&x| differs(x, gapped)).count() as u32;
    let mut out = vec![core];
    let shrink = END_SLACK.min(p.pairs.len().saturating_sub(1) / 2);
    let running = |it: &mut dyn Iterator<Item = &(usize, u8)>, out: &mut Vec<u32>, add: bool| {
        let mut n = core;
        for &x in it {
            if differs(x, gapped) {
                n = if add { n + 1 } else { n - 1 };
            }
            out.push(n);
        }
    };
    running(&mut p.pairs.iter().take(shrink), &mut out, false);
    running(&mut p.head.iter(), &mut out, true);
    running(&mut p.pairs.iter().rev().take(shrink), &mut out, false);
    running(&mut p.tail.iter(), &mut out, true);
    out
}

/// `span_counts` of every same-locus V allele (database index, counts).
fn counts(p: &Projected, v: &SegmentIndex, locus_of: &[Locus], locus: Locus) -> Vec<(usize, Vec<u32>)> {
    (0..v.ids.len()).filter(|&i| locus_of[i] == locus).map(|i| (v.ids[i], span_counts(p, &v.gapped[i]))).collect()
}

/// Assigns V to an oriented read (record 0029).
///
/// 1. **Candidates.** V genes are ranked by their best allele's k-mer count. The best allele of
///    each of the first `TOP` genes within `RATIO` of the leader is scored with
///    `score_semiglobal` (a lone leader is not scored); when the best scores below
///    `MIN_V_SCORE`, up to `FALLBACK` genes are. Below `MIN_V_SCORE` there is no V hit (`None`).
/// 2. **Projection.** The best allele is aligned, the alignment is cut to the read's V-derived
///    bases, and the read is placed on IMGT gapped positions through it.
/// 3. **Refinement.** While that alignment has a gap, or another allele has fewer mismatches
///    than the projecting allele, the fewest-mismatch alleles not yet scored (up to `TOP`) are
///    scored; a higher score projects instead (at most `REFINE` rounds).
/// 4. **Ties.** Every V allele of the projecting allele's locus is counted for mismatches over
///    the V-derived bases at their IMGT positions (inserted bases skipped; a base the allele
///    doesn't cover is a mismatch). The primary is the record 0008 pick among the alleles with
///    the fewest. The tie set adds every allele with the fewest mismatches when either end of
///    the V-derived bases moves by up to `END_SLACK` bases, since the V/N boundary and a
///    5'-truncated start are uncertain by a few bases.
/// 5. The primary's own alignment and projection are reported.
pub fn assign_v(read: &Oriented, alleles: &[Allele], v: &SegmentIndex, scoring: &Scoring) -> Option<VHit> {
    let q = &read.seq;
    let mut hits = v.index.hits(q);
    hits.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
    // The best-ranked allele of each gene, in k-mer order.
    let mut genes: Vec<&str> = Vec::new();
    let gene_best: Vec<(usize, u32)> = hits
        .iter()
        .filter(|&&(k, _)| {
            let gene = alleles[v.ids[k]].gene.as_str();
            let new = genes.len() < FALLBACK && !genes.contains(&gene);
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
        let (k2, sc) = score_best(&mut gene_best.iter().map(|g| g.0), &mut scored)?;
        if sc < MIN_V_SCORE {
            return None;
        }
        k = k2;
        p = project(q, alleles, v, k, scoring);
    }
    let locus = alleles[v.ids[k]].locus;
    let locus_of: Vec<Locus> = v.ids.iter().map(|&id| alleles[id].locus).collect();
    let mut c = counts(&p, v, &locus_of, locus);
    for _ in 0..REFINE {
        let min = c.iter().map(|x| x.1[0]).min().unwrap_or(0);
        if !p.gapped && c.iter().any(|(id, n)| *id == v.ids[k] && n[0] == min) {
            break;
        }
        let mut order: Vec<(u32, usize)> = c.iter().map(|(id, n)| (n[0], *id)).collect();
        order.sort_unstable();
        let mut fresh =
            order.iter().map(|&(_, id)| v.ids.binary_search(&id).expect("V allele")).filter(|k2| !scored.contains(k2)).take(TOP).collect::<Vec<_>>().into_iter();
        match score_best(&mut fresh, &mut scored) {
            Some((k2, sc)) if sc > p.score => {
                k = k2;
                p = project(q, alleles, v, k, scoring);
                c = counts(&p, v, &locus_of, locus);
            }
            _ => break,
        }
    }
    // Per-span counts of every candidate; the primary comes from the estimated span.
    let spans = c;
    let nspan = spans.first().map_or(0, |s| s.1.len());
    let mins: Vec<u32> = (0..nspan).map(|j| spans.iter().map(|s| s.1[j]).min().unwrap_or(0)).collect();
    let best: Vec<usize> = spans.iter().filter(|s| s.1[0] == mins[0]).map(|s| s.0).collect();
    let primary = super::primary_by_tiebreak(&best, alleles);
    let mut ties: Vec<usize> = spans.iter().filter(|s| s.1.iter().zip(&mins).any(|(n, m)| n == m)).map(|s| s.0).collect();
    ties.sort_by(|&x, &y| (x != primary).cmp(&(y != primary)).then(tiebreak_cmp(&alleles[x], &alleles[y])).then(x.cmp(&y)));
    if primary != v.ids[k] {
        p = project(q, alleles, v, v.ids.binary_search(&primary).expect("tie is a V allele"), scoring);
    }
    Some(VHit { primary, ties, aln: p.aln, imgt: p.imgt, mismatches: mins[0] })
}
