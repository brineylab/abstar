use super::ties::tiebreak_cmp;
use super::v::{fewest, trim};
use super::{primary_by_tiebreak, Oriented, SegmentIndex, VHit};
use crate::align::{align_semiglobal, collapse_gap_pairs, place_gaps, score_semiglobal, Alignment, Op, Scoring};
use abstar_germline::{Allele, Locus};

/// The J call for one oriented read.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct JHit {
    /// Database index of the primary allele (record 0008); `ties[0]`.
    pub primary: usize,
    /// Database indices of the V locus's J alleles tied for the fewest mismatches over the
    /// read's J-derived bases (record 0029), primary first, then record 0008 order.
    pub ties: Vec<usize>,
    /// The read against the primary allele's sequence, cut to the J-derived bases.
    pub aln: Alignment,
    /// The tie set's mismatch count over the read's J-derived bases.
    pub mismatches: u32,
}

/// The J search starts this many bases before the V alignment's end (record 0029).
const OVERLAP: usize = 10;
/// The J alleles with the best `score_semiglobal` scores that are cut (`cut`); the best cut
/// score projects (record 0029) ...
const TOP: usize = 3;
/// ... if its score (on `trim`'s ±2 scale: +2 a match, −2 a mismatch) is at least this; else a
/// weak J (below) or a short 3' fragment (`tail_aln`) may project, or there is no J hit.
const MIN_SCORE: i32 = 30;
/// Weak J (record 0029, Wave 2): below `MIN_SCORE`, the locus J alleles are ranked again by
/// `score_semiglobal` with this mismatch score and the `WEAK_TOP` best are cut with it, so a
/// heavily mutated J is not cut short, ...
const WEAK_MISMATCH: i8 = -1;
/// 8 matches cutting every allele on the 12 benchmark sets; 3 cost 0.2 pp of summed J + junction
/// accuracy (Wave 2 fix round 1).
const WEAK_TOP: usize = 8;
/// ... and the best cut that holds the allele's Trp/Phe118 codon is a candidate if it scores at
/// least this (on `trim`'s ±2 scale).
const WEAK_MIN_SCORE: i32 = 12;
/// Low path (record 0029, Wave 3): when neither fallback above calls, the best gapless cut of the
/// `TOP` normal-scoring cuts is called if it scores at least this (on `trim`'s ±2 scale: 3 net
/// matches) and bounds the junction. The lowest floor whose step still adds more right than
/// wrong calls on the benchmark sets.
const LOW_MIN_SCORE: i32 = 6;
/// The later-start count (ruling B4) applies when the projecting cut has at most this many
/// mismatched or gapped bases.
const LATE_MAX_MM: u32 = 3;

/// Aligns `q` to J allele `k` of `j`, collapses spurious insertion/deletion pairs, cuts the
/// 5' end to the read's J-derived bases (`trim`) and places gaps outside the CDR3 side of
/// Trp/Phe118 (the generator's J indel zone starts 5 nt past the anchor's first base).
fn cut(q: &[u8], alleles: &[Allele], j: &SegmentIndex, k: usize, s: &Scoring) -> Alignment {
    let g = &j.seqs[k];
    let raw = align_semiglobal(q, g, s);
    let cut = trim(&collapse_gap_pairs(&raw, q, g, s), s, true);
    place_gaps(&cut, q, g, alleles[j.ids[k]].anchor.map(|w| 0..w + 5), s)
}

/// +1 when read base `r` is ACGT and equals germline base `g`, else −1.
fn hit(r: u8, g: u8) -> i32 {
    if matches!(r, b'A' | b'C' | b'G' | b'T') && r == g { 1 } else { -1 }
}

/// The best gapless alignment of the read's 3' end to `g`: over every diagonal, the suffix of
/// the read with the highest +1/−1 score (the shortest among equal scores). `None` when no
/// suffix scores above 0.
fn tail_aln(q: &[u8], g: &[u8], s: &Scoring) -> Option<(i32, Alignment)> {
    let n = q.len();
    let mut best: Option<(i32, usize, usize)> = None; // (score, length, germline end)
    for ge in 1..=g.len() {
        let (mut sum, mut bs, mut bl) = (0, 0, 0);
        for k in 1..=ge.min(n) {
            sum += hit(q[n - k], g[ge - k]);
            if sum > bs {
                (bs, bl) = (sum, k);
            }
        }
        if bs > 0 && best.is_none_or(|b| bs > b.0) {
            best = Some((bs, bl, ge));
        }
    }
    let (bs, len, ge) = best?;
    let (qs, gs) = (n - len, ge - len);
    let mut ops: Vec<(Op, u32)> = Vec::new();
    let mut score = 0;
    for x in 0..len {
        let op = if hit(q[qs + x], g[gs + x]) > 0 { Op::Match } else { Op::Mismatch };
        score += if op == Op::Match { s.matched as i32 } else { s.mismatch as i32 };
        match ops.last_mut() {
            Some((o, m)) if *o == op => *m += 1,
            _ => ops.push((op, 1)),
        }
    }
    Some((bs, Alignment { score, q_start: qs, q_end: n, g_start: gs, g_end: ge, ops }))
}

/// A weak J or short J fragment is called only when the junction it implies, from V's Cys104 to
/// its Trp/Phe118 codon's end, is 6 to this many nt (bnab truth sets reach 162).
const MAX_JUNCTION: i64 = 180;

/// Score (matches minus mismatches) a 3' J fragment needs for a locus: 3 + ⌈log4(the locus's
/// total J nt)⌉, plus 1 for TR loci: 8 for IGH, IGK and IGL, 9 for TRB, 10 for TRA. Chosen by
/// generated-truth accuracy (record 0029, Waves 2 and 3); the random-tail null is reported, not a
/// bound. 7 would call a J on reads that end in N/D bases resembling a J (`missing_j_is_none`).
fn short_min(j: &SegmentIndex, ks: &[usize], locus: Locus) -> i32 {
    let total: usize = ks.iter().map(|&k| j.seqs[k].len()).sum();
    let tr = matches!(locus, Locus::TRA | Locus::TRB | Locus::TRG | Locus::TRD);
    3 + (usize::BITS - total.saturating_sub(1).leading_zeros()).div_ceil(2) as i32 + i32::from(tr)
}

/// The aligned read bases of `aln` at read positions `>= min_q` as (allele index, read base);
/// inserted bases are skipped.
fn pairs_from(aln: &Alignment, q: &[u8], min_q: usize) -> Vec<(usize, u8)> {
    let (mut out, mut qi, mut gi) = (Vec::new(), aln.q_start, aln.g_start);
    for &(op, n) in &aln.ops {
        let n = n as usize;
        match op {
            Op::Match | Op::Mismatch => {
                out.extend((0..n).filter(|&x| qi + x >= min_q).map(|x| (gi + x, q[qi + x])));
                qi += n;
                gi += n;
            }
            Op::Ins => qi += n,
            Op::Del => gi += n,
        }
    }
    out
}

/// Fewest mismatches of `seq` against the J-derived bases `p` over every shift, a base the
/// shifted allele doesn't cover or a non-ACGT read base counting as a mismatch (truth's
/// sliding rule, record 0025).
fn slide(p: &[(usize, u8)], seq: &[u8]) -> u32 {
    let (Some(lo), Some(hi)) = (p.first().map(|x| x.0 as i64), p.last().map(|x| x.0 as i64)) else { return 0 };
    let hits = |k: i64| {
        p.iter().filter(|&&(u, b)| matches!(b, b'A' | b'C' | b'G' | b'T') && usize::try_from(u as i64 + k).ok().and_then(|x| seq.get(x)) == Some(&b)).count()
    };
    (p.len() - (-hi..seq.len() as i64 - lo).map(hits).max().unwrap_or(0)) as u32
}

/// Mismatch counts of every J allele of `locus` (database index, count) over `p`.
fn counts(p: &[(usize, u8)], j: &SegmentIndex, ks: &[usize]) -> Vec<(usize, u32)> {
    ks.iter().map(|&k| (j.ids[k], slide(p, &j.seqs[k]))).collect()
}

/// Assigns J to an oriented read after its V call (record 0029).
///
/// 1. Every J allele of the V primary's locus is scored (`score_semiglobal`) against the read
///    from `OVERLAP` bases before the V alignment's end.
/// 2. The `TOP` best are aligned, spurious gap pairs are collapsed, and each alignment's 5' end
///    is cut to the read's J-derived bases (the 3' end is the J's or the read's). The best cut
///    score projects. Below `MIN_SCORE` there are two candidates, each needing an implied
///    junction of 6–`MAX_JUNCTION` nt: the weak J (the best `WEAK_MISMATCH` cut of the
///    `WEAK_TOP` alleles by that scoring that holds its allele's Trp/Phe118 codon and scores
///    `WEAK_MIN_SCORE`; a heavily mutated J), and the B1 fragment (the best gapless alignment
///    of the read's 3' end to any locus J allele whose matches minus mismatches reach
///    `short_min`). The fragment goes first when its doubled score is at least the weak cut's
///    (both on the ±2 scale); the other is the fallback. Last comes the low path: the best
///    gapless cut of the `TOP` normal cuts scoring `LOW_MIN_SCORE`, junction bounded. A
///    candidate whose reported alignment (after the pick in step 3) fails its guard is
///    dropped; with none left there is no J hit (`None`).
/// 3. Every J allele of the locus is counted for mismatches over those bases at its best shift
///    (truth's sliding rule: an uncovered base is a mismatch); when the projecting cut has at most
///    `LATE_MAX_MM` differences, only bases from the latest start among the top cuts with its
///    score count (ruling B4). The record 0008 pick among the
///    fewest cuts the J-derived bases itself, and the count is repeated once over its cut: the
///    tie set is the alleles with the fewest mismatches there, the primary their record 0008
///    pick, and the reported alignment the primary's own.
pub fn assign_j(read: &Oriented, v: &VHit, alleles: &[Allele], j: &SegmentIndex, scoring: &Scoring) -> Option<JHit> {
    let from = v.aln.q_end.saturating_sub(OVERLAP);
    let q = &read.seq[from..];
    let locus = alleles[v.primary].locus;
    let ks: Vec<usize> = (0..j.ids.len()).filter(|&k| j.loci[k] == locus).collect();
    let mut ranked: Vec<(usize, i32)> = ks.iter().map(|&k| (k, score_semiglobal(q, &j.seqs[k], scoring).0)).collect();
    ranked.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
    let top: Vec<(usize, Alignment)> = ranked.iter().take(TOP).map(|&(k, _)| (k, cut(q, alleles, j, k, scoring))).collect();
    let best = |c: &mut dyn Iterator<Item = (usize, Alignment)>| c.reduce(|a, b| if b.1.score > a.1.score { b } else { a });
    // the junction a cut of allele k implies: Cys104 on V's last diagonal to Trp/Phe118's end, on
    // the cut's first diagonal; without both anchors it cannot be bounded, so no call
    let bounded = |k: usize, a: &Alignment| {
        let (Some(va), Some(ja)) = (alleles[v.primary].anchor, alleles[j.ids[k]].anchor) else { return false };
        let cys = v.aln.q_end as i64 - (v.aln.g_end as i64 - va as i64);
        let w = (from + a.q_start) as i64 + ja as i64 - a.g_start as i64;
        (6..=MAX_JUNCTION).contains(&(w + 3 - cys))
    };
    let holds_anchor = |k: usize, a: &Alignment| alleles[j.ids[k]].anchor.is_some_and(|w| a.g_start <= w && w + 3 <= a.g_end);
    let weak_s = Scoring { mismatch: WEAK_MISMATCH, ..*scoring };
    // what a non-strong call must satisfy, on the projecting cut and again on the reported one
    let guard = |path: Path, k: usize, a: &Alignment| match path {
        Path::Strong => true,
        Path::Weak => a.score >= WEAK_MIN_SCORE && holds_anchor(k, a) && bounded(k, a),
        Path::Tail => bounded(k, a),
        Path::Low => a.score >= LOW_MIN_SCORE && a.ops.iter().all(|o| matches!(o.0, Op::Match | Op::Mismatch)) && bounded(k, a),
    };
    let finish = |path: Path, mut k: usize, mut aln: Alignment| -> Option<JHit> {
        let recut = |k: usize| match path {
            Path::Strong | Path::Low => cut(q, alleles, j, k, scoring),
            Path::Weak => cut(q, alleles, j, k, &weak_s),
            Path::Tail => tail_aln(q, &j.seqs[k], scoring).map_or_else(|| cut(q, alleles, j, k, scoring), |t| t.1),
        };
        // Among the top cuts with the best score, count from the latest start: N bases that match
        // a sibling's 5' end don't decide the allele (ruling B4). Not on mutated J, whose 5'
        // bases still separate siblings.
        let cut_mm: u32 = aln.ops.iter().filter(|o| o.0 != Op::Match).map(|o| o.1).sum();
        let late = if path != Path::Strong || cut_mm > LATE_MAX_MM {
            0
        } else {
            top.iter().filter(|(_, a)| a.score == aln.score).map(|(_, a)| a.q_start).max().unwrap_or(0)
        };
        // the later start only narrows the counted bases; it never leaves none
        let count = |aln: &Alignment| {
            let p = match pairs_from(aln, q, late) {
                p if p.is_empty() => pairs_from(aln, q, 0),
                p => p,
            };
            counts(&p, j, &ks)
        };
        let mut c = count(&aln);
        let pick = primary_by_tiebreak(&fewest(&c).0, alleles);
        if pick != j.ids[k] {
            k = j.ids.binary_search(&pick).expect("tie is a J allele");
            aln = recut(k);
            c = count(&aln);
        }
        let (mut ties, mismatches) = fewest(&c);
        ties.sort_by(|&x, &y| tiebreak_cmp(&alleles[x], &alleles[y]).then(x.cmp(&y)));
        let primary = ties[0];
        if primary != j.ids[k] {
            k = j.ids.binary_search(&primary).expect("tie is a J allele");
            aln = recut(k);
        }
        // the reported alignment is the primary's re-cut: it must still pass the path's guard
        if !guard(path, k, &aln) {
            return None;
        }
        aln.q_start += from;
        aln.q_end += from;
        Some(JHit { primary, ties, aln, mismatches })
    };

    if let Some((k, a)) = best(&mut top.iter().cloned()).filter(|(_, a)| a.score >= MIN_SCORE) {
        return finish(Path::Strong, k, a);
    }
    // Weak J: a heavily mutated J, realigned with `WEAK_MISMATCH`; the `WEAK_TOP` alleles by
    // that alignment's score are cut.
    let mut wr: Vec<(usize, i32)> = ks.iter().map(|&k| (k, score_semiglobal(q, &j.seqs[k], &weak_s).0)).collect();
    wr.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
    let weak = best(&mut wr.iter().take(WEAK_TOP).map(|&(k, _)| (k, cut(q, alleles, j, k, &weak_s))).filter(|(k, a)| guard(Path::Weak, *k, a)));
    // Short J fragment at the read's 3' end (ruling B1).
    let tail = ks
        .iter()
        .filter_map(|&k| tail_aln(q, &j.seqs[k], scoring).map(|t| (k, t)))
        .reduce(|a, b| if b.1 .0 > a.1 .0 { b } else { a })
        .filter(|(k, (bs, a))| *bs >= short_min(j, &ks, locus) && bounded(*k, a));
    // The better-supported candidate goes first, on the ±2 scale (a fragment's matches minus
    // mismatches, doubled); the other is the fallback when the first fails its guard after the
    // record 0008 pick.
    let tail_first = match (&weak, &tail) {
        (Some((_, w)), Some((_, (bs, _)))) => 2 * bs >= w.score,
        _ => true,
    };
    let weak = weak.map(|(k, a)| (Path::Weak, k, a));
    let tail = tail.map(|(k, (_, a))| (Path::Tail, k, a));
    // Low path: the best gapless normal cut, last (record 0029, Wave 3).
    let low = best(&mut top.iter().filter(|&(k, a)| guard(Path::Low, *k, a)).cloned()).map(|(k, a)| (Path::Low, k, a));
    let [first, second] = if tail_first { [tail, weak] } else { [weak, tail] };
    [first, second, low].into_iter().flatten().find_map(|(path, k, a)| finish(path, k, a))
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Path {
    Strong,
    Weak,
    Tail,
    Low,
}
