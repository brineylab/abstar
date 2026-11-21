use super::ties::tiebreak_cmp;
use super::v::{fewest, trim};
use super::{primary_by_tiebreak, Oriented, SegmentIndex, VHit};
use crate::align::{align_semiglobal, collapse_gap_pairs, place_gaps, score_semiglobal, Alignment, Op, Scoring};
use abstar_germline::Allele;

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
/// ... if its score is at least this; else there is no J hit.
const MIN_SCORE: i32 = 30;

/// Aligns `q` to J allele `k` of `j`, collapses spurious insertion/deletion pairs, cuts the
/// 5' end to the read's J-derived bases (`trim`) and places gaps outside the CDR3 side of
/// Trp/Phe118 (the generator's J indel zone starts 5 nt past the anchor's first base).
fn cut(q: &[u8], alleles: &[Allele], j: &SegmentIndex, k: usize, s: &Scoring) -> Alignment {
    let g = &j.seqs[k];
    let raw = align_semiglobal(q, g, s);
    let cut = trim(&collapse_gap_pairs(&raw, q, g, s), s, true);
    place_gaps(&cut, q, g, alleles[j.ids[k]].anchor.map(|w| 0..w + 5), s)
}

/// The aligned read bases of `aln` as (allele index, read base); inserted bases are skipped.
fn pairs(aln: &Alignment, q: &[u8]) -> Vec<(usize, u8)> {
    let (mut out, mut qi, mut gi) = (Vec::new(), aln.q_start, aln.g_start);
    for &(op, n) in &aln.ops {
        let n = n as usize;
        match op {
            Op::Match | Op::Mismatch => {
                out.extend((0..n).map(|x| (gi + x, q[qi + x])));
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
///    score projects; below `MIN_SCORE` there is no J hit (`None`).
/// 3. Every J allele of the locus is counted for mismatches over those bases at its best shift
///    (truth's sliding rule: an uncovered base is a mismatch). The record 0008 pick among the
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
    let (mut k, mut aln) = ranked
        .iter()
        .take(TOP)
        .map(|&(k, _)| (k, cut(q, alleles, j, k, scoring)))
        .reduce(|a, b| if b.1.score > a.1.score { b } else { a })
        .filter(|(_, a)| a.score >= MIN_SCORE)?;
    let mut c = counts(&pairs(&aln, q), j, &ks);
    let pick = primary_by_tiebreak(&fewest(&c).0, alleles);
    if pick != j.ids[k] {
        k = j.ids.binary_search(&pick).expect("tie is a J allele");
        aln = cut(q, alleles, j, k, scoring);
        c = counts(&pairs(&aln, q), j, &ks);
    }
    let (mut ties, mismatches) = fewest(&c);
    ties.sort_by(|&x, &y| tiebreak_cmp(&alleles[x], &alleles[y]).then(x.cmp(&y)));
    let primary = ties[0];
    if primary != j.ids[k] {
        aln = cut(q, alleles, j, j.ids.binary_search(&primary).expect("tie is a J allele"), scoring);
    }
    aln.q_start += from;
    aln.q_end += from;
    Some(JHit { primary, ties, aln, mismatches })
}
