//! Amino-acid V and J assignment (spec §Assignment).
//!
//! Allele indices index the database's allele slice; coordinates are 0-based half-open on the
//! full query.
use super::germline::AaSegment;
use crate::align::{align_protein, blosum62, Alignment, Op, ProteinScoring};
use crate::assign::tiebreak_cmp;
use abstar_germline::{Allele, Locus};

pub const MIN_V_SCORE_AA: i32 = 50; // calibrated 2026-11-23: aa gate (test_accuracy.py SETS), true v/j/junction accuracy
pub const MIN_J_SCORE_AA: i32 = 15; // calibrated 2026-11-23: aa gate (test_accuracy.py SETS), true v/j/junction accuracy

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AaHit {
    /// Database index of the primary allele (record 0008); `ties[0]`.
    pub primary: usize,
    /// Database indices, primary first, then record 0008 order.
    pub ties: Vec<usize>,
    /// Index into `AaSegment.groups` of the primary's group.
    pub group: usize,
    /// Query vs the group's `seq`, trimmed to the segment-derived residues.
    pub aln: Alignment,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AaCalls {
    pub v: AaHit,
    pub j: Option<AaHit>,
    /// Query residues `..v_end` are V's, `j_start..` J's; between them non-germline.
    pub v_end: usize,
    pub j_start: usize,
}

/// One scored alignment column: a residue pair (BLOSUM62), or a whole gap run (affine cost).
/// `q`/`g` are the query/germline positions it starts at.
#[derive(Clone, Copy)]
struct Col {
    op: Op,
    n: u32,
    q: usize,
    g: usize,
    score: i32,
}

fn columns(aln: &Alignment, q: &[u8], g: &[u8], s: &ProteinScoring) -> Vec<Col> {
    let (mut qi, mut gi) = (aln.q_start, aln.g_start);
    let mut cols = Vec::new();
    for &(op, n) in &aln.ops {
        let k = n as usize;
        match op {
            Op::Match | Op::Mismatch => {
                for d in 0..k {
                    cols.push(Col { op, n: 1, q: qi + d, g: gi + d, score: blosum62(q[qi + d], g[gi + d]) });
                }
                (qi, gi) = (qi + k, gi + k);
            }
            Op::Ins | Op::Del => {
                cols.push(Col { op, n, q: qi, g: gi, score: s.gap_open + (n as i32 - 1) * s.gap_extend });
                match op {
                    Op::Ins => qi += k,
                    _ => gi += k,
                }
            }
        }
    }
    cols
}

/// `aln` cut to its maximum-scoring contiguous stretch of columns: equal scores take the
/// longest, then leftmost. With `keep_end` the stretch must end at the last column and equal
/// scores take the shortest.
fn trim_protein(aln: &Alignment, q: &[u8], g: &[u8], s: &ProteinScoring, keep_end: bool) -> Alignment {
    let cols = columns(aln, q, g, s);
    if cols.is_empty() {
        return aln.clone();
    }
    // Kadane over prefix sums: for each end, the start at the earliest (keep_end: latest)
    // minimal prefix.
    let (mut prefix, mut low, mut low_at) = (0, 0, 0);
    let (mut best, mut start, mut end) = (i32::MIN, 0, 0);
    for (c, col) in cols.iter().enumerate() {
        prefix += col.score;
        let (sum, len) = (prefix - low, c + 1 - low_at);
        if (!keep_end || c + 1 == cols.len()) && (sum > best || (sum == best && len > end - start)) {
            (best, start, end) = (sum, low_at, c + 1);
        }
        if prefix < low || (keep_end && prefix == low) {
            (low, low_at) = (prefix, c + 1);
        }
    }
    let mut ops: Vec<(Op, u32)> = Vec::new();
    for c in &cols[start..end] {
        match ops.last_mut() {
            Some((last, m)) if *last == c.op => *m += c.n,
            _ => ops.push((c.op, c.n)),
        }
    }
    let last = cols[end - 1];
    let (q_end, g_end) = match last.op {
        Op::Match | Op::Mismatch => (last.q + 1, last.g + 1),
        Op::Ins => (last.q + last.n as usize, last.g),
        Op::Del => (last.q, last.g + last.n as usize),
    };
    Alignment { score: best.max(0), q_start: cols[start].q, q_end, g_start: cols[start].g, g_end, ops }
}

/// Mismatched and gapped columns of `aln` inside query `lo..hi` (a deletion counts when it
/// falls between two residues of the range), plus residues of the range `aln` doesn't cover.
fn differences(aln: &Alignment, lo: usize, hi: usize) -> usize {
    let overlap = |a: usize, b: usize| b.min(hi).saturating_sub(a.max(lo));
    let mut d = (hi - lo) - overlap(aln.q_start, aln.q_end);
    let mut qi = aln.q_start;
    for &(op, n) in &aln.ops {
        let k = n as usize;
        match op {
            Op::Match => qi += k,
            Op::Mismatch | Op::Ins => {
                d += overlap(qi, qi + k);
                qi += k;
            }
            Op::Del if lo < qi && qi < hi => d += k,
            Op::Del => {}
        }
    }
    d
}

/// The query residue aligned to germline position `g`, if any.
fn query_at(aln: &Alignment, g: usize) -> Option<usize> {
    let (mut qi, mut gi) = (aln.q_start, aln.g_start);
    for &(op, n) in &aln.ops {
        let k = n as usize;
        match op {
            Op::Match | Op::Mismatch if (gi..gi + k).contains(&g) => return Some(qi + g - gi),
            Op::Match | Op::Mismatch => (qi, gi) = (qi + k, gi + k),
            Op::Ins => qi += k,
            Op::Del => gi += k,
        }
    }
    None
}

/// Aligns `q[from..]` against every group (of `locus`, if given). The best must reach `min`;
/// ties are the alleles of the best's-locus groups with the fewest `differences` over the best's
/// trimmed query range.
#[allow(clippy::too_many_arguments)]
fn hit(q: &[u8], from: usize, seg: &AaSegment, locus: Option<Locus>, alleles: &[Allele], s: &ProteinScoring, min: i32, keep_end: bool) -> Option<AaHit> {
    let alns: Vec<(usize, Alignment)> = seg
        .groups
        .iter()
        .enumerate()
        .filter(|(_, g)| locus.is_none_or(|l| g.locus == l))
        .map(|(k, g)| {
            let mut a = align_protein(&q[from..], &g.seq, s);
            (a.q_start, a.q_end) = (a.q_start + from, a.q_end + from);
            (k, a)
        })
        .collect();
    // Highest score, first group on equal scores.
    let (best, ba) = alns.iter().rev().max_by_key(|(_, a)| a.score)?;
    if ba.score < min {
        return None;
    }
    let locus = seg.groups[*best].locus;
    let t = trim_protein(ba, q, &seg.groups[*best].seq, s, keep_end);
    let diffs: Vec<(usize, usize)> =
        alns.iter().filter(|(k, _)| seg.groups[*k].locus == locus).map(|(k, a)| (*k, differences(a, t.q_start, t.q_end))).collect();
    let fewest = diffs.iter().map(|x| x.1).min()?;
    let mut ties: Vec<usize> = diffs.iter().filter(|x| x.1 == fewest).flat_map(|x| seg.groups[x.0].alleles.iter().copied()).collect();
    // ties[0] is `primary_by_tiebreak(&ties, alleles)`.
    ties.sort_by(|&x, &y| tiebreak_cmp(&alleles[x], &alleles[y]).then(x.cmp(&y)));
    let primary = ties[0];
    let (group, ga) = alns.iter().find(|(k, _)| seg.groups[*k].alleles.contains(&primary))?;
    let aln = trim_protein(ga, q, &seg.groups[*group].seq, s, keep_end);
    Some(AaHit { primary, ties, group: *group, aln })
}

/// V and J for one protein query (spec §Assignment). `None` = no V hit.
pub fn assign_aa(q: &[u8], v: &AaSegment, j: &AaSegment, alleles: &[Allele], s: &ProteinScoring) -> Option<AaCalls> {
    let vh = hit(q, 0, v, None, alleles, s, MIN_V_SCORE_AA, false)?;
    let vg = &v.groups[vh.group];
    // J is searched from V's end or Cys104, whichever is earlier, so V can be cut back.
    let cys = vg.anchor.and_then(|a| query_at(&vh.aln, a));
    let from = cys.map_or(vh.aln.q_end, |c| c.min(vh.aln.q_end));
    let jh = hit(q, from, j, Some(vg.locus), alleles, s, MIN_J_SCORE_AA, true);
    let (mut v_end, mut j_start) = (vh.aln.q_end, jh.as_ref().map_or(q.len(), |h| h.aln.q_start));
    if let Some(jh) = jh.as_ref().filter(|_| j_start < v_end) {
        let vc = columns(&vh.aln, q, &vg.seq, s);
        let jc = columns(&jh.aln, q, &j.groups[jh.group].seq, s);
        let score = |c: usize| -> i32 { vc.iter().filter(|x| x.q < c).map(|x| x.score).sum::<i32>() + jc.iter().filter(|x| x.q >= c).map(|x| x.score).sum::<i32>() };
        // max_by_key keeps the last maximum; reversing makes it the earliest cut.
        let c = (j_start..=v_end).rev().max_by_key(|&c| score(c)).expect("non-empty range");
        (v_end, j_start) = (c, c);
    }
    Some(AaCalls { v: vh, j: jh, v_end, j_start })
}
