use super::ties::tiebreak_cmp;
use super::v::fewest;
use super::{JHit, Oriented, SegmentIndex};
use crate::align::{cut_query_end, Scoring};
use abstar_germline::Allele;

/// The C call for one oriented read.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CHit {
    /// Database index of the primary allele (record 0008); `ties[0]`.
    pub primary: usize,
    /// Database indices of the J locus's C alleles tied for the fewest mismatches over the
    /// C bases at the chosen placement, primary first, then record 0008 order.
    pub ties: Vec<usize>,
    /// Read bases `start..end` are the primary's bases `g_start..g_end` (gapless).
    pub start: usize,
    pub end: usize,
    pub g_start: usize,
    pub g_end: usize,
}

/// Fewer C bases than this give no C call.
const MIN_C_LEN: usize = 10;
/// The placements tried: the C allele's first base this many bases either side of the J
/// alignment's last base (record 0030, 2026-11-21 note).
const WINDOW: usize = 3;
/// No C call below this identity over the C bases. On random tails after J (12,000 generated
/// IGH, IGK, IGL, TRA and TRB reads, tails of 10-100 nt) the best placement reaches it on 5 reads,
/// all with 10-12-nt tails; on the real BCR corpus 28 of 49,932 reads fall below it (record 0030).
pub const MIN_C_IDENTITY: f64 = 0.8;

fn acgt(b: u8) -> bool {
    matches!(b, b'A' | b'C' | b'G' | b'T')
}

/// Assigns C after J (record 0030), before the junction step. Built-in C alleles start with the
/// J's last base (the codon the J and C exons share), so the C allele's base 0 is expected on
/// J's last read base.
///
/// 1. Placement: every C allele of the J primary's locus is placed with its base 0 at each read
///    position within `WINDOW` of J's last base, and compared base for base over the read bases
///    every placement covers; the fewest mismatches wins, ties toward the expected placement,
///    then the earlier one.
/// 2. At that placement the C bases start after the base 0 that J supplies (`g_start` 1), or at
///    base 0 when it lies after J's end (`g_start` 0). When they start before J's alignment end,
///    J is cut back to them (`j.aln`, under `s`), so C and J never overlap. As in truth (record
///    0025), a base the allele doesn't cover or a non-ACGT read base is a mismatch; the alleles
///    with the fewest tie, and the primary is the record 0008 pick.
/// 3. `None`, with J untouched, when fewer than `MIN_C_LEN` C bases remain or the primary's
///    identity over them is below `MIN_C_IDENTITY` (no evidence for C: a random tail, an adapter).
pub fn assign_c(read: &Oriented, j: &mut JHit, alleles: &[Allele], c: &SegmentIndex, s: &Scoring) -> Option<CHit> {
    let (q, e) = (&read.seq, j.aln.q_end);
    if q.len() < e + MIN_C_LEN {
        return None;
    }
    let locus = alleles[j.primary].locus;
    let ks: Vec<usize> = (0..c.ids.len()).filter(|&k| c.loci[k] == locus).collect();
    // mismatches of allele `k` with its base 0 on read position `p`, over read bases `from..`
    let mm = |k: usize, p: usize, from: usize| {
        let g = &c.seqs[k];
        (from..q.len()).filter(|&i| !acgt(q[i]) || g.get(i - p) != Some(&q[i])).count() as u32
    };
    let expected = e - 1;
    let common = e + WINDOW - 1;
    let p = (expected.saturating_sub(WINDOW)..=expected + WINDOW)
        .flat_map(|p| ks.iter().map(move |&k| (p, k)))
        .min_by_key(|&(p, k)| (mm(k, p, common), p.abs_diff(expected), p))
        .map(|(p, _)| p)?;
    let start = if p >= e { p } else { p + 1 };
    let counts: Vec<(usize, u32)> = ks.iter().map(|&k| (c.ids[k], mm(k, p, start))).collect();
    let (mut ties, _) = fewest(&counts);
    ties.sort_by(|&x, &y| tiebreak_cmp(&alleles[x], &alleles[y]).then(x.cmp(&y)));
    let primary = *ties.first()?;
    let g = &c.seqs[c.ids.binary_search(&primary).expect("tie is a C allele")];
    let (g_start, g_end) = (start - p, (q.len() - p).min(g.len()));
    if g_end < g_start + MIN_C_LEN {
        return None;
    }
    let end = p + g_end;
    let matched = (start..end).filter(|&i| acgt(q[i]) && q[i] == g[i - p]).count();
    if (matched as f64) < MIN_C_IDENTITY * (end - start) as f64 {
        return None;
    }
    if start <= j.aln.q_start {
        return None;
    }
    if start < e {
        j.aln = cut_query_end(&j.aln, start, s);
    }
    Some(CHit { primary, ties, start, end, g_start, g_end })
}
