use super::ties::tiebreak_cmp;
use super::v::fewest;
use super::{JHit, Oriented, SegmentIndex};
use abstar_germline::Allele;

/// The C call for one oriented read.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CHit {
    /// Database index of the primary allele (record 0008); `ties[0]`.
    pub primary: usize,
    /// Database indices of the J locus's C alleles tied for the fewest mismatches over the
    /// bases after J, primary first, then record 0008 order.
    pub ties: Vec<usize>,
    /// Read bases `start..end` are the primary's bases `0..g_end`.
    pub start: usize,
    pub end: usize,
    pub g_end: usize,
}

/// Fewer bases than this after J give no C call.
const MIN_C_LEN: usize = 10;

/// Assigns C by a prefix comparison from the J alignment's end: every C allele of the J
/// primary's locus is compared base for base with the read from there (`allele[u]` with read
/// base `start + u`), and the alleles with the fewest mismatches tie. As in truth (record 0025),
/// a base the allele doesn't cover or a non-ACGT read base is a mismatch. `None` when fewer
/// than `MIN_C_LEN` bases follow J.
pub fn assign_c(read: &Oriented, j: &JHit, alleles: &[Allele], c: &SegmentIndex) -> Option<CHit> {
    let start = j.aln.q_end;
    let tail = read.seq.get(start..).filter(|t| t.len() >= MIN_C_LEN)?;
    let locus = alleles[j.primary].locus;
    let counts: Vec<(usize, u32)> = (0..c.ids.len())
        .filter(|&k| c.loci[k] == locus)
        .map(|k| {
            let g = &c.seqs[k];
            let mm = tail.iter().enumerate().filter(|&(u, &b)| !matches!(b, b'A' | b'C' | b'G' | b'T') || g.get(u) != Some(&b)).count();
            (c.ids[k], mm as u32)
        })
        .collect();
    let (mut ties, _) = fewest(&counts);
    if ties.is_empty() {
        return None;
    }
    ties.sort_by(|&x, &y| tiebreak_cmp(&alleles[x], &alleles[y]).then(x.cmp(&y)));
    let primary = ties[0];
    let g_end = tail.len().min(c.seqs[c.ids.binary_search(&primary).expect("tie is a C allele")].len());
    Some(CHit { primary, ties, start, end: start + g_end, g_end })
}
