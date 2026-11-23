//! Orientation and gene assignment (record 0029), junction resolution (records 0030, 0031)
//! and C assignment.
//!
//! Allele indices in results (`primary`, `ties` of `VHit`, `JHit`, `DHit` and `CHit`) index the
//! database's allele slice.
//! Coordinates are 0-based half-open on the oriented read.
mod c;
mod j;
mod junction;
mod orient;
mod ties;
mod v;

pub use c::{assign_c, CHit};
pub use j::{assign_j, JHit};
pub use junction::{resolve_junction, DHit, Junction};
pub use orient::{orient, reverse_complement, Oriented};
pub use ties::primary_by_tiebreak;
pub use v::{assign_v, VHit, MIN_V_SCORE};

use crate::align::{KmerIndex, Scoring};
use abstar_germline::{Allele, Locus, Segment};

/// k-mer length of the prefilter indexes (record 0029).
pub const K: usize = 7;

/// One segment's alleles, prepared once per run: their database indices, uppercase ungapped
/// sequences and their reverse complements, IMGT-gapped sequences, and k-mer index (entry `i` is
/// allele `ids[i]` throughout).
pub struct SegmentIndex {
    /// Database indices, ascending.
    pub ids: Vec<usize>,
    pub seqs: Vec<Vec<u8>>,
    pub rc: Vec<Vec<u8>>,
    pub gapped: Vec<Vec<u8>>,
    pub loci: Vec<Locus>,
    pub index: KmerIndex,
}

impl SegmentIndex {
    /// Every allele of `segment` except known non-functional ones (record 0038).
    pub fn new(alleles: &[Allele], segment: Segment) -> SegmentIndex {
        SegmentIndex::for_loci(alleles, segment, None, false)
    }

    /// `new` restricted to alleles of `loci` (all loci when `None`), with known non-functional
    /// alleles (P, NF) only when `nonfunctional`.
    pub fn for_loci(alleles: &[Allele], segment: Segment, loci: Option<&[Locus]>, nonfunctional: bool) -> SegmentIndex {
        let keep = |a: &Allele| {
            a.segment == segment && loci.is_none_or(|l| l.contains(&a.locus)) && (nonfunctional || !a.functionality.is_nonfunctional())
        };
        let ids: Vec<usize> = (0..alleles.len()).filter(|&i| keep(&alleles[i])).collect();
        let seqs: Vec<Vec<u8>> = ids.iter().map(|&i| alleles[i].sequence().to_ascii_uppercase().into_bytes()).collect();
        let gapped: Vec<Vec<u8>> = ids.iter().map(|&i| alleles[i].sequence_gapped.to_ascii_uppercase().into_bytes()).collect();
        let index = KmerIndex::new(&seqs.iter().map(|s| s.as_slice()).collect::<Vec<_>>(), K);
        let loci = ids.iter().map(|&i| alleles[i].locus).collect();
        let rc = seqs.iter().map(|s| reverse_complement(s)).collect();
        SegmentIndex { ids, seqs, rc, gapped, loci, index }
    }
}

/// Orients the read and assigns V (record 0029). When the chosen strand gives no V hit, the
/// other strand is tried; a read and its reverse complement try the same strands in the same
/// order. `None` means no V hit (`no_v_hit`).
pub fn call_v(read: &[u8], alleles: &[Allele], v: &SegmentIndex, j: &SegmentIndex, scoring: &Scoring) -> Option<(Oriented, VHit)> {
    let o = orient(read, &v.index, &j.index)?;
    if let Some(h) = assign_v(&o, alleles, v, j, scoring) {
        return Some((o, h));
    }
    let flipped = Oriented { seq: reverse_complement(&o.seq), rev_comp: !o.rev_comp };
    assign_v(&flipped, alleles, v, j, scoring).map(|h| (flipped, h))
}
