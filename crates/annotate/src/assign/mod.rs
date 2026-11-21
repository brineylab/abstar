//! Orientation and gene assignment (record 0029).
//!
//! Allele indices in results (`VHit::primary`, `VHit::ties`) index the database's allele slice.
//! Coordinates are 0-based half-open on the oriented read.
mod orient;
mod ties;
mod v;

pub use orient::{orient, Oriented};
pub use ties::primary_by_tiebreak;
pub use v::{assign_v, VHit, MIN_V_SCORE};

use crate::align::KmerIndex;
use abstar_germline::{Allele, Segment};

/// k-mer length of the prefilter indexes (record 0029).
pub const K: usize = 7;

/// One segment's alleles, prepared once per run: their database indices, uppercase ungapped and
/// IMGT-gapped sequences, and k-mer index (entry `i` is allele `ids[i]` throughout).
pub struct SegmentIndex {
    /// Database indices, ascending.
    pub ids: Vec<usize>,
    pub seqs: Vec<Vec<u8>>,
    pub gapped: Vec<Vec<u8>>,
    pub index: KmerIndex,
}

impl SegmentIndex {
    pub fn new(alleles: &[Allele], segment: Segment) -> SegmentIndex {
        let ids: Vec<usize> = (0..alleles.len()).filter(|&i| alleles[i].segment == segment).collect();
        let seqs: Vec<Vec<u8>> = ids.iter().map(|&i| alleles[i].sequence().to_ascii_uppercase().into_bytes()).collect();
        let gapped: Vec<Vec<u8>> = ids.iter().map(|&i| alleles[i].sequence_gapped.to_ascii_uppercase().into_bytes()).collect();
        let index = KmerIndex::new(&seqs.iter().map(|s| s.as_slice()).collect::<Vec<_>>(), K);
        SegmentIndex { ids, seqs, gapped, index }
    }
}
