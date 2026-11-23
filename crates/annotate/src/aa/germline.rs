//! V and J germline alleles translated to protein once per run, grouped by identical translation.
use crate::derive::translate;
use abstar_germline::{Allele, Locus, Segment};
use std::collections::HashMap;

/// Alleles of one locus that share an exact translation.
pub struct AaGroup {
    /// Database indices of every allele with this exact translation, ascending.
    pub alleles: Vec<usize>,
    pub locus: Locus,
    /// Ungapped residues.
    pub seq: Vec<u8>,
    /// V: ungapped residue index -> IMGT gapped residue index (0-based); J: identity.
    pub imgt: Vec<usize>,
    /// V: ungapped index of Cys104 (IMGT 104); J: index of Trp/Phe118. None if the allele has no anchor.
    pub anchor: Option<usize>,
}

pub struct AaSegment {
    pub groups: Vec<AaGroup>,
}

impl AaSegment {
    /// Filtering matches `SegmentIndex::for_loci`. V alleles that are not IMGT-gapped and J
    /// alleles without an anchor are excluded. Groups are in order of first allele.
    pub fn new(alleles: &[Allele], segment: Segment, loci: Option<&[Locus]>, nonfunctional: bool) -> AaSegment {
        let mut groups: Vec<AaGroup> = Vec::new();
        let mut by_key: HashMap<(Locus, Vec<u8>), usize> = HashMap::new();
        for (i, a) in alleles.iter().enumerate() {
            if !(a.segment == segment && loci.is_none_or(|l| l.contains(&a.locus)) && (nonfunctional || !a.functionality.is_nonfunctional())) {
                continue;
            }
            let (gapped, anchor) = match segment {
                Segment::V => {
                    let Some(g) = translate_v(&a.sequence_gapped) else { continue };
                    let anchor = a.anchor.and_then(|_| g.get(..103)).map(|p| p.iter().filter(|&&r| r != b'.').count());
                    (g, anchor)
                }
                Segment::J => {
                    let Some(anchor) = a.anchor else { continue };
                    let (g, idx) = translate_j(&a.sequence(), anchor);
                    (g, Some(idx))
                }
                _ => continue,
            };
            let gi = *by_key.entry((a.locus, gapped.clone())).or_insert_with(|| {
                let mut seq = Vec::new();
                let mut imgt = Vec::new();
                for (k, &r) in gapped.iter().enumerate() {
                    if r != b'.' || segment == Segment::J {
                        seq.push(r);
                        imgt.push(k);
                    }
                }
                groups.push(AaGroup { alleles: Vec::new(), locus: a.locus, seq, imgt, anchor });
                groups.len() - 1
            });
            groups[gi].alleles.push(i);
        }
        AaSegment { groups }
    }
}

/// IMGT-gapped residues, `.` for a fully gapped codon. `None` if there is no gap or a codon
/// mixes `.` and bases.
pub fn translate_v(gapped_nt: &str) -> Option<Vec<u8>> {
    if !gapped_nt.contains('.') {
        return None;
    }
    gapped_nt
        .as_bytes()
        .as_chunks::<3>()
        .0
        .iter()
        .map(|c| match c.iter().filter(|&&b| b == b'.').count() {
            0 => Some(translate(c, 0).as_bytes()[0]),
            3 => Some(b'.'),
            _ => None,
        })
        .collect()
}

/// Residues from `anchor % 3`, and the index of the anchor residue (`anchor / 3`).
pub fn translate_j(ungapped_nt: &str, anchor: usize) -> (Vec<u8>, usize) {
    (translate(ungapped_nt.as_bytes(), anchor % 3).into_bytes(), anchor / 3)
}
