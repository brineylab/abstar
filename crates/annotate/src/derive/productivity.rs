use super::{translate, Regions};
use crate::assign::Oriented;

/// AIRR productivity (record 0032), defined as the generator defines it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Productivity {
    /// `vj_in_frame`, no stop codon, and the junction keeps a conserved Cys104 (TGT/TGC) and
    /// Trp/Phe118 (TGG/TTT/TTC) codon.
    pub productive: bool,
    /// There is a junction, its length is a multiple of 3, and its Trp/Phe118 codon is in the V
    /// reading frame.
    pub vj_in_frame: bool,
    /// A stop codon in the V reading frame over the V(D)J region (`Regions::vdj`, no C).
    pub stop_codon: bool,
}

pub fn productivity(read: &Oriented, regions: &Regions) -> Productivity {
    let q = &read.seq;
    let stop_codon = translate(&q[..regions.vdj.1], regions.frame).contains('*');
    let (vj_in_frame, conserved) = regions.junction.map_or((false, false), |(s, e)| {
        (
            (e - s) % 3 == 0 && (e - 3) % 3 == regions.frame % 3,
            matches!(&q[s..s + 3], b"TGT" | b"TGC") && matches!(&q[e - 3..e], b"TGG" | b"TTT" | b"TTC"),
        )
    });
    Productivity { productive: vj_in_frame && !stop_codon && conserved, vj_in_frame, stop_codon }
}
