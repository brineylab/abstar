use super::mutations::segments;
use super::{translate, Regions};
use crate::assign::{JHit, Junction, Oriented, VHit};
use abstar_germline::Allele;

/// The per-base and per-codon masks over the V(D)J region (record 0034).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Masks {
    /// `0` FWR, `1`/`2`/`3` CDR1/2/3; `None` when CDR3's ends aren't known (CDR3 is null).
    pub cdr: Option<String>,
    pub cdr_aa: Option<String>,
    /// `V`, `N`, `D` or `J`; P bases are `N`, an inserted base takes its segment.
    pub segment: String,
    pub segment_aa: String,
    /// `1` for N/P bases, inserted bases and mismatches, else `0`.
    pub nongermline: String,
    /// `1` when the codon's amino acid differs from the germline's or a base of it has no
    /// germline base (N/P/inserted).
    pub nongermline_aa: String,
}

/// The masks of a read (record 0034). Nt masks have one character per base of `sequence_vdj`
/// (`regions.vdj`), never a C base; aa masks one per codon of `sequence_vdj_aa` (complete V-frame
/// codons from `regions.frame`), each taking its codon's 3rd base's segment and region.
pub fn masks(read: &Oriented, v: &VHit, j: Option<(&JHit, &Junction)>, alleles: &[Allele], regions: &Regions) -> Masks {
    let q = &read.seq;
    let (i0, end) = regions.vdj;
    let (mut seg, mut germ) = (vec![b'N'; end - i0], vec![None; end - i0]);
    for (k, s) in segments(read, v, j, None, alleles).iter().take(3).enumerate() {
        let Some((cols, g)) = s else { continue };
        for &(i, u) in cols {
            if let Some(i) = i.filter(|i| (i0..end).contains(i)) {
                (seg[i - i0], germ[i - i0]) = (b"VDJ"[k], u.map(|u| g[u]));
            }
        }
    }
    let ng: Vec<u8> = (i0..end).map(|i| if germ[i - i0].is_some_and(|g: u8| g.eq_ignore_ascii_case(&q[i])) { b'0' } else { b'1' }).collect();
    let cdr: Option<Vec<u8>> = regions.cdr3_known.then(|| {
        regions.labels[i0..end]
            .iter()
            .map(|x| match x {
                Some(1) => b'1',
                Some(3) => b'2',
                Some(5) => b'3',
                _ => b'0',
            })
            .collect()
    });
    let codons: Vec<usize> = (regions.frame..).step_by(3).take_while(|p| p + 3 <= end).collect();
    let text = |b: &[u8]| String::from_utf8_lossy(b).into_owned();
    let third = |m: &[u8]| codons.iter().map(|p| m[p + 2 - i0] as char).collect::<String>();
    let ng_aa: String = codons
        .iter()
        .map(|&p| match (p..p + 3).map(|i| germ[i - i0]).collect::<Option<Vec<u8>>>() {
            Some(g) if translate(&g, 0) == translate(&q[p..p + 3], 0) => '0',
            _ => '1',
        })
        .collect();
    Masks {
        cdr_aa: cdr.as_deref().map(third),
        cdr: cdr.as_deref().map(text),
        segment_aa: third(&seg),
        segment: text(&seg),
        nongermline_aa: ng_aa,
        nongermline: text(&ng),
    }
}
