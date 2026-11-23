//! IMGT regions, the junction, translation, alignment strings, productivity and the assembled
//! sequences of an annotated read (record 0032), and its per-segment mutations (record 0033) and
//! masks (record 0034), from its V, J, junction and C calls.
//!
//! Coordinates are 0-based half-open on the oriented read. Without a J call (`j` is `None`)
//! V ends at its alignment's end and the V(D)J region runs to the read's end (ruling A8).
mod alignment;
mod masks;
mod mutations;
mod productivity;
mod regions;
mod translate;

pub use alignment::{alignments, Alignments};
pub use masks::{masks, Masks};
pub use mutations::{mutations, segments, SegmentMutations};
pub(crate) use mutations::segment;
pub use productivity::{productivity, Productivity};
pub use regions::{regions, Regions};
pub use translate::{region_aa, translate};

use crate::align::{Alignment, Op};
use crate::assign::{JHit, Junction, Oriented, VHit};
use abstar_germline::Allele;
use std::ops::Range;

/// One alignment column: (read position, ungapped germline index). An inserted base has no
/// germline index, a deleted base no read position.
pub type Col = (Option<usize>, Option<usize>);

/// The columns of `aln` over read positions `span`. Positions of `span` outside the alignment
/// continue its first or last diagonal, as the junction step extends V and J (record 0030). A
/// deletion is kept only between two positions of `span`.
pub(crate) fn columns(aln: &Alignment, span: Range<usize>) -> Vec<Col> {
    let mut out: Vec<Col> = (span.start..aln.q_start.min(span.end))
        .filter_map(|i| (aln.g_start + i).checked_sub(aln.q_start).map(|u| (Some(i), Some(u))))
        .collect();
    let (mut qi, mut gi) = (aln.q_start, aln.g_start);
    for &(op, n) in &aln.ops {
        for _ in 0..n {
            match op {
                Op::Match | Op::Mismatch => {
                    if span.contains(&qi) {
                        out.push((Some(qi), Some(gi)));
                    }
                    (qi, gi) = (qi + 1, gi + 1);
                }
                Op::Ins => {
                    if span.contains(&qi) {
                        out.push((Some(qi), None));
                    }
                    qi += 1;
                }
                Op::Del => {
                    if qi > span.start && span.contains(&qi) {
                        out.push((None, Some(gi)));
                    }
                    gi += 1;
                }
            }
        }
    }
    out.extend((aln.q_end.max(span.start)..span.end).map(|i| (Some(i), Some(aln.g_end + i - aln.q_end))));
    out
}

/// Where each segment lies on the read: V's and J's columns, V's end, J's start and the end
/// of the V(D)J region (J's end, or the read's without J).
struct Layout {
    v: Vec<Col>,
    j: Vec<Col>,
    v_end: usize,
    j_start: usize,
    end: usize,
}

impl Layout {
    fn new(read: &Oriented, v: &VHit, j: Option<(&JHit, &Junction)>) -> Layout {
        let len = read.seq.len();
        match j {
            Some((j, jn)) => Layout {
                v: columns(&v.aln, v.aln.q_start..jn.v_end),
                j: columns(&j.aln, jn.j_start..j.aln.q_end),
                v_end: jn.v_end,
                j_start: jn.j_start,
                end: j.aln.q_end,
            },
            None => Layout { v: columns(&v.aln, v.aln.q_start..v.aln.q_end), j: Vec::new(), v_end: v.aln.q_end, j_start: len, end: len },
        }
    }
}

/// Ungapped index -> IMGT gapped index (0-based) of a V allele.
fn imgt_of(a: &Allele) -> Vec<usize> {
    a.sequence_gapped.bytes().enumerate().filter(|&(_, b)| b != b'.').map(|(g, _)| g).collect()
}
