//! Amino-acid regions, junction, mutations, masks and alignment strings (spec §Regions and
//! junction, §Fields), the residue-level forms of `crate::derive`'s nucleotide rules.
//!
//! Coordinates are 0-based half-open on the query. Without J, V ends at its alignment's end and
//! the V(D)J region runs to the query's end.
use super::assign::AaCalls;
use super::germline::{AaGroup, AaSegment};
use crate::derive::{columns, segment, Col, SegmentMutations};

/// Exclusive IMGT gapped residue ends (0-based) of FWR1, CDR1, FWR2, CDR2 and FWR3, which ends
/// with Cys104.
const V_ENDS: [usize; 5] = [26, 38, 55, 65, 104];
const FWR3: usize = 4;
const CDR3: usize = 5;
const FWR4: usize = 6;
/// Furthest Cys104 is projected past V's last aligned residue (ruling B2's 50 nt, in residues).
const MAX_PROJECTION: usize = 16;
/// Most residues looked back from J's last-diagonal Trp/Phe118 projection (#162).
const MAX_TRP_LOOKBACK: usize = 5;

/// IMGT regions and junction of a protein query. A region is `None` when no residue of it is in
/// the query, and empty (`(p, p)`) when it lies inside the query with every residue deleted.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AaRegions {
    pub fwr1: Option<(usize, usize)>,
    pub cdr1: Option<(usize, usize)>,
    pub fwr2: Option<(usize, usize)>,
    pub cdr2: Option<(usize, usize)>,
    pub fwr3: Option<(usize, usize)>,
    /// `None` unless both Cys104 and Trp/Phe118 are found.
    pub cdr3: Option<(usize, usize)>,
    pub fwr4: Option<(usize, usize)>,
    /// Cys104 through Trp/Phe118; `None` unless both are found.
    pub junction: Option<(usize, usize)>,
    /// `sequence_vdj_aa`: the first V residue through J's end, or the query's end without J.
    pub vdj: (usize, usize),
    /// Region index (0..7, FWR1..FWR4) of each query residue.
    pub labels: Vec<Option<usize>>,
}

/// The masks over `AaRegions::vdj`, one character per residue.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AaMasks {
    /// `0` FWR, `1`/`2`/`3` CDR1/2/3; `None` when CDR3 is.
    pub cdr: Option<String>,
    /// `V`, `N` or `J`; an inserted residue takes its segment.
    pub segment: String,
    /// `1` for N, inserted and mismatched residues, else `0`.
    pub nongermline: String,
}

/// V's and J's columns, cut at `v_end`/`j_start`, their groups and the V(D)J region's end.
struct Layout<'a> {
    v: Vec<Col>,
    j: Vec<Col>,
    vg: &'a AaGroup,
    jg: Option<&'a AaGroup>,
    end: usize,
}

impl<'a> Layout<'a> {
    fn new(len: usize, c: &AaCalls, v: &'a AaSegment, j: &'a AaSegment) -> Layout<'a> {
        let jh = c.j.as_ref();
        Layout {
            v: columns(&c.v.aln, c.v.aln.q_start..c.v_end),
            j: jh.map_or(Vec::new(), |h| columns(&h.aln, c.j_start..h.aln.q_end)),
            vg: &v.groups[c.v.group],
            jg: jh.map(|h| &j.groups[h.group]),
            end: jh.map_or(len, |h| h.aln.q_end).max(c.v_end),
        }
    }
}

/// The query residue aligned to germline residue `u`.
fn at(cols: &[Col], u: usize) -> Option<usize> {
    cols.iter().find(|c| c.1 == Some(u)).and_then(|c| c.0)
}

/// IMGT regions and junction (spec §Regions and junction):
///
/// - A V residue takes the region of its IMGT position (past Cys104: CDR3); an inserted
///   residue takes the region of the residue before it.
/// - Cys104 is the residue aligned to the V anchor; when V's alignment ends before it, the
///   residue its last diagonal projects there (at most `MAX_PROJECTION` on), if `C` and before J.
///   Residues from V's end through it are FWR3, then CDR3.
/// - Trp/Phe118 is the residue aligned to the J anchor; when J's alignment starts past it, the
///   residue J's first diagonal projects there if it is the J allele's own anchor residue, else
///   the nearest such residue up to `MAX_TRP_LOOKBACK` before the last-diagonal projection,
///   after Cys104 either way. FWR4 starts there; J residues before it are CDR3.
pub fn aa_regions(q: &[u8], c: &AaCalls, v: &AaSegment, j: &AaSegment) -> AaRegions {
    let len = q.len();
    let l = Layout::new(len, c, v, j);
    let i0 = c.v.aln.q_start;
    let cys = l.vg.anchor.and_then(|a| {
        at(&l.v, a).or_else(|| {
            let (i, u) = l.v.iter().rev().find_map(|&(i, u)| i.zip(u)).filter(|&(_, u)| u < a && a - u <= MAX_PROJECTION)?;
            Some(i + a - u).filter(|&p| p < c.j_start && q[p] == b'C')
        })
    });
    let trp = (|| {
        let b = l.jg?.anchor?;
        if let Some(t) = at(&l.j, b) {
            return Some(t);
        }
        let own = *l.jg?.seq.get(b)?;
        let aligned = |c: &Col| c.0.zip(c.1);
        let ((i, u), (il, ul)) = (l.j.iter().find_map(aligned)?, l.j.iter().rev().find_map(aligned)?);
        if u <= b {
            return None; // the anchor is deleted
        }
        let ok = |p: &usize| *p > cys.unwrap_or(usize::MAX) && *p < l.end && q[*p] == own;
        i.checked_sub(u - b).filter(ok).or_else(|| {
            let pm = il.checked_sub(ul - b)?;
            (0..=MAX_TRP_LOOKBACK).filter_map(|k| pm.checked_sub(k)).find(ok)
        })
    })();

    let mut vmap = vec![None; len];
    let mut jmap = vec![None; len];
    for (map, cols) in [(&mut vmap, &l.v), (&mut jmap, &l.j)] {
        for &(i, u) in cols {
            if let Some(i) = i {
                map[i] = u;
            }
        }
    }
    let mut reg: Vec<Option<usize>> = vec![None; len];
    for i in i0..l.end {
        let before = i.checked_sub(1).and_then(|p| reg[p]);
        reg[i] = if trp.is_some_and(|t| i >= t) {
            Some(FWR4)
        } else if i < c.v_end {
            vmap[i].map_or(before, |u| Some(V_ENDS.iter().position(|&e| l.vg.imgt[u] < e).unwrap_or(CDR3)))
        } else if cys.is_some_and(|p| i <= p) {
            Some(FWR3)
        } else if trp.is_some() || i < c.j_start {
            Some(CDR3)
        } else {
            jmap[i].map_or(before, |u| Some(if l.jg.and_then(|g| g.anchor).is_some_and(|b| u >= b) { FWR4 } else { CDR3 }))
        };
    }
    let junction = cys.zip(trp).filter(|(s, t)| s < t).map(|(s, t)| (s, t + 1));
    let [fwr1, cdr1, fwr2, cdr2, fwr3, cdr3, fwr4]: [Option<(usize, usize)>; 7] = std::array::from_fn(|k| {
        let (first, last) = (reg.iter().position(|&x| x == Some(k)), reg.iter().rposition(|&x| x == Some(k)));
        if k == CDR3 && junction.is_none() {
            None
        } else if let (Some(a), Some(z)) = (first, last) {
            Some((a, z + 1))
        } else if reg.iter().any(|x| x.is_some_and(|x| x < k)) {
            // every residue deleted, between the regions before and after it
            reg.iter().position(|x| x.is_some_and(|x| x > k)).map(|p| (p, p))
        } else {
            None
        }
    });
    AaRegions { fwr1, cdr1, fwr2, cdr2, fwr3, cdr3, fwr4, junction, vdj: (i0, l.end), labels: reg }
}

/// V and J residue substitutions, insertions and deletions against the primary's group, labelled
/// `IMGT{pos}` (V) or `{pos}` (J). A segment with no aligned residue is `None`.
pub fn aa_mutations(q: &[u8], c: &AaCalls, v: &AaSegment, j: &AaSegment) -> [Option<SegmentMutations>; 2] {
    let l = Layout::new(q.len(), c, v, j);
    [
        segment(q, &l.v, &l.vg.seq, |u| format!("IMGT{}", l.vg.imgt[u] + 1), None),
        l.jg.and_then(|g| segment(q, &l.j, &g.seq, |u| (u + 1).to_string(), None)),
    ]
}

/// The CDR, segment and non-germline masks over `r.vdj`.
pub fn aa_masks(q: &[u8], c: &AaCalls, v: &AaSegment, j: &AaSegment, r: &AaRegions) -> AaMasks {
    let l = Layout::new(q.len(), c, v, j);
    let (i0, end) = r.vdj;
    let (mut seg, mut ng) = (vec![b'N'; end - i0], vec![b'1'; end - i0]);
    for (cols, g, s) in [(&l.v, Some(l.vg), b'V'), (&l.j, l.jg, b'J')] {
        let Some(g) = g else { continue };
        for &(i, u) in cols {
            if let Some(i) = i {
                seg[i - i0] = s;
                ng[i - i0] = if u.is_some_and(|u| g.seq[u].eq_ignore_ascii_case(&q[i])) { b'0' } else { b'1' };
            }
        }
    }
    let cdr = r.cdr3.map(|_| {
        r.labels[i0..end]
            .iter()
            .map(|x| match x {
                Some(1) => '1',
                Some(3) => '2',
                Some(5) => '3',
                _ => '0',
            })
            .collect()
    });
    let text = |b: Vec<u8>| String::from_utf8(b).expect("ASCII");
    AaMasks { cdr, segment: text(seg), nongermline: text(ng) }
}

/// `sequence_alignment_aa` and `germline_alignment_aa` over the V(D)J region: V with its group's
/// IMGT gaps (`.` in both), N residues against `-`, J as aligned; an inserted residue is `-` in
/// the germline, a deleted one `-` in the sequence.
pub fn aa_alignments(q: &[u8], c: &AaCalls, v: &AaSegment, j: &AaSegment) -> (String, String) {
    let l = Layout::new(q.len(), c, v, j);
    let (mut s, mut g) = (String::new(), String::new());
    let mut push = |a: char, b: char| {
        s.push(a);
        g.push(b);
    };
    let res = |i: Option<usize>| i.map_or('-', |i| q[i] as char);
    let mut prev: Option<usize> = None;
    for &(i, u) in &l.v {
        if let Some(u) = u {
            let p = l.vg.imgt[u];
            for _ in 0..prev.map_or(0, |x| p - x - 1) {
                push('.', '.');
            }
            prev = Some(p);
        }
        push(res(i), u.map_or('-', |u| l.vg.seq[u] as char));
    }
    (c.v_end..c.j_start.min(l.end)).for_each(|i| push(q[i] as char, '-'));
    if let Some(jg) = l.jg {
        for &(i, u) in &l.j {
            push(res(i), u.map_or('-', |u| jg.seq[u] as char));
        }
    }
    (s, g)
}
