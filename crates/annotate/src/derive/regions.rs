use super::{imgt_of, Col, Layout};
use crate::assign::{CHit, JHit, Junction, Oriented, VHit};
use abstar_germline::Allele;

/// Exclusive IMGT gapped ends (0-based) of FWR1, CDR1, FWR2, CDR2 and FWR3, which ends with
/// the Cys104 codon.
const V_ENDS: [usize; 5] = [78, 114, 165, 195, 312];
const FWR3: usize = 4;
const CDR3: usize = 5;
const FWR4: usize = 6;
const NAMES: [&str; 7] = ["fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4"];

/// IMGT regions, junction and assembled-sequence extents of a read (record 0032). A region is
/// `None` when no base of it is in the read, and empty (`(p, p)`) when it lies inside the read
/// with every base deleted.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Regions {
    pub fwr1: Option<(usize, usize)>,
    pub cdr1: Option<(usize, usize)>,
    pub fwr2: Option<(usize, usize)>,
    pub cdr2: Option<(usize, usize)>,
    pub fwr3: Option<(usize, usize)>,
    pub cdr3: Option<(usize, usize)>,
    pub fwr4: Option<(usize, usize)>,
    /// From the Cys104 codon's first base through the Trp/Phe118 codon; `None` unless both
    /// codons are in the read.
    pub junction: Option<(usize, usize)>,
    /// The first codon of the V reading frame: codons start in phase with the first V base's
    /// IMGT position.
    pub frame: usize,
    /// `sequence_vdj`: the first V base through J's end, or the read's end without J.
    pub vdj: (usize, usize),
    /// `sequence_vdjc`: `vdj` through C's end; `None` without C.
    pub vdjc: Option<(usize, usize)>,
}

impl Regions {
    /// The seven IMGT regions with their field names, 5' to 3'.
    pub fn all(&self) -> [(&'static str, Option<(usize, usize)>); 7] {
        let r = [self.fwr1, self.cdr1, self.fwr2, self.cdr2, self.fwr3, self.cdr3, self.fwr4];
        std::array::from_fn(|k| (NAMES[k], r[k]))
    }
}

/// Germline index of each read position from `cols`.
fn by_read(cols: &[Col], len: usize) -> Vec<Option<usize>> {
    let mut m = vec![None; len];
    for &(i, u) in cols {
        if let Some(i) = i {
            m[i] = u;
        }
    }
    m
}

/// The read's IMGT regions, junction and assembled-sequence extents (record 0032), as the
/// generator derives them from its origin map, on the called alleles:
///
/// - A V base takes the region of its IMGT gapped position on the primary (V bases past the
///   Cys104 codon are CDR3); an inserted base takes the region of the base before it.
/// - Without the V base, the Cys104 codon is the 3 read bases from the anchor's first base,
///   placed from the read's last V base when the V keeps that base (ruling A7): those are
///   FWR3 whatever segment they came from. N and D bases after them are CDR3; J bases are
///   CDR3 before the J anchor (Trp/Phe118) and FWR4 from it. Bases after J (C) have none.
/// - CDR3 is `None` when the V anchor's first base is not in V, the J has no anchor or its
///   alignment starts past it; it is partial when the read stops inside it (ruling A4).
/// - The junction needs both anchor bases in their segments and the Trp/Phe118 codon in the
///   read.
pub fn regions(read: &Oriented, v: &VHit, j: Option<(&JHit, &Junction)>, c: Option<&CHit>, alleles: &[Allele]) -> Regions {
    let len = read.seq.len();
    let l = Layout::new(read, v, j);
    let va = &alleles[v.primary];
    let imgt = imgt_of(va);
    let (vmap, jmap) = (by_read(&l.v, len), by_read(&l.j, len));
    let i0 = v.aln.q_start;
    let frame = i0 + (3 - imgt[v.aln.g_start] % 3) % 3;
    let jb = j.and_then(|(j, _)| alleles[j.primary].anchor);
    // read position of the Cys104 codon's first base (before the read when truncated away)
    let last_v = l.v.iter().rev().find_map(|&(i, u)| i.zip(u));
    let cys = va.anchor.zip(last_v).filter(|&(a, (_, u))| a <= u).map(|(a, (i, u))| i as i64 - (u - a) as i64);
    let j_first = l.j.iter().find_map(|c| c.1);
    let cdr3_known = cys.is_some() && (j.is_none() || jb.zip(j_first).is_some_and(|(b, u)| u <= b));

    let mut reg: Vec<Option<usize>> = vec![None; len];
    for i in i0..len {
        let before = i.checked_sub(1).and_then(|p| reg[p]);
        reg[i] = if i < l.v_end {
            vmap[i].map_or(before, |u| Some(V_ENDS.iter().position(|&e| imgt[u] < e).unwrap_or(CDR3)))
        } else if i >= l.end {
            None
        } else if cys.is_some_and(|c| (i as i64) < c + 3) {
            Some(FWR3)
        } else if i < l.j_start {
            Some(CDR3)
        } else {
            jmap[i].map_or(before, |u| jb.map(|b| if u < b { CDR3 } else { FWR4 }))
        };
    }
    let [fwr1, cdr1, fwr2, cdr2, fwr3, cdr3, fwr4]: [Option<(usize, usize)>; 7] = std::array::from_fn(|k| {
        let (first, last) = (reg.iter().position(|&x| x == Some(k)), reg.iter().rposition(|&x| x == Some(k)));
        if k == CDR3 && !cdr3_known {
            None
        } else if let (Some(a), Some(z)) = (first, last) {
            Some((a, z + 1))
        } else if reg.iter().any(|x| x.is_some_and(|x| x < k)) {
            // every base deleted, between the regions before and after it
            reg.iter().position(|x| x.is_some_and(|x| x > k)).map(|p| (p, p))
        } else {
            None
        }
    });

    let at = |cols: &[Col], u: usize| cols.iter().find(|c| c.1 == Some(u)).and_then(|c| c.0);
    let start = va.anchor.and_then(|a| at(&l.v, a));
    let end = jb.and_then(|b| at(&l.j, b));
    let junction = match (start, end) {
        (Some(s), Some(e)) if s + 3 <= e && e + 3 <= len => Some((s, e + 3)),
        _ => None,
    };
    Regions { fwr1, cdr1, fwr2, cdr2, fwr3, cdr3, fwr4, junction, frame, vdj: (i0, l.end), vdjc: c.map(|c| (i0, c.end)) }
}
