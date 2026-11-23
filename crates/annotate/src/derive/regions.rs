use super::{imgt_of, Col, Layout};
use crate::align::{align_semiglobal, Op, Scoring};
use crate::assign::{CHit, JHit, Junction, Oriented, VHit};
use abstar_germline::Allele;

/// Exclusive IMGT gapped ends (0-based) of FWR1, CDR1, FWR2, CDR2 and FWR3, which ends with
/// the Cys104 codon.
const V_ENDS: [usize; 5] = [78, 114, 165, 195, 312];
const FWR3: usize = 4;
const CDR3: usize = 5;
const FWR4: usize = 6;
/// Furthest the Cys104 codon is projected past V's last base (ruling B2): the bnab truth set's
/// projections reach 49 nt, all at truth's anchor.
const MAX_PROJECTION: usize = 50;
/// Furthest the realigned Trp/Phe118 codon is placed from J's first diagonal when J starts past
/// it (#162): 15 reaches a 15-nt tandem duplication in a kappa FWR4; on the 1,356 such Jaffe
/// chains the placement finds legacy abstar's codon in 1,342 (record 0032).
const MAX_TRP_SHIFT: usize = 15;
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
    /// Region index (0..7, FWR1..FWR4) of each read position, for the masks.
    pub(crate) labels: Vec<Option<usize>>,
    /// Both CDR3 ends are known (the masks' `cdr_mask` needs it).
    pub(crate) cdr3_known: bool,
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
///   alignment starts past it with no codon placed; it is partial when the read stops inside it (ruling A4).
/// - The junction needs both anchor bases in their segments and the Trp/Phe118 codon in the
///   read. A J starting past its anchor has the codon placed by realigning J from it (#162).
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
    let mut cys = va.anchor.zip(last_v).filter(|&(a, (_, u))| a <= u).map(|(a, (i, u))| i as i64 - (u - a) as i64);
    // V ends before the anchor (a heavily mutated V 3' end): the Cys104 codon is projected along
    // V's last diagonal, at most MAX_PROJECTION nt, when it is TGT or TGC and ends before J (B2)
    let mut projected = None;
    if cys.is_none() && j.is_some() {
        if let Some((a, (i, u))) = va.anchor.zip(last_v).filter(|&(a, (_, u))| a > u && a - u <= MAX_PROJECTION) {
            let p = i + (a - u);
            if p + 3 <= l.j_start && matches!(&read.seq[p..p + 3], b"TGT" | b"TGC") {
                (cys, projected) = (Some(p as i64), Some(p));
            }
        }
    }
    let j_first = l.j.iter().find_map(|c| c.1);
    // J starts past its anchor (a mutated or inserted J 5' end): the Trp/Phe118 codon is the
    // TGG, TTT or TTC that J's first diagonal projects to, else the one J's anchor lands on when J
    // (anchor through alignment end) is realigned to the read after the Cys104 codon, gaps
    // anywhere, if it is in the V frame, at most MAX_TRP_SHIFT nt from the projection and
    // followed by an aligned J base before any gap (a TGG then an insertion is as likely N or D).
    // FWR4 starts there even if V's alignment reached into it by chance (#162)
    let mut trp = None;
    if let (Some(c), Some((jh, _)), Some(b), Some((i, u))) = (cys, j, jb, l.j.iter().find_map(|&(i, u)| i.zip(u))) {
        if u > b && c >= 0 {
            let from = c as usize + 3;
            let codon = |q: usize| q >= from && q + 3 <= l.end && matches!(&read.seq[q..q + 3], b"TGG" | b"TTT" | b"TTC");
            let p = i as i64 - (u - b) as i64;
            trp = usize::try_from(p).ok().filter(|&q| codon(q)).or_else(|| {
                let g = alleles[jh.primary].sequence();
                let a = align_semiglobal(&read.seq[from.min(l.end)..l.end], &g.as_bytes()[b..jh.aln.g_end], &Scoring::default());
                let q = from + a.q_start;
                let run: u32 = a.ops.iter().take_while(|o| matches!(o.0, Op::Match | Op::Mismatch)).map(|o| o.1).sum();
                (a.g_start == 0 && run > 3 && (q - from) % 3 == 0 && (q as i64 - p).unsigned_abs() <= MAX_TRP_SHIFT as u64 && codon(q)).then_some(q)
            });
        }
    }
    let cdr3_known = cys.is_some() && (j.is_none() || trp.is_some() || jb.zip(j_first).is_some_and(|(b, u)| u <= b));

    let mut reg: Vec<Option<usize>> = vec![None; len];
    for i in i0..len {
        let before = i.checked_sub(1).and_then(|p| reg[p]);
        reg[i] = if trp.is_some_and(|t| i >= t && i < l.end) {
            Some(FWR4)
        } else if i < l.v_end {
            vmap[i].map_or(before, |u| Some(V_ENDS.iter().position(|&e| imgt[u] < e).unwrap_or(CDR3)))
        } else if i >= l.end {
            None
        } else if cys.is_some_and(|c| (i as i64) < c + 3) {
            Some(FWR3)
        } else if trp.is_some() || i < l.j_start {
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
    let start = va.anchor.and_then(|a| at(&l.v, a)).or(projected);
    let end = jb.and_then(|b| at(&l.j, b)).or(trp);
    let junction = match (start, end) {
        (Some(s), Some(e)) if s + 3 <= e && e + 3 <= len => Some((s, e + 3)),
        _ => None,
    };
    Regions { fwr1, cdr1, fwr2, cdr2, fwr3, cdr3, fwr4, junction, frame, vdj: (i0, l.end), vdjc: c.map(|c| (i0, c.end)), labels: reg, cdr3_known }
}
