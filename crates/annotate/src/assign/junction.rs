use super::ties::tiebreak_cmp;
use super::{JHit, Oriented, SegmentIndex, VHit};
use crate::align::{Alignment, Op};
use abstar_germline::Allele;
use std::cmp::Reverse;
use std::ops::Range;

/// The D call for one oriented read (records 0030, 0031).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DHit {
    /// Database index of the primary allele (record 0008); `ties[0]`.
    pub primary: usize,
    /// Database indices of the V locus's D alleles reaching the best joint score in the called
    /// orientation, primary first, then record 0008 order.
    pub ties: Vec<usize>,
    /// The D is called in reverse complement (record 0031).
    pub inverted: bool,
    /// Read bases `start..end` are the primary's bases `g_start..g_end`, indexed on the
    /// reverse-complemented allele when `inverted`.
    pub start: usize,
    pub end: usize,
    pub g_start: usize,
    pub g_end: usize,
}

/// The segment boundaries of one oriented read: V covers `..v_end`, J covers `j_start..`
/// (0-based half-open, `v_end <= j_start`), the D (if any) lies between.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Junction {
    pub v_end: usize,
    pub d: Option<DHit>,
    pub j_start: usize,
}

/// A chosen D run shorter than this is no D (record 0025's threshold).
const MIN_D_LEN: usize = 5;

/// A D is called only when its best configuration beats the best no-D configuration by at
/// least this (ruling B3; record 0030): a chance D over N bases rarely gains 6.
const D_MARGIN: i32 = 6;

/// An inverted D is called only when its best score beats the best forward D score by at
/// least this (forward wins ties; record 0031).
const INVERTED_MARGIN: i32 = 2;

/// Score of an impossible cut; sums of three stay far above `i32::MIN`.
const NEG: i32 = -(1 << 28);

/// +1 when read base `r` is ACGT and equals germline base `g`, else −1.
fn hit(r: u8, g: u8) -> i32 {
    if matches!(r, b'A' | b'C' | b'G' | b'T') && r == g { 1 } else { -1 }
}

fn upper(a: &Allele) -> Vec<u8> {
    a.sequence().to_ascii_uppercase().into_bytes()
}

/// V score by V end: `a[e]` sums read positions `aln.q_start..e` (+1/−1 against the allele
/// `g`, 0 for an inserted base); past the alignment's end V continues on its last diagonal
/// while the allele has bases. `NEG` where V cannot end.
fn v_scores(q: &[u8], g: &[u8], aln: &Alignment) -> Vec<i32> {
    let mut a = vec![NEG; q.len() + 1];
    let (mut qi, mut gi, mut sum) = (aln.q_start, aln.g_start, 0);
    for &(op, n) in &aln.ops {
        for _ in 0..n {
            match op {
                Op::Match | Op::Mismatch => {
                    sum += hit(q[qi], g[gi]);
                    (qi, gi) = (qi + 1, gi + 1);
                }
                Op::Ins => qi += 1,
                Op::Del => {
                    gi += 1;
                    continue;
                }
            }
            a[qi] = sum;
        }
    }
    while qi < q.len() && gi < g.len() {
        sum += hit(q[qi], g[gi]);
        (qi, gi) = (qi + 1, gi + 1);
        a[qi] = sum;
    }
    a
}

/// J score by J start: `b[t]` sums read positions `t..aln.q_end` the same way; before the
/// alignment's start J continues on its first diagonal while the allele and the read have
/// bases. `NEG` where J cannot start.
fn j_scores(q: &[u8], g: &[u8], aln: &Alignment) -> Vec<i32> {
    let mut b = vec![NEG; q.len() + 1];
    let (mut qi, mut gi, mut sum) = (aln.q_end, aln.g_end, 0);
    for &(op, n) in aln.ops.iter().rev() {
        for _ in 0..n {
            match op {
                Op::Match | Op::Mismatch => {
                    (qi, gi) = (qi - 1, gi - 1);
                    sum += hit(q[qi], g[gi]);
                }
                Op::Ins => qi -= 1,
                Op::Del => {
                    gi -= 1;
                    continue;
                }
            }
            b[qi] = sum;
        }
    }
    while qi > 0 && gi > 0 {
        (qi, gi) = (qi - 1, gi - 1);
        sum += hit(q[qi], g[gi]);
        b[qi] = sum;
    }
    b
}

/// Every placement of a D of length `lg` overlapping read positions `zone`: (read position of
/// the D's first base, the read positions it covers inside `zone`).
fn placements(lg: usize, zone: Range<usize>) -> impl Iterator<Item = (i64, Range<usize>)> {
    (zone.start as i64 - lg as i64 + 1..zone.end as i64).map(move |s| {
        let lo = s.max(zone.start as i64) as usize;
        (s, lo..((s + lg as i64) as usize).min(zone.end))
    })
}

/// Best V + D + J score with the D at `s` restricted to `w` (V end ≤ D start < D end ≤ J
/// start): D from `i1` to `i2` scores `pa[i1] − P(i1) + P(i2) + sb[i2]`, `P` the running D score.
// ponytail: one serial pass per placement is about 165 µs/read over 62 D orientations; SIMD
// across offsets (independent diagonals, i16 lanes) is the upgrade path if it matters.
fn best_d(q: &[u8], g: &[u8], s: i64, w: Range<usize>, pa: &[i32], sb: &[i32]) -> i32 {
    let gw = &g[(w.start as i64 - s) as usize..(w.end as i64 - s) as usize];
    let (mut p, mut ml, mut best) = (0, NEG, NEG);
    for (((&r, &b), &v), &j) in q[w.clone()].iter().zip(gw).zip(&pa[w.clone()]).zip(&sb[w.start + 1..=w.end]) {
        ml = ml.max(v - p);
        p += hit(r, b);
        best = best.max(ml + p + j);
    }
    best
}

/// Resolves the segment boundaries and the D call of an oriented read after its V and J calls
/// (records 0030, 0031). The V and J tie sets are not changed.
///
/// A configuration gives V the read up to the V end, J the read from the J start, optionally a
/// D placement between them (V end ≤ D start < D end ≤ J start), and scores +1/−1 for each
/// base given to a segment that matches/doesn't match its germline; N/P bases score 0. V and
/// J are the primaries on their alignments' diagonals (`v_scores`, `j_scores`); D is every D
/// allele of the V primary's locus (`d`, uppercase) in both orientations at every offset.
///
/// - **No D** when `locus_has_d` is false, when no D configuration beats the best no-D score by
///   `D_MARGIN`, or when the chosen D configuration's run is shorter than `MIN_D_LEN`: the
///   boundaries are the best no-D configuration's.
/// - The D orientation is inverted only when the best inverted score beats the best forward
///   score by at least `INVERTED_MARGIN` (record 0031). The D tie set is the alleles reaching
///   that orientation's best score, the primary its record 0008 pick, and the boundaries come
///   from the primary's best configuration.
/// - Ties break toward the longest V, then the longest J, then the earliest D start, then the
///   longest D.
pub fn resolve_junction(read: &Oriented, v: &VHit, j: &JHit, alleles: &[Allele], d: &SegmentIndex, locus_has_d: bool) -> Junction {
    let q = &read.seq;
    let len = q.len();
    let (mut pa, mut sb) = (v_scores(q, &upper(&alleles[v.primary]), &v.aln), j_scores(q, &upper(&alleles[j.primary]), &j.aln));
    // pa[c]: best V score with V end <= c, reached last at ea[c]; sb[c]: best J score with J
    // start >= c, reached first at tb[c]
    let (mut ea, mut tb): (Vec<usize>, Vec<usize>) = ((0..=len).collect(), (0..=len).collect());
    for x in 1..=len {
        if pa[x - 1] > pa[x] {
            (pa[x], ea[x]) = (pa[x - 1], ea[x - 1]);
        }
    }
    for y in (0..len).rev() {
        if sb[y + 1] > sb[y] {
            (sb[y], tb[y]) = (sb[y + 1], tb[y + 1]);
        }
    }
    let s0 = (0..=len).map(|c| pa[c] + sb[c]).max().expect("cuts");
    if s0 < NEG / 2 {
        // no cut leaves V and J a base each (J ends where V starts): keep the alignments' ends
        return Junction { v_end: v.aln.q_end, d: None, j_start: j.aln.q_start.max(v.aln.q_end) };
    }
    let (v_end, j_start) =
        (0..=len).filter(|&c| pa[c] + sb[c] == s0).map(|c| (ea[c], tb[c])).max_by_key(|&(e, t)| (e, Reverse(t))).expect("best cut");
    let no_d = Junction { v_end, d: None, j_start };
    let ks: Vec<usize> = (0..d.ids.len()).filter(|&k| locus_has_d && d.loci[k] == alleles[v.primary].locus).collect();
    // A D of length lg starting or ending at cut c caps the score at pa[c] + lg + sb[c], so a D
    // reaching s0 lies between the first and last cuts where that cap reaches s0
    let zone = |lg: usize| {
        let near = |c: &usize| pa[*c] + lg as i32 + sb[*c] >= s0;
        (0..=len).find(near).expect("best cut")..(0..=len).rfind(near).expect("best cut")
    };

    let seqs: Vec<[&[u8]; 2]> = ks.iter().map(|&k| [&d.seqs[k][..], &d.rc[k][..]]).collect();
    let best_of = |g: &[u8]| placements(g.len(), zone(g.len())).map(|(s, w)| best_d(q, g, s, w, &pa, &sb)).max().unwrap_or(NEG);
    let scores: Vec<[i32; 2]> = seqs.iter().map(|[f, r]| [best_of(f), best_of(r)]).collect();
    let best = |o: usize| scores.iter().map(|s| s[o]).max().unwrap_or(NEG);
    let inverted = best(1) - best(0) >= INVERTED_MARGIN;
    let o = usize::from(inverted);
    let sd = best(o);
    if sd < s0 + D_MARGIN {
        return no_d;
    }
    let mut ties: Vec<usize> = (0..ks.len()).filter(|&x| scores[x][o] == sd).collect();
    ties.sort_by(|&x, &y| tiebreak_cmp(&alleles[d.ids[ks[x]]], &alleles[d.ids[ks[y]]]).then(x.cmp(&y)));
    let g = seqs[ties[0]][o];

    // The primary's configurations reaching sd: (V end, J start, D start, D end, D offset).
    let mut pick: Option<(usize, usize, usize, usize, i64)> = None;
    let key = |p: &(usize, usize, usize, usize, i64)| (p.0, Reverse(p.1), Reverse(p.2), p.3);
    for (s, w) in placements(g.len(), zone(g.len())) {
        if best_d(q, g, s, w.clone(), &pa, &sb) != sd {
            continue;
        }
        let mut p = vec![0];
        for i in w.clone() {
            p.push(p[p.len() - 1] + hit(q[i], g[(i as i64 - s) as usize]));
        }
        for k1 in 0..w.len() {
            for k2 in k1 + 1..=w.len() {
                let (i1, i2) = (w.start + k1, w.start + k2);
                let c = (ea[i1], tb[i2], i1, i2, s);
                if pa[i1] - p[k1] + p[k2] + sb[i2] == sd && pick.is_none_or(|b| key(&c) > key(&b)) {
                    pick = Some(c);
                }
            }
        }
    }
    let (v_end, j_start, start, end, s) = pick.expect("a placement reaches sd");
    if end - start < MIN_D_LEN {
        return no_d;
    }
    let ties = ties.iter().map(|&x| d.ids[ks[x]]).collect::<Vec<_>>();
    let (g_start, g_end) = ((start as i64 - s) as usize, (end as i64 - s) as usize);
    Junction { v_end, d: Some(DHit { primary: ties[0], ties, inverted, start, end, g_start, g_end }), j_start }
}
