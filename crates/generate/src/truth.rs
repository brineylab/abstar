//! Truth sets: every equally valid allele call and segment boundary for a simulated read,
//! by the fewest-mismatches rule (decision 0025). Computed from the origin map alone; no
//! read is ever aligned.

use std::collections::BTreeSet;

use abstar_germline::Allele;

use crate::recombine::revcomp;
use crate::{Candidates, Origin, Simulated};

/// A D run shorter than this also accepts "no D" (decision 0025).
const MIN_D_LEN: usize = 5;

/// Score of an impossible cut; sums of up to three stay far below any real score.
const NEG: i32 = -(1 << 28);

/// Allele sets hold sorted allele names. Boundaries are sorted 1-based positions in the
/// final read (AIRR `v_sequence_end`, `d_sequence_start`, `d_sequence_end`,
/// `j_sequence_start`) appearing in any configuration tied for the best score of some
/// combination of valid alleles; a segment with no bases in the read has empty sets and lists.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Truth {
    pub v_valid: Vec<String>,
    pub d_valid: Vec<String>,
    /// The read's D-derived run is shorter than 5 nt (always true without a D).
    pub d_null_valid: bool,
    pub j_valid: Vec<String>,
    pub c_valid: Vec<String>,
    pub v_end_valid: Vec<i64>,
    pub d_start_valid: Vec<i64>,
    pub d_end_valid: Vec<i64>,
    pub j_start_valid: Vec<i64>,
    /// The true allele is not in its (non-empty) set.
    pub unrecoverable_v: bool,
    pub unrecoverable_d: bool,
    pub unrecoverable_j: bool,
}

/// Uppercase candidate sequences, built once per run (see `Candidates::new`).
#[derive(Debug, Clone)]
pub(crate) struct Seqs {
    v_gapped: Vec<Vec<u8>>,
    v: Vec<Vec<u8>>,
    d: Vec<Vec<u8>>,
    d_rc: Vec<Vec<u8>>,
    j: Vec<Vec<u8>>,
    c: Vec<Vec<u8>>,
}

impl Seqs {
    pub(crate) fn new(v: &[Allele], d: &[Allele], j: &[Allele], c: &[Allele]) -> Self {
        let up = |a: &[Allele]| -> Vec<Vec<u8>> { a.iter().map(|a| a.sequence().to_ascii_uppercase().into_bytes()).collect() };
        let d = up(d);
        Self {
            v_gapped: v.iter().map(|a| a.sequence_gapped.to_ascii_uppercase().into_bytes()).collect(),
            v: up(v),
            d_rc: d.iter().map(|s| revcomp(s)).collect(),
            d,
            j: up(j),
            c: up(c),
        }
    }
}

/// +1 if `g` is the read base `r`, −1 otherwise: a different base, a gap or no base (the
/// coverage rule). Both sides are uppercase.
fn score(r: u8, g: Option<&u8>) -> i32 {
    if g == Some(&r) { 1 } else { -1 }
}

fn at(seq: &[u8], x: i64) -> Option<&u8> {
    usize::try_from(x).ok().and_then(|x| seq.get(x))
}

/// Read bases placed at their index (`0` where there is none, which matches nothing), so a
/// candidate is compared with one zip.
fn project(read: &[u8], bases: &[(usize, usize)]) -> Vec<u8> {
    let mut out = vec![0; bases.iter().map(|x| x.1 + 1).max().unwrap_or(0)];
    for &(i, u) in bases {
        out[u] = read[i];
    }
    out
}

fn matches(proj: &[u8], seq: &[u8]) -> usize {
    // u8 lanes vectorize; 255-byte chunks cannot overflow them
    proj.chunks(255)
        .zip(seq.chunks(255))
        .map(|(a, b)| a.iter().zip(b).fold(0u8, |n, (x, y)| n + u8::from(x == y)) as usize)
        .sum()
}

/// Fewest mismatches of `seq` against the `total` projected read bases when index `u` is
/// compared with `seq[u + k]`, over every shift `k` covering at least one of them, and the
/// shifts tied for it.
fn slide(proj: &[u8], total: usize, seq: &[u8]) -> (usize, Vec<i64>) {
    let hi = proj.len() as i64 - 1;
    let lo = proj.iter().position(|&b| b != 0).unwrap_or(0) as i64;
    // m[k + hi]: matches at shift k
    let mut m = vec![0u16; seq.len() + proj.len()];
    for (u, &r) in proj.iter().enumerate().filter(|x| *x.1 != 0) {
        let base = (hi - u as i64) as usize;
        for (c, &b) in m[base..base + seq.len()].iter_mut().zip(seq) {
            *c += u16::from(b == r);
        }
    }
    // pre[x]: read bases among proj[..x]
    let mut pre = vec![0usize; proj.len() + 1];
    for (x, &b) in proj.iter().enumerate() {
        pre[x + 1] = pre[x] + usize::from(b != 0);
    }
    let mut best = (usize::MAX, Vec::new());
    for k in -hi..seq.len() as i64 - lo {
        let (a, b) = ((-k).max(0) as usize, (seq.len() as i64 - k).min(hi + 1) as usize);
        if a >= b || pre[b] == pre[a] {
            continue;
        }
        let mm = total - m[(k + hi) as usize] as usize;
        if mm < best.0 {
            best = (mm, vec![k]);
        } else if mm == best.0 {
            best.1.push(k);
        }
    }
    best
}

/// Indices of the smallest scores, and their sorted, deduplicated names.
fn fewest(alleles: &[Allele], scores: &[usize]) -> (Vec<usize>, Vec<String>) {
    let Some(&min) = scores.iter().min() else { return (vec![], vec![]) };
    let idx: Vec<usize> = (0..scores.len()).filter(|&x| scores[x] == min).collect();
    let mut names: Vec<String> = idx.iter().map(|&x| alleles[x].name.clone()).collect();
    names.sort();
    names.dedup();
    (idx, names)
}

/// Where a segment placement may lie: read positions `start..start + scores.len()`, each
/// with its +1/−1 (0 for inserted bases) against the placed allele.
struct Seg {
    start: usize,
    scores: Vec<i32>,
}

/// V scores by V end (`a[e]`, V covering `..e`) and the best at or before each cut (`pa`);
/// an absent V is pinned to cut 0.
fn v_cut(v: Option<&Seg>, len: usize) -> (Vec<i32>, Vec<i32>) {
    let mut a = vec![NEG; len + 1];
    match v {
        Some(v) => {
            let mut sum = 0;
            for (x, &sc) in v.scores.iter().enumerate() {
                sum += sc;
                a[v.start + x + 1] = sum;
            }
        }
        None => a[0] = 0,
    }
    let mut pa = a.clone();
    for x in 1..=len {
        pa[x] = pa[x].max(pa[x - 1]);
    }
    (a, pa)
}

/// J scores by J start (`b[t]`, J covering `t..`) and the best at or after each cut (`sb`);
/// an absent J is pinned to cut `len`.
fn j_cut(j: Option<&Seg>, len: usize) -> (Vec<i32>, Vec<i32>) {
    let mut b = vec![NEG; len + 1];
    match j {
        Some(j) => {
            let mut sum = 0;
            for (x, &sc) in j.scores.iter().enumerate().rev() {
                sum += sc;
                b[j.start + x] = sum;
            }
        }
        None => b[len] = 0,
    }
    let mut sb = b.clone();
    for y in (0..len).rev() {
        sb[y] = sb[y].max(sb[y + 1]);
    }
    (b, sb)
}

/// A D placement: read positions `start..end`, position `i` holding `q[i + off]`.
struct DPlace<'a> {
    q: &'a [u8],
    off: i64,
    start: usize,
    end: usize,
}

impl DPlace<'_> {
    fn score(&self, read: &[u8], i: usize) -> i32 {
        score(read[i], at(self.q, i as i64 + self.off))
    }
}

/// Best V + D + J with D at `d` (V end ≤ D start < D end ≤ J start).
fn best_d(read: &[u8], pa: &[i32], d: &DPlace, sb: &[i32]) -> i32 {
    // D from i1 to i2 scores pa[i1] − P(i1) + P(i2) + sb[i2], P the running D score
    let (mut p, mut ml, mut best) = (0, NEG, NEG);
    for i in d.start..d.end {
        ml = ml.max(pa[i] - p);
        p += d.score(read, i);
        best = best.max(ml + p + sb[i + 1]);
    }
    best
}

/// The D starts and ends of every configuration reaching `max` with D at `d`.
fn d_ties(read: &[u8], pa: &[i32], d: &DPlace, sb: &[i32], max: i32) -> (Vec<usize>, Vec<usize>) {
    let mut p = vec![0];
    for i in d.start..d.end {
        p.push(p[p.len() - 1] + d.score(read, i));
    }
    let w = d.end - d.start;
    let left = |k: usize| pa[d.start + k] - p[k];
    let right = |k: usize| p[k] + sb[d.start + k];
    let mut ml = vec![NEG; w + 1]; // best left(k1), k1 < k
    for k in 1..=w {
        ml[k] = ml[k - 1].max(left(k - 1));
    }
    let mut mr = vec![NEG; w + 1]; // best right(k2), k2 > k
    for k in (0..w).rev() {
        mr[k] = mr[k + 1].max(right(k + 1));
    }
    let starts = (0..w).filter(|&k| left(k) + mr[k] == max).map(|k| d.start + k).collect();
    let ends = (1..=w).filter(|&k| ml[k] + right(k) == max).map(|k| d.start + k).collect();
    (starts, ends)
}

/// Best V + J with no D (V end ≤ J start), and the cuts `t` tying it (V ends by `t`, J
/// starts from `t`).
fn best_no_d(pa: &[i32], sb: &[i32]) -> Option<(i32, Vec<usize>)> {
    let max = pa.iter().zip(sb).map(|(a, b)| a + b).max()?;
    (max > NEG / 2).then(|| (max, (0..pa.len()).filter(|&t| pa[t] + sb[t] == max).collect()))
}

/// Union over cuts `x` in `xs` of the best-scoring ends in `scores` reached from `x`:
/// indices `..=x` when `forward`, else `x..`.
fn argmax_by_cut(scores: &[i32], xs: &BTreeSet<usize>, forward: bool) -> BTreeSet<usize> {
    let order: Vec<usize> = if forward { (0..scores.len()).collect() } else { (0..scores.len()).rev().collect() };
    let (mut out, mut best, mut arg) = (BTreeSet::new(), NEG / 2, Vec::new());
    for e in order {
        if scores[e] > best {
            best = scores[e];
            arg.clear();
        }
        if scores[e] == best {
            arg.push(e);
        }
        if xs.contains(&e) {
            out.extend(&arg);
        }
    }
    out
}

pub fn truth(s: &Simulated, c: &Candidates) -> Truth {
    let read: Vec<u8> = s.read.to_ascii_uppercase();
    let (read, len, seqs) = (&read[..], read.len(), &c.seqs);
    // (read position, index) of each segment's bases; D bases indexed by their place in the run
    let (mut v_at, mut d_at, mut j_at, mut c_at) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    for (i, &o) in s.origin.iter().enumerate() {
        match o {
            Origin::V(g) => v_at.push((i, g)),
            Origin::D(_) => d_at.push((i, d_at.len())),
            Origin::J(u) => j_at.push((i, u)),
            Origin::C(u) => c_at.push((i, u)),
            _ => {}
        }
    }
    let (v_proj, d_proj, j_proj, c_proj) = (project(read, &v_at), project(read, &d_at), project(read, &j_at), project(read, &c_at));

    // V: IMGT-gapped position by position; C: prefix
    let mm = |seqs: &[Vec<u8>], proj: &[u8], total: usize| -> Vec<usize> {
        if total == 0 { vec![] } else { seqs.iter().map(|q| total - matches(proj, q)).collect() }
    };
    let (v_set, v_valid) = fewest(&c.v, &mm(&seqs.v_gapped, &v_proj, v_at.len()));
    let (_, c_valid) = fewest(&c.c, &mm(&seqs.c, &c_proj, c_at.len()));

    // D: every allele (and its reverse complement when inversion is on) at every shift
    let mut d_scores = Vec::new();
    let mut d_shifts: Vec<Vec<(&[u8], i64)>> = Vec::new();
    if !d_at.is_empty() {
        for (x, fwd) in seqs.d.iter().enumerate() {
            let mut orients = vec![(slide(&d_proj, d_at.len(), fwd), &fwd[..])];
            if c.d_inversion {
                let rc = &seqs.d_rc[x];
                orients.push((slide(&d_proj, d_at.len(), rc), &rc[..]));
            }
            let best = orients.iter().map(|o| o.0 .0).min().unwrap();
            d_scores.push(best);
            d_shifts.push(
                orients.into_iter().filter(|o| o.0 .0 == best).flat_map(|((_, ks), q)| ks.into_iter().map(move |k| (q, k))).collect(),
            );
        }
    }
    let (d_set, d_valid) = fewest(&c.d, &d_scores);

    // J: every allele at every shift of the true J index
    let mut j_scores = Vec::new();
    let mut j_shifts = Vec::new();
    if !j_at.is_empty() {
        for q in &seqs.j {
            let (mm, ks) = slide(&j_proj, j_at.len(), q);
            j_scores.push(mm);
            j_shifts.push(ks);
        }
    }
    let (j_set, j_valid) = fewest(&c.j, &j_scores);

    // Placements. V by its true IMGT mapping, then the allele's own bases past the last
    // V-derived base; D and J at each tied best shift.
    let v_cuts: Vec<(Vec<i32>, Vec<i32>)> = match v_at.last() {
        None => vec![v_cut(None, len)],
        Some(&(pl, gl)) => v_set
            .iter()
            .map(|&x| {
                let g = &seqs.v_gapped[x];
                let mut scores: Vec<i32> = (0..=pl)
                    .map(|i| match s.origin[i] {
                        Origin::V(y) => score(read[i], g.get(y)),
                        _ => 0,
                    })
                    .collect();
                let next = g.iter().take(gl + 1).filter(|&&b| b != b'.').count();
                scores.extend(seqs.v[x].iter().skip(next).zip(&read[pl + 1..]).map(|(b, &r)| score(r, Some(b))));
                v_cut(Some(&Seg { start: 0, scores }), len)
            })
            .collect(),
    };
    let d_groups: Vec<Vec<DPlace>> = d_set
        .iter()
        .map(|&x| {
            d_shifts[x]
                .iter()
                .map(|&(q, k)| {
                    // read position i holds index i − r0 + k
                    let off = k - d_at[0].0 as i64;
                    let start = (-off).max(0) as usize;
                    let end = (q.len() as i64 - off).min(len as i64) as usize;
                    DPlace { q, off, start, end }
                })
                .collect()
        })
        .collect();
    let j_groups: Vec<Vec<(Vec<i32>, Vec<i32>)>> = match (j_at.first(), j_at.last()) {
        (Some(&(pf, uf)), Some(&(pl, _))) => j_set
            .iter()
            .map(|&x| {
                let q = &seqs.j[x];
                j_shifts[x]
                    .iter()
                    .map(|&k| {
                        // J bases before the first J-derived one, as far as the allele and read go
                        let ext = (uf as i64 + k).clamp(0, pf as i64) as usize;
                        let scores = (pf - ext..=pl)
                            .map(|i| match s.origin[i] {
                                Origin::J(u) => score(read[i], at(q, u as i64 + k)),
                                _ if i >= pf => 0,
                                _ => score(read[i], at(q, uf as i64 + k - (pf - i) as i64)),
                            })
                            .collect();
                        j_cut(Some(&Seg { start: pf - ext, scores }), len)
                    })
                    .collect()
            })
            .collect(),
        _ => vec![vec![j_cut(None, len)]],
    };

    // Each combination of valid alleles (V × D-or-no-D × J) keeps the configurations tying
    // its own best over its placements; the boundary lists are the union. No-D
    // combinations run when D is absent or short (d_null_valid).
    let with_d = !d_at.is_empty();
    let no_d = d_at.len() < MIN_D_LEN;
    let mut v_lims: Vec<BTreeSet<usize>> = vec![BTreeSet::new(); v_cuts.len()];
    let mut j_lims: Vec<Vec<BTreeSet<usize>>> = j_groups.iter().map(|g| vec![BTreeSet::new(); g.len()]).collect();
    let (mut d_starts, mut d_ends) = (BTreeSet::new(), BTreeSet::new());
    for (vi, (_, pa)) in v_cuts.iter().enumerate() {
        for (gi, jg) in j_groups.iter().enumerate() {
            if no_d {
                let found: Vec<_> = jg.iter().enumerate().filter_map(|(ji, (_, sb))| Some((best_no_d(pa, sb)?, ji))).collect();
                let max = found.iter().map(|x| x.0 .0).max();
                for ((m, ts), ji) in found {
                    if Some(m) == max {
                        v_lims[vi].extend(&ts);
                        j_lims[gi][ji].extend(ts);
                    }
                }
            }
            for dg in d_groups.iter().filter(|_| with_d) {
                let scored: Vec<(i32, usize, &DPlace)> = jg
                    .iter()
                    .enumerate()
                    .flat_map(|(ji, (_, sb))| dg.iter().map(move |d| (best_d(read, pa, d, sb), ji, d)))
                    .collect();
                let max = scored.iter().map(|x| x.0).max().unwrap_or(NEG);
                if max < NEG / 2 {
                    continue;
                }
                for &(_, ji, d) in scored.iter().filter(|x| x.0 == max) {
                    let (st, en) = d_ties(read, pa, d, &jg[ji].1, max);
                    v_lims[vi].extend(&st);
                    j_lims[gi][ji].extend(&en);
                    d_starts.extend(st);
                    d_ends.extend(en);
                }
            }
        }
    }
    let mut v_ends = BTreeSet::new();
    for ((a, _), lims) in v_cuts.iter().zip(&v_lims) {
        v_ends.extend(argmax_by_cut(a, lims, true));
    }
    let mut j_starts = BTreeSet::new();
    for (g, lims) in j_groups.iter().zip(&j_lims) {
        for ((b, _), l) in g.iter().zip(lims) {
            j_starts.extend(argmax_by_cut(b, l, false));
        }
    }

    let list = |cuts: &BTreeSet<usize>, add: i64, present: bool| -> Vec<i64> {
        if present { cuts.iter().map(|&x| x as i64 + add).collect() } else { vec![] }
    };
    let unrec = |set: &[String], a: Option<&Allele>| a.is_some_and(|a| !set.is_empty() && !set.contains(&a.name));

    Truth {
        v_end_valid: list(&v_ends, 0, !v_at.is_empty()),
        d_start_valid: list(&d_starts, 1, with_d),
        d_end_valid: list(&d_ends, 0, with_d),
        j_start_valid: list(&j_starts, 1, !j_at.is_empty()),
        unrecoverable_v: unrec(&v_valid, Some(&s.v)),
        unrecoverable_d: unrec(&d_valid, s.d.as_ref()),
        unrecoverable_j: unrec(&j_valid, Some(&s.j)),
        d_null_valid: no_d,
        v_valid,
        d_valid,
        j_valid,
        c_valid,
    }
}
