//! Truth sets: every equally valid allele call and segment boundary for a simulated read,
//! by the fewest-mismatches rule (decision 0025). Computed from the origin map alone; no
//! read is ever aligned.

use std::collections::BTreeSet;

use abstar_germline::Allele;

use crate::recombine::revcomp;
use crate::{Candidates, Origin, Simulated};

/// A D run shorter than this also accepts "no D" (decision 0025).
const MIN_D_LEN: usize = 5;

/// Allele sets hold sorted allele names. Boundaries are sorted 1-based positions in the
/// final read (AIRR `v_sequence_end`, `d_sequence_start`, `d_sequence_end`,
/// `j_sequence_start`) appearing in any configuration tied for the best score; a segment
/// with no bases in the read has empty sets and lists.
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

/// +1 if `g` is a base equal to `r` (ignoring case); −1 for a different base, a gap or no
/// base (the coverage rule).
fn score(r: u8, g: Option<&u8>) -> i32 {
    match g {
        Some(&g) if g != b'.' && g.eq_ignore_ascii_case(&r) => 1,
        _ => -1,
    }
}

fn at(seq: &[u8], x: i64) -> Option<&u8> {
    usize::try_from(x).ok().and_then(|x| seq.get(x))
}

/// Fewest mismatches of `seq` against read bases `(read position, index)` when each is
/// compared with `seq[index + k]`, over every shift `k` covering at least one of them, and
/// the shifts tied for it. `bases` is non-empty with ascending indices.
fn slide(read: &[u8], bases: &[(usize, usize)], seq: &[u8]) -> (usize, Vec<i64>) {
    let (lo, hi) = (bases[0].1 as i64, bases[bases.len() - 1].1 as i64);
    let mut best = (usize::MAX, Vec::new());
    for k in -hi..seq.len() as i64 - lo {
        if !bases.iter().any(|&(_, u)| at(seq, u as i64 + k).is_some()) {
            continue;
        }
        let mm = bases.iter().filter(|&&(i, u)| score(read[i], at(seq, u as i64 + k)) < 0).count();
        if mm < best.0 {
            best = (mm, vec![k]);
        } else if mm == best.0 {
            best.1.push(k);
        }
    }
    best
}

/// The sorted names with the fewest mismatches, and the payloads of those entries.
fn fewest<T>(scored: Vec<(&Allele, usize, T)>) -> (Vec<String>, Vec<T>) {
    let Some(min) = scored.iter().map(|x| x.1).min() else { return (vec![], vec![]) };
    let (mut names, mut kept) = (Vec::new(), Vec::new());
    for (a, _, t) in scored.into_iter().filter(|x| x.1 == min) {
        names.push(a.name.clone());
        kept.push(t);
    }
    names.sort();
    names.dedup();
    (names, kept)
}

/// Where a segment placement may lie: read positions `start..start + scores.len()`, each
/// with its +1/−1 (0 for inserted bases) against the placed allele.
struct Seg {
    start: usize,
    scores: Vec<i32>,
}

/// Each placement, or a single "absent" when the segment has no bases in the read.
fn options(segs: &[Seg], present: bool) -> Vec<Option<&Seg>> {
    if present { segs.iter().map(Some).collect() } else { vec![None] }
}

/// Best total of one (V, D, J) placement combination and the cut points of every
/// configuration tying it: `[V end, D start, D end, J start]` as 0-based cuts (V covers
/// `..ve`, D `ds..de`, J `js..`). An absent V or J is pinned to cut 0 or `len`.
fn best(v: Option<&Seg>, d: Option<&Seg>, j: Option<&Seg>, len: usize) -> Option<(i32, [BTreeSet<usize>; 4])> {
    // a[e]: V score for V end e; b[t]: J score for J start t
    let (mut a, mut b) = (vec![None; len + 1], vec![None; len + 1]);
    match v {
        Some(v) => {
            let mut sum = 0;
            for (x, &sc) in v.scores.iter().enumerate() {
                sum += sc;
                a[v.start + x + 1] = Some(sum);
            }
        }
        None => a[0] = Some(0),
    }
    match j {
        Some(j) => {
            let mut sum = 0;
            for (x, &sc) in j.scores.iter().enumerate().rev() {
                sum += sc;
                b[j.start + x] = Some(sum);
            }
        }
        None => b[len] = Some(0),
    }
    // pa[x] = best V ending at or before x; sb[y] = best J starting at or after y
    let mut pa = a.clone();
    for x in 1..=len {
        pa[x] = pa[x].max(pa[x - 1]);
    }
    let mut sb = b.clone();
    for y in (0..len).rev() {
        sb[y] = sb[y].max(sb[y + 1]);
    }

    // middles: (V must end by, D score, J must start from, D cuts)
    let mut mids = Vec::new();
    match d {
        Some(d) => {
            for s in 0..d.scores.len() {
                let mut sum = 0;
                for e in s..d.scores.len() {
                    sum += d.scores[e];
                    mids.push((d.start + s, sum, d.start + e + 1, Some((d.start + s, d.start + e + 1))));
                }
            }
        }
        None => mids.extend((0..=len).map(|t| (t, 0, t, None))),
    }
    let total = |&(x, sc, y, _): &(usize, i32, usize, _)| Some(pa[x]? + sc + sb[y]?);
    let max = mids.iter().filter_map(total).max()?;
    let mut cuts: [BTreeSet<usize>; 4] = Default::default();
    for m in mids.iter().filter(|m| total(m) == Some(max)) {
        let (x, _, y, dc) = *m;
        cuts[0].extend((0..=x).filter(|&e| a[e].is_some() && a[e] == pa[x]));
        cuts[3].extend((y..=len).filter(|&t| b[t].is_some() && b[t] == sb[y]));
        if let Some((s, e)) = dc {
            cuts[1].insert(s);
            cuts[2].insert(e);
        }
    }
    Some((max, cuts))
}

pub fn truth(s: &Simulated, c: &Candidates) -> Truth {
    let (read, len) = (&s.read[..], s.read.len());
    let of = |f: fn(Origin) -> Option<usize>| -> Vec<(usize, usize)> {
        s.origin.iter().enumerate().filter_map(|(i, &o)| f(o).map(|x| (i, x))).collect()
    };
    let v_at = of(|o| if let Origin::V(g) = o { Some(g) } else { None });
    let j_at = of(|o| if let Origin::J(u) = o { Some(u) } else { None });
    let c_at = of(|o| if let Origin::C(u) = o { Some(u) } else { None });
    // D bases indexed by their place in the run
    let d_at: Vec<(usize, usize)> =
        of(|o| if let Origin::D(x) = o { Some(x) } else { None }).iter().enumerate().map(|(x, &(i, _))| (i, x)).collect();

    // V: IMGT-gapped position by position
    let mut v_scored = Vec::new();
    if !v_at.is_empty() {
        for a in &c.v {
            let g = a.sequence_gapped.as_bytes();
            v_scored.push((a, v_at.iter().filter(|&&(i, x)| score(read[i], g.get(x)) < 0).count(), a));
        }
    }
    let (v_valid, v_set) = fewest(v_scored);

    // D: every allele (and its reverse complement when inversion is on) at every shift
    let mut d_scored = Vec::new();
    if !d_at.is_empty() {
        for a in &c.d {
            let fwd = a.sequence().into_bytes();
            let mut orients = vec![(slide(read, &d_at, &fwd), fwd.clone())];
            if c.d_inversion {
                let rc = revcomp(&fwd);
                orients.push((slide(read, &d_at, &rc), rc));
            }
            let mm = orients.iter().map(|x| x.0 .0).min().unwrap();
            let places: Vec<(Vec<u8>, i64)> = orients
                .into_iter()
                .filter(|x| x.0 .0 == mm)
                .flat_map(|((_, ks), seq)| ks.into_iter().map(move |k| (seq.clone(), k)))
                .collect();
            d_scored.push((a, mm, places));
        }
    }
    let (d_valid, d_set) = fewest(d_scored);

    // J: every allele at every shift of the true J index
    let mut j_scored = Vec::new();
    if !j_at.is_empty() {
        for a in &c.j {
            let q = a.sequence().into_bytes();
            let (mm, ks) = slide(read, &j_at, &q);
            j_scored.push((a, mm, ks.into_iter().map(|k| (q.clone(), k)).collect::<Vec<_>>()));
        }
    }
    let (j_valid, j_set) = fewest(j_scored);

    // C: prefix
    let mut c_scored = Vec::new();
    if !c_at.is_empty() {
        for a in &c.c {
            let q = a.sequence().into_bytes();
            c_scored.push((a, c_at.iter().filter(|&&(i, u)| score(read[i], q.get(u)) < 0).count(), ()));
        }
    }
    let (c_valid, _) = fewest(c_scored);

    // Placements. V by its true IMGT mapping, then the allele's own bases past the last
    // V-derived base; D and J at each tied best shift.
    let v_segs: Vec<Seg> = match v_at.last() {
        None => vec![],
        Some(&(pl, gl)) => v_set
            .iter()
            .map(|a| {
                let g = a.sequence_gapped.as_bytes();
                let mut scores: Vec<i32> = (0..=pl)
                    .map(|i| match s.origin[i] {
                        Origin::V(x) => score(read[i], g.get(x)),
                        _ => 0,
                    })
                    .collect();
                let next = g.iter().take(gl + 1).filter(|&&b| b != b'.').count();
                let ung = a.sequence().into_bytes();
                scores.extend(ung.iter().skip(next).zip(&read[pl + 1..]).map(|(b, &r)| score(r, Some(b))));
                Seg { start: 0, scores }
            })
            .collect(),
    };
    let d_segs: Vec<Seg> = d_set
        .iter()
        .flatten()
        .map(|(seq, k)| {
            // read position i holds index i − r0 + k
            let r0 = d_at[0].0 as i64;
            let start = (r0 - k).max(0);
            let end = (r0 - k + seq.len() as i64).min(len as i64);
            Seg { start: start as usize, scores: (start..end).map(|i| score(read[i as usize], at(seq, i - r0 + k))).collect() }
        })
        .collect();
    let j_segs: Vec<Seg> = match (j_at.first(), j_at.last()) {
        (Some(&(pf, uf)), Some(&(pl, _))) => j_set
            .iter()
            .flatten()
            .map(|(q, k)| {
                // J bases before the first J-derived one, as far as the allele and read go
                let ext = (uf as i64 + k).clamp(0, pf as i64) as usize;
                let scores = (pf - ext..=pl)
                    .map(|i| match s.origin[i] {
                        Origin::J(u) => score(read[i], at(q, u as i64 + k)),
                        _ if i >= pf => 0,
                        _ => score(read[i], at(q, uf as i64 + k - (pf - i) as i64)),
                    })
                    .collect();
                Seg { start: pf - ext, scores }
            })
            .collect(),
        _ => vec![],
    };

    // joint maximum over every combination; union the cuts of those reaching it
    let (vs, ds, js) = (options(&v_segs, !v_at.is_empty()), options(&d_segs, !d_at.is_empty()), options(&j_segs, !j_at.is_empty()));
    let mut found: Vec<(i32, [BTreeSet<usize>; 4])> = Vec::new();
    for &v in &vs {
        for &d in &ds {
            for &j in &js {
                found.extend(best(v, d, j, len));
            }
        }
    }
    let max = found.iter().map(|x| x.0).max();
    let mut cuts: [BTreeSet<usize>; 4] = Default::default();
    for (_, c) in found.into_iter().filter(|x| Some(x.0) == max) {
        for (all, these) in cuts.iter_mut().zip(c) {
            all.extend(these);
        }
    }
    let list = |k: usize, add: i64, present: bool| -> Vec<i64> {
        if present { cuts[k].iter().map(|&x| x as i64 + add).collect() } else { vec![] }
    };
    let unrec = |set: &[String], a: Option<&Allele>| a.is_some_and(|a| !set.is_empty() && !set.contains(&a.name));

    Truth {
        v_end_valid: list(0, 0, !v_at.is_empty()),
        d_start_valid: list(1, 1, !d_at.is_empty()),
        d_end_valid: list(2, 0, !d_at.is_empty()),
        j_start_valid: list(3, 1, !j_at.is_empty()),
        unrecoverable_v: unrec(&v_valid, Some(&s.v)),
        unrecoverable_d: unrec(&d_valid, s.d.as_ref()),
        unrecoverable_j: unrec(&j_valid, Some(&s.j)),
        d_null_valid: d_at.len() < MIN_D_LEN,
        v_valid,
        d_valid,
        j_valid,
        c_valid,
    }
}
