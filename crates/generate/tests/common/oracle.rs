//! Brute-force truth reference. Enumerates every candidate × orientation × offset and every
//! ordered (V end, D start, D end, J start) directly; shares no code with `src/truth.rs`.
//! The tuple enumeration skips a branch only when an upper bound on everything it could
//! still add (each remaining segment's unconstrained best) is below the best total
//! already found, so every tuple tying the maximum is visited.
use abstar_generate::{Candidates, Origin, Simulated, Truth};
use abstar_germline::Allele;

fn bytes(a: &Allele) -> Vec<u8> {
    a.sequence().into_bytes()
}

fn rc(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            x => *x,
        })
        .collect()
}

fn hit(read: u8, g: Option<&u8>) -> bool {
    g.is_some_and(|&g| g != b'.' && g.eq_ignore_ascii_case(&read))
}

fn pm(read: u8, g: Option<&u8>) -> i32 {
    if hit(read, g) { 1 } else { -1 }
}

/// Sorted, deduplicated names with the smallest score.
fn fewest(scores: &[(String, usize)]) -> Vec<String> {
    let Some(m) = scores.iter().map(|x| x.1).min() else { return vec![] };
    let mut v: Vec<String> = scores.iter().filter(|x| x.1 == m).map(|x| x.0.clone()).collect();
    v.sort();
    v.dedup();
    v
}

/// Every (name, sequence, shift, mismatches) whose mismatches are that allele's best,
/// for alleles in `set`.
fn tied(all: Vec<(String, Vec<u8>, i64, usize)>, set: &[String]) -> Vec<(String, Vec<u8>, i64)> {
    all.iter()
        .filter(|x| set.contains(&x.0))
        .filter(|x| x.3 == all.iter().filter(|y| y.0 == x.0).map(|y| y.3).min().unwrap())
        .map(|x| (x.0.clone(), x.1.clone(), x.2))
        .collect()
}

/// Items grouped by allele name.
fn by_name<T>(items: Vec<(String, T)>) -> Vec<Vec<T>> {
    let mut m: std::collections::BTreeMap<String, Vec<T>> = Default::default();
    for (n, t) in items {
        m.entry(n).or_default().push(t);
    }
    m.into_values().collect()
}

/// Every ordered (ve, ds, de, js) of one placement of each segment tying its best.
/// A branch is skipped only if what it has plus each remaining segment's
/// unconstrained best is below the best total found.
fn search(vo: &[(usize, i32)], d: Option<&[(usize, usize, i32)]>, jo: &[(usize, i32)]) -> (i32, Vec<[usize; 4]>) {
    let (mut best, mut ties) = (i32::MIN, Vec::<[usize; 4]>::new());
    let mut consider = |t: i32, x: [usize; 4], best: &mut i32| {
        if t > *best {
            *best = t;
            ties.clear();
        }
        if t == *best {
            ties.push(x);
        }
    };
    let maxj = jo.iter().map(|x| x.1).max().unwrap_or(i32::MIN / 4);
    let maxd = d.map_or(0, |d| d.iter().map(|x| x.2).max().unwrap_or(i32::MIN / 4));
    for &(ve, a) in vo.iter().rev() {
        if a + maxd + maxj < best {
            continue;
        }
        match d {
            Some(pairs) => {
                for &(st, en, dd) in pairs {
                    if st < ve || a + dd + maxj < best {
                        continue;
                    }
                    for &(jst, b) in jo {
                        if jst >= en {
                            consider(a + dd + b, [ve, st, en, jst], &mut best);
                        }
                    }
                }
            }
            None => {
                for &(jst, b) in jo {
                    if jst >= ve {
                        consider(a + b, [ve, 0, 0, jst], &mut best);
                    }
                }
            }
        }
    }
    (best, ties)
}

/// The configurations tying the best over all of a combination's placement pairs.
fn combo_ties(found: Vec<(i32, Vec<[usize; 4]>)>) -> Vec<[usize; 4]> {
    let Some(m) = found.iter().filter(|x| !x.1.is_empty()).map(|x| x.0).max() else { return vec![] };
    found.into_iter().filter(|x| x.0 == m).flat_map(|x| x.1).collect()
}

pub fn truth(s: &Simulated, c: &Candidates) -> Truth {
    let (read, org, l) = (&s.read, &s.origin, s.read.len());
    let mut vs = vec![];
    let mut ds = vec![];
    let mut js = vec![];
    let mut cs = vec![];
    for (i, &o) in org.iter().enumerate() {
        match o {
            Origin::V(g) => vs.push((i, g)),
            Origin::D(_) => ds.push(i),
            Origin::J(u) => js.push((i, u)),
            Origin::C(u) => cs.push((i, u)),
            _ => {}
        }
    }

    // allele sets
    let mut v_scores = vec![];
    if !vs.is_empty() {
        for a in &c.v {
            let g = a.sequence_gapped.as_bytes();
            v_scores.push((a.name.clone(), vs.iter().filter(|&&(i, x)| !hit(read[i], g.get(x))).count()));
        }
    }
    let v_valid = fewest(&v_scores);

    // D: shift = candidate index of the run's first base
    let mut d_all = vec![];
    if !ds.is_empty() {
        for a in &c.d {
            let mut orients = vec![bytes(a)];
            if c.d_inversion {
                orients.push(rc(&bytes(a)));
            }
            for o in orients {
                for shift in -(ds.len() as i64 - 1)..o.len() as i64 {
                    let mm = (0..ds.len())
                        .filter(|&x| {
                            let k = shift + x as i64;
                            !(k >= 0 && hit(read[ds[x]], o.get(k as usize)))
                        })
                        .count();
                    d_all.push((a.name.clone(), o.clone(), shift, mm));
                }
            }
        }
    }
    let d_valid = fewest(&d_all.iter().map(|x| (x.0.clone(), x.3)).collect::<Vec<_>>());
    let d_places = tied(d_all, &d_valid);

    // J: shift = candidate index minus true J index
    let mut j_all = vec![];
    if !js.is_empty() {
        let umin = js.iter().map(|x| x.1).min().unwrap() as i64;
        let umax = js.iter().map(|x| x.1).max().unwrap() as i64;
        for b in &c.j {
            let q = bytes(b);
            for sh in -umax..q.len() as i64 - umin {
                let at = |u: usize| u as i64 + sh;
                if !js.iter().any(|&(_, u)| (0..q.len() as i64).contains(&at(u))) {
                    continue;
                }
                let mm = js.iter().filter(|&&(i, u)| !(at(u) >= 0 && hit(read[i], q.get(at(u) as usize)))).count();
                j_all.push((b.name.clone(), q.clone(), sh, mm));
            }
        }
    }
    let j_valid = fewest(&j_all.iter().map(|x| (x.0.clone(), x.3)).collect::<Vec<_>>());
    let j_places = tied(j_all, &j_valid);

    let mut c_scores = vec![];
    if !cs.is_empty() {
        for a in &c.c {
            let q = bytes(a);
            c_scores.push((a.name.clone(), cs.iter().filter(|&&(i, u)| !hit(read[i], q.get(u))).count()));
        }
    }
    let c_valid = fewest(&c_scores);

    // per-position score of each placement; None = the segment cannot cover it
    let v_prof: Vec<Vec<Option<i32>>> = v_valid
        .iter()
        .map(|n| {
            let a = c.v.iter().find(|a| &a.name == n).unwrap();
            let g = a.sequence_gapped.as_bytes();
            let ung = bytes(a);
            let (pl, gl) = *vs.last().unwrap();
            let c0 = g[..(gl + 1).min(g.len())].iter().filter(|&&b| b != b'.').count();
            (0..l)
                .map(|i| {
                    if i <= pl {
                        Some(match org[i] {
                            Origin::V(x) => pm(read[i], g.get(x)),
                            _ => 0,
                        })
                    } else {
                        ung.get(c0 + i - pl - 1).map(|b| pm(read[i], Some(b)))
                    }
                })
                .collect()
        })
        .collect();
    let d_prof: Vec<(String, Vec<Option<i32>>)> = d_places
        .iter()
        .map(|(n, o, shift)| {
            let p = (0..l)
                .map(|i| {
                    let k = shift + i as i64 - ds[0] as i64;
                    (k >= 0 && (k as usize) < o.len()).then(|| pm(read[i], o.get(k as usize)))
                })
                .collect();
            (n.clone(), p)
        })
        .collect();
    let j_prof: Vec<(String, Vec<Option<i32>>)> = j_places
        .iter()
        .map(|(n, q, sh)| {
            let (pf, uf) = js[0];
            let jend = js.last().unwrap().0 + 1;
            let p = (0..l)
                .map(|i| {
                    if i >= jend {
                        None
                    } else if i >= pf {
                        Some(match org[i] {
                            Origin::J(u) => {
                                let k = u as i64 + sh;
                                if k >= 0 { pm(read[i], q.get(k as usize)) } else { -1 }
                            }
                            _ => 0,
                        })
                    } else {
                        let k = uf as i64 + sh - (pf - i) as i64;
                        (k >= 0 && (k as usize) < q.len()).then(|| pm(read[i], q.get(k as usize)))
                    }
                })
                .collect();
            (n.clone(), p)
        })
        .collect();

    // (cut, score) options per placement: V covers [0, ve), J covers [js, J end)
    let v_opts: Vec<Vec<(usize, i32)>> = if vs.is_empty() {
        vec![vec![(0, 0)]]
    } else {
        v_prof
            .iter()
            .map(|p| {
                let (mut out, mut sum) = (vec![], 0);
                for e in 1..=l {
                    let Some(x) = p[e - 1] else { break };
                    sum += x;
                    out.push((e, sum));
                }
                out
            })
            .collect()
    };
    // J and D options grouped by allele
    let j_opts: Vec<Vec<Vec<(usize, i32)>>> = if js.is_empty() {
        vec![vec![vec![(l, 0)]]]
    } else {
        by_name(j_prof)
            .iter()
            .map(|g| g.iter().map(|p| {
                let jend = js.last().unwrap().0 + 1;
                let (mut out, mut sum) = (vec![], 0);
                for st in (0..jend).rev() {
                    let Some(x) = p[st] else { break };
                    sum += x;
                    out.push((st, sum));
                }
                out
            }).collect())
            .collect()
    };
    let d_opts: Vec<Vec<Vec<(usize, usize, i32)>>> = by_name(d_prof)
            .iter()
            .map(|g| g.iter().map(|p| {
                let mut out = vec![];
                for a in 0..l {
                    for b in a + 1..=l {
                        if p[a..b].iter().all(Option::is_some) {
                            out.push((a, b, p[a..b].iter().map(|x| x.unwrap()).sum()));
                        }
                    }
                }
                out
            }).collect())
            .collect();

    // each combination of valid alleles keeps its own best; no-D combinations run when
    // the D run is short and contribute V ends and J starts only
    let (mut ties, mut no_d_ties) = (vec![], vec![]);
    for vo in &v_opts {
        for jg in &j_opts {
            if ds.len() < 5 {
                no_d_ties.extend(combo_ties(jg.iter().map(|jo| search(vo, None, jo)).collect()));
            }
            for dg in &d_opts {
                let found = jg.iter().flat_map(|jo| dg.iter().map(move |d| search(vo, Some(d), jo))).collect();
                ties.extend(combo_ties(found));
            }
        }
    }
    let with_ends: Vec<[usize; 4]> = ties.iter().chain(&no_d_ties).copied().collect();
    let list = |t: &[[usize; 4]], f: &dyn Fn(&[usize; 4]) -> i64, on: bool| {
        let mut v: Vec<i64> = if on { t.iter().map(f).collect() } else { vec![] };
        v.sort();
        v.dedup();
        v
    };
    let unrec = |set: &[String], a: Option<&Allele>| a.is_some_and(|a| !set.is_empty() && !set.contains(&a.name));
    Truth {
        v_end_valid: list(&with_ends, &|x| x[0] as i64, !vs.is_empty()),
        d_start_valid: list(&ties, &|x| x[1] as i64 + 1, !ds.is_empty()),
        d_end_valid: list(&ties, &|x| x[2] as i64, !ds.is_empty()),
        j_start_valid: list(&with_ends, &|x| x[3] as i64 + 1, !js.is_empty()),
        unrecoverable_v: unrec(&v_valid, Some(&s.v)),
        unrecoverable_d: unrec(&d_valid, s.d.as_ref()),
        unrecoverable_j: unrec(&j_valid, Some(&s.j)),
        d_null_valid: ds.len() < 5,
        v_valid,
        d_valid,
        j_valid,
        c_valid,
    }
}
