mod common;
use abstar_generate::{
    bind, load_model, preset, truth, Bound, Candidates, Dist, GenerateParams, Origin, Simulated, Simulator, TrimSet, Truth,
};
use abstar_germline::{Allele, Segment};
use common::*;
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;
use std::ops::Range;
use std::path::Path;

/// `short_d`: D trims up to 8 at each end, so some D runs are shorter than 5 nt.
fn bound(name: &str, alleles: Vec<Allele>, d_inversion: f64, short_d: bool) -> Bound {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join(format!("tests/data/{name}.model.toml"));
    let mut m = load_model(&path).unwrap();
    m.d_inversion = d_inversion;
    m.shm_model = "motif".into();
    if short_d {
        let t = TrimSet { default: Dist { values: (0..=8).collect(), probs: vec![1.0 / 9.0; 9] }, genes: Default::default() };
        (m.trim.d_5, m.trim.d_3) = (Some(t.clone()), Some(t));
    }
    bind(m, &db(alleles)).unwrap()
}

fn cands(alleles: &[Allele], d_inversion: bool) -> Candidates {
    let of = |s: Segment| alleles.iter().filter(|a| a.segment == s).cloned().collect();
    Candidates { v: of(Segment::V), d: of(Segment::D), j: of(Segment::J), c: of(Segment::C), d_inversion }
}

fn get(alleles: &[Allele], name: &str) -> Allele {
    alleles.iter().find(|a| a.name == name).unwrap().clone()
}

fn bytes(a: &Allele) -> Vec<u8> {
    a.sequence().into_bytes()
}

fn names(v: &[&str]) -> Vec<String> {
    v.iter().map(|s| s.to_string()).collect()
}

/// An unmutated read: V's first `v_len` bases, `n1`, then optionally D bases `d.1` and
/// `n2`, then J from index `j_trim`.
fn hand(v: &Allele, v_len: usize, n1: &[u8], d: Option<(&Allele, Range<usize>)>, n2: &[u8], j: &Allele, j_trim: usize) -> Simulated {
    let (mut read, mut origin) = (Vec::new(), Vec::new());
    let gapped: Vec<usize> = v.sequence_gapped.bytes().enumerate().filter(|x| x.1 != b'.').map(|x| x.0).collect();
    let vs = bytes(v);
    for u in 0..v_len {
        read.push(vs[u]);
        origin.push(Origin::V(gapped[u]));
    }
    let mut push = |bases: &[u8], o: Origin| {
        read.extend(bases);
        origin.extend(std::iter::repeat_n(o, bases.len()));
    };
    push(n1, Origin::N1);
    if let Some((da, r)) = &d {
        let ds = bytes(da);
        for x in r.clone() {
            push(&ds[x..x + 1], Origin::D(x));
        }
        push(n2, Origin::N2);
    }
    let js = bytes(j);
    for u in j_trim..js.len() {
        push(&js[u..u + 1], Origin::J(u));
    }
    Simulated {
        read,
        origin,
        v: v.clone(),
        d: d.as_ref().map(|x| x.0.clone()),
        d_inverted: false,
        j: j.clone(),
        c: None,
        c_len: 0,
        v_trim_3: (vs.len() - v_len) as i64,
        d_trim_5: 0,
        d_trim_3: 0,
        j_trim_5: j_trim as i64,
        np1: String::from_utf8(n1.to_vec()).unwrap(),
        np2: String::from_utf8(n2.to_vec()).unwrap(),
        mutations: vec![],
        indels: vec![],
        truncate_5: 0,
        truncate_3: 0,
        v_mutation_frequency: 0.0,
        warnings: vec![],
    }
}

/// A base differing from every base in `avoid`.
fn other(avoid: &[u8]) -> u8 {
    *b"ACGT".iter().find(|b| !avoid.contains(b)).unwrap()
}

/// Brute-force reference. Enumerates every candidate × orientation × offset and every
/// ordered (V end, D start, D end, J start) directly; shares no code with `truth.rs`.
/// The tuple enumeration skips a branch only when an upper bound on everything it could
/// still add (each remaining segment's unconstrained best) is below the best total
/// already found, so every tuple tying the maximum is visited.
mod oracle {
    use super::*;

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
    fn tied(all: Vec<(String, Vec<u8>, i64, usize)>, set: &[String]) -> Vec<(Vec<u8>, i64)> {
        all.iter()
            .filter(|x| set.contains(&x.0))
            .filter(|x| x.3 == all.iter().filter(|y| y.0 == x.0).map(|y| y.3).min().unwrap())
            .map(|x| (x.1.clone(), x.2))
            .collect()
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
        let d_prof: Vec<Vec<Option<i32>>> = d_places
            .iter()
            .map(|(o, shift)| {
                (0..l)
                    .map(|i| {
                        let k = shift + i as i64 - ds[0] as i64;
                        (k >= 0 && (k as usize) < o.len()).then(|| pm(read[i], o.get(k as usize)))
                    })
                    .collect()
            })
            .collect();
        let j_prof: Vec<Vec<Option<i32>>> = j_places
            .iter()
            .map(|(q, sh)| {
                let (pf, uf) = js[0];
                let jend = js.last().unwrap().0 + 1;
                (0..l)
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
                    .collect()
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
        let j_opts: Vec<Vec<(usize, i32)>> = if js.is_empty() {
            vec![vec![(l, 0)]]
        } else {
            j_prof
                .iter()
                .map(|p| {
                    let jend = js.last().unwrap().0 + 1;
                    let (mut out, mut sum) = (vec![], 0);
                    for st in (0..jend).rev() {
                        let Some(x) = p[st] else { break };
                        sum += x;
                        out.push((st, sum));
                    }
                    out
                })
                .collect()
        };
        let d_opts: Vec<Option<Vec<(usize, usize, i32)>>> = if ds.is_empty() {
            vec![None]
        } else {
            d_prof
                .iter()
                .map(|p| {
                    let mut out = vec![];
                    for a in 0..l {
                        for b in a + 1..=l {
                            if p[a..b].iter().all(Option::is_some) {
                                out.push((a, b, p[a..b].iter().map(|x| x.unwrap()).sum()));
                            }
                        }
                    }
                    Some(out)
                })
                .collect()
        };

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
        for vo in &v_opts {
            for d in &d_opts {
                for jo in &j_opts {
                    let maxj = jo.iter().map(|x| x.1).max().unwrap();
                    let maxd = d.as_ref().map_or(0, |d| d.iter().map(|x| x.2).max().unwrap());
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
                }
            }
        }
        let list = |f: &dyn Fn(&[usize; 4]) -> i64, on: bool| {
            let mut v: Vec<i64> = if on { ties.iter().map(f).collect() } else { vec![] };
            v.sort();
            v.dedup();
            v
        };
        let unrec = |set: &[String], a: Option<&Allele>| a.is_some_and(|a| !set.is_empty() && !set.contains(&a.name));
        Truth {
            v_end_valid: list(&|x| x[0] as i64, !vs.is_empty()),
            d_start_valid: list(&|x| x[1] as i64 + 1, !ds.is_empty()),
            d_end_valid: list(&|x| x[2] as i64, !ds.is_empty()),
            j_start_valid: list(&|x| x[3] as i64 + 1, !js.is_empty()),
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
}

#[test]
fn oracle_agreement() {
    let igh_inv = bound("tiny_igh", tricky_igh_alleles(), 0.3, false);
    let igh = bound("tiny_igh", tricky_igh_alleles(), 0.0, true);
    let igk = bound("tiny_igk", tricky_igk_alleles(), 0.0, false);
    let base = GenerateParams {
        mutation: "uniform:0-0.3".parse().unwrap(),
        indel_rate: 0.3,
        truncate_5p: "uniform:0-60".parse().unwrap(),
        truncate_3p: "uniform:0-40".parse().unwrap(),
        ..preset("naive").unwrap()
    };
    let mut n = [0usize; 12];
    for seed in 0..3000u64 {
        let b = [&igh_inv, &igh_inv, &igh, &igk][seed as usize % 4];
        let p = GenerateParams {
            constant_len: (seed % 2 == 0).then_some(30),
            frameshift_indels: seed % 8 < 4,
            ..base.clone()
        };
        let s = Simulator::new(b, &p).unwrap().simulate(&mut ChaCha8Rng::seed_from_u64(seed));
        let t = truth(&s, &b.candidates);
        assert_eq!(t, oracle::truth(&s, &b.candidates), "seed {seed}");
        for (k, x) in [
            t.v_valid.len() > 1,
            t.d_valid.len() > 1,
            t.j_valid.len() > 1,
            t.c_valid.len() > 1,
            t.unrecoverable_v,
            t.unrecoverable_d,
            t.unrecoverable_j,
            t.d_null_valid && s.d.is_some(),
            t.v_end_valid.len() > 1,
            t.d_start_valid.len() > 1 || t.d_end_valid.len() > 1,
            t.j_start_valid.len() > 1,
            !s.indels.is_empty(),
        ]
        .into_iter()
        .enumerate()
        {
            n[k] += usize::from(x);
        }
    }
    eprintln!(
        "ambiguous v/d/j/c {:?}, unrecoverable v/d/j {:?}, d_null {}, multi v_end/d/j_start {:?}, indels {}",
        &n[..4],
        &n[4..7],
        n[7],
        &n[8..11],
        n[11]
    );
    assert!(n.iter().all(|&x| x >= 10), "every case must occur: {n:?}");
}

#[test]
fn n_matching_germline_widens_boundary() {
    let al = tricky_igk_alleles();
    let (v, j) = (get(&al, "IGKV1-1*01"), get(&al, "IGKJ1*01"));
    let (vs, jq) = (bytes(&v), bytes(&j));
    // the N base equals both V's next germline base and the J base before J's first base,
    // so V and J can claim it equally; no further base can be claimed
    let t = (1..jq.len())
        .find(|&t| jq[t - 1] == vs[100] && jq[t] != vs[101] && (t < 2 || jq[t - 2] != vs[99]))
        .unwrap();
    let tr = truth(&hand(&v, 100, &[vs[100]], None, b"", &j, t), &cands(&al, false));
    assert_eq!(tr.v_end_valid, vec![100, 101]);
    assert_eq!(tr.j_start_valid, vec![101, 102]);

    // when only V can claim it, the matching N base moves V's end (+1 beats 0)
    let n = [vs[100], other(&[vs[101]])];
    let tr = truth(&hand(&v, 100, &n, None, b"", &j, 0), &cands(&al, false));
    assert_eq!(tr.v_end_valid, vec![101]);
    assert_eq!(tr.j_start_valid, vec![103]);
}

#[test]
fn mutated_edge_base_goes_to_n() {
    let al = tricky_igk_alleles();
    let (v, j) = (get(&al, "IGKV1-1*01"), get(&al, "IGKJ1*01"));
    let vs = bytes(&v);
    let mut s = hand(&v, 100, &[other(&[vs[100]])], None, b"", &j, 0);
    s.read[99] = other(&[vs[99]]);
    let tr = truth(&s, &cands(&al, false));
    assert_eq!(tr.v_end_valid, vec![99]);
    assert_eq!(tr.v_valid, names(&["IGKV1-1*01", "IGKV1-1*02"]));
}

#[test]
fn closer_allele_makes_unrecoverable() {
    let al = tricky_igh_alleles();
    let (v, d, j) = (get(&al, "IGHV1-1*01"), get(&al, "IGHD2-1*01"), get(&al, "IGHJ2*01"));
    let mut s = hand(&v, bytes(&v).len(), b"GG", Some((&d, 2..20)), b"TC", &j, 0);
    let c = cands(&al, false);
    let tr = truth(&s, &c);
    assert_eq!(tr.v_valid, names(&["IGHV1-1*01"]));
    assert!(!tr.unrecoverable_v);

    // mutate the read toward IGHV1-1*03 one differing position at a time
    let (g1, g3) = (v.sequence_gapped.as_bytes(), get(&al, "IGHV1-1*03").sequence_gapped.into_bytes());
    let diffs: Vec<usize> = (0..g1.len()).filter(|&g| g1[g] != g3[g]).collect();
    assert_eq!(diffs.len(), 2);
    let at = |s: &Simulated, g| s.origin.iter().position(|&o| o == Origin::V(g)).unwrap();
    let i = at(&s, diffs[0]);
    s.read[i] = g3[diffs[0]];
    let tr = truth(&s, &c);
    assert_eq!(tr.v_valid, names(&["IGHV1-1*01", "IGHV1-1*03"]));
    assert!(!tr.unrecoverable_v);
    let i = at(&s, diffs[1]);
    s.read[i] = g3[diffs[1]];
    let tr = truth(&s, &c);
    assert_eq!(tr.v_valid, names(&["IGHV1-1*03"]));
    assert!(tr.unrecoverable_v);
}

#[test]
fn short_d_null_valid() {
    let al = tricky_igh_alleles();
    let (v, d, j) = (get(&al, "IGHV1-1*01"), get(&al, "IGHD2-1*01"), get(&al, "IGHJ2*01"));
    let vl = bytes(&v).len() - 2;
    let c = cands(&al, false);
    for (r, null) in [(5..9, true), (5..10, false)] {
        let tr = truth(&hand(&v, vl, b"CC", Some((&d, r.clone())), b"CC", &j, 0), &c);
        assert_eq!(tr.d_null_valid, null, "D run {r:?}");
        assert!(tr.d_valid.contains(&"IGHD2-1*01".to_string()));
    }
    // a VJ locus has no D: null is the only valid D call
    let al = tricky_igk_alleles();
    let tr = truth(&hand(&get(&al, "IGKV1-1*01"), 200, b"A", None, b"", &get(&al, "IGKJ1*01"), 0), &cands(&al, false));
    assert!(tr.d_null_valid && tr.d_valid.is_empty() && tr.d_start_valid.is_empty() && tr.d_end_valid.is_empty());
}

#[test]
fn truncated_read_coverage() {
    let al = tricky_igh_alleles();
    let (v, d, j) = (get(&al, "IGHV1-1*01"), get(&al, "IGHD1-1*01"), get(&al, "IGHJ1*01"));
    let full = hand(&v, bytes(&v).len(), b"A", Some((&d, 0..17)), b"T", &j, 0);
    let c = cands(&al, false);
    // the 5'-partial *04 lacks bases the full read has: each counts as a mismatch
    assert_eq!(truth(&full, &c).v_valid, names(&["IGHV1-1*01"]));

    let mut tr = full.clone();
    tr.read.drain(..40);
    tr.origin.drain(..40);
    tr.truncate_5 = 40;
    // expected: fewest mismatches over the V positions the truncated read still covers
    let covered: Vec<(usize, usize)> = (0..tr.read.len())
        .filter_map(|i| match tr.origin[i] {
            Origin::V(g) => Some((i, g)),
            _ => None,
        })
        .collect();
    let mm = |a: &Allele| {
        let g = a.sequence_gapped.as_bytes();
        covered.iter().filter(|&&(i, x)| g.get(x).is_none_or(|&b| b != tr.read[i])).count()
    };
    let best = c.v.iter().map(mm).min().unwrap();
    let mut want: Vec<String> = c.v.iter().filter(|a| mm(a) == best).map(|a| a.name.clone()).collect();
    want.sort();
    let got = truth(&tr, &c).v_valid;
    assert_eq!(got, want);
    assert_eq!(got, names(&["IGHV1-1*01", "IGHV1-1*04"]));
}
