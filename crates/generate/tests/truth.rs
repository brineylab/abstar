mod common;
use abstar_generate::{
    bind, load_model, preset, truth, Bound, Candidates, Dist, GenerateParams, Origin, Simulated, Simulator, TrimSet,
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
    Candidates::new(of(Segment::V), of(Segment::D), of(Segment::J), of(Segment::C), d_inversion)
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

/// A V sibling missing the last 3 bases can't take an N base that matches the true V's
/// next germline base, so its own best V end (the true end) is valid too (ruling G7).
#[test]
fn short_v_sibling_keeps_its_own_best_end() {
    let al = tricky_igh_alleles();
    let (v, d, j) = (get(&al, "IGHV1-2*01"), get(&al, "IGHD2-1*01"), get(&al, "IGHJ2*01"));
    let (vs, ds) = (bytes(&v), bytes(&d));
    let l = vs.len();
    let s = hand(&v, l - 3, &[vs[l - 3], other(&[vs[l - 2], ds[1]])], Some((&d, 2..20)), b"TC", &j, 0);
    let c = cands(&al, false);
    let tr = truth(&s, &c);
    assert_eq!(tr.v_valid, names(&["IGHV1-2*01", "IGHV1-2*02"]));
    assert_eq!(tr.v_end_valid, vec![l as i64 - 3, l as i64 - 2]);
    assert_eq!(tr, oracle::truth(&s, &c));
}

/// The J analogue: a sibling missing the first 3 bases can't take an N base matching the
/// true J's previous germline base (ruling G7).
#[test]
fn short_j_sibling_keeps_its_own_best_start() {
    let al = tricky_igh_alleles();
    let (v, d, j) = (get(&al, "IGHV1-1*01"), get(&al, "IGHD2-1*01"), get(&al, "IGHJ1*01"));
    let (ds, jq) = (bytes(&d), bytes(&j));
    let n2 = [other(&[jq[1], ds[20]]), jq[2]];
    let s = hand(&v, bytes(&v).len() - 2, b"CC", Some((&d, 2..20)), &n2, &j, 3);
    let c = cands(&al, false);
    let tr = truth(&s, &c);
    let pf = s.origin.iter().position(|o| matches!(o, Origin::J(_))).unwrap() as i64;
    assert_eq!(tr.j_valid, names(&["IGHJ1*01", "IGHJ1*02", "IGHJ1*03"]));
    assert_eq!(tr.j_start_valid, vec![pf, pf + 1]);
    assert_eq!(tr, oracle::truth(&s, &c));
}

/// With a D run under 5 nt, the V ends and J starts of the best no-D configurations are
/// valid too (ruling G8): here only no-D lets V keep all three D bases.
#[test]
fn short_d_adds_no_d_boundaries() {
    let al = tricky_igh_alleles();
    let (v, j) = (get(&al, "IGHV1-1*01"), get(&al, "IGHJ2*01"));
    let vs = bytes(&v);
    let l = vs.len();
    // a D ending in V's last 3 bases (so it cannot extend 3'), led by a base they lack
    let tail = &vs[l - 3..];
    let b = other(tail) as char;
    let core = String::from_utf8(tail.to_vec()).unwrap();
    let d = Allele {
        sequence_gapped: format!("{b}{b}{b}{b}{core}"),
        ..allele("IGHD9-9*01", abstar_germline::Locus::IGH, Segment::D, "human", false)
    };
    let n2 = [other(&[b as u8]), other(&[b as u8])];
    let s = hand(&v, l - 3, b"", Some((&d, 4..7)), &n2, &j, 0);
    let c = cands(&[v.clone(), d.clone(), j.clone()], false);
    let tr = truth(&s, &c);
    assert!(tr.d_null_valid);
    let l = l as i64;
    // with D: V ends by l − 1 (D needs a base); without D: V takes all three
    assert_eq!(tr.v_end_valid, vec![l - 3, l - 2, l - 1, l]);
    assert_eq!(tr, oracle::truth(&s, &c));
}

#[test]
#[ignore = "timing: cargo test --release -p abstar_generate --test truth -- --ignored --nocapture"]
fn truth_speed_human_igh() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../python/abstar/germline_dbs/human");
    let al: Vec<Allele> = abstar_germline::load(&dir)
        .unwrap()
        .alleles
        .into_iter()
        .filter(|a| a.species == "human" && a.locus == abstar_germline::Locus::IGH)
        .collect();
    let c = cands(&al, true);
    let (d, j) = (get(&al, "IGHD3-10*01"), get(&al, "IGHJ4*02"));
    let cases = [
        ("typical (15 nt D)", "IGHV3-23*01", 3..18, None),
        ("3 nt D", "IGHV3-23*01", 10..13, None),
        ("1 nt D", "IGHV3-23*01", 12..13, None),
        // 5'-truncated to a few V bases: large, ambiguous V sets
        ("1 nt D, IGHV1-69*01 60 V nt", "IGHV1-69*01", 12..13, Some(60)),
        ("1 nt D, IGHV1-69*01 25 V nt", "IGHV1-69*01", 12..13, Some(25)),
        ("1 nt D, IGHV3-30*01 25 V nt", "IGHV3-30*01", 12..13, Some(25)),
        ("1 nt D, IGHV3-30*01 10 V nt", "IGHV3-30*01", 12..13, Some(10)),
    ];
    for (name, vname, r, keep) in cases {
        let v = get(&al, vname);
        let vl = bytes(&v).len() - 2;
        let mut s = hand(&v, vl, b"GGT", Some((&d, r)), b"CCA", &j, 4);
        for i in (7..vl).step_by(20) {
            s.read[i] = other(&[s.read[i]]); // 5% V mutation
        }
        if let Some(k) = keep {
            s.read.drain(..vl - k);
            s.origin.drain(..vl - k);
            s.truncate_5 = vl - k;
        }
        let n = 500;
        let t0 = std::time::Instant::now();
        for _ in 0..n {
            std::hint::black_box(truth(std::hint::black_box(&s), &c));
        }
        let us = t0.elapsed().as_secs_f64() * 1e6 / n as f64;
        eprintln!("{name}: {us:.1} µs/read (V set {})", truth(&s, &c).v_valid.len());
        assert!(us < 500.0, "{name}: {us:.1} µs/read");
    }
}
