mod common;
use abstar_generate::{
    bind, load_model, preset, productivity, replay, Bound, Dist, GenerateParams, IndelKind, Origin, Simulated,
    Simulator, Spec, TrimSet,
};
use abstar_germline::{Allele, Segment};
use common::*;
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;
use std::path::Path;

fn model(name: &str) -> abstar_generate::Model {
    load_model(&Path::new(env!("CARGO_MANIFEST_DIR")).join(format!("tests/data/{name}.model.toml"))).unwrap()
}

fn igh() -> Bound {
    let mut m = model("tiny_igh");
    m.d_inversion = 0.3;
    bind(m, &db(sim_igh_alleles())).unwrap()
}

fn igk() -> Bound {
    let mut m = model("tiny_igk"); // no shm_model field: defaults to "none"
    m.shm_model = "motif".into();
    bind(m, &db(sim_igk_alleles())).unwrap()
}

fn naive() -> GenerateParams {
    preset("naive").unwrap()
}

fn busy() -> GenerateParams {
    GenerateParams {
        mutation: "uniform:0-0.3".parse().unwrap(),
        indel_rate: 0.5,
        frameshift_indels: true,
        truncate_5p: "uniform:0-40".parse().unwrap(),
        truncate_3p: "uniform:0-40".parse().unwrap(),
        constant_len: Some(30),
        ..naive()
    }
}

fn sim(b: &Bound, p: &GenerateParams, seed: u64) -> Simulated {
    Simulator::new(b, p).unwrap().simulate(&mut ChaCha8Rng::seed_from_u64(seed))
}

fn revcomp(s: &[u8]) -> Vec<u8> {
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

fn seq(a: &Allele) -> Vec<u8> {
    a.sequence().into_bytes()
}

/// Gapped index of each ungapped base.
fn gapped_index(a: &Allele) -> Vec<usize> {
    a.sequence_gapped.bytes().enumerate().filter(|&(_, c)| c != b'.').map(|(i, _)| i).collect()
}

fn pos_of(s: &Simulated, o: Origin) -> Option<usize> {
    s.origin.iter().position(|&x| x == o)
}

#[test]
fn replay_rebuilds_read() {
    let (h, k) = (igh(), igk());
    let (mut indels, mut inverted, mut mutated, mut dups) = (0, 0, 0, 0);
    for seed in 0..2000 {
        let p = GenerateParams { indel_duplication: seed % 2 == 1, ..busy() };
        let b = if seed % 4 == 3 { &k } else { &h };
        let s = sim(b, &p, seed);
        assert_eq!(replay(&s), s.read, "seed {seed}");
        assert_eq!(s.origin.len(), s.read.len());
        indels += usize::from(!s.indels.is_empty());
        dups += usize::from(p.indel_duplication && s.indels.iter().any(|i| i.kind == IndelKind::Ins));
        inverted += usize::from(s.d_inverted);
        mutated += usize::from(!s.mutations.is_empty());
    }
    assert!(indels > 800 && dups > 100 && inverted > 300 && mutated > 1900, "{indels} {dups} {inverted} {mutated}");
}

#[test]
fn origin_map_is_exact() {
    let (h, k) = (igh(), igk());
    for seed in 0..2000 {
        let p = GenerateParams { indel_duplication: seed % 2 == 1, ..busy() };
        let s = sim(if seed % 4 == 3 { &k } else { &h }, &p, seed);
        let (v, j) = (s.v.sequence_gapped.as_bytes(), seq(&s.j));
        let d = s.d.as_ref().map(|d| if s.d_inverted { revcomp(&seq(d)) } else { seq(d) });
        let c = s.c.as_ref().map(seq);
        let mutated: Vec<usize> =
            s.mutations.iter().filter_map(|m| m.pos.checked_sub(s.truncate_5)).collect();
        let (mut n1, mut n2) = (Vec::new(), Vec::new());
        for (i, (&o, &base)) in s.origin.iter().zip(&s.read).enumerate() {
            let want = match o {
                Origin::V(g) => v[g],
                Origin::D(x) => d.as_ref().unwrap()[x],
                Origin::J(u) => j[u],
                Origin::C(u) => c.as_ref().unwrap()[u],
                Origin::Inserted => {
                    let ins = s.indels.iter().find(|x| x.kind == IndelKind::Ins).unwrap();
                    ins.bases.as_bytes()[i + s.truncate_5 - ins.read_pos]
                }
                Origin::N1 | Origin::N2 | Origin::P => {
                    if o == Origin::N1 { n1.push(base) } else if o == Origin::N2 { n2.push(base) }
                    continue;
                }
            };
            if !mutated.contains(&i) {
                assert_eq!(base, want, "seed {seed} pos {i} {o:?}");
            }
            if matches!(o, Origin::C(_) | Origin::Inserted) {
                assert!(!mutated.contains(&i), "seed {seed}: {o:?} mutated");
            }
        }
        // untruncated N runs equal the recorded N bases, barring mutation
        let untouched = s.mutations.is_empty() && s.truncate_5 == 0 && s.truncate_3 == 0;
        if untouched {
            assert_eq!(n1, s.np1.as_bytes(), "seed {seed}");
            assert_eq!(n2, s.np2.as_bytes(), "seed {seed}");
        }
        // germline indices strictly increase along the read; jumps only at a deletion
        let idx = |o: &Origin| match *o {
            Origin::V(g) => Some((0, g)),
            Origin::D(x) => Some((1, x)),
            Origin::J(u) => Some((2, u)),
            Origin::C(u) => Some((3, u)),
            _ => None,
        };
        let germ: Vec<(u8, usize)> = s.origin.iter().filter_map(idx).collect();
        assert!(germ.windows(2).all(|w| w[0] < w[1]), "seed {seed}");
    }
}

#[test]
fn p_nucleotides_are_palindromic() {
    let mut m = model("tiny_igh");
    let fixed = |v| TrimSet { default: Dist { values: vec![v], probs: vec![1.0] }, genes: Default::default() };
    m.trim.v_3 = fixed(-2);
    m.trim.j_5 = fixed(-3);
    m.trim.d_5 = Some(fixed(-2));
    m.trim.d_3 = Some(fixed(-2));
    let b = bind(m, &db(sim_igh_alleles())).unwrap();
    for seed in 0..50 {
        let s = sim(&b, &naive(), seed);
        assert_eq!((s.v_trim_3, s.d_trim_5, s.d_trim_3, s.j_trim_5), (-2, -2, -2, -3));
        assert!(!s.d_inverted);
        let d = seq(s.d.as_ref().unwrap());
        let (d0, dl) = (pos_of(&s, Origin::D(0)).unwrap(), pos_of(&s, Origin::D(d.len() - 1)).unwrap());
        assert_eq!(s.read[d0 - 2..d0], revcomp(&d[..2])[..], "seed {seed}");
        assert_eq!(s.read[dl + 1..dl + 3], revcomp(&d[d.len() - 2..])[..], "seed {seed}");
        assert_eq!(s.origin[d0 - 2..d0], [Origin::P; 2]);
        assert_eq!(s.origin[dl + 1..dl + 3], [Origin::P; 2]);
        let v = seq(&s.v);
        let last_v = s.origin.iter().rposition(|o| matches!(o, Origin::V(_))).unwrap();
        assert_eq!(s.read[last_v + 1..last_v + 3], revcomp(&v[v.len() - 2..])[..], "seed {seed}");
        assert_eq!(s.origin[last_v + 1..last_v + 3], [Origin::P, Origin::P]);
        let first_j = pos_of(&s, Origin::J(0)).unwrap();
        assert_eq!(s.read[first_j - 3..first_j], revcomp(&seq(&s.j)[..3])[..], "seed {seed}");
        assert_eq!(s.origin[first_j - 3..first_j], [Origin::P; 3]);
    }
}

#[test]
fn anchors_never_touched_by_indels() {
    let (h, k) = (igh(), igk());
    let p = GenerateParams { indel_rate: 1.0, frameshift_indels: true, ..naive() };
    let mut with_indel = 0;
    for seed in 0..2000 {
        let s = sim(if seed % 4 == 3 { &k } else { &h }, &p, seed);
        with_indel += usize::from(!s.indels.is_empty());
        let (a, b) = (s.v.anchor.unwrap(), s.j.anchor.unwrap());
        let ga = gapped_index(&s.v)[a];
        let pv = pos_of(&s, Origin::V(ga)).unwrap();
        let want: Vec<Origin> = (ga - 3..ga + 3).map(Origin::V).collect();
        assert_eq!(s.origin[pv - 3..pv + 3], want[..], "seed {seed}: {:?}", s.indels);
        assert_eq!(s.read[pv..pv + 3], seq(&s.v)[a..a + 3], "seed {seed}");
        let pj = pos_of(&s, Origin::J(b)).unwrap();
        let want: Vec<Origin> = (b..b + 6).map(Origin::J).collect();
        assert_eq!(s.origin[pj..pj + 6], want[..], "seed {seed}: {:?}", s.indels);
        assert_eq!(s.read[pj..pj + 3], seq(&s.j)[b..b + 3], "seed {seed}");
        for i in &s.indels {
            let (lo, hi) = (i.germline_pos, i.germline_pos + if i.kind == IndelKind::Del { i.len } else { 0 });
            assert!(if i.read_pos < pv { hi + 3 <= a } else { lo >= b + 6 }, "seed {seed}: {i:?}");
        }
    }
    assert!(with_indel > 1900, "{with_indel}");
}

#[test]
fn indel_lengths_and_equivalence() {
    let h = igh();
    let p = GenerateParams { indel_rate: 1.0, ..naive() };
    let (mut n, mut in_cdr) = (0, 0);
    for seed in 0..2000 {
        let s = sim(&h, &p, seed);
        let pv = pos_of(&s, Origin::V(gapped_index(&s.v)[s.v.anchor.unwrap()])).unwrap();
        for i in &s.indels {
            assert!(i.len % 3 == 0 && (3..=36).contains(&i.len), "{i:?}");
            assert_eq!(i.bases.len(), i.len);
            assert!(i.equiv.0 <= i.read_pos && i.read_pos <= i.equiv.1, "{i:?}");
            let g = if i.read_pos < pv { gapped_index(&s.v)[i.germline_pos] } else { 0 };
            n += 1;
            in_cdr += usize::from((69..123).contains(&g) || (156..204).contains(&g));
        }
    }
    // CDR1/CDR2 ± 3 codons hold ~84 of ~300 eligible bases: ~28% unweighted, ~53% at 3x
    let f = in_cdr as f64 / n as f64;
    assert!((0.45..0.62).contains(&f), "{f}");
    // the equivalence range is exactly the set of placements giving the same sequence
    let p = GenerateParams { indel_rate: 1.0, frameshift_indels: true, indel_duplication: true, ..naive() };
    for seed in 0..300 {
        let s = sim(&h, &p, seed);
        let Some(i) = s.indels.first() else { continue };
        let full = replay(&s);
        let mut before = full.clone();
        match i.kind {
            IndelKind::Ins => drop(before.drain(i.read_pos..i.read_pos + i.len)),
            IndelKind::Del => drop(before.splice(i.read_pos..i.read_pos, i.bases.bytes())),
        }
        // placing the indel at r gives the same read
        let same = |r: usize| match i.kind {
            IndelKind::Del => before[..r] == full[..r] && before[r + i.len..] == full[r..],
            IndelKind::Ins => full[..r] == before[..r] && full[r + i.len..] == before[r..],
        };
        let lim = before.len() - if i.kind == IndelKind::Del { i.len } else { 0 };
        let range: Vec<usize> = (0..=lim).filter(|&r| same(r)).collect();
        assert_eq!((range[0], *range.last().unwrap()), i.equiv, "seed {seed}: {i:?}");
        assert_eq!(range.len(), i.equiv.1 - i.equiv.0 + 1, "contiguous");
    }
}

#[test]
fn truncation_longer_than_read() {
    let h = igh();
    for (t5, t3) in [("1000", "0"), ("0", "1000"), ("1000", "1000")] {
        let p = GenerateParams {
            truncate_5p: t5.parse().unwrap(),
            truncate_3p: t3.parse().unwrap(),
            constant_len: Some(20),
            ..busy()
        };
        let s = sim(&h, &p, 1);
        assert!(s.read.is_empty() && s.origin.is_empty());
        assert_eq!(replay(&s), s.read);
        let pr = productivity(&s.read, &s.origin, &s.v, &s.j);
        assert!(pr.junction.is_none() && !pr.productive);
    }
}

#[test]
fn spec_parsing() {
    assert_eq!("0.1".parse::<Spec>().unwrap(), Spec::Fixed(0.1));
    assert_eq!("fixed:0.1".parse::<Spec>().unwrap(), Spec::Fixed(0.1));
    assert_eq!("uniform:0.15-0.35".parse::<Spec>().unwrap(), Spec::Uniform(0.15, 0.35));
    assert_eq!("beta:2,10".parse::<Spec>().unwrap(), Spec::Beta(2.0, 10.0));
    for bad in ["", "abc", "nan", "fixed:", "fixed:x", "uniform:0.1", "uniform:0.3-0.1", "beta:2", "beta:0,1", "gamma:1"] {
        let e = bad.parse::<Spec>().unwrap_err();
        assert_eq!(e.code(), "INVALID_ARGUMENT", "{bad}");
    }
    assert_eq!(preset("naive").unwrap().mutation, Spec::Fixed(0.0));
    assert_eq!(preset("tcr").unwrap().mutation, Spec::Fixed(0.0));
    assert_eq!(preset("memory").unwrap().mutation, Spec::Uniform(0.02, 0.12));
    let b = preset("bnab").unwrap();
    assert_eq!((b.mutation, b.indel_rate), (Spec::Uniform(0.15, 0.35), 0.108));
    assert!(b.long_cdr3_weight > 0.0);
    let e = preset("bnAb2").unwrap_err();
    assert_eq!(e.code(), "INVALID_ARGUMENT");
    assert!(e.to_string().contains("naive, memory, bnab, tcr"), "{e}");
}

/// Read = V up to its Cys codon + `n` + J from its Trp codon on.
fn hand_built(n: &str, cys: &str) -> (Vec<u8>, Vec<Origin>, Allele, Allele) {
    let alleles = sim_igh_alleles();
    let v = alleles.iter().find(|a| a.name == "IGHV1-2*01").unwrap().clone();
    let j = alleles.iter().find(|a| a.name == "IGHJ1*01").unwrap().clone();
    let (a, b) = (v.anchor.unwrap(), j.anchor.unwrap());
    let g = gapped_index(&v);
    let mut read = seq(&v)[..a].to_vec();
    read.extend(cys.bytes());
    read.extend(n.bytes());
    read.extend(&seq(&j)[b..]);
    let mut origin: Vec<Origin> = (0..a + 3).map(|u| Origin::V(g[u])).collect();
    origin.extend(n.bytes().map(|_| Origin::N1));
    origin.extend((b..seq(&j).len()).map(Origin::J));
    (read, origin, v, j)
}

#[test]
fn productive_flags_match_translation() {
    let check = |n: &str, cys: &str, want: (bool, bool, bool)| {
        let (read, origin, v, j) = hand_built(n, cys);
        let p = productivity(&read, &origin, &v, &j);
        let a = v.anchor.unwrap();
        assert_eq!(p.junction, Some(a..a + 6 + n.len()), "{n}");
        assert_eq!((p.vj_in_frame, p.stop_codon, p.productive), want, "{n} {cys}");
    };
    check("A", "TGT", (false, false, false)); // 1 nt junction insert: frameshift
    check("GCT", "TGT", (true, false, true)); // 3 nt: in frame, Ala
    check("TAG", "TGT", (true, true, false)); // in-frame stop
    check("GCT", "TCT", (true, false, false)); // Cys104 lost
}

#[test]
fn frameshift_in_v_is_not_in_frame() {
    // a 1 nt deletion upstream of Cys keeps the junction length a multiple of 3
    // but puts the J anchor out of the V frame
    let (mut read, mut origin, v, j) = hand_built("GCT", "TGT");
    read.remove(30);
    origin.remove(30);
    let p = productivity(&read, &origin, &v, &j);
    assert_eq!(p.junction.as_ref().map(|r| r.len()), Some(9));
    assert!(!p.vj_in_frame && !p.productive);
}

#[test]
fn same_seed_same_sequence() {
    let (h, p) = (igh(), busy());
    for seed in 0..20 {
        assert_eq!(sim(&h, &p, seed), sim(&h, &p, seed));
    }
    assert_ne!(sim(&h, &p, 1).read, sim(&h, &p, 2).read);
}

#[test]
fn realized_v_frequency_and_cap_warning() {
    let h = igh();
    let p = GenerateParams { mutation: Spec::Fixed(0.1), ..naive() };
    for seed in 0..50 {
        let s = sim(&h, &p, seed);
        let nv = s.origin.iter().filter(|o| matches!(o, Origin::V(_))).count() as f64;
        assert!((s.v_mutation_frequency - (0.1 * nv).round() / nv).abs() < 1e-12, "seed {seed}");
        assert!(s.warnings.is_empty());
    }
    let s = sim(&h, &GenerateParams { mutation: Spec::Fixed(1.0), ..naive() }, 1);
    assert!(s.v_mutation_frequency < 1.0);
    assert!(s.warnings.iter().any(|w| w.contains("safety cap")), "{:?}", s.warnings);
}

#[test]
fn mutation_needs_a_mutation_model() {
    let raw = bind(model("tiny_igk"), &db(sim_igk_alleles())).unwrap();
    let p = GenerateParams { mutation: Spec::Fixed(0.1), ..naive() };
    let e = Simulator::new(&raw, &p).err().unwrap();
    assert_eq!(e.code(), "INVALID_ARGUMENT");
    assert!(e.to_string().contains("shm_model = \"motif\""), "{e}");
    for preset_name in ["memory", "bnab"] {
        assert!(Simulator::new(&raw, &preset(preset_name).unwrap()).is_err(), "{preset_name}");
    }
    assert!(Simulator::new(&raw, &naive()).is_ok());
    let p = GenerateParams { shm_model: Some("motif".into()), ..p };
    assert!(!sim(&raw, &p, 1).mutations.is_empty());
    let p = GenerateParams { shm_model: Some("none".into()), ..p };
    assert!(Simulator::new(&igh(), &p).is_err());
}

#[test]
fn mutation_targets_final_read_divergence() {
    let h = igh();
    let p = GenerateParams { mutation: Spec::Fixed(0.2), truncate_5p: Spec::Fixed(60.0), ..naive() };
    for seed in 0..200 {
        let s = sim(&h, &p, seed);
        assert_eq!(s.truncate_5, 60);
        let nv = s.origin.iter().filter(|o| matches!(o, Origin::V(_))).count() as f64;
        // exactly the target, to mutate's rounding: round(0.2 nV) / nV
        let f = s.v_mutation_frequency;
        assert!((f - (0.2 * nv).round() / nv).abs() < 1e-12, "seed {seed}: {f}");
        assert!(f >= 0.2 - 0.5 / nv && f < 0.2 + 2.0 / nv, "seed {seed}: {f}");
        assert!(s.mutations.iter().all(|m| m.pos >= 60 && m.pos < 60 + s.read.len()));
    }
}

#[test]
fn isotype_weights_choose_constant() {
    let p = GenerateParams {
        constant_len: Some(20),
        isotype_weights: [("IGHG1".to_string(), 1.0)].into(),
        ..naive()
    };
    for seed in 0..20 {
        let s = sim(&igh(), &p, seed);
        assert_eq!(s.c.as_ref().unwrap().gene, "IGHG1");
        let want: Vec<Origin> = (0..20).map(Origin::C).collect();
        assert_eq!(s.origin[s.origin.len() - 20..], want[..]);
    }
}

#[test]
fn long_cdr3_weight_uses_long_tail() {
    let p = GenerateParams { long_cdr3_weight: 1.0, ..naive() };
    for seed in 0..50 {
        let s = sim(&igh(), &p, seed);
        assert!([10, 15, 20].contains(&s.np1.len()) && [10, 15, 20].contains(&s.np2.len()));
    }
}

#[test]
fn simulator_rejects_bad_params() {
    let h = igh();
    let bad = [
        GenerateParams { shm_model: Some("s5f".into()), ..naive() },
        GenerateParams { constant_len: Some(10), isotype_weights: [("IGHA1".to_string(), 1.0)].into(), ..naive() },
        GenerateParams { indel_rate: 1.5, ..naive() },
        GenerateParams { long_cdr3_weight: -0.1, ..naive() },
        GenerateParams { mutation: Spec::Fixed(1.5), ..naive() },
        GenerateParams { truncate_5p: Spec::Fixed(-1.0), ..naive() },
        GenerateParams { shm_table: Some("/nonexistent/s5f".into()), mutation: Spec::Fixed(0.1), ..naive() },
    ];
    for p in bad {
        let e = Simulator::new(&h, &p).err().unwrap_or_else(|| panic!("accepted {p:?}"));
        assert_eq!(e.code(), "INVALID_ARGUMENT", "{e}");
    }
    let no_c = bind(model("tiny_igh"), &db(sim_igh_alleles().into_iter().filter(|a| a.segment != Segment::C).collect()))
        .unwrap();
    let e = Simulator::new(&no_c, &GenerateParams { constant_len: Some(10), ..naive() }).err().unwrap();
    assert_eq!(e.code(), "INVALID_ARGUMENT");
}
