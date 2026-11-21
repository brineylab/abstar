//! Validation of the generator itself (spec: *Validation of the generator itself*) on the
//! built-in human models and database. Every test uses fixed seeds, so it cannot flake.
//! Too slow for debug builds; CI runs
//! `cargo test -p abstar_generate --release --test validation -- --include-ignored`.
mod common;
use std::collections::{BTreeSet, HashMap};
use std::path::Path;

use abstar_generate::{
    bind, generate_to_files, load_model, preset, productivity, replay, seed_for, truth, Bound, Dist, GenerateMeta,
    GenerateParams, Simulated, Simulator,
};
use common::oracle;
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;
use rayon::prelude::*;

fn bound(name: &str) -> Bound {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../python/abstar");
    let db = abstar_germline::load(&root.join("germline_dbs/human")).unwrap();
    bind(load_model(&root.join("models").join(name).join("model.toml")).unwrap(), &db).unwrap()
}

/// Read `i` of a run, seeded exactly as the pipeline seeds it.
fn sim(s: &Simulator, seed: u64, i: usize) -> Simulated {
    s.simulate(&mut ChaCha8Rng::seed_from_u64(seed_for(seed, i as u64)))
}

/// P(χ²_k ≥ x) = Q(k/2, x/2), the regularized upper incomplete gamma function
/// (series / continued fraction, Numerical Recipes §6.2). k/2 is a half-integer, so
/// ln Γ(k/2) is exact by Γ(a+1) = aΓ(a) from Γ(1) = 1 or Γ(1/2) = √π.
fn chi2_sf(x: f64, k: usize) -> f64 {
    let (a, x) = (k as f64 / 2.0, x / 2.0);
    let (mut lg, mut t) = if k.is_multiple_of(2) { (0.0, 1.0) } else { (0.5 * std::f64::consts::PI.ln(), 0.5) };
    while t < a {
        lg += f64::ln(t);
        t += 1.0;
    }
    let front = (a * x.ln() - x - lg).exp();
    if x < a + 1.0 {
        let (mut ap, mut del, mut sum) = (a, 1.0 / a, 1.0 / a);
        while del.abs() > sum.abs() * 1e-15 {
            ap += 1.0;
            del *= x / ap;
            sum += del;
        }
        1.0 - sum * front
    } else {
        let tiny = 1e-300;
        let (mut b, mut c) = (x + 1.0 - a, 1.0 / tiny);
        let mut d = 1.0 / b;
        let mut h = d;
        for i in 1..10_000 {
            let an = -(i as f64) * (i as f64 - a);
            b += 2.0;
            d = 1.0 / Some(an * d + b).filter(|v| v.abs() >= tiny).unwrap_or(tiny);
            c = Some(b + an / c).filter(|v| v.abs() >= tiny).unwrap_or(tiny);
            h *= d * c;
            if (d * c - 1.0).abs() < 1e-15 {
                break;
            }
        }
        front * h
    }
}

/// (table, category) → count or probability.
type Table = HashMap<(String, String), f64>;

fn add(t: &mut Table, table: &str, cat: impl ToString, x: f64) {
    *t.entry((table.to_string(), cat.to_string())).or_default() += x;
}

fn merge(mut a: Table, b: Table) -> Table {
    for (k, v) in b {
        *a.entry(k).or_default() += v;
    }
    a
}

/// Trim value `d` maps to after the sampler clamps it to a segment of `len` nt.
fn clamped(d: &Dist, len: usize) -> impl Iterator<Item = (i64, f64)> + '_ {
    d.values.iter().zip(&d.probs).map(move |(&v, &p)| (v.clamp(-(len as i64), len as i64), p))
}

/// Every table's distribution implied by the bound model, with the sampler's clamping
/// (a D 3' trim is clamped to what the 5' trim leaves).
fn expected(b: &Bound) -> Table {
    let (m, mut t) = (&b.model, Table::new());
    let gene_dist = |ts: &abstar_generate::TrimSet, g: &str| ts.genes.get(g).unwrap_or(&ts.default).clone();
    for (a, p) in &b.v {
        add(&mut t, "v", &a.name, *p);
        for (x, q) in clamped(&gene_dist(&m.trim.v_3, &a.gene), a.sequence().len()) {
            add(&mut t, "v_3", x, p * q);
        }
    }
    for (a, p) in &b.j {
        add(&mut t, "j", &a.name, *p);
        for (x, q) in clamped(&gene_dist(&m.trim.j_5, &a.gene), a.sequence().len()) {
            add(&mut t, "j_5", x, p * q);
        }
        let Some(ds) = b.d_given_j.get(&a.gene) else { continue };
        for (d, pd) in ds {
            add(&mut t, &format!("d|{}", a.gene), &d.name, *pd * p); // rescaled per table below
            let len = d.sequence().len();
            let (d5, d3) = (m.trim.d_5.as_ref().unwrap(), m.trim.d_3.as_ref().unwrap());
            for (x5, q5) in clamped(&gene_dist(d5, &d.gene), len) {
                add(&mut t, "d_5", x5, p * pd * q5);
                for (x3, q3) in clamped(&gene_dist(d3, &d.gene), len - x5.max(0) as usize) {
                    add(&mut t, "d_3", x3, p * pd * q5 * q3);
                }
            }
        }
    }
    for (key, d) in &m.insertions {
        for (&v, &p) in d.values.iter().zip(&d.probs) {
            add(&mut t, &format!("n_{key}"), v.max(0), p);
        }
    }
    for (prev, row) in &m.n_markov {
        for (base, p) in "ACGT".chars().zip(row) {
            add(&mut t, &format!("n_after_{prev}"), base, *p);
        }
    }
    t
}

fn by_table(t: &Table) -> HashMap<&str, Vec<(String, f64)>> {
    let mut out: HashMap<&str, Vec<(String, f64)>> = HashMap::new();
    for ((name, cat), x) in t {
        out.entry(name.as_str()).or_default().push((cat.clone(), *x));
    }
    out
}

fn observe(t: &mut Table, s: &Simulated) {
    add(t, "v", &s.v.name, 1.0);
    add(t, "j", &s.j.name, 1.0);
    add(t, "v_3", s.v_trim_3, 1.0);
    add(t, "j_5", s.j_trim_5, 1.0);
    if let Some(d) = &s.d {
        add(t, &format!("d|{}", s.j.gene), &d.name, 1.0);
        add(t, "d_5", s.d_trim_5, 1.0);
        add(t, "d_3", s.d_trim_3, 1.0);
    }
    let ns = if s.d.is_some() { vec![("vd", &s.np1), ("dj", &s.np2)] } else { vec![("vj", &s.np1)] };
    for (key, n) in ns {
        add(t, &format!("n_{key}"), n.len(), 1.0);
        let mut prev = "start".to_string();
        for c in n.chars() {
            add(t, &format!("n_after_{prev}"), c, 1.0);
            prev = c.to_string();
        }
    }
}

/// Chi-square of one observed table against its expected distribution. Categories with an
/// expected count under 5 are pooled into one bin (merged into the smallest bin if the
/// pool is itself under 5). Returns (n, χ², df, p); with under two bins, (n, TV, 0, NaN).
fn compare(obs: &[(String, f64)], exp: &[(String, f64)]) -> (f64, f64, usize, f64) {
    let n: f64 = obs.iter().map(|x| x.1).sum();
    let total: f64 = exp.iter().map(|x| x.1).sum();
    let get = |v: &[(String, f64)], k: &str| v.iter().find(|x| x.0 == k).map_or(0.0, |x| x.1);
    let cats: BTreeSet<&str> = obs.iter().chain(exp).map(|x| x.0.as_str()).collect();
    let (mut bins, mut pool) = (Vec::new(), (0.0, 0.0));
    let mut tv = 0.0;
    for k in cats {
        let (o, e) = (get(obs, k), n * get(exp, k) / total);
        assert!(o == 0.0 || e > 0.0, "category {k}: {o} draws of a zero-probability outcome");
        tv += (o - e).abs() / n / 2.0;
        if e >= 5.0 {
            bins.push((o, e));
        } else {
            pool = (pool.0 + o, pool.1 + e);
        }
    }
    if pool.1 > 0.0 {
        match bins.iter_mut().min_by(|a, b| a.1.total_cmp(&b.1)) {
            Some(small) if pool.1 < 5.0 => *small = (small.0 + pool.0, small.1 + pool.1),
            _ => bins.push(pool),
        }
    }
    if bins.len() < 2 {
        return (n, tv, 0, f64::NAN);
    }
    let chi: f64 = bins.iter().map(|(o, e)| (o - e).powi(2) / e).sum();
    (n, chi, bins.len() - 1, chi2_sf(chi, bins.len() - 1))
}

#[test]
fn chi2_sf_matches_tables() {
    // χ² critical values at p = 0.001 for df 1, 3, 10, 100
    for (x, k) in [(10.828, 1), (16.266, 3), (29.588, 10), (149.449, 100)] {
        let p = chi2_sf(x, k);
        assert!((p - 0.001).abs() < 2e-6, "df {k}: {p}");
    }
    assert!((chi2_sf(2.0, 2) - (-1.0f64).exp()).abs() < 1e-12);
}

#[test]
#[ignore = "release only: cargo test -p abstar_generate --release --test validation -- --include-ignored"]
fn distribution_recovery() {
    let b = bound("human-igh");
    let p = GenerateParams { n: 200_000, seed: 2026, ..preset("naive").unwrap() };
    let s = Simulator::new(&b, &p).unwrap();
    let obs = (0..p.n)
        .into_par_iter()
        .fold(Table::new, |mut t, i| {
            let r = sim(&s, p.seed, i);
            assert!(!r.d_inverted, "read {i}: inverted D under d_inversion = 0");
            observe(&mut t, &r);
            t
        })
        .reduce(Table::new, merge);
    let exp = expected(&b);
    let (obs, exp) = (by_table(&obs), by_table(&exp));
    let names: BTreeSet<&str> = obs.keys().chain(exp.keys()).copied().collect();
    let mut failed = vec![];
    for name in names {
        let (o, e) = (obs.get(name).cloned().unwrap_or_default(), exp.get(name).cloned().unwrap_or_default());
        if o.is_empty() {
            assert!(name.starts_with("n_after_"), "table {name}: no observations");
            continue;
        }
        let (n, stat, df, pv) = compare(&o, &e);
        let ok = if df == 0 { stat < 0.01 } else { pv > 0.001 };
        eprintln!("{name}: n {n} {} {stat:.2} df {df} p {pv:.4}", if df == 0 { "TV" } else { "chi2" });
        if !ok {
            failed.push(name);
        }
    }
    assert!(failed.is_empty(), "tables that do not match the bound model: {failed:?}");

    // the built-in model has d_inversion = 0 (checked above); test the rate where it is not
    let mut inv = b.clone();
    inv.model.d_inversion = 0.1;
    let s = Simulator::new(&inv, &p).unwrap();
    let k = (0..p.n).into_par_iter().filter(|&i| sim(&s, p.seed, i).d_inverted).count() as f64;
    let (rate, se) = (k / p.n as f64, (0.1f64 * 0.9 / p.n as f64).sqrt());
    eprintln!("d_inversion 0.1: rate {rate:.5}, {:.2} SE", (rate - 0.1) / se);
    assert!((rate - 0.1).abs() <= 3.0 * se, "inversion rate {rate}");
}

/// Read `i` under `productive_only`, drawn as the pipeline draws it (attempt `j` is stream
/// `i × 1000 + j`, at most 1000 attempts), and the number of attempts it took.
fn productive_read(s: &Simulator, seed: u64, i: usize) -> (Simulated, usize) {
    (0..1000)
        .map(|j| s.simulate(&mut ChaCha8Rng::seed_from_u64(seed_for(seed, i as u64 * 1000 + j))))
        .zip(1..)
        .find(|(r, _)| productivity(&r.read, &r.origin, &r.v, &r.j).productive)
        .unwrap_or_else(|| panic!("read {i}: no productive read in 1000 attempts"))
}

/// The bnab targets hold on productive output (`--preset bnab --productive-only`), the
/// preset's intended use: real bnAbs are productive. Denominator: productive reads.
#[test]
#[ignore = "release only: cargo test -p abstar_generate --release --test validation -- --include-ignored"]
fn bnab_targets() {
    let b = bound("human-igh");
    let p = GenerateParams { n: 50_000, seed: 2026, productive_only: true, ..preset("bnab").unwrap() };
    let s = Simulator::new(&b, &p).unwrap();
    // (V mutation frequency, V-derived base count, CDRH3 aa, has_indel, attempts)
    let reads: Vec<(f64, usize, usize, bool, usize)> = (0..p.n)
        .into_par_iter()
        .map(|i| {
            let (r, attempts) = productive_read(&s, p.seed, i);
            let junction = productivity(&r.read, &r.origin, &r.v, &r.j).junction.unwrap();
            let nv = r.origin.iter().filter(|o| matches!(o, abstar_generate::Origin::V(_))).count();
            // CDRH3 = junction_aa length − 2
            (r.v_mutation_frequency, nv, junction.len() / 3 - 2, !r.indels.is_empty(), attempts)
        })
        .collect();

    // Mutation hits only V bases that reach the final read (deleted and truncated bases are
    // never mutable), so realized divergence is round(target × nV) / nV: within 0.5/nV of
    // [0.15, 0.35] (measured: 0.40/nV below, 0.45/nV above). The bound allows a symmetric
    // 2/nV either side, so a future change that lets an indel move mutated V bases (e.g. a
    // deletion removing them) does not trip it spuriously.
    for (i, &(f, nv, ..)) in reads.iter().enumerate() {
        let (lo, hi) = (0.15 - 2.0 / nv as f64, 0.35 + 2.0 / nv as f64);
        assert!(f >= lo - 1e-12 && f < hi, "read {i}: V mutation {f} outside [{lo}, {hi}) (nV {nv})");
    }

    let n = p.n as f64;
    let frac = |min: usize| reads.iter().filter(|r| r.2 >= min).count() as f64 / n;
    let (f24, f35) = (frac(24), frac(35));
    let indel = reads.iter().filter(|r| r.3).count() as f64 / n;
    let attempts = reads.iter().map(|r| r.4).sum::<usize>() as f64 / n;
    eprintln!(
        "long_cdr3_weight {}: productive reads {}: CDRH3 >= 24 aa {f24:.4}, >= 35 aa {f35:.4}; \
         has_indel {indel:.4}; mean attempts {attempts:.2}",
        p.long_cdr3_weight, p.n
    );
    assert!(f24 >= 0.25 && f35 >= 0.02, "CDRH3 >= 24 aa {f24}, >= 35 aa {f35}");
    assert!((indel - 0.10).abs() <= 0.01, "has_indel {indel}");
}

/// `--preset bnab --productive-only` through the pipeline: no read hits the attempt limit,
/// and attempts per read stay near 1 / (in-frame junction rate) ≈ 3.1, since selection
/// keeps mutation from adding stops or breaking the anchor codons. Measured 3.98; the rest
/// is stops in N/D and inserted bases.
#[test]
#[ignore = "release only: cargo test -p abstar_generate --release --test validation -- --include-ignored"]
fn bnab_productive_only_attempts() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../python/abstar");
    let meta = GenerateMeta {
        abstar_version: "test".into(),
        germline_db: "human".into(),
        db_path: root.join("germline_dbs/human"),
        model_name: "human-igh".into(),
        model_path: root.join("models/human-igh/model.toml"),
        preset: Some("bnab".into()),
    };
    let p = GenerateParams { n: 5_000, seed: 3, productive_only: true, ..preset("bnab").unwrap() };
    let dir = tempfile::tempdir().unwrap();
    let sum = generate_to_files(&bound("human-igh"), &p, dir.path(), "bnab", &[], &meta).unwrap();
    let mean = sum.attempts as f64 / sum.n as f64;
    eprintln!("bnab productive_only: {} reads, mean attempts {mean:.2}; warnings {:?}", sum.n, sum.warnings);
    assert_eq!((sum.n, sum.productive), (5_000, 5_000));
    assert!(mean < 4.5, "mean attempts {mean}");
}

#[test]
#[ignore = "release only: cargo test -p abstar_generate --release --test validation -- --include-ignored"]
fn event_replay_builtin() {
    for name in ["human-igh", "human-igk", "human-igl", "human-tra", "human-trb"] {
        let b = bound(name);
        let base = GenerateParams {
            n: 20_000,
            seed: 7,
            mutation: "uniform:0-0.35".parse().unwrap(),
            long_cdr3_weight: 0.5,
            indel_rate: 0.5,
            truncate_5p: "uniform:0-40".parse().unwrap(),
            truncate_3p: "uniform:0-40".parse().unwrap(),
            constant_len: Some(40),
            shm_model: Some("motif".into()), // TCR models do not mutate by default
            ..preset("naive").unwrap()
        };
        // half the reads with frameshifting, duplicated insertions and stop-codon selection
        // (productive_only), half without
        let on =
            GenerateParams { frameshift_indels: true, indel_duplication: true, productive_only: true, ..base.clone() };
        let sims = [Simulator::new(&b, &base).unwrap(), Simulator::new(&b, &on).unwrap()];
        let indels = (0..base.n)
            .into_par_iter()
            .map(|i| {
                let r = sim(&sims[i % 2], base.seed, i);
                assert_eq!(replay(&r), r.read, "{name} read {i}");
                usize::from(!r.indels.is_empty())
            })
            .sum::<usize>();
        eprintln!("{name}: {} reads replayed, {indels} with an indel", base.n);
        assert!(indels > base.n / 5, "{name}: {indels} indels");
    }
}

#[test]
#[ignore = "release only: cargo test -p abstar_generate --release --test validation -- --include-ignored"]
fn oracle_agreement_builtin() {
    for name in ["human-igh", "human-igl"] {
        let b = bound(name);
        let base = GenerateParams {
            mutation: "uniform:0-0.3".parse().unwrap(),
            indel_rate: 0.3,
            truncate_5p: "uniform:0-60".parse().unwrap(),
            truncate_3p: "uniform:0-40".parse().unwrap(),
            ..preset("naive").unwrap()
        };
        let variants = [
            GenerateParams { constant_len: Some(30), frameshift_indels: true, ..base.clone() },
            base.clone(),
        ];
        let sims: Vec<Simulator> = variants.iter().map(|p| Simulator::new(&b, p).unwrap()).collect();
        let ambiguous = (0..2000usize)
            .into_par_iter()
            .map(|i| {
                let r = sim(&sims[i % 2], 11, i);
                let t = truth(&r, &b.candidates);
                assert_eq!(t, oracle::truth(&r, &b.candidates), "{name} read {i}");
                usize::from(t.v_valid.len() > 1 || t.unrecoverable_v || t.v_end_valid.len() > 1)
            })
            .sum::<usize>();
        eprintln!("{name}: 2000 reads agree; {ambiguous} with an ambiguous or unrecoverable V call or V end");
    }
}


/// Pins neutral-mode output (mutation on, `productive_only` off): any change to what the
/// generator draws or writes in that mode changes these hashes. Update them only for an
/// intended change, and say why in the commit.
#[test]
fn neutral_output_pinned() {
    use sha2::{Digest, Sha256};
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../python/abstar");
    let meta = GenerateMeta {
        abstar_version: "test".into(),
        germline_db: "human".into(),
        db_path: root.join("germline_dbs/human"),
        model_name: "human-igh".into(),
        model_path: root.join("models/human-igh/model.toml"),
        preset: Some("memory".into()),
    };
    let p = GenerateParams { n: 300, seed: 42, indel_rate: 0.2, ..preset("memory").unwrap() };
    let dir = tempfile::tempdir().unwrap();
    generate_to_files(&bound("human-igh"), &p, dir.path(), "pin", &[abstar_table::writers::Format::Airr], &meta).unwrap();
    let hash = |f: &str| {
        let bytes = std::fs::read(dir.path().join(f)).unwrap();
        Sha256::digest(&bytes).iter().map(|b| format!("{b:02x}")).collect::<String>()
    };
    assert_eq!(hash("pin.fasta"), "1e60db4508fc6450f63399da818ee460cb64ecadd31c696dcb5b315157585ef9");
    assert_eq!(hash("pin.truth.tsv"), "745e6f134bddb99aa69205be87a8e3ac91bed44b9c5018d0ef53ecdbb1def80a");
}
