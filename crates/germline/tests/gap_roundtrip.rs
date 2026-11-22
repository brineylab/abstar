//! Round trips: strip shipped IMGT-gapped V alleles, regap them, and compare with IMGT's gaps.
//! Run with `cargo test -p abstar_germline --release --test gap_roundtrip -- --nocapture`.
use abstar_germline::gap::gap_v;
use abstar_germline::{load, Allele, Locus, Segment};
use rayon::prelude::*;
use std::collections::BTreeMap;
use std::path::PathBuf;

fn alleles(db: &str) -> Vec<Allele> {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../python/abstar/germline_dbs")
        .join(db);
    load(&dir).unwrap().alleles
}

fn anchored_v(db: &[Allele]) -> Vec<&Allele> {
    db.iter()
        .filter(|a| a.segment == Segment::V && a.anchor.is_some())
        .collect()
}

#[derive(Default, Clone, Copy)]
struct Tally {
    n: usize,
    exact: usize,
    cys104: usize,
    failed: usize,
}

/// Regap each query with `templates(query)`; tally per locus.
fn measure(
    label: &str,
    queries: &[&Allele],
    templates: impl Fn(&Allele) -> Vec<Allele> + Sync,
) -> BTreeMap<Locus, Tally> {
    let rows: Vec<(Locus, bool, bool, bool)> = queries
        .par_iter()
        .map(|a| {
            let r = gap_v(&[(a.id.clone(), a.sequence())], &templates(a))
                .pop()
                .unwrap();
            let g = r.gapped.as_deref().unwrap_or("");
            // Cys104 placed where IMGT places it: the same input codon at gapped nt 309..312.
            let cys = g.get(309..312).is_some_and(|c| !c.contains('.'))
                && g[..309].bytes().filter(|&b| b != b'.').count() == a.anchor.unwrap();
            let exact = g == a.sequence_gapped;
            if !exact {
                eprintln!("  {label} wrong: {} {:?}", a.id, r.failure.or(r.template));
            }
            (a.locus, exact, cys, r.gapped.is_none())
        })
        .collect();
    let mut t: BTreeMap<Locus, Tally> = BTreeMap::new();
    for (locus, exact, cys, failed) in rows {
        let e = t.entry(locus).or_default();
        e.n += 1;
        e.exact += exact as usize;
        e.cys104 += cys as usize;
        e.failed += failed as usize;
    }
    eprintln!("{label}\n  locus      n  exact  failed  cys104 (of gapped)");
    for (l, x) in &t {
        eprintln!(
            "  {l:<5} {:>6} {:>6} {:>7} {:>7}",
            x.n, x.exact, x.failed, x.cys104
        );
    }
    t
}

fn total(t: &BTreeMap<Locus, Tally>) -> Tally {
    t.values().fold(Tally::default(), |a, x| Tally {
        n: a.n + x.n,
        exact: a.exact + x.exact,
        cys104: a.cys104 + x.cys104,
        failed: a.failed + x.failed,
    })
}

fn pct(k: usize, n: usize) -> f64 {
    100.0 * k as f64 / n as f64
}

#[test]
fn self_round_trip_leave_one_out() {
    for db in ["human", "c57bl6", "balbc"] {
        let all = alleles(db);
        let t = measure(db, &anchored_v(&all), |q| {
            all.iter().filter(|a| a.id != q.id).cloned().collect()
        });
        let x = total(&t);
        assert!(pct(x.exact, x.n) >= 99.0, "{db}: {}/{} exact", x.exact, x.n);
    }
}

/// Exact-match floors per locus, in whole percent of all alleles (a failure counts against it).
/// Cys104 must be where IMGT puts it in at least 99% of the alleles written; a refusal places no
/// Cys104, right or wrong, so failures are left out of that rate.
fn assert_floors(label: &str, t: &BTreeMap<Locus, Tally>, floors: &[(Locus, f64)]) {
    for &(locus, floor) in floors {
        let x = t[&locus];
        let written = x.n - x.failed;
        assert!(
            pct(x.exact, x.n) >= floor,
            "{label} {locus}: {}/{} exact",
            x.exact,
            x.n
        );
        assert!(
            pct(x.cys104, written) >= 99.0,
            "{label} {locus}: {}/{written} Cys104",
            x.cys104
        );
    }
}

#[test]
fn mouse_from_human() {
    let human = alleles("human");
    let floors = [(Locus::IGH, 100.0), (Locus::IGK, 98.0), (Locus::IGL, 100.0)];
    for db in ["c57bl6", "balbc"] {
        let mouse = alleles(db);
        let t = measure(&format!("{db} from human"), &anchored_v(&mouse), |_| {
            human.clone()
        });
        assert_floors(db, &t, &floors);
    }
}

#[test]
fn human_ig_from_mouse() {
    let human = alleles("human");
    let mouse: Vec<Allele> = [alleles("c57bl6"), alleles("balbc")].concat();
    let queries: Vec<&Allele> = anchored_v(&human)
        .into_iter()
        .filter(|a| matches!(a.locus, Locus::IGH | Locus::IGK | Locus::IGL))
        .collect();
    let t = measure("human IG from mouse", &queries, |_| mouse.clone());
    let floors = [(Locus::IGH, 100.0), (Locus::IGK, 100.0), (Locus::IGL, 78.0)];
    assert_floors("human", &t, &floors);
}
