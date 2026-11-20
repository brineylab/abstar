mod common;
use abstar_generate::{bind, load_model, Model};
use abstar_germline::{Locus, Segment};
use common::*;
use std::path::Path;

fn tiny() -> Model {
    load_model(&Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/tiny_igh.model.toml")).unwrap()
}

fn p(v: &[(abstar_germline::Allele, f64)], name: &str) -> f64 {
    v.iter().find(|(a, _)| a.name == name).unwrap().1
}

#[test]
fn binds_tiny_model() {
    let b = bind(tiny(), &tiny_igh_db()).unwrap();
    assert!(b.warnings.is_empty());
    assert!((p(&b.v, "IGHV1-1*01") - 0.3).abs() < 1e-9);
    assert!((p(&b.v, "IGHV1-1*02") - 0.2).abs() < 1e-9);
    assert!((p(&b.v, "IGHV1-2*01") - 0.3).abs() < 1e-9);
    assert_eq!(b.d_given_j["IGHJ1"].len(), 2);
    assert_eq!(b.candidates.c.len(), 1);
    assert_eq!(b.candidates.v.len(), 4);
}

#[test]
fn drops_missing_gene_and_renormalizes() {
    let mut m = tiny();
    m.allele_usage.clear();
    m.v_usage = [("IGHV1-1".into(), 0.5), ("IGHV1-2".into(), 0.45), ("IGHVZ".into(), 0.05)].into();
    let b = bind(m, &tiny_igh_db()).unwrap();
    assert_eq!(b.warnings.len(), 1);
    assert!(b.warnings[0].contains("IGHVZ"));
    // IGHV1-1 has two alleles splitting 0.5/0.95 equally.
    assert!((p(&b.v, "IGHV1-2*01") - 0.4737).abs() < 1e-3);
    assert!((b.v.iter().map(|x| x.1).sum::<f64>() - 1.0).abs() < 1e-9);
}

#[test]
fn low_coverage_fails() {
    let mut m = tiny();
    m.allele_usage.clear();
    m.v_usage = [("IGHV1-1".into(), 0.5), ("IGHV1-2".into(), 0.3), ("IGHVZ".into(), 0.2)].into();
    let e = bind(m, &tiny_igh_db()).unwrap_err();
    assert_eq!(e.code(), "MODEL_INVALID");
    assert!(e.to_string().contains("IGHVZ"), "{e}");
}

#[test]
fn anchorless_allele_not_generated_but_is_candidate() {
    let mut a = tiny_igh_alleles();
    a.iter_mut().find(|x| x.name == "IGHV1-1*02").unwrap().anchor = None;
    let b = bind(tiny(), &db(a)).unwrap();
    assert!(b.v.iter().all(|(x, _)| x.name != "IGHV1-1*02"));
    assert!((p(&b.v, "IGHV1-1*01") - 0.5).abs() < 1e-9);
    assert!(b.candidates.v.iter().any(|x| x.name == "IGHV1-1*02"));
}

#[test]
fn gene_with_only_anchorless_alleles_is_dropped() {
    let mut a = tiny_igh_alleles();
    a.iter_mut().find(|x| x.name == "IGHV2-1*01").unwrap().anchor = None;
    let mut m = tiny();
    m.v_usage = [("IGHV1-1".into(), 0.5), ("IGHV1-2".into(), 0.45), ("IGHV2-1".into(), 0.05)].into();
    let b = bind(m, &db(a)).unwrap();
    assert!(b.warnings[0].contains("IGHV2-1"));
}

#[test]
fn species_filter_in_mixed_db() {
    let mut a = tiny_igh_alleles();
    let mut mm = allele("IGHV1-1*01", Locus::IGH, Segment::V, "mouse", true);
    mm.id = "mm|IGHV1-1*01".into();
    a.push(mm);
    let b = bind(tiny(), &db(a)).unwrap();
    assert!(b.v.iter().all(|(x, _)| x.species == "human"));
    assert_eq!(b.candidates.v.len(), 4);
}

#[test]
fn ungapped_v_candidates_rejected() {
    let mut a = tiny_igh_alleles();
    a[2].sequence_gapped = a[2].sequence();
    let e = bind(tiny(), &db(a)).unwrap_err();
    assert_eq!(e.code(), "MODEL_INVALID");
    assert!(e.to_string().contains("IGHV1-2*01"), "{e}");
}

#[test]
fn d_rows_renormalized_per_j() {
    let mut m = tiny();
    m.d_given_j.get_mut("IGHJ1").unwrap().insert("IGHD9".into(), 0.0);
    // drop IGHD2-1 from the db: row mass 0.7 of 1.0 for J1 (0.2/1 for J2) -> coverage < 0.9
    let mut a = tiny_igh_alleles();
    a.retain(|x| x.name != "IGHD2-1*01");
    assert!(bind(tiny(), &db(a)).is_err());
    let b = bind(m, &tiny_igh_db()).unwrap();
    assert_eq!(b.d_given_j["IGHJ1"].len(), 2);
}
