use abstar_generate::{load_model, GenerateError};
use abstar_germline::Locus;
use std::path::{Path, PathBuf};

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data").join(name)
}

#[test]
fn tiny_model_loads() {
    let m = load_model(&data("tiny_igh.model.toml")).unwrap();
    assert_eq!(m.locus, Locus::IGH);
    assert_eq!(m.v_usage.len(), 3);
    assert_eq!(m.j_usage.len(), 2);
    assert_eq!(m.d_given_j["IGHJ1"].len(), 2);
    assert_eq!(m.trim.v_3.genes.len(), 1);
    assert_eq!(m.sources.len(), 1);
    assert!(m.trim.v_3.default.values.contains(&-2));
    assert_eq!(m.n_markov.len(), 5);
    assert_eq!(load_model(&data("tiny_igk.model.toml")).unwrap().locus, Locus::IGK);
}

#[test]
fn missing_file_is_model_not_found() {
    let e = load_model(&data("nope.model.toml")).unwrap_err();
    assert!(matches!(e, GenerateError::ModelNotFound(_)));
    assert_eq!(e.code(), "MODEL_NOT_FOUND");
}

#[test]
fn probabilities_must_sum_to_one() {
    let e = load_model(&data("bad_sum.model.toml")).unwrap_err();
    assert_eq!(e.code(), "MODEL_INVALID");
    assert!(e.to_string().contains("v_usage"), "{e}");
}

#[test]
fn d_tables_required_for_igh_forbidden_for_igk() {
    let e = load_model(&data("bad_igh_no_d.model.toml")).unwrap_err();
    assert_eq!(e.code(), "MODEL_INVALID");
    assert!(e.to_string().contains("d_given_j"), "{e}");
    let e = load_model(&data("bad_igk_with_d.model.toml")).unwrap_err();
    assert_eq!(e.code(), "MODEL_INVALID");
    assert!(e.to_string().contains("d_given_j"), "{e}");
}
