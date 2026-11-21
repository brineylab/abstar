use abstar_generate::{bind, load_model};
use std::collections::BTreeMap;
use std::path::Path;

fn kept(usage: &BTreeMap<String, f64>, bound: &[String]) -> f64 {
    usage.iter().filter(|(g, _)| bound.contains(g)).map(|(_, p)| p).sum()
}

#[test]
fn each_builtin_model_binds() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../python/abstar");
    let db = abstar_germline::load(&root.join("germline_dbs/human")).unwrap();
    for name in ["human-igh", "human-igk", "human-igl", "human-tra", "human-trb"] {
        let m = load_model(&root.join("models").join(name).join("model.toml")).unwrap();
        assert_eq!(m.name, name);
        assert_eq!(m.shm_model, if name.contains("-ig") { "motif" } else { "none" }, "{name}");
        assert!(!m.sources.is_empty(), "{name}: no [[source]] entries");
        let b = bind(m.clone(), &db).unwrap_or_else(|e| panic!("{name}: {e}"));
        eprintln!("{name}: warnings {:?}", b.warnings);
        assert!(b.warnings.is_empty(), "{name}: {:?}", b.warnings);

        let genes = |v: &[(abstar_germline::Allele, f64)]| v.iter().map(|(a, _)| a.gene.clone()).collect::<Vec<_>>();
        let cov_v = kept(&m.v_usage, &genes(&b.v));
        let cov_j = kept(&m.j_usage, &genes(&b.j));
        assert!(cov_v >= 0.90 && cov_j >= 0.90, "{name}: V {cov_v} J {cov_j}");
        for (j, row) in &m.d_given_j {
            let cov_d = kept(row, &genes(&b.d_given_j[j]));
            assert!(cov_d >= 0.90, "{name}: D|{j} {cov_d}");
        }
        if name == "human-igh" {
            assert!(m.insertions_long_tail.contains_key("vd") && m.insertions_long_tail.contains_key("dj"));
        }
    }
}

/// Fraction of reads without a junction (e.g. the J anchor trimmed away) over 5,000 unmutated reads.
fn no_junction_rate(name: &str) -> f64 {
    use arrow_array::Array;
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../python/abstar");
    let db = abstar_germline::load(&root.join("germline_dbs/human")).unwrap();
    let b = bind(load_model(&root.join("models").join(name).join("model.toml")).unwrap(), &db).unwrap();
    let mut p = abstar_generate::preset("naive").unwrap();
    p.n = 5000;
    let (mut n, mut none) = (0, 0);
    for rb in abstar_generate::generate_in_memory(&b, &p).unwrap() {
        let c = rb.column_by_name("junction_aa").unwrap();
        n += c.len();
        none += c.null_count();
    }
    none as f64 / n as f64
}

#[test]
fn light_chain_trims_keep_the_junction() {
    for name in ["human-igk", "human-igl"] {
        let r = no_junction_rate(name);
        eprintln!("{name}: no-junction rate {r:.4}");
        assert!(r < 0.08, "{name}: {r}");
    }
}
