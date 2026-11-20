use std::collections::BTreeMap;
use std::path::Path;

use abstar_germline::Locus;
use serde::Deserialize;

use crate::GenerateError;

const TOL: f64 = 1e-6;

/// Discrete distribution; negative values in trim tables are P-nucleotides.
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Dist {
    pub values: Vec<i64>,
    pub probs: Vec<f64>,
}

#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Source {
    pub citation: String,
    pub used_for: String,
}

/// One trim: a segment-wide default plus per-gene overrides.
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TrimSet {
    pub default: Dist,
    #[serde(default)]
    pub genes: BTreeMap<String, Dist>,
}

#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Trims {
    pub v_3: TrimSet,
    pub d_5: Option<TrimSet>,
    pub d_3: Option<TrimSet>,
    pub j_5: TrimSet,
}

#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Model {
    pub format_version: u32,
    pub name: String,
    pub species: String,
    pub locus: Locus,
    #[serde(default)]
    pub description: String,
    #[serde(default, rename = "source")]
    pub sources: Vec<Source>,
    pub v_usage: BTreeMap<String, f64>,
    #[serde(default)]
    pub d_given_j: BTreeMap<String, BTreeMap<String, f64>>,
    pub j_usage: BTreeMap<String, f64>,
    #[serde(default)]
    pub allele_usage: BTreeMap<String, f64>,
    pub trim: Trims,
    pub insertions: BTreeMap<String, Dist>,
    #[serde(default)]
    pub insertions_long_tail: BTreeMap<String, Dist>,
    pub n_markov: BTreeMap<String, [f64; 4]>,
    #[serde(default)]
    pub d_inversion: f64,
    #[serde(default = "none")]
    pub shm_model: String,
}

fn none() -> String {
    "none".into()
}

fn bad<T>(msg: String) -> Result<T, GenerateError> {
    Err(GenerateError::ModelInvalid(msg))
}

fn check_sum(key: &str, probs: impl Iterator<Item = f64>) -> Result<(), GenerateError> {
    let (mut sum, mut n) = (0.0, 0);
    for p in probs {
        if !(0.0..=1.0).contains(&p) {
            return bad(format!("{key}: probability {p} outside [0, 1]"));
        }
        sum += p;
        n += 1;
    }
    if n == 0 || (sum - 1.0).abs() > TOL {
        return bad(format!("{key}: probabilities sum to {sum}, expected 1"));
    }
    Ok(())
}

fn check_dist(key: &str, d: &Dist) -> Result<(), GenerateError> {
    if d.values.len() != d.probs.len() {
        return bad(format!(
            "{key}: values has {} entries but probs has {}",
            d.values.len(),
            d.probs.len()
        ));
    }
    check_sum(key, d.probs.iter().copied())
}

fn check_trim(key: &str, t: &TrimSet) -> Result<(), GenerateError> {
    check_dist(&format!("{key}.default"), &t.default)?;
    for (g, d) in &t.genes {
        check_dist(&format!("{key}.genes.{g}"), d)?;
    }
    Ok(())
}

fn validate(m: &Model) -> Result<(), GenerateError> {
    if m.format_version != 1 {
        return bad(format!("format_version: unsupported value {}, expected 1", m.format_version));
    }
    let has_d = matches!(m.locus, Locus::IGH | Locus::TRB | Locus::TRD);
    if !(0.0..=1.0).contains(&m.d_inversion) {
        return bad(format!("d_inversion: {} outside [0, 1]", m.d_inversion));
    }
    if !matches!(m.shm_model.as_str(), "HH_S5F" | "HKL_S5F" | "none") {
        return bad(format!("shm_model: unknown value '{}'", m.shm_model));
    }
    check_sum("v_usage", m.v_usage.values().copied())?;
    check_sum("j_usage", m.j_usage.values().copied())?;
    check_trim("trim.v_3", &m.trim.v_3)?;
    check_trim("trim.j_5", &m.trim.j_5)?;

    let keys: &[&str] = if has_d { &["vd", "dj"] } else { &["vj"] };
    for (table, name) in [(&m.insertions, "insertions"), (&m.insertions_long_tail, "insertions_long_tail")] {
        if let Some(k) = table.keys().find(|k| !keys.contains(&k.as_str())) {
            return bad(format!("{name}.{k}: not valid for locus {}", m.locus));
        }
        for (k, d) in table {
            check_dist(&format!("{name}.{k}"), d)?;
        }
    }
    if let Some(k) = keys.iter().find(|k| !m.insertions.contains_key(**k)) {
        return bad(format!("insertions.{k}: required for locus {}", m.locus));
    }

    if has_d {
        if m.d_given_j.is_empty() {
            return bad(format!("d_given_j: required for locus {}", m.locus));
        }
        for (j, row) in &m.d_given_j {
            check_sum(&format!("d_given_j.{j}"), row.values().copied())?;
        }
        for (t, name) in [(&m.trim.d_5, "trim.d_5"), (&m.trim.d_3, "trim.d_3")] {
            match t {
                Some(t) => check_trim(name, t)?,
                None => return bad(format!("{name}: required for locus {}", m.locus)),
            }
        }
    } else {
        if !m.d_given_j.is_empty() {
            return bad(format!("d_given_j: not allowed for locus {}", m.locus));
        }
        if m.trim.d_5.is_some() || m.trim.d_3.is_some() {
            return bad(format!("trim.d_5/trim.d_3: not allowed for locus {}", m.locus));
        }
        if m.d_inversion != 0.0 {
            return bad(format!("d_inversion: must be 0 for locus {}", m.locus));
        }
    }

    if let Some(k) = m.n_markov.keys().find(|k| !["start", "A", "C", "G", "T"].contains(&k.as_str())) {
        return bad(format!("n_markov.{k}: unknown key"));
    }
    check_allele_usage(m)?;
    for k in ["start", "A", "C", "G", "T"] {
        match m.n_markov.get(k) {
            Some(r) => check_sum(&format!("n_markov.{k}"), r.iter().copied())?,
            None => return bad(format!("n_markov.{k}: missing")),
        }
    }
    Ok(())
}

fn check_allele_usage(m: &Model) -> Result<(), GenerateError> {
    let mut sums: BTreeMap<&str, f64> = BTreeMap::new();
    for (a, &p) in &m.allele_usage {
        if !(0.0..=1.0).contains(&p) {
            return bad(format!("allele_usage.{a}: probability {p} outside [0, 1]"));
        }
        let gene = a.split('*').next().unwrap_or(a);
        let known = m.v_usage.contains_key(gene)
            || m.j_usage.contains_key(gene)
            || m.d_given_j.values().any(|r| r.contains_key(gene));
        if !known {
            return bad(format!("allele_usage.{a}: gene {gene} is not in v_usage, j_usage or d_given_j"));
        }
        *sums.entry(gene).or_default() += p;
    }
    for (g, s) in sums {
        if (s - 1.0).abs() > TOL {
            return bad(format!("allele_usage (gene {g}): probabilities sum to {s}, expected 1"));
        }
    }
    Ok(())
}

pub fn load_model(path: &Path) -> Result<Model, GenerateError> {
    let text = match std::fs::read_to_string(path) {
        Ok(t) => t,
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => {
            return Err(GenerateError::ModelNotFound(path.to_path_buf()))
        }
        Err(e) => return Err(e.into()),
    };
    let m: Model = toml::from_str(&text)
        .map_err(|e| GenerateError::ModelInvalid(format!("{}: {e}", path.display())))?;
    validate(&m)?;
    Ok(m)
}
