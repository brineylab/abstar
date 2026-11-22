use crate::{GermlineError, Segment};
use serde::Deserialize;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Recipe {
    pub name: String,
    #[serde(default)]
    pub description: String,
    pub keep: Option<PathBuf>,
    #[serde(default)]
    pub prefixes: BTreeMap<String, String>,
    #[serde(default, rename = "source")]
    pub sources: Vec<Source>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Source {
    pub species: Option<String>,
    pub database: Option<String>,
    #[serde(default)]
    pub description: String,
    #[serde(default)]
    pub ogrdb_json: Vec<PathBuf>,
    #[serde(default)]
    pub gap_v: bool,
    pub gap_templates: Option<Vec<String>>,
    #[serde(default)]
    v_fasta: Vec<PathBuf>,
    #[serde(default)]
    d_fasta: Vec<PathBuf>,
    #[serde(default)]
    j_fasta: Vec<PathBuf>,
    #[serde(default)]
    c_fasta: Vec<PathBuf>,
}

impl Source {
    /// FASTA files with the segment their key implies.
    pub fn fasta(&self) -> impl Iterator<Item = (Segment, &PathBuf)> {
        [
            (Segment::V, &self.v_fasta),
            (Segment::D, &self.d_fasta),
            (Segment::J, &self.j_fasta),
            (Segment::C, &self.c_fasta),
        ]
        .into_iter()
        .flat_map(|(s, v)| v.iter().map(move |p| (s, p)))
    }
}

pub(crate) fn read(path: &Path) -> Result<Recipe, GermlineError> {
    let bad = |m: String| GermlineError::RecipeInvalid(format!("{}: {m}", path.display()));
    let text = std::fs::read_to_string(path).map_err(|e| bad(e.to_string()))?;
    let r: Recipe = toml::from_str(&text).map_err(|e| bad(e.to_string()))?;
    // the name becomes a directory name: one plain path component, never a traversal
    let n = r.name.as_str();
    if n.is_empty() || n == "." || n == ".." || n.chars().any(|c| c == '/' || c == '\\' || c.is_control()) {
        return Err(bad(format!("name {n:?} must be a single plain path component")));
    }
    if r.sources.is_empty() {
        return Err(bad("no [[source]] entries".into()));
    }
    for (i, s) in r.sources.iter().enumerate() {
        let n = i + 1;
        let files = !s.ogrdb_json.is_empty() || s.fasta().next().is_some();
        if s.gap_v && s.database.is_some() {
            return Err(bad(format!("source {n}: gap_v does not apply to database sources")));
        }
        if s.gap_templates.is_some() && !s.gap_v {
            return Err(bad(format!("source {n}: gap_templates requires gap_v = true")));
        }
        match (&s.database, &s.species) {
            (Some(_), Some(_)) => return Err(bad(format!("source {n}: database sources take no species"))),
            (Some(_), None) if files => return Err(bad(format!("source {n}: database sources take no files"))),
            (None, None) => return Err(bad(format!("source {n}: species is required unless database is set"))),
            (None, Some(_)) if !files => return Err(bad(format!("source {n}: no files listed"))),
            _ => {}
        }
    }
    Ok(r)
}
