//! What a `generate.json` records beyond the parameters, and the run summary.
use std::collections::BTreeMap;
use std::path::PathBuf;

/// Provenance supplied by the caller (CLI or Python).
#[derive(Debug, Clone)]
pub struct GenerateMeta {
    pub abstar_version: String,
    pub germline_db: String,
    /// Database directory; `generate.json` takes its name and `alleles_sha256` from the manifest here.
    pub db_path: PathBuf,
    pub model_name: String,
    /// Model file; `generate.json` records its sha256.
    pub model_path: PathBuf,
    /// Preset the parameters started from, if any.
    pub preset: Option<String>,
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct GenerateSummary {
    pub n: u64,
    pub productive: u64,
    /// Reads simulated, including those rejected by `productive_only`.
    pub attempts: u64,
    pub fasta: Option<PathBuf>,
    /// Format name (`airr`, `parquet`) → truth file.
    pub truth: BTreeMap<String, PathBuf>,
    pub generate_json: Option<PathBuf>,
    pub seed: u64,
    pub warnings: Vec<String>,
}
