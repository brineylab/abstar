//! `run.json` support.
use std::path::PathBuf;
pub use abstar_table::now_rfc3339;

pub struct RunMeta {
    pub abstar_version: String,
    pub germline_db: String,
    pub db_path: PathBuf,
}
