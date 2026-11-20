//! Generates `TruthRecord`, `arrow_schema()` and `to_record_batch()` from `truth.toml`.
use std::path::Path;

fn main() {
    let out = Path::new(&std::env::var("OUT_DIR").unwrap()).join("truth_generated.rs");
    abstar_table::codegen::generate(Path::new("truth.toml"), &out, "TruthRecord", "TRUTH_TOML");
}
