//! Generates `Record`, `arrow_schema()` and `to_record_batch()` from `fields.toml`.
use std::path::Path;

fn main() {
    let out = Path::new(&std::env::var("OUT_DIR").unwrap()).join("schema_generated.rs");
    abstar_table::codegen::generate(Path::new("fields.toml"), &out, "Record", "FIELDS_TOML");
}
