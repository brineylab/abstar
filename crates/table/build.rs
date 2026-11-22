//! Compiles the codegen module into the build script to generate a fixture for `tests/codegen.rs`.
#[path = "src/codegen.rs"]
mod codegen;
use std::path::Path;

fn main() {
    let out = Path::new(&std::env::var("OUT_DIR").unwrap()).join("fixture_generated.rs");
    codegen::generate(Path::new("tests/fixture.toml"), &out, "Fixture", "FIXTURE_TOML");
}
