fn main() {
    // macOS: let the cdylib link without libpython (-undefined dynamic_lookup),
    // so `cargo test --workspace` builds it outside maturin.
    pyo3_build_config::add_extension_module_link_args();
}
