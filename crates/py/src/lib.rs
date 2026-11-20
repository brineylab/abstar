use pyo3::prelude::*;

/// The crate version.
#[pyfunction]
fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(version, m)?)
}

#[cfg(test)]
mod tests {
    #[test]
    fn version_is_dev() {
        assert_eq!(super::version(), "1.0.0-dev.0");
    }
}
