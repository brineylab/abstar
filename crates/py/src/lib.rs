//! `abstar._core`: a thin translation layer over the Rust crates. Arguments are
//! converted with the GIL held, the Rust work runs detached, and every Rust error
//! becomes `CoreError(code, message)`.
use abstar_annotate::input::Input;
use abstar_annotate::pipeline::{RunOptions, RunSummary, StatusCounts};
use abstar_annotate::runmeta::RunMeta;
use abstar_annotate::schema::{arrow_schema, FIELDS_TOML};
use abstar_annotate::writers::Format;
use abstar_germline::{BuildOptions, BuildWarning, Database, Locus};
use pyo3::create_exception;
use pyo3::exceptions::PyException;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use pyo3_arrow::PyTable;
use std::collections::BTreeMap;
use std::fmt::Display;
use std::path::PathBuf;

create_exception!(_core, CoreError, PyException, "A Rust error: args == (code, message).");

const VERSION: &str = env!("CARGO_PKG_VERSION");

fn err(code: &str, e: impl Display) -> PyErr {
    CoreError::new_err((code.to_string(), e.to_string()))
}

fn germline_err(e: abstar_germline::GermlineError) -> PyErr {
    err(e.code(), e)
}

fn run_err(e: abstar_annotate::error::RunError) -> PyErr {
    err(e.code(), e)
}

fn load(path: &str) -> PyResult<Database> {
    abstar_germline::load(&PathBuf::from(path)).map_err(germline_err)
}

fn inputs(paths: Vec<String>) -> Vec<Input> {
    paths.into_iter().map(|p| if p == "-" { Input::Stdin } else { Input::File(p.into()) }).collect()
}

fn options(threads: Option<usize>, strict: bool, loci: Option<Vec<String>>) -> PyResult<RunOptions> {
    let loci = loci
        .map(|l| l.iter().map(|s| s.parse::<Locus>()).collect::<Result<Vec<_>, _>>())
        .transpose()
        .map_err(|e| err("INVALID_ARGUMENT", e))?;
    Ok(RunOptions { threads, strict, loci, ..Default::default() })
}

fn counts_dict<'py>(py: Python<'py>, c: &StatusCounts) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    d.set_item("annotated", c.annotated)?;
    d.set_item("unassigned", c.unassigned)?;
    d.set_item("invalid_input", c.invalid_input)?;
    d.set_item("error", c.error)?;
    Ok(d)
}

fn warnings(w: &[BuildWarning]) -> Vec<BTreeMap<&str, &str>> {
    w.iter().map(|w| BTreeMap::from([("allele", w.allele.as_str()), ("message", w.message.as_str())])).collect()
}

/// The crate version.
#[pyfunction]
fn version() -> &'static str {
    VERSION
}

/// The output field definitions (fields.toml).
#[pyfunction]
fn fields_toml() -> &'static str {
    FIELDS_TOML
}

#[pyfunction]
fn germline_build<'py>(
    py: Python<'py>,
    recipe: PathBuf,
    output_dir: PathBuf,
    overwrite: bool,
    database_paths: BTreeMap<String, PathBuf>,
) -> PyResult<Bound<'py, PyDict>> {
    let opts = BuildOptions { overwrite, abstar_version: VERSION.into() };
    let r = py
        .detach(|| abstar_germline::build(&recipe, &output_dir, &database_paths, &opts))
        .map_err(germline_err)?;
    let d = PyDict::new(py);
    d.set_item("name", &r.name)?;
    d.set_item("path", r.path.to_string_lossy())?;
    d.set_item("n_alleles", r.n_alleles)?;
    d.set_item("counts", &r.counts)?;
    d.set_item("restricted", &r.restricted)?;
    d.set_item("passed_through", &r.passed_through)?;
    d.set_item("warnings", warnings(&r.warnings))?;
    Ok(d)
}

/// Manifest fields plus "counts" and "n_alleles".
#[pyfunction]
fn germline_summary<'py>(py: Python<'py>, path: String) -> PyResult<Bound<'py, PyAny>> {
    let (json, counts, n) = py.detach(|| -> PyResult<_> {
        let db = load(&path)?;
        let json = serde_json::to_string(&db.manifest).map_err(|e| err("INTERNAL", e))?;
        Ok((json, db.counts(), db.alleles.len()))
    })?;
    let d = py.import("json")?.call_method1("loads", (json,))?;
    d.set_item("counts", counts)?;
    d.set_item("n_alleles", n)?;
    Ok(d)
}

fn summary_dict<'py>(py: Python<'py>, s: &RunSummary) -> PyResult<Bound<'py, PyDict>> {
    let ins = s
        .inputs
        .iter()
        .map(|i| {
            let d = PyDict::new(py);
            d.set_item("source", &i.source)?;
            let outs: BTreeMap<&str, String> =
                i.outputs.iter().map(|(k, v)| (k.as_str(), v.to_string_lossy().into_owned())).collect();
            d.set_item("outputs", outs)?;
            d.set_item("counts", counts_dict(py, &i.counts)?)?;
            Ok(d)
        })
        .collect::<PyResult<Vec<_>>>()?;
    let d = PyDict::new(py);
    d.set_item("inputs", ins)?;
    d.set_item("counts", counts_dict(py, &s.counts)?)?;
    d.set_item("run_json", s.run_json.as_ref().map(|p| p.to_string_lossy().into_owned()))?;
    Ok(d)
}

#[pyfunction]
#[allow(clippy::too_many_arguments)]
fn run_to_files<'py>(
    py: Python<'py>,
    inputs: Vec<String>,
    records: Vec<(String, String)>,
    output_dir: PathBuf,
    db_path: String,
    formats: Vec<String>,
    threads: Option<usize>,
    strict: bool,
    loci: Option<Vec<String>>,
) -> PyResult<Bound<'py, PyDict>> {
    let mut ins = self::inputs(inputs);
    if !records.is_empty() {
        ins.push(Input::Records(records));
    }
    let opts = options(threads, strict, loci)?;
    let formats =
        formats.iter().map(|f| f.parse::<Format>()).collect::<Result<Vec<_>, _>>().map_err(|e| err("INVALID_ARGUMENT", e))?;
    let s = py.detach(|| {
        let db = load(&db_path)?;
        let meta = RunMeta { abstar_version: VERSION.into(), germline_db: db.manifest.name.clone(), db_path: db_path.clone().into() };
        abstar_annotate::pipeline::run_to_files(&ins, &output_dir, &db, &formats, &opts, &meta).map_err(run_err)
    })?;
    summary_dict(py, &s)
}

/// TSV to the process's stdout (fd 1), bypassing `sys.stdout`. A closed pipe returns normally.
#[pyfunction]
#[allow(clippy::too_many_arguments)]
fn run_to_stdout<'py>(
    py: Python<'py>,
    inputs: Vec<String>,
    records: Vec<(String, String)>,
    db_path: String,
    threads: Option<usize>,
    strict: bool,
    loci: Option<Vec<String>>,
) -> PyResult<Bound<'py, PyDict>> {
    let mut ins = self::inputs(inputs);
    if !records.is_empty() {
        ins.push(Input::Records(records));
    }
    let opts = options(threads, strict, loci)?;
    let c = py.detach(|| {
        let db = load(&db_path)?;
        abstar_annotate::pipeline::run_to_writer(&ins, &db, &opts, std::io::stdout().lock()).map_err(run_err)
    })?;
    counts_dict(py, &c)
}

/// Annotations as an Arrow table (exposes `__arrow_c_stream__`).
#[pyfunction]
fn run_in_memory(
    py: Python<'_>,
    inputs: Vec<String>,
    records: Vec<(String, String)>,
    db_path: String,
    threads: Option<usize>,
    strict: bool,
    loci: Option<Vec<String>>,
) -> PyResult<PyTable> {
    let mut ins = self::inputs(inputs);
    if !records.is_empty() {
        ins.push(Input::Records(records));
    }
    let opts = options(threads, strict, loci)?;
    let batches = py.detach(|| {
        let db = load(&db_path)?;
        abstar_annotate::pipeline::run_in_memory(&ins, &db, &opts).map_err(run_err)
    })?;
    PyTable::try_new(batches, arrow_schema())
}

#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("CoreError", m.py().get_type::<CoreError>())?;
    m.add_function(wrap_pyfunction!(version, m)?)?;
    m.add_function(wrap_pyfunction!(fields_toml, m)?)?;
    m.add_function(wrap_pyfunction!(germline_build, m)?)?;
    m.add_function(wrap_pyfunction!(germline_summary, m)?)?;
    m.add_function(wrap_pyfunction!(run_to_files, m)?)?;
    m.add_function(wrap_pyfunction!(run_to_stdout, m)?)?;
    m.add_function(wrap_pyfunction!(run_in_memory, m)?)
}

#[cfg(test)]
mod tests {
    #[test]
    fn version_is_dev() {
        assert_eq!(super::version(), "1.0.0-dev.0");
    }
}
