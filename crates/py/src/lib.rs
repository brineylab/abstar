//! `abstar._core`: a thin translation layer over the Rust crates. Arguments are
//! converted with the GIL held, the Rust work runs detached, and every Rust error
//! becomes `CoreError(code, message)`.
use abstar_annotate::input::Input;
use abstar_annotate::pipeline::{RunOptions, RunSummary, StatusCounts};
use abstar_annotate::runmeta::RunMeta;
use abstar_annotate::schema::{arrow_schema, FIELDS_TOML};
use abstar_annotate::writers::Format;
use abstar_generate as gen;
use abstar_germline::{BuildOptions, BuildWarning, Database, Locus};
use pyo3::conversion::FromPyObjectOwned;
use pyo3::create_exception;
use pyo3::exceptions::PyException;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};
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
    let parse = |s: &String| {
        s.parse::<Locus>().map_err(|_| {
            err("INVALID_ARGUMENT", format!("unknown locus {s:?}; valid loci: {}", Locus::NAMES.join(", ")))
        })
    };
    let loci = loci.map(|l| l.iter().map(parse).collect::<PyResult<Vec<_>>>()).transpose()?;
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

/// Gap V sequences against the V alleles of the template databases.
#[pyfunction]
fn germline_gap<'py>(
    py: Python<'py>,
    queries: Vec<(String, String)>,
    template_paths: Vec<String>,
) -> PyResult<Bound<'py, PyList>> {
    let results = py.detach(|| -> PyResult<_> {
        let mut templates = Vec::new();
        for p in &template_paths {
            templates.extend(load(p)?.alleles);
        }
        Ok(abstar_germline::gap::gap_v(&queries, &templates))
    })?;
    let out = PyList::empty(py);
    for r in results {
        let d = PyDict::new(py);
        d.set_item("name", r.name)?;
        d.set_item("gapped", r.gapped)?;
        d.set_item("template", r.template)?;
        d.set_item("template_species", r.template_species)?;
        d.set_item("identity", r.identity)?;
        d.set_item("cdr1_length", r.cdr1_length)?;
        d.set_item("cdr2_length", r.cdr2_length)?;
        d.set_item("warnings", r.warnings)?;
        d.set_item("failure", r.failure)?;
        out.append(d)?;
    }
    Ok(out)
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

fn gen_err(e: gen::GenerateError) -> PyErr {
    err(e.code(), e.message())
}

fn opt<'py, T: FromPyObjectOwned<'py>>(key: &str, v: &Bound<'py, PyAny>) -> PyResult<T> {
    v.extract::<T>().map_err(|e| err("INVALID_ARGUMENT", format!("{key}: {}", Into::<PyErr>::into(e))))
}

/// A number (fixed) or a spec string (`fixed:X`, `uniform:LO-HI`, `beta:A,B`).
fn spec(key: &str, v: &Bound<'_, PyAny>) -> PyResult<gen::Spec> {
    match v.extract::<String>() {
        Ok(s) => s.parse().map_err(gen_err),
        Err(_) => Ok(gen::Spec::Fixed(opt(key, v)?)),
    }
}

/// Generation parameters: the preset (default `naive`), then every other key of `d`
/// overrides it; a `None` value leaves the preset's. `locus` is checked against the model.
fn gen_params(d: &Bound<'_, PyDict>) -> PyResult<(gen::GenerateParams, Option<String>, Option<Locus>)> {
    let preset: Option<String> = match d.get_item("preset")? {
        Some(v) if !v.is_none() => Some(opt("preset", &v)?),
        _ => None,
    };
    let mut p = gen::preset(preset.as_deref().unwrap_or("naive")).map_err(gen_err)?;
    let mut locus = None;
    for (k, v) in d.iter() {
        let key: String = opt("option name", &k)?;
        if v.is_none() {
            continue;
        }
        let k = key.as_str();
        match k {
            "preset" => {}
            "locus" => {
                let s: String = opt(k, &v)?;
                locus = Some(s.parse::<Locus>().map_err(|_| {
                    err("INVALID_ARGUMENT", format!("unknown locus {s:?}; valid loci: {}", Locus::NAMES.join(", ")))
                })?);
            }
            "n" => p.n = opt(k, &v)?,
            "seed" => p.seed = opt(k, &v)?,
            "threads" => p.threads = Some(opt(k, &v)?),
            "mutation" => p.mutation = spec(k, &v)?,
            "long_cdr3_weight" => p.long_cdr3_weight = opt(k, &v)?,
            "indels" => p.indel_rate = opt(k, &v)?,
            "frameshift_indels" => p.frameshift_indels = opt(k, &v)?,
            "indel_duplication" => p.indel_duplication = opt(k, &v)?,
            "truncate_5p" => p.truncate_5p = spec(k, &v)?,
            "truncate_3p" => p.truncate_3p = spec(k, &v)?,
            "constant" => p.constant_len = Some(opt(k, &v)?),
            "isotype_weights" => p.isotype_weights = opt(k, &v)?,
            "productive_only" => p.productive_only = opt(k, &v)?,
            "shm_table" => p.shm_table = Some(opt(k, &v)?),
            "shm_model" => p.shm_model = Some(opt(k, &v)?),
            _ => {
                return Err(err(
                    "INVALID_ARGUMENT",
                    format!(
                        "unknown option {k:?}; valid: mutation, long_cdr3_weight, indels, frameshift_indels, \
                         indel_duplication, truncate_5p, truncate_3p, constant, isotype_weights, productive_only, \
                         shm_table, shm_model"
                    ),
                ))
            }
        }
    }
    Ok((p, preset, locus))
}

/// Load the model, check it is for `locus`, and bind it to the database.
fn gen_bind(db_path: &str, model_path: &str, locus: Option<Locus>) -> PyResult<(gen::Bound, Database)> {
    let model = gen::load_model(model_path.as_ref()).map_err(gen_err)?;
    if let Some(l) = locus.filter(|l| *l != model.locus) {
        return Err(err("INVALID_ARGUMENT", format!("model '{}' is for {}, not {l}", model.name, model.locus)));
    }
    let db = load(db_path)?;
    Ok((gen::bind(model, &db).map_err(gen_err)?, db))
}

fn gen_summary<'py>(py: Python<'py>, s: &gen::GenerateSummary) -> PyResult<Bound<'py, PyDict>> {
    let path = |p: &Option<PathBuf>| p.as_ref().map(|p| p.to_string_lossy().into_owned());
    let d = PyDict::new(py);
    d.set_item("n", s.n)?;
    d.set_item("productive", s.productive)?;
    d.set_item("attempts", s.attempts)?;
    d.set_item("seed", s.seed)?;
    d.set_item("fasta", path(&s.fasta))?;
    let truth: BTreeMap<&str, String> = s.truth.iter().map(|(k, v)| (k.as_str(), v.to_string_lossy().into_owned())).collect();
    d.set_item("truth", truth)?;
    d.set_item("generate_json", path(&s.generate_json))?;
    d.set_item("warnings", &s.warnings)?;
    Ok(d)
}

/// The truth table field definitions (truth.toml).
#[pyfunction]
fn truth_toml() -> &'static str {
    gen::TRUTH_TOML
}

/// name, species, locus, description and sources of a recombination model file.
#[pyfunction]
fn model_summary<'py>(py: Python<'py>, path: PathBuf) -> PyResult<Bound<'py, PyDict>> {
    let m = gen::load_model(&path).map_err(gen_err)?;
    let d = PyDict::new(py);
    d.set_item("name", &m.name)?;
    d.set_item("species", &m.species)?;
    d.set_item("locus", m.locus.to_string())?;
    d.set_item("description", &m.description)?;
    let sources: Vec<BTreeMap<&str, &str>> = m
        .sources
        .iter()
        .map(|s| BTreeMap::from([("citation", s.citation.as_str()), ("used_for", s.used_for.as_str())]))
        .collect();
    d.set_item("sources", sources)?;
    Ok(d)
}

#[pyfunction]
fn generate_to_files<'py>(
    py: Python<'py>,
    db_path: String,
    model_path: String,
    params: Bound<'py, PyDict>,
    output_dir: PathBuf,
    name: String,
    formats: Vec<String>,
) -> PyResult<Bound<'py, PyDict>> {
    let (p, preset, locus) = gen_params(&params)?;
    let formats =
        formats.iter().map(|f| f.parse::<Format>()).collect::<Result<Vec<_>, _>>().map_err(|e| err("INVALID_ARGUMENT", e))?;
    let s = py.detach(|| {
        let (bound, db) = gen_bind(&db_path, &model_path, locus)?;
        let meta = gen::GenerateMeta {
            abstar_version: VERSION.into(),
            germline_db: db.manifest.name.clone(),
            db_path: db_path.clone().into(),
            model_name: bound.model.name.clone(),
            model_path: model_path.clone().into(),
            preset,
        };
        gen::generate_to_files(&bound, &p, &output_dir, &name, &formats, &meta).map_err(gen_err)
    })?;
    gen_summary(py, &s)
}

/// FASTA to the process's stdout (fd 1), bypassing `sys.stdout`. A closed pipe returns normally.
#[pyfunction]
fn generate_fasta_to_stdout<'py>(
    py: Python<'py>,
    db_path: String,
    model_path: String,
    params: Bound<'py, PyDict>,
) -> PyResult<Bound<'py, PyDict>> {
    let (p, _, locus) = gen_params(&params)?;
    let s = py.detach(|| {
        let (bound, _) = gen_bind(&db_path, &model_path, locus)?;
        gen::generate_fasta_to_writer(&bound, &p, std::io::stdout().lock()).map_err(gen_err)
    })?;
    gen_summary(py, &s)
}

/// Truth rows as an Arrow table (exposes `__arrow_c_stream__`).
#[pyfunction]
fn generate_in_memory<'py>(
    py: Python<'py>,
    db_path: String,
    model_path: String,
    params: Bound<'py, PyDict>,
) -> PyResult<PyTable> {
    let (p, _, locus) = gen_params(&params)?;
    let batches = py.detach(|| {
        let (bound, _) = gen_bind(&db_path, &model_path, locus)?;
        gen::generate_in_memory(&bound, &p).map_err(gen_err)
    })?;
    PyTable::try_new(batches, gen::arrow_schema())
}

#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(truth_toml, m)?)?;
    m.add_function(wrap_pyfunction!(model_summary, m)?)?;
    m.add_function(wrap_pyfunction!(generate_to_files, m)?)?;
    m.add_function(wrap_pyfunction!(generate_fasta_to_stdout, m)?)?;
    m.add_function(wrap_pyfunction!(generate_in_memory, m)?)?;
    m.add("CoreError", m.py().get_type::<CoreError>())?;
    m.add_function(wrap_pyfunction!(version, m)?)?;
    m.add_function(wrap_pyfunction!(fields_toml, m)?)?;
    m.add_function(wrap_pyfunction!(germline_build, m)?)?;
    m.add_function(wrap_pyfunction!(germline_summary, m)?)?;
    m.add_function(wrap_pyfunction!(germline_gap, m)?)?;
    m.add("GAP_DEFAULT_TEMPLATES", abstar_germline::gap::DEFAULT_TEMPLATES.to_vec())?;
    m.add_function(wrap_pyfunction!(run_to_files, m)?)?;
    m.add_function(wrap_pyfunction!(run_to_stdout, m)?)?;
    m.add_function(wrap_pyfunction!(run_in_memory, m)?)
}
