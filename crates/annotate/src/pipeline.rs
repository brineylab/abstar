//! Streaming run pipeline: read batches, annotate in parallel, write in input order.
use crate::error::RunError;
use crate::input::{Input, Source};
use crate::runmeta::{now_rfc3339, RunMeta};
use crate::schema::{arrow_schema, to_record_batch, Record};
use crate::stub::annotate_stub;
use crate::writers::{Format, ParquetSink, TsvWriter};
use abstar_germline::{Database, Locus};
use arrow_array::RecordBatch;
use rayon::prelude::*;
use std::cell::Cell;
use std::collections::BTreeMap;
use std::io::{BufWriter, ErrorKind, Write};
use std::panic::{self, AssertUnwindSafe};
use std::path::{Path, PathBuf};
use std::sync::Once;

#[derive(Clone, Debug)]
pub struct RunOptions {
    pub threads: Option<usize>,
    pub strict: bool,
    pub loci: Option<Vec<Locus>>,
    pub batch_size: usize,
}

impl Default for RunOptions {
    fn default() -> Self {
        Self { threads: None, strict: false, loci: None, batch_size: 4096 }
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct StatusCounts {
    pub annotated: u64,
    pub unassigned: u64,
    pub invalid_input: u64,
    pub error: u64,
}

impl StatusCounts {
    fn add(&mut self, o: &StatusCounts) {
        self.annotated += o.annotated;
        self.unassigned += o.unassigned;
        self.invalid_input += o.invalid_input;
        self.error += o.error;
    }
    fn to_json(&self) -> serde_json::Value {
        serde_json::json!({"annotated": self.annotated, "unassigned": self.unassigned,
            "invalid_input": self.invalid_input, "error": self.error})
    }
}

#[derive(Clone, Debug)]
pub struct InputSummary {
    pub source: String,
    pub outputs: BTreeMap<String, PathBuf>,
    pub counts: StatusCounts,
}

#[derive(Clone, Debug)]
pub struct RunSummary {
    pub inputs: Vec<InputSummary>,
    pub counts: StatusCounts,
    pub run_json: Option<PathBuf>,
}

type Annotator<'a> = &'a (dyn Fn(&str, &str, &Database) -> Record + Sync);
type Sink<'a> = &'a mut dyn FnMut(&RecordBatch) -> Result<(), RunError>;

thread_local! { static QUIET: Cell<bool> = const { Cell::new(false) }; }

/// Per-record panics are reported in the error row, so the default hook's stderr
/// message is suppressed only while a record is being annotated (thread-local flag);
/// every other panic still goes to the previously installed hook.
fn install_quiet_hook() {
    static ONCE: Once = Once::new();
    ONCE.call_once(|| {
        let prev = panic::take_hook();
        panic::set_hook(Box::new(move |info| {
            if !QUIET.with(Cell::get) {
                prev(info);
            }
        }));
    });
}

fn annotate_guarded(f: Annotator, id: &str, seq: &str, db: &Database) -> Record {
    QUIET.with(|q| q.set(true));
    let r = panic::catch_unwind(AssertUnwindSafe(|| f(id, seq, db)));
    QUIET.with(|q| q.set(false));
    r.unwrap_or_else(|p| {
        let msg = p
            .downcast_ref::<&str>()
            .map(|s| s.to_string())
            .or_else(|| p.downcast_ref::<String>().cloned())
            .unwrap_or_else(|| "unknown panic".into());
        Record {
            sequence_id: id.into(),
            sequence: seq.into(),
            annotation_status: "error".into(),
            failure_reason: Some(format!("internal_error: {msg}")),
            ..Default::default()
        }
    })
}

fn validate(db: &Database, opts: &RunOptions) -> Result<(), RunError> {
    if opts.batch_size == 0 {
        return Err(RunError::InvalidArgument("batch_size must be at least 1".into()));
    }
    if opts.threads == Some(0) {
        return Err(RunError::InvalidArgument("threads must be at least 1".into()));
    }
    for l in opts.loci.iter().flatten() {
        if !db.manifest.loci.contains(l) {
            let have: Vec<String> = db.manifest.loci.iter().map(|l| l.to_string()).collect();
            return Err(RunError::InvalidArgument(format!(
                "locus {l} is not in database `{}` (has: {})",
                db.manifest.name,
                have.join(", ")
            )));
        }
    }
    Ok(())
}

fn pool(opts: &RunOptions) -> Result<rayon::ThreadPool, RunError> {
    install_quiet_hook();
    rayon::ThreadPoolBuilder::new()
        .num_threads(opts.threads.unwrap_or(0))
        .build()
        .map_err(|e| RunError::Internal(e.to_string()))
}

/// Tabs and newlines cannot be represented in TSV; replace them so one record stays one row.
fn tsv_safe(s: &mut String) {
    if s.contains(['\t', '\n', '\r']) {
        *s = s.replace(['\t', '\n', '\r'], " ");
    }
}

/// Runs one opened input; `counts` is updated as batches complete so callers keep partial counts on error.
#[allow(clippy::too_many_arguments)]
fn run_input(
    input: &Input,
    mut src: Source,
    db: &Database,
    opts: &RunOptions,
    pool: &rayon::ThreadPool,
    ann: Annotator,
    counts: &mut StatusCounts,
    sink: Sink,
) -> Result<(), RunError> {
    let source_file = match input {
        Input::Records(_) => None,
        i => Some(i.source()),
    };
    loop {
        // ponytail: sequential batches; overlap read/annotate/write if benches show I/O stalls
        let batch = src.next_batch(opts.batch_size)?;
        if batch.is_empty() {
            return Ok(());
        }
        let mut recs: Vec<Record> = pool.install(|| {
            batch
                .par_iter()
                .map(|(id, seq, utf8)| {
                    if *utf8 {
                        annotate_guarded(ann, id, seq, db)
                    } else {
                        Record {
                            sequence_id: id.clone(),
                            sequence: seq.clone(),
                            annotation_status: "invalid_input".into(),
                            failure_reason: Some("invalid_encoding".into()),
                            ..Default::default()
                        }
                    }
                })
                .collect()
        });
        for r in &mut recs {
            r.source_file = source_file.clone();
            tsv_safe(&mut r.sequence_id);
            tsv_safe(&mut r.sequence);
            match r.annotation_status.as_str() {
                "annotated" => counts.annotated += 1,
                "unassigned" => counts.unassigned += 1,
                "invalid_input" => counts.invalid_input += 1,
                _ => {
                    counts.error += 1;
                    if opts.strict {
                        return Err(RunError::StrictRecordFailure(format!(
                            "{}: {}",
                            r.sequence_id,
                            r.failure_reason.as_deref().unwrap_or("error")
                        )));
                    }
                }
            }
        }
        sink(&to_record_batch(&recs))?;
    }
}

fn io_err(what: &str, e: impl std::fmt::Display) -> RunError {
    RunError::Io(format!("{what}: {e}"))
}

fn stem_conflicts(inputs: &[Input]) -> Result<Vec<String>, RunError> {
    let stems: Vec<String> = inputs.iter().map(Input::stem).collect();
    let mut by: BTreeMap<&str, Vec<String>> = BTreeMap::new();
    for (s, i) in stems.iter().zip(inputs) {
        by.entry(s).or_default().push(i.source());
    }
    let clash: Vec<String> = by
        .iter()
        .filter(|(_, v)| v.len() > 1)
        .map(|(s, v)| format!("`{s}` from {}", v.join(", ")))
        .collect();
    if clash.is_empty() {
        Ok(stems)
    } else {
        Err(RunError::OutputNameConflict(clash.join("; ")))
    }
}

fn fmt_name(f: Format) -> &'static str {
    match f {
        Format::Airr => "airr",
        Format::Parquet => "parquet",
    }
}

pub fn run_to_files(
    inputs: &[Input],
    output_dir: &Path,
    db: &Database,
    formats: &[Format],
    opts: &RunOptions,
    meta: &RunMeta,
) -> Result<RunSummary, RunError> {
    run_to_files_with(inputs, output_dir, db, formats, opts, meta, &annotate_stub)
}

fn run_to_files_with(
    inputs: &[Input],
    output_dir: &Path,
    db: &Database,
    formats: &[Format],
    opts: &RunOptions,
    meta: &RunMeta,
    ann: Annotator,
) -> Result<RunSummary, RunError> {
    validate(db, opts)?;
    let mut formats = formats.to_vec();
    formats.sort_by_key(|f| fmt_name(*f));
    formats.dedup();
    if formats.is_empty() {
        return Err(RunError::InvalidArgument("at least one output format is required".into()));
    }
    let stems = stem_conflicts(inputs)?;
    for i in inputs {
        if let Input::File(p) = i {
            if !p.exists() {
                return Err(RunError::InputNotFound(i.source()));
            }
        }
    }
    let pool = pool(opts)?;
    std::fs::create_dir_all(output_dir).map_err(|e| io_err(&output_dir.display().to_string(), e))?;
    // A stale run.json would mark a failed re-run as finished; it is rewritten last.
    let run_json = output_dir.join("run.json");
    match std::fs::remove_file(&run_json) {
        Err(e) if e.kind() != ErrorKind::NotFound => return Err(io_err(&run_json.display().to_string(), e)),
        _ => {}
    }
    let started = now_rfc3339();
    let mut total = StatusCounts::default();
    let mut summaries = vec![];
    for (input, stem) in inputs.iter().zip(&stems) {
        let src = Source::open(input)?; // before truncating this input's previous outputs
        let mut outputs = BTreeMap::new();
        let (mut tsv, mut pq) = (None, None);
        for f in &formats {
            let ext = if *f == Format::Airr { "tsv" } else { "parquet" };
            let path = output_dir.join(format!("{stem}.{ext}"));
            let what = path.display().to_string();
            match f {
                Format::Airr => {
                    let file = std::fs::File::create(&path).map_err(|e| io_err(&what, e))?;
                    tsv = Some((TsvWriter::new(BufWriter::new(file)).map_err(|e| io_err(&what, e))?, what));
                }
                Format::Parquet => pq = Some((ParquetSink::create(&path).map_err(|e| io_err(&what, e))?, what)),
            }
            outputs.insert(fmt_name(*f).to_string(), path);
        }
        let mut counts = StatusCounts::default();
        run_input(input, src, db, opts, &pool, ann, &mut counts, &mut |b| {
            if let Some((w, what)) = tsv.as_mut() {
                w.write(b).map_err(|e| io_err(what, e))?;
            }
            if let Some((w, what)) = pq.as_mut() {
                w.write(b).map_err(|e| io_err(what, e))?;
            }
            Ok(())
        })?;
        if let Some((w, what)) = tsv {
            w.finish().map_err(|e| io_err(&what, e))?;
        }
        if let Some((w, what)) = pq {
            w.finish().map_err(|e| io_err(&what, e))?;
        }
        total.add(&counts);
        summaries.push(InputSummary { source: input.source(), outputs, counts });
    }
    let doc = serde_json::json!({
        "abstar_version": meta.abstar_version,
        "airr_version": arrow_schema().metadata().get("airr_version"),
        "started": started,
        "finished": now_rfc3339(),
        "parameters": {
            "germline_db": meta.germline_db,
            "db_path": meta.db_path,
            "formats": formats.iter().map(|f| fmt_name(*f)).collect::<Vec<_>>(),
            "threads": opts.threads,
            "strict": opts.strict,
            "loci": opts.loci.as_ref().map(|l| l.iter().map(|l| l.to_string()).collect::<Vec<_>>()),
        },
        "database": {"name": db.manifest.name, "alleles_sha256": db.manifest.alleles_sha256},
        "inputs": summaries.iter().map(|s| serde_json::json!({
            "source": s.source, "outputs": s.outputs, "counts": s.counts.to_json(),
        })).collect::<Vec<_>>(),
        "counts": total.to_json(),
    });
    let text = serde_json::to_string_pretty(&doc).map_err(|e| RunError::Internal(e.to_string()))?;
    std::fs::write(&run_json, text + "\n").map_err(|e| io_err(&run_json.display().to_string(), e))?;
    Ok(RunSummary { inputs: summaries, counts: total, run_json: Some(run_json) })
}

/// TSV to `out`. A closed pipe (e.g. `| head`) ends the output and returns normally.
pub fn run_to_writer<W: Write>(
    inputs: &[Input],
    db: &Database,
    opts: &RunOptions,
    out: W,
) -> Result<StatusCounts, RunError> {
    validate(db, opts)?;
    let pool = pool(opts)?;
    let pipe = |e: &std::io::Error| e.kind() == ErrorKind::BrokenPipe;
    let mut w = match TsvWriter::new(BufWriter::with_capacity(1 << 16, out)) {
        Ok(w) => w,
        Err(e) if pipe(&e) => return Ok(StatusCounts::default()),
        Err(e) => return Err(io_err("output", e)),
    };
    let mut counts = StatusCounts::default();
    let mut closed = false;
    for input in inputs {
        let mut c = StatusCounts::default();
        let r = run_input(input, Source::open(input)?, db, opts, &pool, &annotate_stub, &mut c, &mut |b| match w.write(b) {
            Err(e) if pipe(&e) => {
                closed = true;
                Err(RunError::Io("broken pipe".into()))
            }
            r => r.map_err(|e| io_err("output", e)),
        });
        counts.add(&c);
        match r {
            Err(_) if closed => return Ok(counts),
            r => r?,
        }
    }
    match w.finish() {
        Err(e) if !pipe(&e) => Err(io_err("output", e)),
        _ => Ok(counts),
    }
}

pub fn run_in_memory(inputs: &[Input], db: &Database, opts: &RunOptions) -> Result<Vec<RecordBatch>, RunError> {
    run_in_memory_with(inputs, db, opts, &annotate_stub)
}

/// Same as `run_in_memory` with an injected annotator (tests).
#[doc(hidden)]
pub fn run_in_memory_with(
    inputs: &[Input],
    db: &Database,
    opts: &RunOptions,
    ann: &(dyn Fn(&str, &str, &Database) -> Record + Sync),
) -> Result<Vec<RecordBatch>, RunError> {
    validate(db, opts)?;
    let pool = pool(opts)?;
    let mut out = vec![];
    for input in inputs {
        run_input(input, Source::open(input)?, db, opts, &pool, ann, &mut StatusCounts::default(), &mut |b| {
            out.push(b.clone());
            Ok(())
        })?;
    }
    if out.is_empty() {
        out.push(to_record_batch(&[]));
    }
    Ok(out)
}
