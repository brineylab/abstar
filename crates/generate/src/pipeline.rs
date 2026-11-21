//! Deterministic, parallel generation: reads, truth rows (decision 0026) and `generate.json`.
//!
//! Read `i` (0-based) is simulated with a `ChaCha8Rng` seeded by `seed_for(seed, i)`, or
//! `seed_for(seed, i * 1000 + j)` for attempt `j` under `productive_only`, so output depends
//! only on the seed, never on threads or batch size.

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fs::File;
use std::io::{self, BufWriter, ErrorKind, Write};
use std::path::Path;

use abstar_germline::{Allele, Locus, Manifest};
use abstar_table::now_rfc3339;
use abstar_table::writers::{Format, ParquetSink, TsvWriter};
use arrow_array::RecordBatch;
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;
use rayon::prelude::*;
use sha2::{Digest, Sha256};

use crate::schema::{arrow_schema, to_record_batch, TruthRecord};
use crate::{
    productivity, truth, Bound, GenerateError, GenerateMeta, GenerateParams, GenerateSummary, IndelKind, Origin,
    Productivity, Simulated, Simulator,
};

/// Attempts per read under `productive_only` before the run fails.
const MAX_ATTEMPTS: u64 = 1000;

fn splitmix64(mut z: u64) -> u64 {
    z = z.wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// RNG seed of stream `i` of a run with seed `global`.
pub fn seed_for(global: u64, i: u64) -> u64 {
    splitmix64(global ^ splitmix64(i))
}

fn invalid(msg: String) -> GenerateError {
    GenerateError::InvalidArgument(msg)
}

fn io_at(path: &Path, e: impl std::fmt::Display, kind: ErrorKind) -> GenerateError {
    GenerateError::Io(io::Error::new(kind, format!("{}: {e}", path.display())))
}

fn io_err(path: &Path, e: io::Error) -> GenerateError {
    io_at(path, &e, e.kind())
}

/// Standard genetic code; anything but ACGT translates to `X`.
pub(crate) fn translate(nt: &[u8]) -> String {
    const AA: &[u8; 64] = b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";
    let code = |b: u8| match b.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    };
    nt.as_chunks::<3>()
        .0
        .iter()
        .map(|c| match (code(c[0]), code(c[1]), code(c[2])) {
            (Some(a), Some(b), Some(d)) => AA[a * 16 + b * 4 + d] as char,
            _ => 'X',
        })
        .collect()
}

/// Ungapped index of IMGT-gapped position `g` of `a`.
fn ungapped(a: &Allele, g: usize) -> usize {
    a.sequence_gapped.as_bytes()[..g].iter().filter(|&&c| c != b'.').count()
}

/// Run totals; updated per batch, so a broken pipe still reports what was written.
#[derive(Default)]
struct Counts {
    n: u64,
    productive: u64,
    attempts: u64,
    /// Per-read warning (text before the first `:`) → reads that raised it.
    warnings: BTreeMap<String, u64>,
}

struct Row {
    rec: TruthRecord,
    productive: bool,
    attempts: u64,
    warnings: Vec<String>,
}

type Sink<'a> = &'a mut dyn FnMut(&[TruthRecord]) -> Result<(), GenerateError>;

/// A validated run: everything built once and shared by every read.
struct Run<'a> {
    bound: &'a Bound,
    params: &'a GenerateParams,
    sim: Simulator<'a>,
    pool: rayon::ThreadPool,
    /// Candidate allele name → gene.
    genes: HashMap<&'a str, &'a str>,
    width: usize,
}

impl<'a> Run<'a> {
    fn new(bound: &'a Bound, params: &'a GenerateParams) -> Result<Self, GenerateError> {
        let sim = Simulator::new(bound, params)?;
        if params.batch_size == 0 {
            return Err(invalid("batch_size must be at least 1".into()));
        }
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(params.threads.unwrap_or(0))
            .build()
            .map_err(|e| GenerateError::Internal(e.to_string()))?;
        let c = &bound.candidates;
        let genes = c.v.iter().chain(&c.d).chain(&c.j).chain(&c.c).map(|a| (a.name.as_str(), a.gene.as_str())).collect();
        Ok(Self { bound, params, sim, pool, genes, width: params.n.to_string().len() })
    }

    /// Simulates every read in index order, handing each batch to `sink`. `full = false`
    /// fills only `sequence_id` and `sequence` (FASTA output needs no truth sets).
    fn rows(&self, full: bool, counts: &mut Counts, sink: Sink) -> Result<(), GenerateError> {
        let bs = self.params.batch_size;
        for start in (0..self.params.n).step_by(bs) {
            let end = (start + bs).min(self.params.n);
            let rows: Vec<Row> =
                self.pool.install(|| (start..end).into_par_iter().map(|i| self.one(i, full)).collect::<Result<_, _>>())?;
            let mut recs = Vec::with_capacity(rows.len());
            for r in rows {
                counts.n += 1;
                counts.productive += u64::from(r.productive);
                counts.attempts += r.attempts;
                for w in r.warnings {
                    let key = w.split(':').next().unwrap_or(&w).to_string();
                    *counts.warnings.entry(key).or_default() += 1;
                }
                recs.push(r.rec);
            }
            sink(&recs)?;
        }
        Ok(())
    }

    fn one(&self, i: usize, full: bool) -> Result<Row, GenerateError> {
        let p = self.params;
        let id = format!("g{:0w$}", i + 1, w = self.width);
        let draw = |k: u64| {
            let s = self.sim.simulate(&mut ChaCha8Rng::seed_from_u64(seed_for(p.seed, k)));
            let pr = productivity(&s.read, &s.origin, &s.v, &s.j);
            (s, pr)
        };
        let (s, pr, attempts) = if p.productive_only {
            let base = i as u64 * MAX_ATTEMPTS;
            (0..MAX_ATTEMPTS)
                .map(|j| (draw(base + j), j + 1))
                .find(|((_, pr), _)| pr.productive)
                .map(|((s, pr), k)| (s, pr, k))
                .ok_or_else(|| {
                    invalid(format!(
                        "productive rate too low: no productive read for {id} in {MAX_ATTEMPTS} attempts; \
                         check the model's trims and insertions, or turn off productive_only"
                    ))
                })?
        } else {
            let (s, pr) = draw(i as u64);
            (s, pr, 1)
        };
        let rec = if full {
            self.record(id, &s, &pr)
        } else {
            TruthRecord { sequence_id: id, sequence: String::from_utf8_lossy(&s.read).into_owned(), ..Default::default() }
        };
        Ok(Row { rec, productive: pr.productive, attempts, warnings: s.warnings })
    }

    /// The truth row of one read. All read positions are 1-based in the final read.
    fn record(&self, id: String, s: &Simulated, pr: &Productivity) -> TruthRecord {
        let t = truth(s, &self.bound.candidates);
        let (read, o) = (&s.read, &s.origin);
        let text = |b: &[u8]| String::from_utf8_lossy(b).into_owned();
        // 0-based first and last read index of each segment
        let span = |f: fn(&Origin) -> bool| Some((o.iter().position(f)?, o.iter().rposition(f)?));
        let v = span(|x| matches!(x, Origin::V(_)));
        let d = span(|x| matches!(x, Origin::D(_)));
        let j = span(|x| matches!(x, Origin::J(_)));
        let c = span(|x| matches!(x, Origin::C(_)));
        let start = |x: Option<(usize, usize)>| x.map(|x| x.0 as i64 + 1);
        let end = |x: Option<(usize, usize)>| x.map(|x| x.1 as i64 + 1);
        // germline index (1-based) of read index `k`
        let germ = |k: usize| {
            Some(match o[k] {
                Origin::V(g) => ungapped(&s.v, g),
                Origin::D(u) | Origin::J(u) | Origin::C(u) => u,
                _ => return None,
            } as i64
                + 1)
        };
        let gstart = |x: Option<(usize, usize)>| x.and_then(|x| germ(x.0));
        let gend = |x: Option<(usize, usize)>| x.and_then(|x| germ(x.1));
        // N/P bases between adjacent segments, as they are in the read
        let between = |l: Option<(usize, usize)>, r: Option<(usize, usize)>| match (l, r) {
            (Some(l), Some(r)) => Some(text(&read[l.1 + 1..r.0])),
            _ => None,
        };
        let np1 = between(v, d.or(j));
        let np2 = d.and(between(d, j));

        let junction = pr.junction.clone().map(|r| &read[r]);
        let junction_aa = junction.map(translate);
        let cdr3_min = if self.bound.model.locus == Locus::IGH { 24 } else { 16 };

        let (nv, same) = read.iter().zip(o).fold((0usize, 0usize), |(n, m), (&b, x)| match *x {
            Origin::V(g) => (n + 1, m + usize::from(b.eq_ignore_ascii_case(&s.v.sequence_gapped.as_bytes()[g]))),
            _ => (n, m),
        });

        // net substitution per read position
        let t5 = s.truncate_5;
        let mut net: BTreeMap<usize, (u8, u8)> = BTreeMap::new();
        for m in &s.mutations {
            net.entry(m.pos - t5).and_modify(|e| e.1 = m.to).or_insert((m.from, m.to));
        }
        let mutations: Vec<String> = net
            .into_iter()
            .filter(|(_, (f, to))| f != to)
            .map(|(p, (f, to))| format!("{}:{}>{}", p + 1, f as char, to as char))
            .collect();

        // indels with a read base on both sides, equivalence range clipped to the same
        let read_end = t5 + read.len(); // untruncated index just past the read
        let indels: Vec<String> = s
            .indels
            .iter()
            .filter_map(|x| {
                let (kind, inner) = match x.kind {
                    IndelKind::Ins => ("ins", x.len),
                    IndelKind::Del => ("del", 0),
                };
                let (lo_ok, hi_ok) = (t5 + 1, read_end.checked_sub(1 + inner)?);
                if x.read_pos < lo_ok || x.read_pos > hi_ok {
                    return None;
                }
                let pos = |u: usize| u - t5 + 1;
                Some(format!(
                    "{kind}:{}:{}:{}:{}-{}",
                    pos(x.read_pos),
                    x.len,
                    x.bases,
                    pos(x.equiv.0.max(lo_ok)),
                    pos(x.equiv.1.min(hi_ok))
                ))
            })
            .collect();

        let genes = |set: &[String]| -> Vec<String> {
            set.iter().map(|a| self.genes[a.as_str()].to_string()).collect::<BTreeSet<_>>().into_iter().collect()
        };
        let name = |a: &Allele| a.name.clone();
        let gene = |a: &Allele| a.gene.clone();
        let vmf = s.v_mutation_frequency;

        let mut rec = TruthRecord {
            sequence_id: id,
            sequence: text(read),
            rev_comp: Some(false),
            productive: Some(pr.productive),
            vj_in_frame: Some(pr.vj_in_frame),
            stop_codon: Some(pr.stop_codon),
            locus: Some(self.bound.model.locus.to_string()),
            v_call: Some(name(&s.v)),
            d_call: s.d.as_ref().map(name),
            j_call: Some(name(&s.j)),
            c_call: s.c.as_ref().map(name),
            junction: junction.map(text),
            junction_length: junction.map(|x| x.len() as i64),
            long_cdr3: junction_aa.as_ref().map(|aa| aa.len() as i64 - 2 >= cdr3_min),
            junction_aa,
            v_sequence_start: start(v),
            v_sequence_end: end(v),
            d_sequence_start: start(d),
            d_sequence_end: end(d),
            j_sequence_start: start(j),
            j_sequence_end: end(j),
            c_sequence_start: start(c),
            c_sequence_end: end(c),
            v_germline_start: gstart(v),
            v_germline_end: gend(v),
            d_germline_start: gstart(d),
            d_germline_end: gend(d),
            j_germline_start: gstart(j),
            j_germline_end: gend(j),
            np1_length: np1.as_ref().map(|x| x.len() as i64),
            np1,
            np2_length: np2.as_ref().map(|x| x.len() as i64),
            np2,
            v_identity: (nv > 0).then(|| same as f64 / nv as f64),
            v_gene: Some(gene(&s.v)),
            d_gene: s.d.as_ref().map(gene),
            j_gene: Some(gene(&s.j)),
            c_gene: s.c.as_ref().map(gene),
            v_gene_valid: Some(genes(&t.v_valid)),
            d_gene_valid: Some(genes(&t.d_valid)),
            j_gene_valid: Some(genes(&t.j_valid)),
            ambiguous_v: Some(t.v_valid.len() > 1),
            ambiguous_d: Some(t.d_valid.len() + usize::from(t.d_null_valid) > 1),
            ambiguous_j: Some(t.j_valid.len() > 1),
            v_call_valid: Some(t.v_valid),
            d_call_valid: Some(t.d_valid),
            j_call_valid: Some(t.j_valid),
            c_call_valid: Some(t.c_valid),
            d_call_null_valid: Some(t.d_null_valid),
            v_sequence_end_valid: Some(t.v_end_valid),
            d_sequence_start_valid: Some(t.d_start_valid),
            d_sequence_end_valid: Some(t.d_end_valid),
            j_sequence_start_valid: Some(t.j_start_valid),
            v_trim_3: Some(s.v_trim_3),
            d_trim_5: s.d.as_ref().map(|_| s.d_trim_5),
            d_trim_3: s.d.as_ref().map(|_| s.d_trim_3),
            j_trim_5: Some(s.j_trim_5),
            d_inverted: s.d.as_ref().map(|_| s.d_inverted),
            mutations: Some(mutations),
            has_indel: Some(!indels.is_empty()),
            indels: Some(indels),
            truncate_5p: Some(s.truncate_5 as i64),
            truncate_3p: Some(s.truncate_3 as i64),
            truncated_5p: Some(s.truncate_5 > 0),
            truncated_3p: Some(s.truncate_3 > 0),
            v_mutation_frequency: Some(vmf),
            high_mutation: Some(vmf >= 0.20),
            unrecoverable_v: Some(t.unrecoverable_v),
            unrecoverable_d: Some(t.unrecoverable_d),
            unrecoverable_j: Some(t.unrecoverable_j),
            ..Default::default()
        };
        crate::derive::fill(&mut rec, s, pr.junction.is_some());
        rec
    }
}

/// Bind warnings, then each per-read warning once with the number of reads that raised it.
fn warnings(bound: &Bound, counts: &Counts) -> Vec<String> {
    let per_read = counts.warnings.iter().map(|(w, &k)| format!("{w} ({k} sequence{})", if k == 1 { "" } else { "s" }));
    bound.warnings.iter().cloned().chain(per_read).collect()
}

fn summary(params: &GenerateParams, bound: &Bound, counts: &Counts) -> GenerateSummary {
    GenerateSummary {
        n: counts.n,
        productive: counts.productive,
        attempts: counts.attempts,
        seed: params.seed,
        warnings: warnings(bound, counts),
        ..Default::default()
    }
}

fn fmt_name(f: Format) -> &'static str {
    match f {
        Format::Airr => "airr",
        Format::Parquet => "parquet",
    }
}

fn sha256_file(path: &Path) -> Result<String, GenerateError> {
    let bytes = std::fs::read(path).map_err(|e| io_err(path, e))?;
    Ok(Sha256::digest(&bytes).iter().map(|b| format!("{b:02x}")).collect())
}

fn manifest(db_path: &Path) -> Result<Manifest, GenerateError> {
    let path = db_path.join("manifest.toml");
    let text = std::fs::read_to_string(&path).map_err(|e| io_err(&path, e))?;
    toml::from_str(&text).map_err(|e| io_at(&path, e, ErrorKind::InvalidData))
}

/// Writes `<name>.fasta`, `<name>.truth.tsv` and/or `<name>.truth.parquet`, then
/// `generate.json`, into `out_dir`. An existing `generate.json` is deleted before any
/// output is written, so a failed rerun cannot look finished.
pub fn generate_to_files(
    bound: &Bound,
    params: &GenerateParams,
    out_dir: &Path,
    name: &str,
    formats: &[Format],
    meta: &GenerateMeta,
) -> Result<GenerateSummary, GenerateError> {
    if name.is_empty() || name == "." || name == ".." || name.contains(['/', '\\']) {
        return Err(invalid(format!("output name '{name}' must be a plain file name")));
    }
    if out_dir.exists() && !out_dir.is_dir() {
        return Err(invalid(format!("output directory {} exists and is not a directory", out_dir.display())));
    }
    let mut formats = formats.to_vec();
    formats.sort_by_key(|f| fmt_name(*f));
    formats.dedup();
    let run = Run::new(bound, params)?;
    let model_sha256 = sha256_file(&meta.model_path)?;
    let db = manifest(&meta.db_path)?;

    std::fs::create_dir_all(out_dir).map_err(|e| io_err(out_dir, e))?;
    let json_path = out_dir.join("generate.json");
    match std::fs::remove_file(&json_path) {
        Err(e) if e.kind() != ErrorKind::NotFound => return Err(io_err(&json_path, e)),
        _ => {}
    }
    let started = now_rfc3339();

    let fasta_path = out_dir.join(format!("{name}.fasta"));
    let mut fasta = BufWriter::new(File::create(&fasta_path).map_err(|e| io_err(&fasta_path, e))?);
    let mut truth_paths = BTreeMap::new();
    let (mut tsv, mut pq) = (None, None);
    for f in &formats {
        let path = out_dir.join(format!("{name}.truth.{}", if *f == Format::Airr { "tsv" } else { "parquet" }));
        match f {
            Format::Airr => {
                let file = File::create(&path).map_err(|e| io_err(&path, e))?;
                let w = TsvWriter::new(BufWriter::new(file), &arrow_schema()).map_err(|e| io_err(&path, e))?;
                tsv = Some((w, path.clone()));
            }
            Format::Parquet => {
                let w = ParquetSink::create(&path, arrow_schema()).map_err(|e| io_at(&path, e, ErrorKind::Other))?;
                pq = Some((w, path.clone()));
            }
        }
        truth_paths.insert(fmt_name(*f).to_string(), path);
    }

    let mut counts = Counts::default();
    run.rows(!formats.is_empty(), &mut counts, &mut |recs| {
        write_fasta(&mut fasta, recs).map_err(|e| io_err(&fasta_path, e))?;
        if tsv.is_some() || pq.is_some() {
            let b = to_record_batch(recs);
            if let Some((w, path)) = tsv.as_mut() {
                w.write(&b).map_err(|e| io_err(path, e))?;
            }
            if let Some((w, path)) = pq.as_mut() {
                w.write(&b).map_err(|e| io_at(path, e, ErrorKind::Other))?;
            }
        }
        Ok(())
    })?;
    fasta.flush().map_err(|e| io_err(&fasta_path, e))?;
    if let Some((w, path)) = tsv {
        w.finish().map_err(|e| io_err(&path, e))?;
    }
    if let Some((w, path)) = pq {
        w.finish().map_err(|e| io_at(&path, e, ErrorKind::Other))?;
    }

    let mut sum = summary(params, bound, &counts);
    // destructured so a new parameter cannot be left out of generate.json
    let GenerateParams {
        n,
        seed,
        mutation,
        long_cdr3_weight,
        indel_rate,
        frameshift_indels,
        indel_duplication,
        truncate_5p,
        truncate_3p,
        constant_len,
        isotype_weights,
        productive_only,
        shm_model,
        shm_table,
        threads,
        batch_size,
    } = params;
    let doc = serde_json::json!({
        "abstar_version": meta.abstar_version,
        "started": started,
        "finished": now_rfc3339(),
        "seed": seed,
        "parameters": {
            "preset": meta.preset,
            "n": n,
            "seed": seed,
            "mutation": mutation.to_string(),
            "long_cdr3_weight": long_cdr3_weight,
            "indel_rate": indel_rate,
            "frameshift_indels": frameshift_indels,
            "indel_duplication": indel_duplication,
            "truncate_5p": truncate_5p.to_string(),
            "truncate_3p": truncate_3p.to_string(),
            "constant_len": constant_len,
            "isotype_weights": isotype_weights,
            "productive_only": productive_only,
            "shm_model": shm_model,
            "shm_table": shm_table,
            "threads": threads,
            "batch_size": batch_size,
            "formats": formats.iter().map(|f| fmt_name(*f)).collect::<Vec<_>>(),
            "germline_db": meta.germline_db,
            "db_path": meta.db_path,
            "model_path": meta.model_path,
        },
        "model": {"name": meta.model_name, "sha256": model_sha256},
        "database": {"name": db.name, "alleles_sha256": db.alleles_sha256},
        "counts": {"n": sum.n, "productive": sum.productive, "attempts": sum.attempts},
        "warnings": sum.warnings,
    });
    let text = serde_json::to_string_pretty(&doc).map_err(|e| GenerateError::Internal(e.to_string()))?;
    std::fs::write(&json_path, text + "\n").map_err(|e| io_err(&json_path, e))?;

    sum.fasta = Some(fasta_path);
    sum.truth = truth_paths;
    sum.generate_json = Some(json_path);
    Ok(sum)
}

fn write_fasta(w: &mut impl Write, recs: &[TruthRecord]) -> io::Result<()> {
    for r in recs {
        writeln!(w, ">{}\n{}", r.sequence_id, r.sequence)?;
    }
    Ok(())
}

/// FASTA only, to `out`. A closed pipe (e.g. `| head`) ends the output and returns normally,
/// with counts for the reads generated so far.
pub fn generate_fasta_to_writer<W: Write>(
    bound: &Bound,
    params: &GenerateParams,
    out: W,
) -> Result<GenerateSummary, GenerateError> {
    let run = Run::new(bound, params)?;
    let mut w = BufWriter::with_capacity(1 << 16, out);
    let mut counts = Counts::default();
    let mut closed = false;
    let pipe = |e: &io::Error| e.kind() == ErrorKind::BrokenPipe;
    let r = run.rows(false, &mut counts, &mut |recs| match write_fasta(&mut w, recs) {
        Err(e) if pipe(&e) => {
            closed = true;
            Err(GenerateError::Io(e))
        }
        r => r.map_err(GenerateError::Io),
    });
    match r {
        Err(_) if closed => {}
        r => r?,
    }
    match w.flush() {
        Err(e) if !pipe(&e) => Err(GenerateError::Io(e)),
        _ => Ok(summary(params, bound, &counts)),
    }
}

/// Truth rows as Arrow batches of at most `batch_size` rows; always at least one batch.
pub fn generate_in_memory(bound: &Bound, params: &GenerateParams) -> Result<Vec<RecordBatch>, GenerateError> {
    let run = Run::new(bound, params)?;
    let mut out = Vec::new();
    run.rows(true, &mut Counts::default(), &mut |recs| {
        out.push(to_record_batch(recs));
        Ok(())
    })?;
    if out.is_empty() {
        out.push(to_record_batch(&[]));
    }
    Ok(out)
}
