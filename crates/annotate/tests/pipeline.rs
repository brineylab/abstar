use abstar_annotate::error::RunError;
use abstar_annotate::input::Input;
use abstar_annotate::pipeline::*;
use abstar_annotate::runmeta::RunMeta;
use abstar_annotate::schema::{arrow_schema, Record};
use abstar_annotate::stub::annotate_stub;
use abstar_annotate::writers::Format;
use abstar_germline::{Database, Locus, Manifest};
use arrow_array::cast::AsArray;
use arrow_array::Array;
use std::collections::BTreeMap;
use std::io::Write;
use std::path::{Path, PathBuf};

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data").join(name)
}

fn db() -> Database {
    Database {
        manifest: Manifest {
            format_version: 1,
            name: "t".into(),
            description: String::new(),
            abstar_version: "1.0.0".into(),
            built: "2026-11-20".into(),
            alleles_sha256: "abc".into(),
            species: vec!["human".into()],
            loci: vec![Locus::IGH, Locus::IGK],
            prefixes: BTreeMap::new(),
            provenance: vec![],
            warnings: vec![],
        },
        alleles: vec![],
    }
}

fn meta() -> RunMeta {
    RunMeta { abstar_version: "1.0.0".into(), germline_db: "t".into(), db_path: "/db".into() }
}

fn opts() -> RunOptions {
    RunOptions { threads: Some(2), ..Default::default() }
}

fn rows(inputs: &[Input]) -> Vec<(String, String, String, Option<String>)> {
    let mut out = vec![];
    for b in run_in_memory(inputs, &db(), &opts()).unwrap() {
        let s = b.schema();
        let col = |n: &str| b.column(s.index_of(n).unwrap()).as_string::<i32>().clone();
        let (id, seq, st, fr) = (col("sequence_id"), col("sequence"), col("annotation_status"), col("failure_reason"));
        for i in 0..b.num_rows() {
            let f = (!fr.is_null(i)).then(|| fr.value(i).to_string());
            out.push((id.value(i).into(), seq.value(i).into(), st.value(i).into(), f));
        }
    }
    out
}

fn tsv(inputs: &[Input], o: &RunOptions) -> Vec<u8> {
    let mut v = vec![];
    run_to_writer(inputs, &db(), o, &mut v).unwrap();
    v
}

#[test]
fn messy_fasta_one_row_per_record() {
    let r = rows(&[Input::File(data("messy.fasta"))]);
    let ids: Vec<&str> = r.iter().map(|x| x.0.as_str()).collect();
    assert_eq!(ids, ["10E8", "seq2", "lower", "empty", "dup", "dup", "0012"]);
    assert_eq!(r[0].1, "ACGTACGTACGTAC");
    assert_eq!(r[1].1, "ACGTACGTACGT");
    assert_eq!(r[2].1, "acgtnacgt");
    for (i, x) in r.iter().enumerate() {
        if i == 3 {
            assert_eq!((x.2.as_str(), x.3.as_deref()), ("invalid_input", Some("empty_sequence")));
        } else {
            assert_eq!((x.2.as_str(), x.3.as_deref()), ("unassigned", Some("engine_not_implemented")));
        }
    }
}

#[test]
fn gzip_and_fastq_inputs() {
    let d = tempfile::tempdir().unwrap();
    let p = d.path().join("messy.fasta");
    std::fs::copy(data("messy.fasta"), &p).unwrap();
    let plain = tsv(&[Input::File(p.clone())], &opts());
    std::fs::copy(data("messy.fasta.gz"), &p).unwrap(); // gzip content, detected by magic bytes
    assert_eq!(plain, tsv(&[Input::File(p)], &opts()));
    let r = rows(&[Input::File(data("reads.fastq"))]);
    assert_eq!(r.len(), 2);
    assert_eq!((r[0].0.as_str(), r[0].1.as_str()), ("r1", "ACGT"));
    assert_eq!((r[1].0.as_str(), r[1].1.as_str()), ("r2", "GGNN"));
}

#[test]
fn invalid_characters() {
    let r = rows(&[Input::Records(vec![("a".into(), "ACGTX".into()), ("b".into(), "ACGTRYN".into()), ("c".into(), "acgtrywskmbdhvn".into())])]);
    assert_eq!(r[0].3.as_deref(), Some("invalid_characters"));
    assert_eq!(r[0].2, "invalid_input");
    assert_eq!(r[1].2, "unassigned");
    assert_eq!(r[2].2, "unassigned");
}

#[test]
fn deterministic_across_threads_and_batches() {
    let one = std::fs::read_to_string(data("messy.fasta")).unwrap();
    let d = tempfile::tempdir().unwrap();
    let p = d.path().join("big.fasta");
    std::fs::write(&p, one.repeat(10_000 / 7 + 1)).unwrap();
    let a = RunOptions { threads: Some(1), batch_size: 7, ..Default::default() };
    let b = RunOptions { threads: Some(4), batch_size: 4096, ..Default::default() };
    let (ta, tb) = (tsv(&[Input::File(p.clone())], &a), tsv(&[Input::File(p)], &b));
    assert!(ta.iter().filter(|&&c| c == b'\n').count() > 10_000);
    assert_eq!(ta, tb);
}

#[test]
fn output_stem_conflict() {
    let d = tempfile::tempdir().unwrap();
    for s in ["a", "b"] {
        std::fs::create_dir(d.path().join(s)).unwrap();
        std::fs::copy(data("messy.fasta"), d.path().join(s).join("s1.fasta")).unwrap();
    }
    let out = d.path().join("out");
    let inputs = [Input::File(d.path().join("a/s1.fasta")), Input::File(d.path().join("b/s1.fasta"))];
    let e = run_to_files(&inputs, &out, &db(), &[Format::Airr], &opts(), &meta()).unwrap_err();
    assert_eq!(e.code(), "OUTPUT_NAME_CONFLICT");
    assert!(!out.exists() || std::fs::read_dir(&out).unwrap().next().is_none());
}

#[test]
fn stems_and_run_json() {
    let d = tempfile::tempdir().unwrap();
    let gz = d.path().join("x.fq.gz");
    std::fs::write(&gz, std::fs::read(data("messy.fasta.gz")).unwrap()).unwrap();
    let out = d.path().join("out");
    let inputs = [Input::File(gz), Input::File(data("reads.fastq")), Input::Records(vec![("q".into(), "AC".into())])];
    let s = run_to_files(&inputs, &out, &db(), &[Format::Airr, Format::Parquet], &opts(), &meta()).unwrap();
    for f in ["x.tsv", "x.parquet", "reads.tsv", "records.parquet", "run.json"] {
        assert!(out.join(f).exists(), "{f}");
    }
    assert_eq!(s.run_json, Some(out.join("run.json")));
    let j: serde_json::Value = serde_json::from_slice(&std::fs::read(out.join("run.json")).unwrap()).unwrap();
    let c = &j["counts"];
    assert_eq!(c["unassigned"].as_u64().unwrap() + c["invalid_input"].as_u64().unwrap(), 7 + 2 + 1);
    assert_eq!(j["airr_version"], "2.0");
    assert_eq!(j["database"]["alleles_sha256"], "abc");
    assert_eq!(j["inputs"].as_array().unwrap().len(), 3);
    assert!(j["started"].as_str().unwrap().ends_with('Z'));
    // re-run overwrites
    run_to_files(&inputs, &out, &db(), &[Format::Airr], &opts(), &meta()).unwrap();
    assert_eq!(std::fs::read_to_string(out.join("reads.tsv")).unwrap().lines().count(), 3);
}

fn boom(id: &str, seq: &str, db: &Database) -> Record {
    if id == "boom" {
        panic!("kaboom");
    }
    annotate_stub(id, seq, db)
}

fn recs() -> Vec<Input> {
    vec![Input::Records(vec![("a".into(), "ACGT".into()), ("boom".into(), "ACGT".into()), ("c".into(), "ACGT".into())])]
}

#[test]
fn panic_becomes_error_row() {
    let b = run_in_memory_with(&recs(), &db(), &opts(), &boom).unwrap();
    let col = b[0].column(arrow_schema().index_of("failure_reason").unwrap()).as_string::<i32>().clone();
    assert_eq!(b[0].num_rows(), 3);
    assert_eq!(col.value(1), "internal_error: kaboom");
    let st = b[0].column(arrow_schema().index_of("annotation_status").unwrap()).as_string::<i32>().clone();
    assert_eq!((st.value(0), st.value(1), st.value(2)), ("unassigned", "error", "unassigned"));
}

#[test]
fn strict_aborts() {
    let o = RunOptions { strict: true, ..opts() };
    let e = run_in_memory_with(&recs(), &db(), &o, &boom).unwrap_err();
    assert_eq!(e.code(), "STRICT_RECORD_FAILURE");
}

#[test]
fn zero_records() {
    let inputs = [Input::File(data("empty.fasta"))];
    let t = String::from_utf8(tsv(&inputs, &opts())).unwrap();
    assert_eq!(t.lines().count(), 1);
    let b = run_in_memory(&inputs, &db(), &opts()).unwrap();
    assert_eq!(b.len(), 1);
    assert_eq!(b[0].num_rows(), 0);
    assert_eq!(b[0].schema(), arrow_schema());
}

#[test]
fn argument_and_input_errors() {
    let o = RunOptions { loci: Some(vec![Locus::TRG]), ..opts() };
    let e = run_in_memory(&recs(), &db(), &o).unwrap_err();
    assert_eq!(e.code(), "INVALID_ARGUMENT");
    let o = RunOptions { loci: Some(vec![Locus::IGH]), ..opts() };
    assert!(run_in_memory(&recs(), &db(), &o).is_ok());
    let e = run_in_memory(&[Input::File("/nonexistent/x.fa".into())], &db(), &opts()).unwrap_err();
    assert!(matches!(e, RunError::InputNotFound(_)));
    assert_eq!(e.code(), "INPUT_NOT_FOUND");
    let d = tempfile::tempdir().unwrap();
    let p = d.path().join("bad.fasta");
    std::fs::write(&p, "not a fasta\n").unwrap();
    let e = run_in_memory(&[Input::File(p)], &db(), &opts()).unwrap_err();
    assert_eq!(e.code(), "INPUT_UNREADABLE");
}

struct PipeAfter(usize);
impl Write for PipeAfter {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        if self.0 == 0 {
            return Err(std::io::ErrorKind::BrokenPipe.into());
        }
        let n = b.len().min(self.0);
        self.0 -= n;
        Ok(n)
    }
    fn flush(&mut self) -> std::io::Result<()> {
        if self.0 == 0 { Err(std::io::ErrorKind::BrokenPipe.into()) } else { Ok(()) }
    }
}

#[test]
fn broken_pipe_is_end_of_output() {
    let one = std::fs::read_to_string(data("messy.fasta")).unwrap();
    let d = tempfile::tempdir().unwrap();
    let p = d.path().join("big.fasta");
    std::fs::write(&p, one.repeat(2000)).unwrap();
    let o = RunOptions { batch_size: 100, ..opts() };
    assert!(run_to_writer(&[Input::File(p.clone())], &db(), &o, PipeAfter(50)).is_ok());
    assert!(run_to_writer(&[Input::File(p)], &db(), &o, PipeAfter(200_000)).is_ok());
}

#[test]
fn tsv_hostile_characters_do_not_panic() {
    let inputs = [Input::Records(vec![("a\tb".into(), "AC\tGT".into())])];
    let t = String::from_utf8(tsv(&inputs, &opts())).unwrap();
    assert_eq!(t.lines().count(), 2);
}

#[test]
fn failed_rerun_removes_run_json_and_keeps_unopened_outputs() {
    let d = tempfile::tempdir().unwrap();
    let (a, b) = (d.path().join("a.fasta"), d.path().join("b.fasta"));
    std::fs::copy(data("messy.fasta"), &a).unwrap();
    std::fs::copy(data("messy.fasta"), &b).unwrap();
    let out = d.path().join("out");
    let inputs = [Input::File(a), Input::File(b.clone())];
    run_to_files(&inputs, &out, &db(), &[Format::Airr], &opts(), &meta()).unwrap();
    let good = std::fs::read(out.join("b.tsv")).unwrap();
    assert!(out.join("run.json").exists());
    std::fs::write(&b, "not a fasta\n").unwrap();
    let e = run_to_files(&inputs, &out, &db(), &[Format::Airr], &opts(), &meta()).unwrap_err();
    assert_eq!(e.code(), "INPUT_UNREADABLE");
    assert!(!out.join("run.json").exists(), "a failed run must not look finished");
    assert_eq!(std::fs::read(out.join("b.tsv")).unwrap(), good, "unopenable input truncated its output");
}
