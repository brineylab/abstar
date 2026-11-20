mod common;
use std::io::{self, Write};
use std::path::{Path, PathBuf};

use abstar_generate::{
    bind, generate_fasta_to_writer, generate_in_memory, generate_to_files, load_model, preset, seed_for, Bound, Dist,
    GenerateMeta, GenerateParams, TrimSet, TRUTH_TOML,
};
use abstar_germline::Allele;
use abstar_table::writers::Format;
use arrow_array::cast::AsArray;
use arrow_array::types::{Float64Type, Int64Type};
use arrow_array::{Array, RecordBatch};
use common::*;

fn model(name: &str) -> abstar_generate::Model {
    load_model(&model_path(name)).unwrap()
}

fn model_path(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(format!("tests/data/{name}.model.toml"))
}

fn igh() -> Bound {
    let mut m = model("tiny_igh");
    m.d_inversion = 0.3;
    bind(m, &db(sim_igh_alleles())).unwrap()
}

fn fixed(v: i64) -> Dist {
    Dist { values: vec![v], probs: vec![1.0] }
}

/// IGK with no trims and `n` N bases, all `A`: the read is V + A^n + J exactly.
fn igk_fixed(n: i64) -> Bound {
    let mut m = model("tiny_igk");
    m.trim.v_3 = TrimSet { default: fixed(0), genes: Default::default() };
    m.trim.j_5 = TrimSet { default: fixed(0), genes: Default::default() };
    m.insertions.insert("vj".into(), fixed(n));
    for k in ["start", "A", "C", "G", "T"] {
        m.n_markov.insert(k.into(), [1.0, 0.0, 0.0, 0.0]);
    }
    bind(m, &db(sim_igk_alleles())).unwrap()
}

fn params(n: usize, seed: u64) -> GenerateParams {
    GenerateParams { n, seed, ..preset("naive").unwrap() }
}

/// Mutation, indels (frameshifting too), truncation and a C prefix.
fn busy(n: usize, seed: u64) -> GenerateParams {
    GenerateParams {
        mutation: "uniform:0-0.3".parse().unwrap(),
        indel_rate: 0.5,
        frameshift_indels: true,
        truncate_5p: "uniform:0-40".parse().unwrap(),
        truncate_3p: "uniform:0-40".parse().unwrap(),
        constant_len: Some(30),
        ..params(n, seed)
    }
}

/// A germline database on disk plus run metadata pointing at it and the IGH model file.
fn meta(dir: &Path, alleles: Vec<Allele>) -> GenerateMeta {
    let db_path = dir.join("db");
    abstar_germline::write(&db(alleles), &db_path, false).unwrap();
    GenerateMeta {
        abstar_version: "1.0.0-test".into(),
        germline_db: "tiny".into(),
        db_path,
        model_name: "tiny_igh".into(),
        model_path: model_path("tiny_igh"),
        preset: Some("naive".into()),
    }
}

fn truth_names() -> Vec<String> {
    let doc: toml::Table = TRUTH_TOML.parse().unwrap();
    doc["field"].as_array().unwrap().iter().map(|f| f["name"].as_str().unwrap().to_string()).collect()
}

fn read(p: &Path) -> String {
    std::fs::read_to_string(p).unwrap()
}

fn col<'a>(b: &'a RecordBatch, name: &str) -> &'a dyn Array {
    b.column_by_name(name).unwrap_or_else(|| panic!("no column {name}")).as_ref()
}

fn s(b: &RecordBatch, name: &str, row: usize) -> Option<String> {
    let c = col(b, name).as_string::<i32>();
    c.is_valid(row).then(|| c.value(row).to_string())
}

fn i(b: &RecordBatch, name: &str, row: usize) -> Option<i64> {
    let c = col(b, name).as_primitive::<Int64Type>();
    c.is_valid(row).then(|| c.value(row))
}

fn f(b: &RecordBatch, name: &str, row: usize) -> Option<f64> {
    let c = col(b, name).as_primitive::<Float64Type>();
    c.is_valid(row).then(|| c.value(row))
}

fn bl(b: &RecordBatch, name: &str, row: usize) -> Option<bool> {
    let c = col(b, name).as_boolean();
    c.is_valid(row).then(|| c.value(row))
}

fn ls(b: &RecordBatch, name: &str, row: usize) -> Vec<String> {
    let v = col(b, name).as_list::<i32>().value(row);
    let v = v.as_string::<i32>();
    (0..v.len()).map(|k| v.value(k).to_string()).collect()
}

fn li(b: &RecordBatch, name: &str, row: usize) -> Vec<i64> {
    let v = col(b, name).as_list::<i32>().value(row);
    v.as_primitive::<Int64Type>().values().to_vec()
}

fn one(batches: &[RecordBatch]) -> RecordBatch {
    assert_eq!(batches.len(), 1, "expected one batch");
    batches[0].clone()
}

#[test]
fn seed_for_is_splitmix() {
    fn sm(mut z: u64) -> u64 {
        z = z.wrapping_add(0x9E3779B97F4A7C15);
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }
    for (g, k) in [(0, 0), (7, 1), (u64::MAX, 12345)] {
        assert_eq!(seed_for(g, k), sm(g ^ sm(k)));
    }
}

#[test]
fn byte_identical_across_threads_and_batches() {
    let b = igh();
    let dir = tempfile::tempdir().unwrap();
    let m = meta(dir.path(), sim_igh_alleles());
    let run = |threads, batch_size, sub: &str| {
        let p = GenerateParams { threads: Some(threads), batch_size, ..busy(5000, 7) };
        let out = dir.path().join(sub);
        generate_to_files(&b, &p, &out, "x", &[Format::Airr], &m).unwrap();
        (read(&out.join("x.fasta")), read(&out.join("x.truth.tsv")))
    };
    let (fa1, tsv1) = run(1, 3, "a");
    let (fa2, tsv2) = run(4, 4096, "b");
    assert_eq!(fa1.lines().count(), 10000);
    assert_eq!(tsv1.lines().count(), 5001);
    assert!(fa1 == fa2, "FASTA differs across threads/batches");
    assert!(tsv1 == tsv2, "TSV differs across threads/batches");
}

#[test]
fn files_and_meta() {
    let b = igh();
    let dir = tempfile::tempdir().unwrap();
    let m = meta(dir.path(), sim_igh_alleles());
    let out = dir.path().join("out");
    let p = busy(50, 11);
    let sum = generate_to_files(&b, &p, &out, "run1", &[Format::Airr, Format::Parquet], &m).unwrap();
    for f in ["run1.fasta", "run1.truth.tsv", "run1.truth.parquet", "generate.json"] {
        assert!(out.join(f).is_file(), "{f} missing");
    }
    assert_eq!((sum.n, sum.seed), (50, 11));
    assert_eq!(sum.attempts, 50);
    assert_eq!(sum.fasta, Some(out.join("run1.fasta")));
    assert_eq!(sum.truth["airr"], out.join("run1.truth.tsv"));
    assert_eq!(sum.truth["parquet"], out.join("run1.truth.parquet"));
    assert_eq!(sum.generate_json, Some(out.join("generate.json")));

    let j: serde_json::Value = serde_json::from_str(&read(&out.join("generate.json"))).unwrap();
    let keys: Vec<&str> = j.as_object().unwrap().keys().map(String::as_str).collect();
    let mut want = vec![
        "abstar_version", "started", "finished", "seed", "parameters", "model", "database", "counts", "warnings",
    ];
    want.sort();
    let mut got = keys.clone();
    got.sort();
    assert_eq!(got, want);
    assert_eq!(j["seed"], 11);
    assert_eq!(j["abstar_version"], "1.0.0-test");
    assert_eq!(j["counts"]["n"], 50);
    assert_eq!(j["counts"]["attempts"], 50);
    assert_eq!(j["counts"]["productive"], sum.productive);
    assert_eq!(j["parameters"]["preset"], "naive");
    assert_eq!(j["parameters"]["indel_rate"], 0.5);
    assert_eq!(j["parameters"]["truncate_5p"], "uniform:0-40");
    assert_eq!(j["parameters"]["constant_len"], 30);
    assert_eq!(j["model"]["name"], "tiny_igh");
    assert_eq!(j["model"]["sha256"].as_str().unwrap().len(), 64);
    assert_eq!(j["database"]["name"], "tiny");
    assert_eq!(j["database"]["alleles_sha256"].as_str().unwrap().len(), 64);

    let tsv = read(&out.join("run1.truth.tsv"));
    let header: Vec<String> = tsv.lines().next().unwrap().split('\t').map(String::from).collect();
    assert_eq!(header, truth_names());
    assert_eq!(tsv.lines().count(), 51);
    let fa = read(&out.join("run1.fasta"));
    assert!(fa.starts_with(">g01\n"), "{}", &fa[..20]);
    assert!(fa.contains("\n>g50\n"));

    // invalid parameters fail before any file is touched: the previous run stays whole
    let bad = GenerateParams { indel_rate: 2.0, ..p };
    assert_eq!(generate_to_files(&b, &bad, &out, "run1", &[Format::Airr], &m).unwrap_err().code(), "INVALID_ARGUMENT");
    assert!(out.join("generate.json").exists());
    assert_eq!(read(&out.join("run1.truth.tsv")), tsv);
}

#[test]
fn zero_n() {
    let b = igh();
    let dir = tempfile::tempdir().unwrap();
    let m = meta(dir.path(), sim_igh_alleles());
    let out = dir.path().join("out");
    let sum = generate_to_files(&b, &params(0, 1), &out, "z", &[Format::Airr], &m).unwrap();
    assert_eq!((sum.n, sum.productive, sum.attempts), (0, 0, 0));
    assert_eq!(read(&out.join("z.fasta")), "");
    assert_eq!(read(&out.join("z.truth.tsv")), truth_names().join("\t") + "\n");
    let batches = generate_in_memory(&b, &params(0, 1)).unwrap();
    assert_eq!(batches.len(), 1);
    assert_eq!(batches[0].num_rows(), 0);
}

#[test]
fn productive_only_limit() {
    // V and J joined with exactly one N base: always out of frame
    let b = igk_fixed(1);
    let p = GenerateParams { productive_only: true, ..params(3, 1) };
    let e = generate_in_memory(&b, &p).unwrap_err();
    assert_eq!(e.code(), "INVALID_ARGUMENT");
    assert!(e.to_string().contains("productive rate too low"), "{e}");

    // the same failure through generate_to_files removes a stale generate.json
    let dir = tempfile::tempdir().unwrap();
    let m = meta(dir.path(), sim_igk_alleles());
    let out = dir.path().join("out");
    std::fs::create_dir_all(&out).unwrap();
    std::fs::write(out.join("generate.json"), "{}").unwrap();
    assert_eq!(generate_to_files(&b, &p, &out, "x", &[Format::Airr], &m).unwrap_err().code(), "INVALID_ARGUMENT");
    assert!(!out.join("generate.json").exists());
}

#[test]
fn productive_only_counts_attempts() {
    // three N bases keep the frame; zero or one N breaks it
    let mut m = model("tiny_igk");
    m.insertions.insert("vj".into(), Dist { values: vec![1, 3], probs: vec![0.5, 0.5] });
    m.trim.v_3 = TrimSet { default: fixed(0), genes: Default::default() };
    m.trim.j_5 = TrimSet { default: fixed(0), genes: Default::default() };
    let b = bind(m, &db(sim_igk_alleles())).unwrap();
    let p = GenerateParams { productive_only: true, ..params(200, 5) };
    let mut buf = Vec::new();
    let sum = generate_fasta_to_writer(&b, &p, &mut buf).unwrap();
    assert_eq!((sum.n, sum.productive), (200, 200));
    assert!(sum.attempts > 250 && sum.attempts < 600, "{}", sum.attempts);
    let rows = one(&generate_in_memory(&b, &p).unwrap());
    assert!((0..200).all(|r| bl(&rows, "productive", r) == Some(true)));
}

struct Failing(usize);

impl Write for Failing {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        if self.0 + buf.len() > 100 {
            return Err(io::ErrorKind::BrokenPipe.into());
        }
        self.0 += buf.len();
        Ok(buf.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

#[test]
fn fasta_broken_pipe() {
    assert!(generate_fasta_to_writer(&igh(), &params(500, 3), Failing(0)).is_ok());
}

#[test]
fn junction_removed_by_truncation() {
    let b = igk_fixed(3);
    let p = GenerateParams { truncate_3p: "fixed:40".parse().unwrap(), ..params(20, 2) };
    let rows = one(&generate_in_memory(&b, &p).unwrap());
    assert_eq!(rows.num_rows(), 20);
    for r in 0..20 {
        assert!(s(&rows, "sequence", r).unwrap().len() > 300);
        assert_eq!(s(&rows, "junction", r), None);
        assert_eq!(s(&rows, "junction_aa", r), None);
        assert_eq!(i(&rows, "junction_length", r), None);
        assert_eq!(bl(&rows, "long_cdr3", r), None);
        assert_eq!(bl(&rows, "productive", r), Some(false));
        assert_eq!(bl(&rows, "truncated_3p", r), Some(true));
        assert_eq!(i(&rows, "truncate_3p", r), Some(40));
    }
}

/// Every coordinate of one fully determined read: IGKV1-1*01 (311 nt, Cys104 at 289),
/// `AAA`, IGKJ1*01 (38 nt, Phe118 at 8), the first 10 nt of IGKC*01; 5 nt cut from the 5' end.
#[test]
fn hand_traced_row() {
    let b = igk_fixed(3);
    let p = GenerateParams {
        truncate_5p: "fixed:5".parse().unwrap(),
        constant_len: Some(10),
        ..params(1, 9)
    };
    let rows = one(&generate_in_memory(&b, &p).unwrap());
    let a = sim_igk_alleles();
    let (v, j, c) = (a[0].sequence(), a[1].sequence(), a[2].sequence());
    assert_eq!((v.len(), j.len()), (311, 38));
    let read = format!("{}AAA{}{}", &v[5..], j, &c[..10]);
    assert_eq!(s(&rows, "sequence_id", 0).as_deref(), Some("g1"));
    assert_eq!(s(&rows, "sequence", 0), Some(read.clone()));
    assert_eq!(bl(&rows, "rev_comp", 0), Some(false));
    assert_eq!(s(&rows, "locus", 0).as_deref(), Some("IGK"));
    assert_eq!(s(&rows, "v_call", 0).as_deref(), Some("IGKV1-1*01"));
    assert_eq!(s(&rows, "d_call", 0), None);
    assert_eq!(s(&rows, "j_call", 0).as_deref(), Some("IGKJ1*01"));
    assert_eq!(s(&rows, "c_call", 0).as_deref(), Some("IGKC*01"));
    assert_eq!(s(&rows, "v_gene", 0).as_deref(), Some("IGKV1-1"));
    assert_eq!(s(&rows, "c_gene", 0).as_deref(), Some("IGKC"));
    // junction: Cys104 (V 289-291) .. Phe118 (J 8-10), 23 + 3 + 10 nt
    let junction = format!("{}AAA{}", &v[288..], &j[..10]);
    assert_eq!(junction, "TGTCAACAGTATAATAGTTACCCAAAGTGGACGTTC");
    assert_eq!(s(&rows, "junction", 0), Some(junction));
    assert_eq!(s(&rows, "junction_aa", 0).as_deref(), Some("CQQYNSYPKWTF"));
    assert_eq!(i(&rows, "junction_length", 0), Some(36));
    assert_eq!(bl(&rows, "productive", 0), Some(true));
    assert_eq!(bl(&rows, "vj_in_frame", 0), Some(true));
    assert_eq!(bl(&rows, "stop_codon", 0), Some(false));
    // V: read 1-306 = germline 6-311; N: 307-309; J: 310-347 = germline 1-38; C: 348-357
    let get = |n: &str| i(&rows, n, 0);
    assert_eq!((get("v_sequence_start"), get("v_sequence_end")), (Some(1), Some(306)));
    assert_eq!((get("v_germline_start"), get("v_germline_end")), (Some(6), Some(311)));
    assert_eq!((get("d_sequence_start"), get("d_germline_start")), (None, None));
    assert_eq!((get("j_sequence_start"), get("j_sequence_end")), (Some(310), Some(347)));
    assert_eq!((get("j_germline_start"), get("j_germline_end")), (Some(1), Some(38)));
    assert_eq!((get("c_sequence_start"), get("c_sequence_end")), (Some(348), Some(357)));
    assert_eq!(s(&rows, "np1", 0).as_deref(), Some("AAA"));
    assert_eq!(get("np1_length"), Some(3));
    assert_eq!((s(&rows, "np2", 0), get("np2_length")), (None, None));
    assert_eq!(f(&rows, "v_identity", 0), Some(1.0));
    assert_eq!(ls(&rows, "v_call_valid", 0), ["IGKV1-1*01"]);
    assert_eq!(ls(&rows, "j_gene_valid", 0), ["IGKJ1"]);
    assert_eq!(ls(&rows, "d_call_valid", 0), Vec::<String>::new());
    assert_eq!(bl(&rows, "d_call_null_valid", 0), Some(true));
    assert_eq!((get("v_trim_3"), get("j_trim_5"), get("d_trim_5")), (Some(0), Some(0), None));
    assert_eq!(bl(&rows, "d_inverted", 0), None);
    assert_eq!(ls(&rows, "mutations", 0), Vec::<String>::new());
    assert_eq!(ls(&rows, "indels", 0), Vec::<String>::new());
    assert_eq!((get("truncate_5p"), get("truncate_3p")), (Some(5), Some(0)));
    assert_eq!(f(&rows, "v_mutation_frequency", 0), Some(0.0));
    for (flag, want) in [
        ("high_mutation", false),
        ("long_cdr3", false),
        ("has_indel", false),
        ("truncated_5p", true),
        ("truncated_3p", false),
        ("unrecoverable_v", false),
        ("ambiguous_v", false),
        ("ambiguous_d", false),
    ] {
        assert_eq!(bl(&rows, flag, 0), Some(want), "{flag}");
    }
}

/// Event strings point at the read: each mutation's `to` base is at its position, each
/// indel's position lies in its equivalence range, and coordinates stay inside the read.
#[test]
fn events_are_consistent_with_the_read() {
    let rows = one(&generate_in_memory(&igh(), &GenerateParams { batch_size: 1000, ..busy(400, 4) }).unwrap());
    let (mut muts, mut indels) = (0, 0);
    for r in 0..rows.num_rows() {
        let seq = s(&rows, "sequence", r).unwrap().into_bytes();
        let mut last = 0;
        for m in ls(&rows, "mutations", r) {
            let (pos, sub) = m.split_once(':').unwrap();
            let (from, to) = sub.split_once('>').unwrap();
            let pos: usize = pos.parse().unwrap();
            assert!(pos > last, "{m}: positions sorted and distinct");
            last = pos;
            assert_ne!(from, to);
            assert_eq!(seq[pos - 1] as char, to.chars().next().unwrap(), "{m}");
            muts += 1;
        }
        let ind = ls(&rows, "indels", r);
        assert_eq!(bl(&rows, "has_indel", r), Some(!ind.is_empty()));
        for x in ind {
            let p: Vec<&str> = x.split(':').collect();
            let (pos, len): (i64, usize) = (p[1].parse().unwrap(), p[2].parse().unwrap());
            let (lo, hi) = p[4].split_once('-').unwrap();
            let (lo, hi): (i64, i64) = (lo.parse().unwrap(), hi.parse().unwrap());
            assert!(lo <= pos && pos <= hi && lo >= 2, "{x}");
            assert_eq!(p[3].len(), len);
            if p[0] == "ins" {
                assert_eq!(&seq[pos as usize - 1..pos as usize - 1 + len], p[3].as_bytes(), "{x}");
            } else {
                assert_eq!(p[0], "del");
            }
            assert!((hi as usize) <= seq.len());
            indels += 1;
        }
        let ends = ["v_sequence_end", "d_sequence_end", "j_sequence_end", "c_sequence_end"];
        for e in ends {
            if let Some(x) = i(&rows, e, r) {
                assert!(x as usize <= seq.len());
            }
        }
        let vmf = f(&rows, "v_mutation_frequency", r).unwrap();
        if let Some(id) = f(&rows, "v_identity", r) {
            assert!((id + vmf - 1.0).abs() < 1e-12);
        }
        assert_eq!(bl(&rows, "high_mutation", r), Some(vmf >= 0.2));
        assert_eq!(bl(&rows, "ambiguous_v", r), Some(ls(&rows, "v_call_valid", r).len() > 1));
        let _ = li(&rows, "v_sequence_end_valid", r);
    }
    assert!(muts > 1000 && indels > 50, "{muts} {indels}");
}
