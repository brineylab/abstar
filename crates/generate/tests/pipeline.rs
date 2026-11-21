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
        // no J base: sequence_vdj runs to the last V/N base, the masks index it, and CDR3
        // holds the V tail and one N base (read 292-312; codons from read 1)
        let read = s(&rows, "sequence", r).unwrap();
        assert_eq!(read.len(), 312);
        assert_eq!(s(&rows, "sequence_vdj", r), Some(read.clone()));
        assert_eq!(s(&rows, "sequence_vdjc", r), None);
        assert_eq!((i(&rows, "cdr3_start", r), i(&rows, "cdr3_end", r), i(&rows, "cdr3_length", r)), (Some(292), Some(312), Some(7)));
        let mask = s(&rows, "cdr_mask", r).unwrap();
        assert!(mask.len() == 312 && mask.ends_with(&"3".repeat(21)) && mask.as_bytes()[290] == b'0', "{mask}");
        assert_eq!(s(&rows, "gene_segment_mask", r), Some("V".repeat(311) + "N"));
        assert_eq!(s(&rows, "nongermline_mask_aa", r).map(|x| x.len()), Some(104));
        assert_eq!(s(&rows, "fwr4", r), None);
    }
}

/// V trimmed into the Cys104 codon, keeping its first base: the junction exists, and the
/// regions follow read offsets from the Cys codon start. Read: V 6-289 (Cys T at read 284),
/// `AAA` (285-287), IGKJ1*01 (288-325, Phe118 at 295-297).
#[test]
fn cdr3_by_cys_offset_when_v_trimmed_into_cys() {
    let mut m = model("tiny_igk");
    m.trim.v_3 = TrimSet { default: fixed(22), genes: Default::default() };
    m.trim.j_5 = TrimSet { default: fixed(0), genes: Default::default() };
    m.insertions.insert("vj".into(), fixed(3));
    for k in ["start", "A", "C", "G", "T"] {
        m.n_markov.insert(k.into(), [1.0, 0.0, 0.0, 0.0]);
    }
    let b = bind(m, &db(sim_igk_alleles())).unwrap();
    let p = GenerateParams { truncate_5p: "fixed:5".parse().unwrap(), ..params(1, 9) };
    let rows = one(&generate_in_memory(&b, &p).unwrap());
    let read = s(&rows, "sequence", 0).unwrap();
    let get = |n: &str| i(&rows, n, 0);
    let junction = s(&rows, "junction", 0).unwrap();
    assert_eq!(junction, read[283..297]);
    assert_eq!(junction, "TAAAGTGGACGTTC");
    // FWR3 keeps the Cys codon position (T + two N bases); CDR3 = junction[3..len-3]
    assert_eq!((get("fwr3_start"), get("fwr3_end")), (Some(170), Some(286)));
    assert!(s(&rows, "fwr3_aa", 0).unwrap().ends_with('*'));
    assert_eq!((get("cdr3_start"), get("cdr3_end")), (Some(287), Some(294)));
    assert_eq!(s(&rows, "cdr3", 0).as_deref(), Some(&junction[3..junction.len() - 3]));
    assert_eq!((s(&rows, "cdr3_aa", 0).as_deref(), get("cdr3_length")), (Some("SG"), Some(2)));
    assert_eq!((get("fwr4_start"), get("fwr4_end")), (Some(295), Some(325)));
    let rep = |c: &str, n: usize| c.repeat(n);
    let mask = [rep("0", 70), rep("1", 24), rep("0", 51), rep("2", 24), rep("0", 117), rep("3", 8), rep("0", 31)].concat();
    assert_eq!(s(&rows, "cdr_mask", 0), Some(mask));
    assert_eq!(s(&rows, "gene_segment_mask", 0), Some(["V".repeat(284), "NNN".into(), "J".repeat(38)].concat()));
}

/// No V base in the read: only FWR4 is valued, phased by the J anchor.
#[test]
fn fwr4_without_v() {
    let b = igk_fixed(3);
    // V (311) + AAA + 6 J bases cut: the read is IGKJ1*01 7-38, Phe118 at read 2-4
    let p = GenerateParams { truncate_5p: "fixed:320".parse().unwrap(), ..params(1, 9) };
    let rows = one(&generate_in_memory(&b, &p).unwrap());
    assert_eq!(s(&rows, "sequence", 0).unwrap().len(), 32);
    assert_eq!((i(&rows, "fwr4_start", 0), i(&rows, "fwr4_end", 0)), (Some(2), Some(32)));
    assert_eq!((s(&rows, "fwr4_aa", 0).as_deref(), i(&rows, "fwr4_length", 0)), (Some("FGQGTKVEIK"), Some(10)));
    for k in ["cdr3", "fwr3", "sequence_aa", "sequence_vdj", "cdr_mask", "gene_segment_mask"] {
        assert_eq!(s(&rows, k, 0), None, "{k}");
    }
}

/// A V allele that is not IMGT-gapped: every V mutation field is null, frequencies included.
/// `bind` refuses such alleles, so the bound V is ungapped after binding.
#[test]
fn ungapped_v_has_null_v_fields() {
    let mut m = model("tiny_igk");
    m.shm_model = "motif".into();
    let mut b = bind(m, &db(sim_igk_alleles())).unwrap();
    for (a, _) in &mut b.v {
        a.sequence_gapped.retain(|c| c != '.');
        a.anchor = None;
    }
    let p = GenerateParams { mutation: "fixed:0.1".parse().unwrap(), ..params(5, 3) };
    let rows = one(&generate_in_memory(&b, &p).unwrap());
    for r in 0..5 {
        assert!(col(&rows, "v_mutations").is_null(r));
        assert_eq!((f(&rows, "v_mutation_frequency", r), f(&rows, "v_mutation_frequency_aa", r)), (None, None));
        assert!(!ls(&rows, "j_mutations", r).is_empty() || i(&rows, "j_mutation_count", r) == Some(0));
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

/// Regions, amino-acid fields and masks of the hand-traced IGK row. Read 1 is IMGT nt 6, so
/// codons start at read positions ≡ 2 (mod 3). Gaps: IMGT codons 10, 31-34 and 60-61.
#[test]
fn region_truth_hand_traced() {
    let b = igk_fixed(3);
    let p = GenerateParams { truncate_5p: "fixed:5".parse().unwrap(), constant_len: Some(10), ..params(1, 9) };
    let rows = one(&generate_in_memory(&b, &p).unwrap());
    let read = s(&rows, "sequence", 0).unwrap();
    let get = |n: &str| i(&rows, n, 0);
    let st = |n: &str| s(&rows, n, 0);
    // (region, start, end, aa length): FWR1 IMGT 6-78 less one gap codon = 70 nt, CDR1 24,
    // FWR2 51, CDR2 24, FWR3 117 (Cys104 at read 284-286), CDR3 30, FWR4 J 8-38 = 31 nt
    let want = [
        ("fwr1", 1, 70, 23),
        ("cdr1", 71, 94, 8),
        ("fwr2", 95, 145, 17),
        ("cdr2", 146, 169, 8),
        ("fwr3", 170, 286, 39),
        ("cdr3", 287, 316, 10),
        ("fwr4", 317, 347, 10),
    ];
    let mut aa = String::new();
    for (r, start, end, len) in want {
        assert_eq!((get(&format!("{r}_start")), get(&format!("{r}_end"))), (Some(start), Some(end)), "{r}");
        assert_eq!(st(r).as_deref(), Some(&read[start as usize - 1..end as usize]), "{r}");
        let r_aa = st(&format!("{r}_aa")).unwrap();
        assert_eq!((r_aa.len() as i64, get(&format!("{r}_length"))), (len, Some(len)), "{r}");
        aa += &r_aa;
    }
    let junction = st("junction").unwrap();
    assert_eq!(st("cdr3").as_deref(), Some(&junction[3..junction.len() - 3]));
    assert_eq!(st("cdr3_aa").as_deref(), Some("QQYNSYPKWT"));
    assert_eq!(st("fwr4_aa").as_deref(), Some("FGQGTKVEIK"));
    assert!(aa.contains("CQQYNSYPKWTF"));
    // V(D)J = read 1-347; the C codons (J's last base + C) translate to PNC
    assert_eq!(st("sequence_vdj").as_deref(), Some(&read[..347]));
    assert_eq!(st("sequence_vdjc"), Some(read.clone()));
    assert_eq!(st("sequence_vdj_aa"), Some(aa.clone()));
    assert_eq!(st("sequence_vdjc_aa"), Some(format!("{aa}PNC")));
    assert_eq!(st("sequence_aa"), Some(format!("{aa}PNC")));
    let rep = |c: &str, n: usize| c.repeat(n);
    let cdr = [rep("0", 70), rep("1", 24), rep("0", 51), rep("2", 24), rep("0", 117), rep("3", 30), rep("0", 31)];
    assert_eq!(st("cdr_mask"), Some(cdr.concat()));
    assert_eq!(st("gene_segment_mask"), Some([rep("V", 306), rep("N", 3), rep("J", 38)].concat()));
    assert_eq!(st("nongermline_mask"), Some([rep("0", 306), rep("1", 3), rep("0", 38)].concat()));
    let cdr_aa = [rep("0", 23), rep("1", 8), rep("0", 17), rep("2", 8), rep("0", 39), rep("3", 10), rep("0", 10)];
    assert_eq!(st("cdr_mask_aa"), Some(cdr_aa.concat()));
    // codons 102 (read 305-307) and 103 (308-310) hold the N bases
    assert_eq!(st("gene_segment_mask_aa"), Some([rep("V", 101), rep("N", 1), rep("J", 13)].concat()));
    assert_eq!(st("nongermline_mask_aa"), Some([rep("0", 101), rep("1", 2), rep("0", 12)].concat()));
    // unmutated: empty lists, zero counts; no D; C codons are J-frame codons 3-5, 6-8 of C
    for seg in ["v", "j", "c"] {
        for k in ["mutations", "mutations_aa", "insertions", "deletions"] {
            assert_eq!(ls(&rows, &format!("{seg}_{k}"), 0), Vec::<String>::new(), "{seg}_{k}");
        }
        assert_eq!(get(&format!("{seg}_mutation_count")), Some(0));
        assert_eq!(get(&format!("{seg}_mutation_count_aa")), Some(0));
        assert_eq!(f(&rows, &format!("{seg}_mutation_frequency_aa"), 0), Some(0.0), "{seg}");
    }
    assert_eq!(f(&rows, "j_mutation_frequency", 0), Some(0.0));
    for k in ["d_mutations", "d_mutations_aa", "d_insertions", "d_deletions"] {
        assert!(col(&rows, k).is_null(0), "{k}");
    }
    assert_eq!((get("d_mutation_count"), f(&rows, "d_mutation_frequency", 0)), (None, None));
}

/// Read `i` of a run, seeded exactly as the pipeline seeds it (no `productive_only`).
fn simulated(b: &Bound, p: &GenerateParams, i: usize) -> abstar_generate::Simulated {
    use rand::SeedableRng;
    let sim = abstar_generate::Simulator::new(b, p).unwrap();
    sim.simulate(&mut rand_chacha::ChaCha8Rng::seed_from_u64(seed_for(p.seed, i as u64)))
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter().rev().map(|b| match b { b'A' => b'T', b'C' => b'G', b'G' => b'C', _ => b'A' }).collect()
}

/// The germline base behind a read base, or `None` for N, P and inserted bases.
fn germline_base(s: &abstar_generate::Simulated, o: abstar_generate::Origin) -> Option<u8> {
    use abstar_generate::Origin::*;
    let d = s.d.as_ref().map(|d| if s.d_inverted { revcomp(d.sequence().as_bytes()) } else { d.sequence().into_bytes() });
    match o {
        V(g) => Some(s.v.sequence_gapped.as_bytes()[g]),
        D(u) => Some(d.unwrap()[u]),
        J(u) => Some(s.j.sequence().as_bytes()[u]),
        C(u) => Some(s.c.as_ref().unwrap().sequence().as_bytes()[u]),
        _ => None,
    }
}

/// Every option on: masks agree with the origin map and have the right lengths.
#[test]
fn masks_consistent() {
    let b = igh();
    let p = GenerateParams { indel_duplication: true, batch_size: 2000, ..busy(2000, 21) };
    let rows = one(&generate_in_memory(&b, &p).unwrap());
    let (mut checked, mut empty) = (0, 0);
    for r in 0..rows.num_rows() {
        let sim = simulated(&b, &p, r);
        assert_eq!(s(&rows, "sequence", r).unwrap().as_bytes(), &sim.read[..]);
        if let Some(j) = s(&rows, "junction", r) {
            assert_eq!(s(&rows, "cdr3", r).as_deref(), Some(&j[3..j.len() - 3]));
        }
        if s(&rows, "cdr3", r).is_some() {
            assert_eq!(i(&rows, "cdr3_length", r), Some(s(&rows, "cdr3_aa", r).unwrap().len() as i64));
            // the regions in the read tile sequence_vdj, in nt and in aa
            let regions = ["fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4"];
            let cat = |sfx: &str| regions.iter().filter_map(|x| s(&rows, &format!("{x}{sfx}"), r)).collect::<String>();
            if let Some(vdj) = s(&rows, "sequence_vdj", r) {
                assert_eq!(cat(""), vdj, "row {r}");
                assert_eq!(cat("_aa"), s(&rows, "sequence_vdj_aa", r).unwrap(), "row {r}");
            }
            empty += regions.iter().filter(|x| s(&rows, x, r).is_some_and(|v| v.is_empty())).count();
        } else {
            assert_eq!(s(&rows, "cdr_mask", r), None);
        }
        // the masks cover the V(D)J region only, never the constant region
        let (nt, aa) = (s(&rows, "sequence_vdj", r), s(&rows, "sequence_vdj_aa", r));
        let Some(nt) = nt else {
            assert_eq!(s(&rows, "gene_segment_mask", r), None);
            continue;
        };
        let aa = aa.unwrap();
        let off = i(&rows, "v_sequence_start", r).unwrap() as usize - 1;
        assert_eq!(&sim.read[off..off + nt.len()], nt.as_bytes());
        for m in ["gene_segment_mask", "nongermline_mask", "cdr_mask"] {
            if let Some(x) = s(&rows, m, r) {
                assert_eq!(x.len(), nt.len(), "{m}");
            }
            if let Some(x) = s(&rows, &format!("{m}_aa"), r) {
                assert_eq!(x.len(), aa.len(), "{m}_aa");
            }
        }
        let seg = s(&rows, "gene_segment_mask", r).unwrap().into_bytes();
        let ng = s(&rows, "nongermline_mask", r).unwrap().into_bytes();
        let mut prev = b'V';
        for k in 0..nt.len() {
            use abstar_generate::Origin::*;
            let o = sim.origin[off + k];
            let want = match o {
                V(_) => b'V',
                D(_) => b'D',
                J(_) => b'J',
                C(_) => b'C',
                N1 | N2 | P => b'N',
                Inserted => prev,
            };
            prev = want;
            assert_eq!(seg[k], want, "row {r} base {k}: {o:?}");
            let differs = germline_base(&sim, o).is_none_or(|g| g != sim.read[off + k]);
            assert_eq!(ng[k], if differs { b'1' } else { b'0' }, "row {r} base {k}: {o:?}");
        }
        checked += 1;
    }
    assert!(checked > 1500 && empty > 0, "{checked} {empty}");
}

/// Per-segment mutations are the recorded events in germline numbering; indels are the
/// recorded indels; frequencies use the segment's bases in the read as denominator.
#[test]
fn segment_mutations_from_events() {
    use abstar_generate::{IndelKind, Origin::*};
    let b = igh();
    let p = GenerateParams {
        mutation: "uniform:0.05-0.3".parse().unwrap(),
        indel_duplication: true,
        batch_size: 2000,
        ..busy(2000, 22)
    };
    let rows = one(&generate_in_memory(&b, &p).unwrap());
    let (mut inverted, mut v_aa, mut seg_indels, mut shifted, mut aa_without_nt) = (0, 0, 0, 0, 0);
    for r in 0..rows.num_rows() {
        let sim = simulated(&b, &p, r);
        let t5 = sim.truncate_5;
        let mut net = std::collections::BTreeMap::new();
        for m in &sim.mutations {
            net.entry(m.pos - t5).and_modify(|e: &mut (u8, u8)| e.1 = m.to).or_insert((m.from, m.to));
        }
        let (mut want, mut np): ([Vec<String>; 3], usize) = (Default::default(), 0);
        for (&pos, &(from, to)) in net.iter().filter(|(_, (f, t))| f != t) {
            let (k, label) = match sim.origin[pos] {
                V(g) => (0, format!("IMGT{}", g + 1)),
                D(u) => (1, (u + 1).to_string()),
                J(u) => (2, (u + 1).to_string()),
                N1 | N2 | P => {
                    np += 1;
                    continue;
                }
                o => panic!("mutation at {o:?}"),
            };
            want[k].push(format!("{label}:{}>{}", from as char, to as char));
        }
        let mut total = np;
        for (k, seg) in ["v", "d", "j"].iter().enumerate() {
            let n = sim.origin.iter().filter(|o| matches!((k, o), (0, V(_)) | (1, D(_)) | (2, J(_)))).count();
            if n == 0 {
                assert!(col(&rows, &format!("{seg}_mutations")).is_null(r), "{seg} row {r}");
                assert_eq!(f(&rows, &format!("{seg}_mutation_frequency"), r), None, "{seg} row {r}");
                continue;
            }
            let got = ls(&rows, &format!("{seg}_mutations"), r);
            assert_eq!(got, want[k], "{seg} row {r}");
            assert_eq!(i(&rows, &format!("{seg}_mutation_count"), r), Some(got.len() as i64));
            let freq = f(&rows, &format!("{seg}_mutation_frequency"), r).unwrap();
            assert!((freq - got.len() as f64 / n as f64).abs() < 1e-12, "{seg} row {r}");
            total += got.len();
            if k == 1 && sim.d_inverted && !got.is_empty() {
                inverted += 1;
            }
        }
        assert_eq!(total, ls(&rows, "mutations", r).len(), "row {r}");
        // each V aa change sits in an IMGT codon holding a V nt change, except next to a deletion: the
        // left-normalized alignment can change a codon there with no substitution (decision 0026)
        let gapped = sim.v.sequence_gapped.as_bytes();
        let gi: Vec<usize> = (0..gapped.len()).filter(|&g| gapped[g] != b'.').collect();
        let near_deletion: Vec<usize> = ls(&rows, "v_deletions", r)
            .iter()
            .flat_map(|d| {
                let pos: usize = d[4..d.find(':').unwrap()].parse().unwrap();
                let u = gi.iter().position(|&g| g + 1 == pos).unwrap();
                (u..u + d.rsplit(':').next().unwrap().len()).map(|x| gi[x] / 3 + 1).collect::<Vec<_>>()
            })
            .collect();
        for m in ls(&rows, "v_mutations_aa", r) {
            let c: usize = m[4..m.find(':').unwrap()].parse().unwrap();
            if !want[0].iter().any(|x| (3 * c - 2..=3 * c).contains(&x[4..x.find(':').unwrap()].parse().unwrap())) {
                assert!(near_deletion.iter().any(|&n| n.abs_diff(c) <= 1), "{m}");
                aa_without_nt += 1;
            }
            v_aa += 1;
        }
        // the observable indel, in germline numbering of its segment
        let mut got: Vec<String> = Vec::new();
        for seg in ["v", "d", "j", "c"] {
            if !col(&rows, &format!("{seg}_insertions")).is_null(r) {
                got.extend(ls(&rows, &format!("{seg}_insertions"), r).into_iter().map(|x| format!("{seg} {x}")));
                got.extend(ls(&rows, &format!("{seg}_deletions"), r).into_iter().map(|x| format!("{seg} {x}")));
            }
        }
        let want: Vec<String> = sim
            .indels
            .iter()
            .filter(|_| !ls(&rows, "indels", r).is_empty())
            .map(|x| {
                let ins = x.kind == IndelKind::Ins;
                let after = x.read_pos - t5 + if ins { x.len } else { 0 };
                let gapped = |u: usize| sim.v.sequence_gapped.bytes().enumerate().filter(|c| c.1 != b'.').nth(u).unwrap().0;
                // leftmost equivalent placement whose left neighbour is in the read and segment
                let is_v = matches!(sim.origin[after], V(_));
                let first = sim.origin.iter().position(|o| if is_v { matches!(o, V(_)) } else { matches!(o, J(_)) });
                let shift = x.read_pos - x.equiv.0.max(t5 + first.unwrap() + 1);
                shifted += usize::from(shift > 0);
                let germ = if is_v { sim.v.sequence() } else { sim.j.sequence() };
                let mut bases = x.bases.clone();
                let u = if ins { x.germline_pos - 1 - shift } else { x.germline_pos - shift };
                if ins {
                    let mut b = bases.into_bytes();
                    b.rotate_right(shift % x.len);
                    bases = String::from_utf8(b).unwrap();
                } else {
                    bases = germ[u..u + x.len].to_string();
                }
                let (seg, pos) = match (is_v, sim.origin[after]) {
                    (true, _) => ("v", format!("IMGT{}", gapped(u) + 1)),
                    (false, J(_)) => ("j", (u + 1).to_string()),
                    (_, o) => panic!("indel before {o:?}"),
                };
                format!("{seg} {pos}:{}:{bases}", if ins { "ins" } else { "del" })
            })
            .collect();
        assert_eq!(got, want, "row {r}");
        seg_indels += got.len();
    }
    let counts = (inverted, v_aa, seg_indels, shifted);
    eprintln!("aa changes beside a deletion without an nt change: {aa_without_nt}");
    assert!(inverted > 20 && v_aa > 1000 && seg_indels > 300 && shifted > 50, "{counts:?}");
}

/// A homopolymer indel generated right of its leftmost equivalent placement is reported at
/// the run's leftmost germline position (VCF-style left normalization).
#[test]
fn segment_indels_left_normalized() {
    use abstar_generate::{IndelKind, Origin::V};
    let b = igh();
    let p = GenerateParams { indel_rate: 1.0, frameshift_indels: true, ..params(3000, 23) };
    let rows = one(&generate_in_memory(&b, &p).unwrap());
    let mut tested = 0;
    for r in 0..rows.num_rows() {
        let sim = simulated(&b, &p, r);
        let Some(x) = sim.indels.first() else { continue };
        let c = x.bases.as_bytes()[0];
        let ins = x.kind == IndelKind::Ins;
        let after = x.read_pos + if ins { x.len } else { 0 };
        if x.bases.bytes().any(|y| y != c) || !matches!(sim.origin[after], V(_)) {
            continue;
        }
        // walk left through the run of `c` on the germline V, keeping a V base on the left
        let g = sim.v.sequence().into_bytes();
        let mut u = x.germline_pos; // first deleted base, or the base after the insertion
        while u > 1 && g[u - 1] == c {
            u -= 1;
        }
        if u == x.germline_pos {
            continue; // already leftmost
        }
        let gapped: Vec<usize> = sim.v.sequence_gapped.bytes().enumerate().filter(|y| y.1 != b'.').map(|y| y.0).collect();
        let bases = x.bases.clone();
        let want = if ins {
            format!("IMGT{}:ins:{bases}", gapped[u - 1] + 1)
        } else {
            format!("IMGT{}:del:{bases}", gapped[u] + 1)
        };
        let got = ls(&rows, if ins { "v_insertions" } else { "v_deletions" }, r);
        assert_eq!(got, [want], "row {r}: generated at germline {}", x.germline_pos);
        tested += 1;
    }
    assert!(tested > 50, "{tested}");
}

/// The aa denominator drops the codons holding a deleted base at the left-normalized placement
/// (the one the `v_deletions` strings give), not at the generated placement.
#[test]
fn aa_denominator_follows_left_normalized_deletions() {
    use abstar_generate::IndelKind;
    let b = igh();
    let p = GenerateParams { indel_rate: 1.0, mutation: "uniform:0.05-0.2".parse().unwrap(), ..params(4000, 25) };
    let rows = one(&generate_in_memory(&b, &p).unwrap());
    let (mut tested, mut moved) = (0, 0);
    for r in 0..rows.num_rows() {
        let dels = ls(&rows, "v_deletions", r);
        let (Some(freq), Some(count)) = (f(&rows, "v_mutation_frequency_aa", r), i(&rows, "v_mutation_count_aa", r)) else { continue };
        if dels.is_empty() {
            continue;
        }
        let sim = simulated(&b, &p, r);
        let gapped = sim.v.sequence_gapped.as_bytes();
        let gi: Vec<usize> = (0..gapped.len()).filter(|&g| gapped[g] != b'.').collect();
        // IMGT positions (1-based) of the normalized deleted bases
        let mut deleted = Vec::new();
        for d in &dels {
            let (pos, bases) = (d[4..d.find(':').unwrap()].parse::<usize>().unwrap(), d.rsplit(':').next().unwrap().len());
            let u = gi.iter().position(|&g| g + 1 == pos).unwrap();
            deleted.extend((u..u + bases).map(|x| gi[x] + 1));
        }
        let (lo, hi) = (gi[i(&rows, "v_germline_start", r).unwrap() as usize - 1] + 1, gi[i(&rows, "v_germline_end", r).unwrap() as usize - 1] + 1);
        let den = (1..=gapped.len() / 3)
            .filter(|n| 3 * n - 2 >= lo && 3 * n <= hi && !gapped[3 * n - 3..3 * n].contains(&b'.'))
            .filter(|n| !deleted.iter().any(|&x| x.div_ceil(3) == *n))
            .count();
        assert!(den > 0 && (freq - count as f64 / den as f64).abs() < 1e-12, "row {r}: {freq} {count} {den}");
        // the generated placement falls in other codons than the normalized one
        if let Some(x) = sim.indels.iter().find(|x| x.kind == IndelKind::Del) {
            let generated: Vec<usize> = (x.germline_pos..x.germline_pos + x.len).map(|u| gi[u] / 3 + 1).collect();
            let normalized: Vec<usize> = deleted.iter().map(|x| x.div_ceil(3)).collect();
            moved += usize::from(generated != normalized);
        }
        tested += 1;
    }
    assert!(tested > 100 && moved > 5, "{tested} {moved}");
}

/// Rebuilding the V aa mutation list from the read and the normalized deletion placement gives truth's list.
#[test]
fn v_aa_mutations_from_normalized_alignment() {
    // standard genetic code, TCAG order
    const AA: &[u8; 64] = b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
    let tr = |c: &[u8]| {
        let k = |x: u8| "TCAG".find(x.to_ascii_uppercase() as char).unwrap();
        AA[16 * k(c[0]) + 4 * k(c[1]) + k(c[2])] as char
    };
    let b = igh();
    let p = GenerateParams { indel_rate: 1.0, mutation: "uniform:0.05-0.25".parse().unwrap(), ..params(4000, 27) };
    let rows = one(&generate_in_memory(&b, &p).unwrap());
    let mut tested = 0;
    for r in 0..rows.num_rows() {
        let dels = ls(&rows, "v_deletions", r);
        if dels.is_empty() || col(&rows, "v_mutations_aa").is_null(r) {
            continue;
        }
        let sim = simulated(&b, &p, r);
        let gapped = sim.v.sequence_gapped.as_bytes();
        let gi: Vec<usize> = (0..gapped.len()).filter(|&g| gapped[g] != b'.').collect();
        let g = sim.v.sequence().into_bytes();
        let mut deleted = vec![false; g.len()];
        for d in &dels {
            let pos: usize = d[4..d.find(':').unwrap()].parse().unwrap();
            let u = gi.iter().position(|&x| x + 1 == pos).unwrap();
            (u..u + d.rsplit(':').next().unwrap().len()).for_each(|x| deleted[x] = true);
        }
        // the segment's read bases, in order, laid onto the germline positions not deleted
        let mut bases = sim.origin.iter().enumerate().filter(|(_, o)| matches!(o, abstar_generate::Origin::V(_))).map(|(k, _)| sim.read[k]);
        let first = sim.origin.iter().find_map(|o| if let abstar_generate::Origin::V(x) = o { gi.iter().position(|&y| y == *x) } else { None }).unwrap();
        let mut aligned: Vec<Option<u8>> = vec![None; g.len()];
        for u in (first..g.len()).filter(|&u| !deleted[u]) {
            aligned[u] = bases.next();
        }
        let mut want = Vec::new();
        for n in 1..=gapped.len() / 3 {
            if gapped[3 * n - 3..3 * n].contains(&b'.') {
                continue;
            }
            let u = gi.iter().position(|&x| x == 3 * n - 3).unwrap();
            let Some(read3) = (u..u + 3).map(|x| aligned[x]).collect::<Option<Vec<u8>>>() else { continue };
            let (from, to) = (tr(&g[u..u + 3]), tr(&read3));
            if from != to {
                want.push(format!("IMGT{n}:{from}>{to}"));
            }
        }
        assert_eq!(ls(&rows, "v_mutations_aa", r), want, "row {r}");
        tested += 1;
    }
    assert!(tested > 150, "{tested}");
}

/// The constant region's germline span is truth too, so its mutation fields can be graded.
#[test]
fn c_germline_coordinates() {
    let b = igh();
    let p = busy(300, 26);
    let rows = one(&generate_in_memory(&b, &p).unwrap());
    let mut with_c = 0;
    for r in 0..rows.num_rows() {
        let (a, z) = (i(&rows, "c_germline_start", r), i(&rows, "c_germline_end", r));
        assert_eq!(a.is_some(), i(&rows, "c_sequence_start", r).is_some());
        if let (Some(a), Some(z), Some(rs), Some(re)) = (a, z, i(&rows, "c_sequence_start", r), i(&rows, "c_sequence_end", r)) {
            assert_eq!((a, z - a), (1, re - rs), "row {r}");
            with_c += 1;
        }
    }
    assert!(with_c > 100, "{with_c}");
}

/// A read cut inside CDR3, before the J's Phe118 codon: CDR3 holds the covered bases, the
/// masks cover it, and the junction and FWR4 are null.
#[test]
fn partial_cdr3_when_truncated() {
    let b = igk_fixed(3);
    // V (306 nt after the 5' cut) + AAA + the first 5 J bases; Phe118 is J 8-10
    let p = GenerateParams {
        truncate_5p: "fixed:5".parse().unwrap(),
        truncate_3p: "fixed:33".parse().unwrap(),
        ..params(1, 9)
    };
    let rows = one(&generate_in_memory(&b, &p).unwrap());
    let read = s(&rows, "sequence", 0).unwrap();
    assert_eq!(read.len(), 314);
    assert_eq!((s(&rows, "junction", 0), s(&rows, "fwr4", 0)), (None, None));
    assert_eq!((i(&rows, "cdr3_start", 0), i(&rows, "cdr3_end", 0)), (Some(287), Some(314)));
    assert_eq!(s(&rows, "cdr3", 0).as_deref(), Some(&read[286..]));
    // codons start at read 2 (mod 3): 287-289 … 311-313 are whole, 314 is a lone base
    assert_eq!(s(&rows, "cdr3_aa", 0).map(|x| x.len()), Some(9));
    assert_eq!(i(&rows, "cdr3_length", 0), Some(9));
    assert_eq!(s(&rows, "sequence_vdj", 0), Some(read.clone()));
    let rep = |c: &str, n: usize| c.repeat(n);
    let cdr = [rep("0", 70), rep("1", 24), rep("0", 51), rep("2", 24), rep("0", 117), rep("3", 28)];
    assert_eq!(s(&rows, "cdr_mask", 0), Some(cdr.concat()));
    let cdr_aa = [rep("0", 23), rep("1", 8), rep("0", 17), rep("2", 8), rep("0", 39), rep("3", 9)];
    assert_eq!(s(&rows, "cdr_mask_aa", 0), Some(cdr_aa.concat()));
    assert_eq!(s(&rows, "sequence_vdj_aa", 0).map(|x| x.len()), Some(104));
}
