//! Pipeline throughput (parse, stub-annotate, write TSV) over the first 10k corpus records.
use abstar_annotate::input::Input;
use abstar_annotate::pipeline::{run_to_writer, RunOptions};
use criterion::{criterion_group, criterion_main, Criterion};
use std::io::{sink, Read};

fn pipeline_10k(c: &mut Criterion) {
    let root = concat!(env!("CARGO_MANIFEST_DIR"), "/../..");
    let gz = std::fs::File::open(format!("{root}/legacy/test_data/bcr_corpus/sequences.fasta.gz")).unwrap();
    let mut text = String::new();
    flate2::read::MultiGzDecoder::new(gz).read_to_string(&mut text).unwrap();
    // 2 lines per record (single-line FASTA).
    let head: String = text.lines().take(20_000).flat_map(|l| [l, "\n"]).collect();
    let path = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join("pipeline_10k.fasta");
    std::fs::write(&path, head).unwrap();

    let db = abstar_germline::load(&std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../python/abstar/germline_dbs/human")).unwrap();
    let inputs = [Input::File(path)];
    let opts = RunOptions::default();
    let mut g = c.benchmark_group("pipeline");
    g.sample_size(20);
    g.bench_function("pipeline_10k", |b| b.iter(|| run_to_writer(&inputs, &db, &opts, sink()).unwrap()));
    g.finish();
}

criterion_group!(benches, pipeline_10k);
criterion_main!(benches);
