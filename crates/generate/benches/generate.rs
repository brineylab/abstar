//! Speed target (spec: *Validation of the generator itself*, item 6): 1M `human-igh` reads,
//! default (`naive`) parameters, FASTA plus truth TSV, in ≤ 10 s on the M-series laptop.
use abstar_generate::{bind, generate_to_files, load_model, preset, GenerateMeta, GenerateParams};
use abstar_table::writers::Format;
use criterion::{criterion_group, criterion_main, Criterion};
use std::path::Path;

fn generate_1m(c: &mut Criterion) {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../python/abstar");
    let (db_path, model_path) = (root.join("germline_dbs/human"), root.join("models/human-igh/model.toml"));
    let b = bind(load_model(&model_path).unwrap(), &abstar_germline::load(&db_path).unwrap()).unwrap();
    let p = GenerateParams { n: 1_000_000, ..preset("naive").unwrap() };
    let meta = GenerateMeta {
        abstar_version: "bench".into(),
        germline_db: "human".into(),
        db_path,
        model_name: "human-igh".into(),
        model_path,
        preset: Some("naive".into()),
    };
    let dir = tempfile::tempdir().unwrap();
    let mut g = c.benchmark_group("generate");
    g.sample_size(10);
    g.bench_function("human_igh_1m_truth_tsv", |x| {
        x.iter(|| generate_to_files(&b, &p, dir.path(), "bench", &[Format::Airr], &meta).unwrap())
    });
    g.finish();
}

criterion_group!(benches, generate_1m);
criterion_main!(benches);
