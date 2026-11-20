use abstar_annotate::schema::{arrow_schema, to_record_batch, Record, FIELDS_TOML};
use abstar_annotate::writers::{ParquetSink, TsvWriter};
use arrow_schema::DataType;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

fn toml_names() -> Vec<String> {
    FIELDS_TOML
        .lines()
        .filter_map(|l| l.strip_prefix("name = \""))
        .map(|l| l.trim_end_matches('"').to_string())
        .collect()
}

fn rec(id: &str) -> Record {
    Record {
        sequence_id: id.into(),
        sequence: "ACGT".into(),
        annotation_status: "ok".into(),
        ..Default::default()
    }
}

fn tsv(recs: &[Record]) -> String {
    let mut w = TsvWriter::new(Vec::new()).unwrap();
    if !recs.is_empty() {
        w.write(&to_record_batch(recs)).unwrap();
    }
    String::from_utf8(w.finish().unwrap()).unwrap()
}

#[test]
fn schema_matches_toml() {
    let s = arrow_schema();
    let names: Vec<_> = s.fields().iter().map(|f| f.name().clone()).collect();
    assert_eq!(names, toml_names());
    assert!(!s.field_with_name("annotation_status").unwrap().is_nullable());
    assert!(!s.field_with_name("sequence_id").unwrap().is_nullable());
    assert!(!s.field_with_name("sequence").unwrap().is_nullable());
    assert_eq!(s.fields().iter().filter(|f| !f.is_nullable()).count(), 3);
    let t = s.field_with_name("v_call_ties").unwrap().data_type().clone();
    match t {
        DataType::List(f) => assert_eq!(f.data_type(), &DataType::Utf8),
        o => panic!("{o:?}"),
    }
    assert_eq!(s.metadata()["airr_version"], "2.0");
    assert_eq!(s.metadata()["abstar_version"], env!("CARGO_PKG_VERSION"));
}

#[test]
fn tsv_encoding() {
    let mut r = rec("a");
    r.rev_comp = Some(true);
    r.v_call_ties = Some(vec!["IGHV1-2*02".into(), "IGHV1-2*04".into()]);
    let out = tsv(&[r]);
    let mut lines = out.lines();
    let h: Vec<_> = lines.next().unwrap().split('\t').collect();
    let row: Vec<_> = lines.next().unwrap().split('\t').collect();
    assert_eq!(h.len(), row.len());
    let cell = |n: &str| row[h.iter().position(|x| *x == n).unwrap()];
    assert_eq!(cell("rev_comp"), "T");
    assert_eq!(cell("productive"), "");
    assert_eq!(cell("v_call_ties"), "IGHV1-2*02,IGHV1-2*04");
    let mut r = rec("b");
    r.productive = Some(false);
    assert!(tsv(&[r]).lines().nth(1).unwrap().split('\t').any(|c| c == "F"));
}

#[test]
fn tsv_header_only_for_zero_rows() {
    assert_eq!(tsv(&[]), format!("{}\n", toml_names().join("\t")));
}

#[test]
#[should_panic(expected = "sequence")]
fn tsv_rejects_tab() {
    let mut r = rec("a");
    r.sequence = "AC\tGT".into();
    tsv(&[r]);
}

fn parquet_round_trip_of(tag: &str, recs: &[Record]) -> arrow_array::RecordBatch {
    let dir = std::env::temp_dir().join(format!("abstar_pq_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join(format!("{tag}.parquet"));
    let batch = to_record_batch(recs);
    let mut sink = ParquetSink::create(&path).unwrap();
    sink.write(&batch).unwrap();
    sink.finish().unwrap();
    let b = ParquetRecordBatchReaderBuilder::try_new(std::fs::File::open(&path).unwrap()).unwrap();
    let schema = b.schema().clone();
    let md = schema.metadata().clone();
    assert_eq!(md["airr_version"], "2.0");
    assert_eq!(md["abstar_version"], env!("CARGO_PKG_VERSION"));
    let got = b.build().unwrap().next().unwrap().unwrap();
    // The reader's per-batch schema drops file metadata; it is checked above on the builder schema.
    let got = arrow_array::RecordBatch::try_new(schema, got.columns().to_vec()).unwrap();
    assert_eq!(got, batch);
    got
}

#[test]
fn parquet_round_trip() {
    let mut r = rec("a");
    r.rev_comp = Some(true);
    r.v_call_ties = Some(vec!["x".into(), "y".into()]);
    parquet_round_trip_of("rt", &[r, rec("b")]);
}

#[test]
fn id_strings_preserved() {
    let recs = [rec("10E8"), rec("0012")];
    let out = tsv(&recs);
    let ids: Vec<_> = out.lines().skip(1).map(|l| l.split('\t').next().unwrap()).collect();
    assert_eq!(ids, ["10E8", "0012"]);
    let b = parquet_round_trip_of("ids", &recs);
    let col = b.column(0).as_any().downcast_ref::<arrow_array::StringArray>().unwrap();
    assert_eq!((col.value(0), col.value(1)), ("10E8", "0012"));
}
