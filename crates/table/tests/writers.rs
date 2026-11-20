use abstar_table::writers::{ParquetSink, TsvWriter};
use arrow_array::builder::{Int64Builder, ListBuilder};
use arrow_array::{ArrayRef, RecordBatch};
use arrow_schema::{DataType, Field, Schema, SchemaRef};
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use std::sync::Arc;

fn schema() -> SchemaRef {
    Arc::new(Schema::new(vec![Field::new(
        "xs",
        DataType::List(Arc::new(Field::new("item", DataType::Int64, true))),
        true,
    )]))
}

fn batch() -> RecordBatch {
    let mut b = ListBuilder::new(Int64Builder::new());
    b.values().append_value(3);
    b.values().append_value(4);
    b.append(true);
    b.append(false);
    RecordBatch::try_new(schema(), vec![Arc::new(b.finish()) as ArrayRef]).unwrap()
}

#[test]
fn list_int64_round_trip() {
    let mut w = TsvWriter::new(Vec::new(), &schema()).unwrap();
    w.write(&batch()).unwrap();
    let out = String::from_utf8(w.finish().unwrap()).unwrap();
    assert_eq!(out, "xs\n3,4\n\n");

    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("t.parquet");
    let mut s = ParquetSink::create(&path, schema()).unwrap();
    s.write(&batch()).unwrap();
    s.finish().unwrap();
    let r = ParquetRecordBatchReaderBuilder::try_new(std::fs::File::open(&path).unwrap()).unwrap();
    assert_eq!(r.build().unwrap().next().unwrap().unwrap(), batch());
}

#[test]
fn header_follows_given_schema() {
    let s: SchemaRef = Arc::new(Schema::new(vec![
        Field::new("b", DataType::Utf8, true),
        Field::new("a", DataType::Int64, true),
    ]));
    let out = String::from_utf8(TsvWriter::new(Vec::new(), &s).unwrap().finish().unwrap()).unwrap();
    assert_eq!(out, "b\ta\n");
}
