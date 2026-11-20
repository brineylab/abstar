//! Compiles generated code for every supported type, then round-trips a row.
#![cfg(feature = "writers")]
use abstar_table::writers::TsvWriter;
#[allow(unused_imports)]
use arrow_array::builder::{BooleanBuilder, Float64Builder, Int64Builder, ListBuilder, StringBuilder};
use arrow_array::{ArrayRef, RecordBatch};
use arrow_schema::{DataType, Field, Schema, SchemaRef};
use std::collections::HashMap;
use std::sync::Arc;

include!(concat!(env!("OUT_DIR"), "/fixture_generated.rs"));

#[test]
fn every_type_compiles_and_round_trips() {
    let rows = [
        Fixture { id: "a".into(), s: Some("x".into()), li: Some(vec![3, 4]), ls: Some(vec!["p".into(), "q".into()]), b: Some(true), i: Some(7), f: Some(1.5) },
        Fixture { id: "b".into(), ..Default::default() },
    ];
    let batch = to_record_batch(&rows);
    let mut w = TsvWriter::new(Vec::new(), &arrow_schema()).unwrap();
    w.write(&batch).unwrap();
    let out = String::from_utf8(w.finish().unwrap()).unwrap();
    assert_eq!(out, "id\ts\tb\ti\tf\tls\tli\na\tx\tT\t7\t1.5\tp,q\t3,4\nb\t\t\t\t\t\t\n");
    assert!(FIXTURE_TOML.contains("list<int64>"));
}
