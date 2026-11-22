//! Output schema, generated from `fields.toml` by `build.rs`.
#[allow(unused_imports)]
use arrow_array::builder::{BooleanBuilder, Float64Builder, Int64Builder, ListBuilder, StringBuilder};
use arrow_array::{ArrayRef, RecordBatch};
use arrow_schema::{DataType, Field, Schema, SchemaRef};
use std::collections::HashMap;
use std::sync::Arc;

include!(concat!(env!("OUT_DIR"), "/schema_generated.rs"));
