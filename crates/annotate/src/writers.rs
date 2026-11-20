//! TSV and Parquet output writers.
use crate::schema::arrow_schema;
use arrow_array::cast::AsArray;
use arrow_array::{Array, RecordBatch};
use arrow_schema::DataType;
use parquet::arrow::ArrowWriter;
use parquet::basic::{Compression, ZstdLevel};
use parquet::errors::ParquetError;
use parquet::file::properties::WriterProperties;
use std::io::{self, Write};
use std::path::Path;
use std::str::FromStr;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Format {
    Airr,
    Parquet,
}

impl FromStr for Format {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, String> {
        match s {
            "airr" => Ok(Format::Airr),
            "parquet" => Ok(Format::Parquet),
            _ => Err(format!("unknown format `{s}` (expected `airr` or `parquet`)")),
        }
    }
}

pub struct TsvWriter<W: Write> {
    w: W,
}

fn cell(out: &mut String, name: &str, s: &str) {
    assert!(
        !s.contains(['\t', '\n', '\r']),
        "field `{name}` contains a tab or newline"
    );
    out.push_str(s);
}

impl<W: Write> TsvWriter<W> {
    pub fn new(mut w: W) -> io::Result<Self> {
        let schema = arrow_schema();
        let names: Vec<&str> = schema.fields().iter().map(|f| f.name().as_str()).collect();
        writeln!(w, "{}", names.join("\t"))?;
        Ok(Self { w })
    }

    pub fn write(&mut self, b: &RecordBatch) -> io::Result<()> {
        let schema = b.schema();
        let mut line = String::new();
        for row in 0..b.num_rows() {
            line.clear();
            for (i, f) in schema.fields().iter().enumerate() {
                if i > 0 {
                    line.push('\t');
                }
                let col = b.column(i);
                if col.is_null(row) {
                    continue;
                }
                match f.data_type() {
                    DataType::Utf8 => cell(&mut line, f.name(), col.as_string::<i32>().value(row)),
                    DataType::Boolean => line.push(if col.as_boolean().value(row) { 'T' } else { 'F' }),
                    DataType::Int64 => line.push_str(&col.as_primitive::<arrow_array::types::Int64Type>().value(row).to_string()),
                    DataType::Float64 => line.push_str(&col.as_primitive::<arrow_array::types::Float64Type>().value(row).to_string()),
                    DataType::List(_) => {
                        let items = col.as_list::<i32>().value(row);
                        for (j, s) in items.as_string::<i32>().iter().enumerate() {
                            if j > 0 {
                                line.push(',');
                            }
                            cell(&mut line, f.name(), s.unwrap_or(""));
                        }
                    }
                    t => panic!("unsupported TSV type {t:?} for field `{}`", f.name()),
                }
            }
            line.push('\n');
            self.w.write_all(line.as_bytes())?;
        }
        Ok(())
    }

    pub fn finish(mut self) -> io::Result<W> {
        self.w.flush()?;
        Ok(self.w)
    }
}

pub struct ParquetSink {
    w: ArrowWriter<std::fs::File>,
}

impl ParquetSink {
    pub fn create(path: &Path) -> Result<Self, ParquetError> {
        let schema = arrow_schema();
        let props = WriterProperties::builder()
            .set_compression(Compression::ZSTD(ZstdLevel::default()))
            .set_key_value_metadata(Some(
                schema
                    .metadata()
                    .iter()
                    .map(|(k, v)| parquet::file::metadata::KeyValue::new(k.clone(), v.clone()))
                    .collect(),
            ))
            .build();
        let file = std::fs::File::create(path).map_err(|e| ParquetError::External(Box::new(e)))?;
        Ok(Self { w: ArrowWriter::try_new(file, schema, Some(props))? })
    }

    pub fn write(&mut self, b: &RecordBatch) -> Result<(), ParquetError> {
        self.w.write(b)
    }

    pub fn finish(self) -> Result<(), ParquetError> {
        self.w.close().map(|_| ())
    }
}
