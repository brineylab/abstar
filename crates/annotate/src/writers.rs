//! TSV and Parquet output writers bound to the annotation schema.
use crate::schema::arrow_schema;
use arrow_array::RecordBatch;
use parquet::errors::ParquetError;
use std::io::{self, Write};
use std::path::Path;

pub use abstar_table::writers::Format;

pub struct TsvWriter<W: Write>(abstar_table::writers::TsvWriter<W>);

impl<W: Write> TsvWriter<W> {
    pub fn new(w: W) -> io::Result<Self> {
        abstar_table::writers::TsvWriter::new(w, &arrow_schema()).map(Self)
    }
    pub fn write(&mut self, b: &RecordBatch) -> io::Result<()> {
        self.0.write(b)
    }
    pub fn finish(self) -> io::Result<W> {
        self.0.finish()
    }
}

pub struct ParquetSink(abstar_table::writers::ParquetSink);

impl ParquetSink {
    pub fn create(path: &Path) -> Result<Self, ParquetError> {
        abstar_table::writers::ParquetSink::create(path, arrow_schema()).map(Self)
    }
    pub fn write(&mut self, b: &RecordBatch) -> Result<(), ParquetError> {
        self.0.write(b)
    }
    pub fn finish(self) -> Result<(), ParquetError> {
        self.0.finish()
    }
}
