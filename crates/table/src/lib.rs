//! Schema code generation and TSV/Parquet writers shared by the annotator and the generator.
#[cfg(feature = "codegen")]
pub mod codegen;
#[cfg(feature = "writers")]
pub mod writers;
