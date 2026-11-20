//! Germline database types, on-disk format, and validation.
mod error;
mod io;
mod model;
mod validate;

pub use error::GermlineError;
pub use io::{load, write};
pub use model::*;
pub use validate::validate;
