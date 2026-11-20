mod error;
mod model;

pub use error::GenerateError;
pub use model::{load_model, Dist, Model, Source, TrimSet, Trims};
