use std::path::PathBuf;

#[derive(Debug, thiserror::Error)]
pub enum GermlineError {
    #[error("germline database not found: {}", .0.display())]
    NotFound(PathBuf),
    #[error("invalid germline database: {0}")]
    Invalid(String),
    #[error("germline database already exists: {}", .0.display())]
    Exists(PathBuf),
    #[error("invalid recipe: {0}")]
    RecipeInvalid(String),
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
}

impl GermlineError {
    pub fn code(&self) -> &'static str {
        match self {
            Self::NotFound(_) => "DATABASE_NOT_FOUND",
            Self::Invalid(_) => "DATABASE_INVALID",
            Self::Exists(_) => "DATABASE_EXISTS",
            Self::RecipeInvalid(_) => "RECIPE_INVALID",
            Self::Io(_) => "IO_ERROR",
        }
    }
}
