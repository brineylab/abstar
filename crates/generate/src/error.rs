use std::path::PathBuf;

#[derive(Debug, thiserror::Error)]
pub enum GenerateError {
    #[error("recombination model not found: {}", .0.display())]
    ModelNotFound(PathBuf),
    #[error("invalid recombination model: {0}")]
    ModelInvalid(String),
    #[error("invalid argument: {0}")]
    InvalidArgument(String),
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    #[error("internal error: {0}")]
    Internal(String),
}

impl GenerateError {
    /// The message without the Display prefix, for callers that show `code()` beside it.
    pub fn message(&self) -> String {
        match self {
            Self::ModelNotFound(p) => format!("{} does not exist", p.display()),
            Self::ModelInvalid(m) | Self::InvalidArgument(m) | Self::Internal(m) => m.clone(),
            Self::Io(e) => e.to_string(),
        }
    }

    pub fn code(&self) -> &'static str {
        match self {
            Self::ModelNotFound(_) => "MODEL_NOT_FOUND",
            Self::ModelInvalid(_) => "MODEL_INVALID",
            Self::InvalidArgument(_) => "INVALID_ARGUMENT",
            Self::Io(_) => "IO_ERROR",
            Self::Internal(_) => "INTERNAL",
        }
    }
}
