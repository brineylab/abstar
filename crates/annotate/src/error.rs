//! Errors for the run path. `code()` is the stable machine-readable code.
use thiserror::Error;

#[derive(Debug, Error)]
pub enum RunError {
    #[error("input not found: {0}")]
    InputNotFound(String),
    #[error("input unreadable: {0}")]
    InputUnreadable(String),
    #[error("invalid argument: {0}")]
    InvalidArgument(String),
    #[error("output name conflict: {0}")]
    OutputNameConflict(String),
    #[error("I/O error: {0}")]
    Io(String),
    #[error("internal error: {0}")]
    Internal(String),
    #[error("strict mode: record failed: {0}")]
    StrictRecordFailure(String),
}

impl RunError {
    pub fn code(&self) -> &'static str {
        match self {
            Self::InputNotFound(_) => "INPUT_NOT_FOUND",
            Self::InputUnreadable(_) => "INPUT_UNREADABLE",
            Self::InvalidArgument(_) => "INVALID_ARGUMENT",
            Self::OutputNameConflict(_) => "OUTPUT_NAME_CONFLICT",
            Self::Io(_) => "IO_ERROR",
            Self::Internal(_) => "INTERNAL",
            Self::StrictRecordFailure(_) => "STRICT_RECORD_FAILURE",
        }
    }
}
