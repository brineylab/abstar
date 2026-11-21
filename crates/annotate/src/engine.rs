//! Run-level annotation engine: built once per run, shared by every worker thread.
use crate::error::RunError;
use crate::schema::Record;
use crate::stub::annotate_stub;
use abstar_germline::Database;

/// Annotates one sequence. Implementations are shared across threads.
pub trait Annotate: Sync {
    fn annotate(&self, id: &str, seq: &str) -> Record;
}

#[derive(Clone, Debug, Default)]
pub struct EngineOptions {
    pub threads: Option<usize>,
}

/// Indexes and germline data built once per run.
pub struct Engine {
    db: Database, // ponytail: only the stub reads this; Task 10 replaces it with indexes
}

impl Engine {
    pub fn new(db: &Database, _opts: &EngineOptions) -> Result<Engine, RunError> {
        Ok(Engine { db: db.clone() })
    }
}

impl Annotate for Engine {
    // Until the real annotator lands (Task 10) the engine delegates to the stub.
    fn annotate(&self, id: &str, seq: &str) -> Record {
        annotate_stub(id, seq, &self.db)
    }
}

/// The placeholder annotator as an `Annotate`.
pub struct Stub<'a>(pub &'a Database);

impl Annotate for Stub<'_> {
    fn annotate(&self, id: &str, seq: &str) -> Record {
        annotate_stub(id, seq, self.0)
    }
}
