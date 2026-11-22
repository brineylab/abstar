//! Placeholder annotator, kept for the engine seam tests.
use crate::record::{failed, invalid_input};
use crate::schema::Record;
use abstar_germline::Database;

pub fn annotate_stub(id: &str, seq: &str, _db: &Database) -> Record {
    match invalid_input(seq) {
        Some(reason) => failed(id, seq, "invalid_input", reason),
        None => failed(id, seq, "unassigned", "engine_not_implemented"),
    }
}
