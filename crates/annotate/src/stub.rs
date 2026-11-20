//! Placeholder annotator; the annotator project replaces it.
use crate::schema::Record;
use abstar_germline::Database;

pub fn annotate_stub(id: &str, seq: &str, _db: &Database) -> Record {
    let (status, reason) = if seq.is_empty() {
        ("invalid_input", "empty_sequence")
    } else if !seq.bytes().all(|b| b"ACGTNRYSWKMBDHV".contains(&b.to_ascii_uppercase())) {
        ("invalid_input", "invalid_characters")
    } else {
        ("unassigned", "engine_not_implemented")
    };
    Record {
        sequence_id: id.into(),
        sequence: seq.into(),
        annotation_status: status.into(),
        failure_reason: Some(reason.into()),
        ..Default::default()
    }
}
