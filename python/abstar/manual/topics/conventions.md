Conventions every abstar output follows. Decision records 0004, 0016, 0018 and 0020 hold the reasoning (`abstar manual decisions 0016`).

## Coordinate spaces

- Output coordinates (TSV and Parquet) follow AIRR: 1-based, closed intervals.
- Germline anchor positions (the first base of the V Cys104 codon, or of the J Trp118/Phe118 codon) are 0-based offsets into the ungapped allele sequence. They are stored in the germline database, not in the output.
- Conversion between the two happens at named boundaries only. The `coordinates` attribute of a field in `abstar manual fields` names its space where it matters.

## TSV and Parquet encoding

- TSV: booleans are `T` and `F`, nulls are empty cells, list fields (such as `v_call_ties`) are comma-joined.
- Parquet carries identical content with native types: bool, int, null, and real lists.
- Both writers read the same in-memory record.
- A tab, CR or LF inside a sequence ID or sequence (possible only from in-memory input or malformed files) is replaced by a space, so every record stays one line.
- `sequence` is written exactly as read; it is not uppercased.

## Record statuses

Every input record produces exactly one output row. Records are never dropped, so record conservation is a row-count check. `annotation_status` and `failure_reason` say what happened:

| `annotation_status` | Meaning |
|---|---|
| `annotated` | The record was annotated. |
| `unassigned` | The input is valid but biologically unassignable (for example, no V gene found). `failure_reason` says why. |
| `invalid_input` | The input is not a usable sequence (non-nucleotide characters, or empty). |
| `error` | An internal per-record error. `failure_reason` carries the message, for example `internal_error: <message>`. |

`--strict` turns any `error` row into a fatal run error (exit 3, `STRICT_RECORD_FAILURE`). Fatal run errors (unreadable input, bad database, I/O failure) stop the run with the codes in `abstar manual errors`.

Current state: the foundation ships a stub annotator. Until the real annotator lands, every valid record is reported as `unassigned` with `failure_reason` `engine_not_implemented`. The statuses above are the designed behaviour.

When writing to disk, `run.json` is written last, beside the outputs: its presence marks a finished run. It records parameters, abstar version, database name and content hash, AIRR schema version, and counts by status.

## `sequence_id` rules

- `sequence_id` is the FASTA or FASTQ header's first whitespace-delimited token, byte-for-byte, and is always a string.
- A bare sequence string is one record with ID `"1"`. A list of sequence strings gets IDs `"1"` to `"N"` by position. A list of `(id, sequence)` pairs keeps its IDs.
- Valid sequence characters, case-insensitive: `ACGTN` plus the IUPAC codes `RYSWKMBDHV`.
- Functionality values are `F`, `ORF`, `P`, `NF` and `unknown`.
