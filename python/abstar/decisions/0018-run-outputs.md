# 0018: Run outputs: stems, conflicts, overwrite, run.json, panics as error rows

## Status

Accepted (2026-11-20)

## Context

A run turns FASTA/FASTQ inputs into TSV and/or Parquet files plus a record of how the run was made. Users rerun into the same directory, pass several inputs, and feed millions of records, so naming, collisions, and per-record failures need fixed rules. Every input record must yield exactly one output row, in input order (0016).

## Decision

**Output stems.** The stem is the input file name minus a trailing `.gz`, then minus one of `.fasta`, `.fa`, `.fna`, `.fastq`, `.fq`. Standard input is `stdin`; in-memory records are `records`. Outputs are `<output_dir>/<stem>.tsv` and `<stem>.parquet`.

**Reading input.** Compression is detected by content, not by extension: a file (or stdin) starting with the gzip magic bytes is decompressed, including multi-member gzip. bzip2, xz and zstd are not supported; a file ending in `.bz2`, `.xz` or `.zst` fails with `INPUT_UNREADABLE` ("compression not supported: use gzip"), and directory expansion includes such files so the run reports them instead of silently skipping them. A leading UTF-8 byte-order mark and leading blank lines (LF or CRLF) are skipped. ASCII spaces and tabs inside sequence lines (such as trailing whitespace) are removed; `sequence` is otherwise exactly as read, not uppercased. A record whose ID or sequence bytes are not valid UTF-8 is not annotated: it becomes an `invalid_input` row with `failure_reason` `invalid_encoding`, and its `sequence_id` and `sequence` are written lossily decoded (invalid bytes become U+FFFD), because Arrow and TSV strings must be UTF-8. This is the one exception to byte-for-byte IDs (0016); it keeps two IDs that differ only in invalid bytes from silently merging into one annotated ID.

**Conflicts.** If two inputs share a stem, the run fails with `OUTPUT_NAME_CONFLICT`, naming the colliding sources, before any file is created. Missing input files fail with `INPUT_NOT_FOUND` at the same point.

**Overwrite.** Existing outputs with the same names are overwritten without prompting; reruns are the normal workflow. Other files in the directory are left alone.

**run.json.** Written last to `<output_dir>/run.json`, so its presence marks a finished run. A run deletes any existing `run.json` before it creates or truncates an output file, so a failed re-run leaves no `run.json` and cannot look finished. Each input is opened before its output files are created, so an input that cannot be opened or parsed fails without truncating that input's previous outputs. Keys: `abstar_version`, `airr_version`, `started` and `finished` (RFC 3339 UTC), `parameters` (`germline_db`, `db_path`, `formats`, `threads`, `strict`, `loci`), `database` (`name`, `alleles_sha256`), `inputs` (`source`, `outputs`, `counts` per input) and the total `counts` (`annotated`, `unassigned`, `invalid_input`, `error`). `source` is the path as given, `-` for stdin, `records` for in-memory input. `source_file` in each row carries the same value, and is null for in-memory records.

**Panics become error rows.** Each record is annotated under `catch_unwind`. A panic yields a row with status `error` and `failure_reason` of `internal_error: <message>`; the run continues. Under `strict`, the first such row aborts the run with `STRICT_RECORD_FAILURE`. The default panic message is suppressed only while a record is annotated, so the failure appears once, in the row.

**Rows stay single lines.** A tab, CR or LF in a sequence ID or sequence (possible only from in-memory input or malformed files) is replaced by a space, because TSV cannot hold it. This keeps one row per record instead of aborting the run.

**Other.** The `loci` option must be a subset of the database's loci (`INVALID_ARGUMENT` otherwise) and is otherwise only recorded. A closed output pipe ends TSV streaming without an error. An empty input yields zero rows, not an error.

## Alternatives considered

- Refuse to overwrite existing outputs: makes the common rerun fail; users can choose a new directory.
- Suffix colliding stems (`s1`, `s1_2`): names would depend on input order and hide the problem.
- Abort the run on a panic: one bad record would lose a long run; strict mode gives that behavior on request.
- Reject records containing tabs: breaks the one-row-per-record guarantee.
- Decode non-UTF-8 bytes lossily and annotate anyway: distinct IDs could merge without notice.
- Fail the whole input on a non-UTF-8 record, or on a leading BOM or blank line: common FASTA variants would abort long runs.
- Detect compression by extension: a gzip file without `.gz`, or stdin, would not be read.

## Consequences

- Output names are predictable from input names, and conflicts surface before work is done.
- `run.json` is enough to see which database and parameters produced a directory.
- A directory without `run.json` holds an unfinished or failed run; outputs from an input that failed mid-stream are partial.
- A custom panic hook is installed once per process; it forwards all other panics to the previous hook.
