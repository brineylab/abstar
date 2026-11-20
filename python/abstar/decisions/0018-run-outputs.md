# 0018: Run outputs: stems, conflicts, overwrite, run.json, panics as error rows

## Status

Accepted (2026-11-20)

## Context

A run turns FASTA/FASTQ inputs into TSV and/or Parquet files plus a record of how the run was made. Users rerun into the same directory, pass several inputs, and feed millions of records, so naming, collisions, and per-record failures need fixed rules. Every input record must yield exactly one output row, in input order (0016).

## Decision

**Output stems.** The stem is the input file name minus a trailing `.gz`, then minus one of `.fasta`, `.fa`, `.fna`, `.fastq`, `.fq`. Standard input is `stdin`; in-memory records are `records`. Outputs are `<output_dir>/<stem>.tsv` and `<stem>.parquet`.

**Conflicts.** If two inputs share a stem, the run fails with `OUTPUT_NAME_CONFLICT`, naming the colliding sources, before any file is created. Missing input files fail with `INPUT_NOT_FOUND` at the same point.

**Overwrite.** Existing outputs with the same names are overwritten without prompting; reruns are the normal workflow. Other files in the directory are left alone.

**run.json.** Written last to `<output_dir>/run.json`, so its presence marks a finished run. Keys: `abstar_version`, `airr_version`, `started` and `finished` (RFC 3339 UTC), `parameters` (`germline_db`, `db_path`, `formats`, `threads`, `strict`, `loci`), `database` (`name`, `alleles_sha256`), `inputs` (`source`, `outputs`, `counts` per input) and the total `counts` (`annotated`, `unassigned`, `invalid_input`, `error`). `source` is the path as given, `-` for stdin, `records` for in-memory input. `source_file` in each row carries the same value, and is null for in-memory records.

**Panics become error rows.** Each record is annotated under `catch_unwind`. A panic yields a row with status `error` and `failure_reason` of `internal_error: <message>`; the run continues. Under `strict`, the first such row aborts the run with `STRICT_RECORD_FAILURE`. The default panic message is suppressed only while a record is annotated, so the failure appears once, in the row.

**Rows stay single lines.** A tab, CR or LF in a sequence ID or sequence (possible only from in-memory input or malformed files) is replaced by a space, because TSV cannot hold it. This keeps one row per record instead of aborting the run.

**Other.** The `loci` option must be a subset of the database's loci (`INVALID_ARGUMENT` otherwise) and is otherwise only recorded. A closed output pipe ends TSV streaming without an error. An empty input yields zero rows, not an error.

## Alternatives considered

- Refuse to overwrite existing outputs: makes the common rerun fail; users can choose a new directory.
- Suffix colliding stems (`s1`, `s1_2`): names would depend on input order and hide the problem.
- Abort the run on a panic: one bad record would lose a long run; strict mode gives that behavior on request.
- Reject records containing tabs: breaks the one-row-per-record guarantee.

## Consequences

- Output names are predictable from input names, and conflicts surface before work is done.
- `run.json` is enough to see which database and parameters produced a directory.
- A custom panic hook is installed once per process; it forwards all other panics to the previous hook.
