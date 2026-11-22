# 0020: abstar.run() input normalization

## Status

Accepted (2026-11-20)

## Context

`abstar.run()` is the one entry point for Python and the CLI. Callers pass a file, a directory, a bare sequence, a list of sequences, or `(id, sequence)` pairs, and a string is ambiguous between a path and a sequence. A typo in a path must not silently be annotated as a sequence.

## Decision

**Path or sequence.** A `str` or `Path` that exists on disk is a path. A `Path` that does not exist raises `INPUT_NOT_FOUND`. A string that does not exist is a sequence only if it consists solely of valid nucleotide characters (`ACGTN` plus IUPAC `RYSWKMBDHV`, case-insensitive, per 0016; the empty string counts and becomes an `invalid_input` row); any other string raises `INPUT_NOT_FOUND` ("not an existing file and not a nucleotide sequence"). This is the rule in 0003. An earlier version of this record treated every non-existing string without `/`, `\` or a FASTA/FASTQ extension as a sequence, so `run("sample1")` annotated the word as a sequence and succeeded; that was reversed.

**Existence checks.** Only a "name too long" or "invalid name" error from the file system means a string cannot be a path (a long sequence). Any other error, such as permission denied on a parent directory, raises `INPUT_UNREADABLE` instead of being taken as "does not exist".

**CLI.** The SEQUENCES positional of `abstar run` is always a path (or `-`, 0021), never a sequence: `abstar run ACGT` raises `INPUT_NOT_FOUND`. Use `--seq` for a sequence.

**Directories** expand recursively to files with those extensions, sorted by relative path. `source_file` is the path as joined from the directory given (`<dir>/<relative path>`), not resolved.

**Default IDs.** A bare string is one sequence with ID `"1"`. An iterable of sequence strings gets IDs `"1"` to `"N"` by position. An iterable of `(id, sequence)` string pairs keeps its IDs. IDs are always strings.

**Mixing** paths with sequences or records raises `INVALID_ARGUMENT`. Sequences and records may also be written to files: they form one in-memory input with output stem `records` (`records.tsv`, `records.parquet`) and a null `source_file` (0018).

**Iterables** are consumed exactly once (generators work). `formats` must be a non-empty subset of `airr`, `parquet`; `threads` is None or an integer of at least 1; both are checked in Python before the core is called. Without `output_dir`, `formats` is ignored.

## Alternatives considered

- Treat every string as a sequence unless a flag says otherwise: a mistyped filename would be annotated as garbage.
- Treat a string as a path only if it contains `/` or `\` or ends in a FASTA/FASTQ extension (the earlier rule): `sample1` or `reads.txt` was annotated as a sequence and exited 0.
- Require explicit `Path` objects for files: unfriendly for the CLI and notebooks.
- Zero-based or per-file default IDs: `"1"`-based positions match the single-sequence case and need no file context.

## Consequences

- Typos in paths fail fast with a clear code, whatever the file name looks like.
- A string with a non-nucleotide character (such as `-` gaps or `X`) cannot be passed as a bare sequence; pass it as an `(id, sequence)` pair to get an `invalid_input` row instead.
- A sequence string that happens to be an existing file name is read as a file; a bare nucleotide string cannot realistically collide.
- `abstar run --seq ACGT OUT_DIR` works through the same path as file inputs.
