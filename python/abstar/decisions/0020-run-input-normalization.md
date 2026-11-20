# 0020: abstar.run() input normalization

## Status

Accepted (2026-11-20)

## Context

`abstar.run()` is the one entry point for Python and the CLI. Callers pass a file, a directory, a bare sequence, a list of sequences, or `(id, sequence)` pairs, and a string is ambiguous between a path and a sequence. A typo in a path must not silently be annotated as a sequence.

## Decision

**Path or sequence.** A `str` or `Path` that exists on disk is a path. A string that does not exist but contains `/` or `\`, or ends in `.fasta`, `.fa`, `.fna`, `.fastq` or `.fq` (optionally `.gz`), raises `INPUT_NOT_FOUND`. Any other string (including the empty string) is a sequence.

**Directories** expand recursively to files with those extensions, sorted by relative path. `source_file` is the path as joined from the directory given (`<dir>/<relative path>`), not resolved.

**Default IDs.** A bare string is one sequence with ID `"1"`. An iterable of sequence strings gets IDs `"1"` to `"N"` by position. An iterable of `(id, sequence)` string pairs keeps its IDs. IDs are always strings.

**Mixing** paths with sequences or records raises `INVALID_ARGUMENT`. Sequences and records need no files, so combining them with `output_dir` also raises `INVALID_ARGUMENT`.

**Iterables** are consumed exactly once (generators work). `formats` must be a non-empty subset of `airr`, `parquet`; `threads` is None or an integer of at least 1; both are checked in Python before the core is called. Without `output_dir`, `formats` is ignored.

## Alternatives considered

- Treat every string as a sequence unless a flag says otherwise: a mistyped filename would be annotated as garbage.
- Require explicit `Path` objects for files: unfriendly for the CLI and notebooks.
- Zero-based or per-file default IDs: `"1"`-based positions match the single-sequence case and need no file context.

## Consequences

- Typos in paths fail fast with a clear code.
- A sequence string that happens to be an existing file name is read as a file; a bare nucleotide string cannot realistically collide.
- Writing in-memory sequences to files is not supported in this version.
