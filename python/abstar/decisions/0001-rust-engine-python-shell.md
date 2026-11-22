# 0001: Rust engine behind a Python shell; Typer CLI

## Status

Accepted (2026-11-20)

## Context

abstar's accuracy problems and bugs trace back to its structure: Python orchestration around MMseqs and Parasail subprocesses, temp-file chunking, and several hand-synchronized schema, serializer and doc surfaces. The rewrite needs a fast core and a CLI and API that cannot drift apart.

## Decision

A Python package (`abstar`) provides the Typer CLI and the Python API. A Rust engine, exposed through PyO3 and built with maturin into a single wheel, owns every hot path: FASTA/FASTQ parsing, batching, multithreading, annotation, ordered reassembly and TSV/Parquet writing. Python passes paths and options and never touches per-sequence data on the bulk path. Rust releases the GIL. Typer, the extension module and Rich are imported lazily so `abstar --help` and small runs start fast.

The CLI is built with Typer. This was a user decision.

## Alternatives considered

- **Rust all the way down** (clap CLI in Rust, Python only an API binding): millisecond startup and a Python-free binary, but the CLI and API would be defined twice in two languages, with a higher bar for lab contributors and awkward Python-ecosystem commands. Rejected; it can be added later over the same crates.
- **Rust binary called as a subprocess**: this is the MMseqs pattern being abandoned (IPC and serialization cost, no in-memory DataFrames, startup and temp-file overhead). Dropped.

## Consequences

- One language defines the CLI and the API, avoiding today's sync problems.
- Python contributors can add commands.
- About 100 ms of interpreter startup is accepted; it is small next to the MMseqs startup it replaces, and speed is unchanged because all heavy work is in Rust.
- A standalone Rust binary remains possible later on the same crates.
