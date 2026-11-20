# 0003: `run()` returns a DataFrame or a RunSummary, chosen by the caller

## Status

Accepted (2026-11-20)

Partly superseded: 0020 gives the exact path-or-sequence rule and how the CLI treats its positional, and 0021 makes `-` mean stdin.

## Context

Return types that vary with the data (for example a single object for one sequence and a table for many, or a different type when some records fail) make API use brittle, especially for agents. This is an explicit anti-goal.

## Decision

`abstar.run(...)` has exactly two return shapes, chosen by the caller and never by the data:

- `output_dir=None`: always a Polars DataFrame, including for one sequence or zero. A single record is `df.row(0, named=True)`.
- `output_dir` given: results stream to disk and a small `RunSummary` is returned (counts annotated/failed/unassigned, output paths, run metadata path).

A string argument is treated as a file path if it exists; otherwise it is validated as a nucleotide sequence, and anything else is an error. 0020 refines this rule and 0021 adds `-` for stdin.

## Alternatives considered

- Return an `abutils.Sequence`-like object for single inputs: not carried forward, since abutils compatibility is a non-goal.
- Shape depending on count or on how many records succeeded: rejected as data-dependent.

## Consequences

- Callers can write code without checking what came back.
- In-memory results cross from Rust as Arrow data and become a Polars DataFrame without copying.
- Streaming to disk keeps memory flat for very large inputs.
