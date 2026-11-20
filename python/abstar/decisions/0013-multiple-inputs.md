# 0013: Multiple inputs: per-input output files, flat layout, `source_file` column

## Status

Accepted (2026-11-20)

## Context

`abstar run` accepts a file, a directory or several inputs. Output naming must be predictable, and in-memory results must not vary in shape with the number of inputs.

## Decision

Each input file produces its own output file per format, named after the input file's stem (as today). There is one `run.json` per run, with counts per input. In memory, multiple inputs come back as one DataFrame in input order, always with a `source_file` extension column, null for sequences passed directly, so the columns never vary. This was a call made in the design.

## Alternatives considered

- Columns that appear only when there are multiple inputs: rejected, since the schema would depend on the data.

## Consequences

- Output files map to inputs by name.
- Downstream code sees the same columns regardless of input count.
