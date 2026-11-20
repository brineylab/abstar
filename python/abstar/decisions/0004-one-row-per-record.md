# 0004: One output row per input record; status columns replace failure logs

## Status

Accepted (2026-11-20)

## Context

Today failed records are omitted from the output and described in per-record diagnostic files under `logs/`. That is hard to audit and invisible to agents. A silent drop is worse than an explicit failure.

## Decision

Every input record produces exactly one output row; records are never silently dropped. The columns `annotation_status` and `failure_reason` distinguish annotated, biologically unassigned (for example no V found), invalid input (for example non-nucleotide characters, or empty), and per-record internal error (with message). Record conservation becomes a row-count check. `--strict` turns any per-record internal error into a fatal run error. Fatal run errors (unreadable input, bad database, I/O failure) stop the run with a documented exit code and error code. When writing to disk, `run.json` records parameters, abstar version, database name and content hash, AIRR schema version, and counts by status.

This was a call made in the design.

## Alternatives considered

- Keep omitting failed rows and writing diagnostic files: more complex, and the `logs/` per-record file tree has to be found and read separately.

## Consequences

- Failures are visible directly in the data, to people and agents alike.
- The `logs/` file tree goes away.
- Consumers filtering to annotated records must filter on `annotation_status`.
