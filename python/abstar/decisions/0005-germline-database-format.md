# 0005: Germline database format: `manifest.toml` + `alleles.tsv`; a database is a set of loci

## Status

Accepted (2026-11-20)

## Context

The engine needs a database format that is diffable, readable by agents, and fast to load. Today's `gl.py` has an open TODO about combined databases, and a required receptor argument prevents mixed BCR/TCR inputs.

## Decision

A database is a directory with `manifest.toml` (name, format version, receptor(s), species, provenance with source/release/date, content hash) and `alleles.tsv`, one row per allele. Required columns: opaque allele ID (never parsed), locus and segment (V/D/J/C), species, functionality, IMGT-gapped sequence, V and J anchor positions, source, and `gene` (see below). The ungapped sequence is derived by stripping gaps. Exact columns are the builder's call.

A database is a set of loci, not a receptor: `human` holds both IG and TR loci. `--receptor`/`--locus` becomes an optional filter rather than a required selector. IG and TR V genes are far too divergent for cross-receptor confusion to be a practical risk.

The `gene` column is filled at build time from IMGT-style names when they parse, otherwise from a user-supplied mapping, otherwise each allele is its own gene. IDs are parsed only at build time, never during annotation.

Builds validate unique IDs, consistent gapped/ungapped sequences, locus and segment presence, and detectable anchors (every allele whose anchor cannot be found is reported by name). Builds are staged and appear atomically. Search indexes are built from the text at load time.

## Alternatives considered

- A compiled or cached binary format: not added; add only if load time is measured to matter.
- One database per receptor: rejected, for the reasons above.

## Consequences

- Plain text is diffable and agent-readable, and loads in milliseconds.
- Mixed BCR/TCR inputs work and one required argument disappears.
- Commands: `abstar germline build`, `list` and `info`.
