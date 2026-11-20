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

**Commands.** A build is staged in a sibling `<dir>.staging-<pid>` and swapped in. `overwrite` replaces only an existing database (a directory holding `manifest.toml`); any other existing path is refused with `DATABASE_EXISTS` and left untouched. `list` shows built-in, then user databases. A database that fails to load is still listed, with `n_alleles` 0 and its error in `error` (load errors name the database directory), so one corrupt database cannot break the listing. Leftover `*.staging-*` and `*.old-*` directories from an interrupted build are never listed or resolved. `info` shows the whole manifest (built date, abstar version, `alleles_sha256`, prefixes, provenance, warnings) plus counts by locus/segment/species; `--json` carries the same fields. A human-readable `build` prints, when a keep-list was applied, the `restricted:` and `passed through:` locus/segment lists, and each warning as `warning: <allele>: <message>`.

## Alternatives considered

- A compiled or cached binary format: not added; add only if load time is measured to matter.
- One database per receptor: rejected, for the reasons above.
- Fail `list` on the first unreadable database: one corrupt donor database would hide all the others.
- Let `overwrite` replace any existing path: a mistyped `--output` could delete an unrelated directory.

## Consequences

- Plain text is diffable and agent-readable, and loads in milliseconds.
- Mixed BCR/TCR inputs work and one required argument disappears.
- Commands: `abstar germline build`, `list` and `info`.
