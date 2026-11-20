# 0022: Shared table crate for schema codegen and writers

## Status

Accepted (2026-11-20)

## Context

`abstar generate` writes tables in the same AIRR TSV and Parquet encodings as `abstar run`, and builds its record types from a `fields.toml` in the same way. The schema generator and writers lived in `abstar_annotate`. 0014 requires the generator to be independent of the annotator, so it cannot depend on that crate.

## Decision

A new crate, `abstar_table`, depends on no other abstar crate. It holds the `fields.toml` code generator (`codegen`, used from `build.rs` as a build-dependency) and the TSV and Parquet writers (`writers`), which take the Arrow schema as an argument instead of reading a fixed one. Cargo features `codegen` and `writers` keep Arrow and Parquet out of the build-script dependency graph. `abstar_annotate::writers` stays as thin wrappers bound to its own schema, so `abstar run` output is unchanged. `list<int64>` is added as a field type (TSV: integers joined by `,`). A Python test enforces that `abstar_table` has no abstar dependencies and that neither it nor `abstar_generate` depends on `abstar_annotate`.

## Alternatives considered

- Have the generator depend on `abstar_annotate` for the writers: violates 0014 and pulls annotation code into the generator.
- Copy the writers and codegen into the generator crate: two encoders could drift, breaking the shared TSV and Parquet conventions.
- Put them in `abstar_germline`: unrelated responsibility, and the annotator would gain no benefit.

## Consequences

- One implementation of the output encodings serves both tools.
- Generated code still expects the Arrow imports to be in scope in the including module.
- Changes to the writers or codegen affect both tools and must keep `abstar run` output stable.
