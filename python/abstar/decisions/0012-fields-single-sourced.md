# 0012: Output fields single-sourced in `fields.toml`

## Status

Accepted (2026-11-20)

## Context

Today the output schema, serializers and docs are separate, hand-synchronized surfaces.

## Decision

`schema/fields.toml` defines every output field: name, Arrow type, nullability, description, coordinate space where relevant, origin (`airr` required, `airr` optional, or `abstar` extension) and the targeted AIRR schema version. That one file drives the column set and types in the Rust TSV and Parquet writers, `abstar manual fields`, and the Sphinx field reference. Whether it drives Rust by code generation or load-time validation is the builder's call. Both writers read the same in-memory record. Every `abstar`-origin field names the decision record that defines it (`decision = "NNNN"`), enforced by a test. The foundation supplies the mechanism and AIRR-required fields; the full list is decided in the annotator design.

## Alternatives considered

- Separate writer, manual and docs definitions: the drift problem being removed.

## Consequences

- One edit changes writers, manual and docs together.
- Every extension field has a traceable rationale.
