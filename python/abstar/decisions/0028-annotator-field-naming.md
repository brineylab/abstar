# 0028: Annotator output field naming

## Status

Accepted (2026-11-20)

## Context

The annotator's output schema (`crates/annotate/fields.toml`) grows from the foundation columns to the full v1 list. The generator's truth table (record 0026) already names many of the same quantities, and the grader joins truth and annotation on column names, so the two sides must agree wherever both exist. Legacy abstar also had established column names that existing users rely on.

## Decision

- **AIRR names where AIRR defines them.** A field AIRR 2.0 defines (calls, coordinates, regions, `junction`, `locus`, `vj_in_frame`, `stop_codon`, `np1`, ...) keeps the AIRR name and meaning, with origin `airr-optional`.
- **Legacy names are reused** where abstar already had the quantity (`fwr1_aa`, `cdr3_length`, `sequence_vdj`, the masks).
- **New-field pattern.** Anything else is `<segment>_<quantity>` (`v_mutations`, `j_insertions`, `c_mutation_frequency_aa`) or `<region>_<quantity>` (`cdr3_length`), with an `_aa` suffix for the amino-acid variant. Segments are `v`, `d`, `j`, `c`.
- **Truth parity.** Where the generator's truth table has a column, the annotator uses the identical name, type and encoding.
- **Coordinates.** Each field states its coordinate space; `*_start`/`*_end` are 1-based closed in `sequence`, germline positions are 1-based in the ungapped called allele, V mutation positions are IMGT-gapped.
- **Masks** cover only the V(D)J region: nt masks index `sequence_vdj`, aa masks `sequence_vdj_aa`, and `gene_segment_mask` uses `V`, `N`, `D`, `J`.

## Alternatives considered

- Keep legacy-only names: rejected, they lack segments beyond V and diverge from the truth columns.
- Prefix every abstar-specific field: rejected, it would break the AIRR-style names users already read.

## Consequences

- Output columns are stable for graders and downstream tools.
- Every new abstar-origin field cites this record until a later record covers it more specifically (inverted D, region lengths, mutation encodings and masks are re-pointed as those records land).
