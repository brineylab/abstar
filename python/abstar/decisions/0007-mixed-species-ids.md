# 0007: Mixed-species IDs: internal prefix, plain names in output, per-segment species

## Status

Accepted (2026-11-20)

## Context

Mixed-species databases (for example `human+c57bl6`) can contain identical allele names from different species, such as `IGHV1-2*01` in both human and mouse. Output must stay compatible with AIRR tools, and humanized-mouse sequences mix species within a single rearrangement (human V/D/J with mouse C).

## Decision

Inside a mixed-species database, allele IDs are made unique with a short species prefix (for example `hs|IGHV1-2*01`, `mm|IGHV1-2*01`). Single-species databases use plain names. Output never shows the prefix: call fields contain standard allele names only. Species is reported per segment (`v_species`, `v_species_ties`, and so on). A cross-species tie with identical names reads `v_call_ties = [IGHV1-2*01, IGHV1-2*01]` and `v_species_ties = [human, mouse]`, with nothing deduplicated. This was a user decision; per-segment species rather than one `species` column was a call made in the design.

## Alternatives considered

Single-species databases always keep plain names. For mixed-species databases:

- **A. Suffix every allele** (the legacy behavior), e.g. `IGHV1-2*01__homo_sapiens`. Consistent, but every call carries the suffix, which tools like Change-O do not expect.
- **B. Suffix only colliding names.** Shorter, but one database mixes plain and suffixed names, and adding a species later can rename existing alleles.
- **C. Prefix every allele with a short species tag**, e.g. `hs|IGHV1-2*01`. Like A but shorter and easier to strip, still non-standard in `v_call`.
- **D (chosen, proposed by the user).** Use C's prefix internally (the allele `id` in the germline database, so names stay unique), but output only the plain gene/allele name in `v_call` and the other call fields, plus per-segment `{s}_species` columns. `v_call` stays standard and compatible with other AIRR tools, while users can filter on gene and species together.

## Consequences

- Call fields are AIRR-compatible.
- Ties across species remain visible and are aligned entry by entry with `*_call_ties`.
- Allele IDs are opaque strings internally and the prefix is stripped at output.
