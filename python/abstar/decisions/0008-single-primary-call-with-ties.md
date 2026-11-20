# 0008: Single primary gene call plus `*_call_ties`; tie-break order

## Status

Accepted (2026-11-20)

## Context

Many alleles can tie for the top score. Comma-joined tie strings in the call field are hard to search, count and sort, but dropping ties loses information.

## Decision

Applied identically to V, D, J and C (shown for V). `v_call` is exactly one allele, always. `v_gene` and `v_species` describe it. `v_call_ties` lists every allele with the top score, including `v_call` (a native list in Parquet, comma-joined in TSV). `v_species_ties` is aligned entry by entry with it. Tied alleles share the top score, so `v_score` covers them all and there is no score list. There is no ambiguity boolean: it is `len(v_call_ties) > 1`. The name `v_calls` is avoided because it is one letter from `v_call`.

`v_call` is chosen among the ties deterministically, in order: (1) best alignment score; (2) functional before ORF before pseudogene; (3) population usage frequency, if the database provides an optional `usage` column; (4) a fixed natural sort order on allele name. Same input and database always give the same call. This was a user decision.

Evaluation: a call is correct if `v_call` is in the generator's set of equally valid answers; coverage of that set by `v_call_ties` is scored separately. Because tie-break 3 shares usage data with the generator, evaluation also includes generated sets with uniform allele usage.

## Alternatives considered

- **Comma-separated ties in the call field**: rejected; the field would no longer be searchable, countable or sortable.
- **A score list**: unnecessary, since tied alleles share the top score.
- **An ambiguity boolean**: not added; add one only if TSV users ask.

## Consequences

- `v_call` is always usable directly.
- The known bias is small and in the least harmful direction: allele counts within tie sets lean toward common, functional alleles. Rigorous analyses use `v_call_ties`.
- Ties are retained, not discarded.
