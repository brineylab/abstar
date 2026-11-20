# 0023: Recombination model format and binding to a germline database

## Status

Accepted (2026-11-20)

## Context

`abstar generate` draws V(D)J rearrangements from a recombination model: usage, trimming and insertion distributions. A model is written against gene names, while generation needs alleles with anchors, and users may supply any germline database, so some model genes will not exist in it.

## Decision

A model is a versioned TOML file (`format_version = 1`) with strict keys; unknown keys and probability tables that do not sum to 1 are `MODEL_INVALID`. Genes are matched to database alleles by plain `gene` name, restricted to the model's `species` and `locus`, so mixed-species databases (prefixed ids, plain names) work. D usage is conditioned on the J gene (`d_given_j`) and bound per J row.

Allele probability is the gene probability times the model's `allele_usage` split, or an equal split across the gene's generatable alleles. Only V and J alleles with an anchor can be generated; an anchorless allele's share goes to its gene's anchored alleles. D alleles need no anchor. A gene with no generatable allele is dropped with a warning naming it, and the remaining mass is renormalized (per J row for D).

Binding fails with `MODEL_INVALID`, listing the dropped genes, if less than 0.90 of the model's V, D or J mass (D weighted by J usage) is retained. Every allele of the locus and species, anchored or not, is a candidate for later equal-answer comparison; V candidates must have a gapped sequence containing `.`, otherwise binding fails (this rejects databases without IMGT gaps, such as macaque).

## Alternatives considered

- Fail on any missing gene: makes models unusable with slightly different database releases.
- Silently drop missing genes: hides a model that does not fit the database. The 0.90 threshold keeps drift visible.
- Generate anchorless alleles: junction positions would be undefined.
- Restrict candidates to generatable alleles: reads could then be wrongly credited with a unique answer when an anchorless allele is equally good.

## Consequences

- Warnings name every dropped gene; the 0.90 threshold is a fixed constant.
- Candidate sets include anchorless alleles, so equal-answer sets can contain alleles that are never generated.
- Changing the format requires a new `format_version`.
