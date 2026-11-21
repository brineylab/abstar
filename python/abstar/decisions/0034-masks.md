# 0034: Masks

## Status

Accepted (2026-11-21)

## Context

The six masks (`cdr_mask`, `gene_segment_mask`, `nongermline_mask` and their `_aa` forms) label each base or codon of the V(D)J region. The generator derives them from each base's origin (`crates/generate/src/derive.rs`, record 0026). The user ruled that masks cover the V(D)J region only (ruling A12). The grader compares `gene_segment_mask` with slack where the `_valid` cut lists leave a base's segment open. It compares the other masks where both gene masks agree (`python/abstar/_grade.py`). The annotator reimplements the rules on its own calls (record 0014).

## Decision

**Extent** (`derive::masks`).
- Nt masks have one character per base of `sequence_vdj` (`Regions::vdj`, record 0032). That is the first V base through J's end, or through the read's end without J. No constant-region base is ever in a mask.
- Aa masks have one character per codon of `sequence_vdj_aa`: the complete codons in the V reading frame (`Regions::frame`).

**Per base.**
- `gene_segment_mask`: `V`, `D` or `J` for bases in that segment's alignment (record 0033's alignments), else `N`. P nucleotides are `N`; the annotator cannot tell them from N bases. An inserted base takes the segment whose alignment holds it (ruling A5).
- `nongermline_mask`: `0` when the base is aligned to a germline base equal to it (case-insensitive), else `1`. So N/P bases, inserted bases and mismatches (an `N` read base included) are `1`.
- `cdr_mask`: the base's IMGT region from the region labelling (`Regions`, record 0032): `1`/`2`/`3` in CDR1/2/3, else `0`. An inserted base takes the region of the base before it. It is null when CDR3's ends are unknown, as in truth: the V keeps no Cys104 base, or the J allele has no anchor or starts past it. A CDR3 that is only partial (ruling A4) still gives a mask.

**Per codon.**
- `gene_segment_mask_aa` and `cdr_mask_aa` take the character of the codon's 3rd base. So a codon crossing a boundary belongs to the downstream segment or region.
- `nongermline_mask_aa` is `1` when any of the codon's bases has no germline base (N, P or inserted). It is also `1` when the translation of its germline bases differs from the read codon's. A codon's germline bases may come from two segments.

## Alternatives considered

- **Masks over `sequence_vdjc`:** superseded by the user's ruling A12.
- **The codon's majority segment for `_aa`:** the 3rd-base rule matches the regions' downstream rule (record 0032) and truth.
- **`nongermline_mask_aa` from substitutions only (ignoring N bases that match germline by chance):** truth marks every N/P base, so a codon with one is non-germline whatever its amino acid.

## Consequences

Measured in `tests/mutations.rs` by a Rust port of the grader's mask rule, on reads whose V, J (and D, when truth has one) calls are the true alleles. With truth's calls, V start, cuts and J end, all six masks equal truth's exactly: 334/334 naive reads and 120/120 at 0–30% mutation. On the indel and bnab sets, 266/281 and 17/18. Those misses are indel placements (below).

| Set (seed) | `cdr` | `gene_segment` | `nongermline` | `cdr_aa` | `gene_segment_aa` | `nongermline_aa` |
|---|---|---|---|---|---|---|
| memory + 40 nt C (12) | 98.11% | 99.39% | 99.46% | 98.24% | 99.53% | 82.62% |
| naive (1) | 98.70% | 99.71% | 99.76% | 98.70% | 99.71% | 81.27% |
| mutation 0–30% (1) | 96.85% | 98.52% | 98.24% | 97.13% | 98.98% | 82.98% |
| memory + 20% indels (11) | 97.59% | 99.09% | 95.69% | 98.05% | 99.14% | 78.93% |
| bnab (1) | 94.69% | 95.36% | 89.55% | 94.86% | 95.52% | 78.61% |

- **`nongermline_mask_aa` (about 80%)** is low because of the grading rule, not the derivation. Take a valid alternative cut that moves a boundary inside a codon but leaves its 3rd base's segment unchanged. Truth marks an N base there that matches the germline by chance as `1`. The annotator counts that base as V/D/J and marks the codon `0`. The grader compares that codon because the 3rd-base gene labels agree. With truth's cuts the codons agree exactly (above). Comparing aa masks only where all three bases' gene labels agree would grade the derivation.
- **`cdr_mask`** misses are mostly length mismatches. They come from a different V start or J end, or a CDR3 known on one side only, when the junction step's valid cut decides an anchor (record 0032).
- **`nongermline_mask` on indel reads.** Truth's masks follow the generated indel placement (its origin map). Its mutation lists follow the normalized placement (rulings A16, A17). On 36 of 577 indel reads (memory + 20% indels, seed 11), truth's V nongermline count differs from its own `v_mutation_count` plus inserted bases. The annotator's masks follow its one alignment, the one its mutation lists use. A follow-up ruling should move truth's masks to the normalized alignment, as ruling A16 did for nt lists.
