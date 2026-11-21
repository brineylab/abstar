# 0034: Masks

## Status

Accepted (2026-11-21)

## Context

The six masks (`cdr_mask`, `gene_segment_mask`, `nongermline_mask` and their `_aa` forms) label each base or codon of the V(D)J region. The generator derives them from each base's origin (`crates/generate/src/derive.rs`, record 0026). The user ruled that masks cover the V(D)J region only (ruling A12). The grader compares `gene_segment_mask` with slack where the `_valid` cut lists leave a base's segment open. It compares the other masks where both nt gene masks agree, and the aa masks at codons whose three bases agree (`python/abstar/_grade.py`, ruling A30). The annotator reimplements the rules on its own calls (record 0014).

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

**Amended 2026-11-21 (rulings A29, A30).**
- Truth's nongermline masks now follow the same normalized indel placement as its mutation lists (record 0026 amendment). So on both sides one alignment per segment feeds the lists, counts, frequencies and masks.
- The grader now compares the aa masks other than `gene_segment_mask_aa` only at codons whose three bases have the same `gene_segment_mask` label in truth and call. Boundary ambiguity stays graded by `gene_segment_mask` and the boundary fields. The tables below are measured after both changes.

Measured in `tests/mutations.rs` by a Rust port of the grader's mask rule. The graded reads are those whose V, J (and D, when truth has one) calls are the true alleles. With truth's calls, V start, cuts and J end, all six masks equal truth's exactly on naive reads (334/334) and at 0–30% mutation (120/120). On the indel set it is 277/281 and on bnab 17/18; the remaining misses are indels the aligner represents differently from truth (record 0033).

| Set (seed) | `cdr` | `gene_segment` | `nongermline` | `cdr_aa` | `gene_segment_aa` | `nongermline_aa` |
|---|---|---|---|---|---|---|
| memory + 40 nt C (12) | 98.11% | 99.39% | 99.46% | 98.11% | 99.53% | 99.59% |
| naive (1) | 98.70% | 99.71% | 99.76% | 98.70% | 99.71% | 99.82% |
| mutation 0–30% (1) | 96.85% | 98.52% | 98.24% | 96.95% | 98.98% | 98.24% |
| memory + 20% indels (11) | 97.59% | 99.09% | 98.91% | 98.05% | 99.14% | 99.05% |
| bnab (1) | 94.69% | 95.36% | 92.54% | 94.86% | 95.52% | 93.20% |

- **`cdr_mask` misses** are mostly length mismatches. They come from a different V start or J end, or from a CDR3 known on one side only, when the junction step's valid cut decides an anchor (record 0032).
- **The other masks' misses** come from upstream calls: the V start, cuts outside the valid lists, D placement, and indel representation (record 0033).
- **Before ruling A30,** `nongermline_mask_aa` agreed on only about 80% of reads. The cause was a valid cut inside a codon that left its 3rd base's segment unchanged: truth marked an N base that matched the germline by chance as `1`, the annotator counted it as V/D/J, and the grader compared that codon anyway.
- **Before ruling A29,** `nongermline_mask` agreed on 95.69% of the indel set. Truth's masks followed the generated placement while its lists followed the normalized one.
- **2026-11-21 (ruling A33):** truth's `cdr_mask` now follows the regions on each segment's normalized indel placement, as the annotator's does. `cdr_mask` on igh_bnab rose from 96.55% to 96.91% on the accuracy gate.
