# 0033: Per-segment mutation numbering, indels and frequencies

## Status

Accepted (2026-11-21)

## Context

Every segment (V, D, J, C) gets substitution lists (nt and aa), insertion and deletion lists, counts and frequencies (spec *Pipeline* step 7). The generator derives the same fields from each base's origin and its normalized indel placement (`crates/generate/src/derive.rs`, record 0026, rulings A13, A15–A18). The grader compares them by the shared-span rule (`python/abstar/_grade.py`): over the germline span that the annotator and truth both cover, when the annotator's span edges are truth's or valid alternative cuts. The annotator reimplements truth's rules on its own alignment (record 0014). Coordinates are 0-based half-open on the oriented read (record 0016); the strings below are 1-based.

## Decision

**One alignment per segment** (`derive::mutations`). Every field of a segment comes from the same alignment, the one the regions use (record 0032):
- V and J: the aligner's alignment to the primary allele, cut at the junction step's V end and J start (record 0030). Its indels are already at the leftmost minimum-mismatch placement within the indel zone (`place_gaps`, rulings A17–A19). V includes its 5′ start extension (ruling A28): the extended bases are aligned V positions.
- D: the junction step's gapless placement on the primary, on the reverse-complemented allele when the D is inverted (record 0031).
- C: the gapless prefix comparison from J's end (`assign_c`), over the covered C bases only.

**Numbering.**
- V nt: `IMGT{pos}:{ref}>{alt}`, `pos` the IMGT-gapped position. D, J, C nt: `{pos}:{ref}>{alt}`, `pos` 1-based in the allele as used. For an inverted D, that is the reverse-complemented allele, the same coordinates as `d_germline_start`/`d_germline_end`.
- A read base other than the germline base (case-insensitive) is a substitution, `N` included.
- Lists are in germline order.

**Codons (aa lists).**
- V: IMGT codons (`IMGT{n}`, `n` the IMGT codon number). Codons that hold an IMGT gap are skipped.
- D: the V reading frame. The codon phase of the D's first base is counted from the read's last V base and its IMGT position. D codon 1 starts at the first D germline position in that phase.
- J: the frame set by Trp/Phe118. Codon 1 starts at `anchor % 3`.
- C: the J frame continued past the J allele's end, `(anchor − len(J)) mod 3`.
- A codon counts only when its three germline bases are all aligned. A codon with a deleted or unaligned base is not counted. It is listed as `{n}:{ref aa}>{alt aa}` when the translations differ. A deletion edge can change an amino acid without a nt substitution (ruling A15).
- The aa list is null when the segment has no frame: a J allele without an anchor, and C after it.

**Indels.**
- An indel is a run of inserted or deleted alignment columns with an aligned column on each side.
- Insertion: `{pos}:ins:{bases}`, `pos` the germline position the inserted bases follow, `bases` the inserted read bases.
- Deletion: `{pos}:del:{bases}`, `pos` the first deleted germline position, `bases` the deleted germline bases.
- V positions carry the `IMGT` prefix. Indels are not substitutions and are not counted in them.

**Counts and frequencies.**
- `count` is the nt list's length, `count_aa` the aa list's.
- `frequency` is `count` over the aligned germline positions (deleted positions excluded). `frequency_aa` is `count_aa` over the counted codons. It is null when no codon is counted.
- A segment with no aligned base has every field null (ruling A6). So does V when its allele isn't IMGT-gapped (ruling A10).

`SegmentMutations::germline` carries each segment's aligned germline span (0-based half-open on the allele as used), for `{seg}_germline_start`/`_end` in the record.

## Alternatives considered

- **D and C positions on the forward allele for an inverted D:** they would disagree with `d_germline_start`/`_end` and with truth (record 0031).
- **The D frame from the first V base:** an out-of-frame V indel would shift it. Truth phases D from the last V base.
- **Re-placing indels in the mutation step:** the aligner already reports truth's placement (rulings A17–A19). A second placement could disagree with the regions and masks.
- **C over the whole allele:** the read covers only a prefix, and truth counts covered positions only.

## Consequences

Measured in `tests/mutations.rs` by a Rust port of the grader's shared-span rule, on reads whose call for the segment is the true allele. With the true calls, cuts and germline span ("exact"), the derivation equals truth on every D and C read, on every V read outside the indel and bnab sets, and on every naive J read.

| Set (seed) | V nt | D nt | J nt | C nt | V exact | J exact |
|---|---|---|---|---|---|---|
| memory + 40 nt C (12) | 99.95% | 99.81% | 99.34% | 100% (215) | 1,331/1,331 | 1,101/1,106 |
| naive (1) | 99.90% | 99.71% | 99.95% | — | 1,457/1,457 | 1,274/1,274 |
| mutation 0–30% (1) | 99.02% | 99.29% | 98.57% | — | 1,034/1,034 | 886/889 |
| memory + 20% indels (11) | 99.12% | 99.49% | 99.36% | — | 1,884/1,887 | 1,673/1,675 |
| bnab (1) | 93.18% | 97.40% | 97.03% | — | 709/726 | 558/561 |

The aa lists, indel lists, counts and frequencies agree within 1 point of the nt list on each set. Every miss traces to an upstream call:
- **The V start** (no valid list in truth): the 5′ extension (ruling A28) stops short on heavily mutated starts. It can also follow a shifted diagonal past a 5′ indel.
- **Invalid V-end or J-start cuts** (records 0029, 0030).
- **Mutated last J bases:** the J aligner places a 1 nt insertion near J's 3′ end instead of two end mismatches. The C bases after J make that alignment score as well. This causes most J exact misses.
- **Indel representation:** the aligner splits one indel into two, or adds a compensating 1 nt gap, where truth has one event. These are the V exact misses (indel and bnab sets).
- C is graded on 215 of 2,000 reads: 40 nt rarely separate C alleles, so the C primary is the true allele only there.

`mutation_accuracy_memory` asserts at least 99% nt-list agreement for each segment on the memory set.
