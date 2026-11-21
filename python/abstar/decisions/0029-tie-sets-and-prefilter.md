# 0029: Orientation, V search and tie sets by IMGT projection

## Status

Accepted (2026-11-20)

## Context

The annotator's first steps (spec *Pipeline* steps 2–3) orient the read, find the V gene and report a tie set by the same fewest-mismatches rule the generator's truth uses (record 0025). Truth sets are computed from each read base's known origin; an annotator only sees the read.

Aligning every V allele of a locus is too slow: one 300-nt V alignment costs about 104 µs score-only and about 148 µs with traceback on the development machine. The spec's first idea was a k-mer prefilter with a margin that keeps every possible tie, but no useful margin exists: at high mutation the true allele's k-mer count is within the noise of unrelated alleles.

Ruling A20 (2026-11-20) replaced it. Rank alleles with k-mers, align only a few, project the read onto IMGT gapped positions through the best alignment, and count mismatches for every allele of the locus at those positions. Rulings A21–A23 (2026-11-20) settle the exact tie rule, the recovery of reads the first pass misses, and spurious gap pairs.

## Decision

**Alignment conventions** (`crates/annotate/src/align`)
- k-mer prefilter with k = 7. A k-mer containing a non-ACGT base is skipped.
  - `hits` counts the read k-mer positions whose k-mer an allele contains.
  - `diagonal_hits` counts only those on the allele's best diagonal band: two adjacent bins of 8 diagonals.
- Scoring: match +2, mismatch −3, gap open −10, gap extend −1; a gap of length n costs −10 − (n − 1). These values keep a 36-nt deletion and a 27-nt insertion as single gaps.
- A non-ACGT read base always scores as a mismatch, including a read N against a germline N.
- Overlap (semi-global) alignment: both read ends and both germline ends are free.
  - Terminal gaps are penalized when choosing the alignment, then dropped, so an alignment starts and ends on aligned bases.
  - `score` is the DP optimum, equal to `score_semiglobal`'s.
  - An optimum ≤ 0 gives an empty alignment.
- `place_gaps` moves each gap to the placement of its kind and length with the fewest mismatches, leftmost on ties (rulings A17–A19).
  - For V it is called with `protect` = the germline from 2 nt before the Cys104 codon to its end, as in the generator.
  - It counts a read N against a germline N as a mismatch.
- `collapse_gap_pairs` (A23) replaces each balanced group of consecutive gaps by the gapless diagonal over the same stretch.
  - A balanced group has equal inserted and deleted lengths, typically an equal-length insertion+deletion pair, with at most 40 aligned columns between its first and last gap.
  - The replacement happens when the diagonal scores within one gap open (10) of the gapped stretch. Ties go to the diagonal.
  - Such groups are compensations for clustered mismatches. The generator places at most one indel per read, so they are never true events there.

**Orientation.**
- For each strand, sum the best V allele's k-mer hits and the best J allele's.
- The larger total wins. Equal totals choose the lexicographically smaller strand, so a read and its reverse complement orient identically. This replaces the spec's "ties go forward", which made a read and its reverse complement differ.
- No V k-mer hit on either strand gives no V hit.
- The best allele is used rather than the sum over alleles, because that sum is dominated by k-mers shared by many alleles. On 4,000 reads with 0–30% mutation, summing misoriented 15 reads and the best-allele rule none.
- If the chosen strand gives no V hit, the other strand is tried (`call_v`).

**Candidates and the no-hit threshold.**
1. V genes are ranked by their best allele's k-mer count. The best allele of each of the first 3 genes within 90% of the leader's count is scored with `score_semiglobal`. A lone leader is aligned directly.
2. If the best scores below 50, genes are re-ranked by `diagonal_hits` and the first 12 are scored. If no allele has 8 diagonal-banded hits, this step is skipped. Random reads have a median of 7.
3. If that also fails, every V gene of one locus is scored, in diagonal order. The locus comes from either of two gates:
   - the locus of a J allele (best of the top 3 by k-mers) whose alignment score is at least 30; random reads reach at most 27;
   - otherwise, the locus of the best V so far, if that V scored at least 20; random reads, both strands and 13 genes, reach at most 23.
4. Below 50 there is no V hit (`no_v_hit`). The floor keeps clear of every non-antibody score measured:
   - 20,000 random reads of 30–1,030 nt, best of 12 candidates: at most 31;
   - every human TR V allele read against the IG-only database: at most 12.

Diagnosis of the no-hit losses on bnab-preset reads (A22):
- Misoriented heavily mutated reads are rescued by the strand retry.
- Candidate recall: the true gene often ranked 12th to 100th by plain 7-mer counts, because chance k-mers of the large IGHV3 family outnumber the few intact true k-mers. Diagonal banding and the locus-wide search fix most of these.
- About 2 reads are limited by the floor: their true allele scores 40–46.

**Projection.**
- The best allele is aligned with `align_semiglobal`, and gap pairs are collapsed.
- The alignment is cut to the read's V-derived bases: its maximum-scoring contiguous run of columns (Kadane; longest, then leftmost, among equals). The cut scores match +2, mismatch −2 and gaps as above. The symmetric mismatch avoids biasing the cut toward the aligned allele, and the gap costs keep a true indel near the V end inside the cut.
- The cut also ends before any gap whose right flank reaches the protected zone. V has no indels there, so such a gap marks where V ended.
- `place_gaps` then runs. Each aligned read base takes the allele's IMGT gapped index; inserted bases take none.

**Refinement.** Refinement runs when the cut alignment has a gap, or when another allele has fewer mismatches than the projecting allele.
- The (up to 3) fewest-mismatch alleles not yet scored are scored.
- One with a higher DP score projects instead.
- There are at most 2 rounds.

**Tie set (A21: exact rule).**
- Every V allele of the projecting allele's locus is counted, across all species in a mixed-species database.
- Its count is the mismatches over the V-derived bases at their IMGT positions:
  - inserted bases are skipped, and deleted positions have no read base;
  - a gap `.`, a position past the allele's end, or a non-ACGT read base is a mismatch;
  - allele positions beyond the read are ignored.
- The record 0008 pick among the fewest-mismatch alleles is aligned. Its own cut defines the V-derived bases, and every allele is counted once more over that cut.
- The tie set is exactly the alleles with the fewest mismatches there; nothing is added for boundary uncertainty.
- The primary is their record 0008 pick: F > ORF > P > NF > unknown, then higher `usage` (a missing usage after any given one), then natural sort of names. `mismatches` is the count.
- The reported alignment and IMGT projection are the primary's own.

## Alternatives considered

- **Prefilter margin** (the brief's first rule): no margin both bounds the work and guarantees completeness. At 25–30% mutation the true allele's 7-mer count (about 24) is within random-read noise (30 for 600 nt).
- **Align every allele:** about 200 × 104 µs per IGH read; rejected on cost (A20).
- **Summed k-mer orientation:** misorients heavily mutated reads (above).
- **Plain 7-mer ranking for the fallback:** the true gene of missed bnab reads ranked 12th to 100th, so 24 of 1,729 recoverable bnab reads (1.39%) stayed unassigned. Diagonal banding brought that to 5 (0.29%).
- **A J-only gate for the locus-wide search:** J is as heavily mutated as V in these reads, so its score rarely reaches 30. The gate on the best V score adds the rest.
- **An early reject on plain k-mer evidence:** best-allele 7-mer counts overlap between random reads (99th percentile 40) and bnab reads (1st percentile 23). The diagonal count separates them better (random median 7), so its threshold of 8 is used.
- **End slack** (union of fewest-mismatch sets with either end of the V-derived bases moved by up to 3 nt): it reported alleles with strictly more mismatches, which contradicts records 0008 and 0025, so A21 rejected it.
- **Cutting with the projecting allele instead of the primary:** measured within noise of the primary re-cut. Coverage on 5'-truncated memory reads was 98.41% against 98.34%. The primary re-cut is kept because the V-derived bases then belong to the reported allele.
- **Other cut scores:**
  - +1/−1 with alignment gaps loses true indels near V ends: single-indel placement agreed with truth on 468 of 517 reads, against 494 of 518.
  - +2/−3 biases the cut toward the aligned allele: bnab primary accuracy 96.8%, against 98.9%.

## Consequences

Measured with the built-in `human` database (release build). Rates are over recoverable reads.

| Set | Primary in `v_call_valid` | `v_call_valid` ⊆ ties | No V hit | Mean ties |
|---|---|---|---|---|
| naive, seed 1 (2,000) | 100.00% | 99.85% | 0 | 1.00 |
| mutation 0–30%, seed 1 (1,871) | 99.36% | 99.68% | 0 | 1.15 |
| memory, 5' truncation 0–250 nt, seed 13 (2,886) | 99.34% | 98.34% | 2 | 2.06 |
| memory + 20% indels, seed 11 (2,913) | 99.86% | 99.83% | 0 | 1.05 |
| bnab preset, seed 1 (1,729) | 98.90% | 98.44% | 5 | 1.27 |
| mutation 30–40%, seed 1 (1,641) | 64.5% | 64.5% | 542 | 1.37 |

Spurious gaps, measured as V alignments with a gap on indel-free reads, with `collapse_gap_pairs` off and on:

| Set | Off | On |
|---|---|---|
| mutation 0–30% | 1.87% | 0.43% |
| bnab | 8.14% | 1.11% |

On reads with a single true V indel where the primary is the true allele, the placement matches truth on 494 of 518 with collapse off and on. The remaining mismatches are indels within a few bases of a V end that the cut removes.

- **Ties are exact but not always truth's set.** The read alone cannot show where V ends. N bases that continue an allele's 3' tail count for that allele, while truth counts only V-origin bases. This dominates the misses on 5'-truncated reads, whose short V leaves few positions. Above about 25% mutation, the remaining compensating gaps and alignments to the wrong CDR length misplace bases.
- **Junk-read cost.** No-hit input costs about 1.3 ms per read, both strands included.
- **Typical cost:** 0.3 ms per naive read.
- **Gapping limits.** An allele whose IMGT gapping is shifted (#138), or that lacks gaps (macaque, #137), is compared at the wrong positions, as in the truth sets.
- **Constants.** Changing any constant here changes calls, so it needs a superseding record: k, scores, 90%, 3, 12, bin 8, 8 diagonal hits, J 30, near V 20, floor 50, the cut scores, or the pair span of 40.
