# 0029: Orientation, V search and tie sets by IMGT projection

## Status

Accepted (2026-11-20)

## Context

The annotator's first steps (spec *Pipeline* steps 2–3) orient the read, find the V gene and report a tie set by the same fewest-mismatches rule the generator's truth uses (record 0025). Truth sets are computed from each read base's known origin; an annotator only sees the read.

Aligning every V allele of a locus is too slow. On the development machine one 300-nt V alignment costs about 104 µs score-only and about 148 µs with traceback. The spec's first idea was a k-mer prefilter whose margin keeps every possible tie, but no useful margin is derivable: at high mutation the k-mer counts of the true allele and of unrelated alleles overlap. Ruling A20 (2026-11-20) replaced it: rank alleles with k-mers, align only a few, project the read onto IMGT gapped positions through the best alignment, and count mismatches for every allele of the locus at those positions.

## Decision

**Alignment conventions** (`crates/annotate/src/align`, from Task 4)
- k-mer prefilter with k = 7. A k-mer containing a non-ACGT base is skipped. An allele's hit count is the number of read k-mer positions whose k-mer the allele contains.
- Scoring: match +2, mismatch −3, gap open −10, gap extend −1; a gap of length n costs −10 − (n − 1). These values keep a 36-nt deletion and a 27-nt insertion as single gaps.
- A non-ACGT read base always scores as a mismatch, including a read N against a germline N.
- Overlap (semi-global) alignment: both read ends and both germline ends are free. Terminal gaps are penalized when choosing the alignment, then dropped, so an alignment starts and ends on aligned bases; `score` is the DP optimum, equal to the score-only pass's.
- An optimum ≤ 0 gives an empty alignment (no hit).
- `place_gaps` moves each gap to the placement of its kind and length with the fewest mismatches, leftmost on ties (rulings A17–A19). For V it is called with `protect` = the germline from 2 nt before the Cys104 codon to the end, so no gap lies within 3 nt upstream of Cys104 or on its CDR3 side, as in the generator.

**Orientation.** For each strand: the best V allele's k-mer hits plus the best J allele's k-mer hits. The larger total wins; a tie goes forward. No V k-mer hit on either strand gives no V hit (`no_v_hit`). The best allele, not the sum over alleles, is used because the sum is dominated by k-mers shared by many alleles: on 4,000 reads with 0–30% mutation, summing misoriented 15 reads and the best-allele rule none.

**Candidates and the no-hit threshold.**
- V genes are ranked by their best allele's k-mer count. The best allele of each of the first 3 genes whose count is at least 90% of the leader's is scored with `score_semiglobal`. A lone leader is not scored; it is aligned directly.
- If the best scores below 50, up to 12 genes are scored. Below 50 there is no V hit.
- The threshold 50 sits above every score seen from non-antibody input: the best of 12 candidates over 20,000 random reads of 30–1,030 nt reached 31, and every human TR V allele read against the IG-only database reached at most 12.
- Real reads with 0–30% mutation scored at least 132 in 99% of cases.

**Projection.**
- The best allele is aligned with `align_semiglobal`.
- The alignment is cut to the read's V-derived bases: its best-scoring stretch under the boundary score (match +1, mismatch −1, gaps as above), the 3' end at the last running-score peak and the 5' end at the first. This is the same ±1 score record 0025 uses for boundaries. It removes the alignment's forced extension through N/D bases to the germline's end.
- Gaps are then placed with `place_gaps` as above.
- Each aligned read base takes the allele's IMGT gapped index. Inserted bases take none.

**Refinement.** The projection may come from an allele whose CDR lengths differ from the true allele's, which misplaces bases. A self-consistency test catches this: if the cut alignment has a gap, or another allele has fewer mismatches than the projecting allele, the (up to 3) fewest-mismatch alleles not yet scored are scored. One with a higher DP score projects instead. There are at most 2 rounds. Typical reads never trigger it.

**Tie set.**
- Every V allele of the projecting allele's locus is counted, across all species in a mixed-species database.
- Its count is the mismatches over the V-derived bases at their IMGT positions:
  - inserted bases are skipped, and deleted positions have no read base;
  - a read base the allele doesn't cover is a mismatch: a gap `.`, or a position past its end;
  - a non-ACGT read base is a mismatch for every allele;
  - allele positions beyond the read are ignored.
- The primary is chosen among the alleles with the fewest mismatches by record 0008: F > ORF > P > NF > unknown, then higher `usage` (a missing usage after any given one), then natural sort of names. `mismatches` is that count.
- The tie set also adds every allele that has the fewest mismatches when one end of the V-derived bases moves by up to 3 bases. The end may move inward, or outward along the alignment's diagonal within the read and allele.
  - The read cannot show where V ends and N begins, nor whether a 5'-truncated read's first bases are V.
  - Cutting by the projecting allele alone drops end bases where it mismatches and siblings match, and keeps N bases that happen to match it.
  - On 5'-truncated memory reads this slack raised tie coverage from 98.4% to 99.8%, while the mean tie set grew from 1.40 to about 1.47 alleles.
- `ties` lists the primary first, then record 0008 order.
- The reported alignment and IMGT projection are the primary's own, cut and gap-placed as above.

## Alternatives considered

- **Prefilter margin** (the brief's first rule): align every allele of every gene within a k-mer margin of the best gene. No margin both bounds the work and guarantees completeness, because at 25–30% mutation the true allele's k-mer count (about 24) is within the noise of unrelated alleles (random 600-nt reads reach 30).
- **Align every allele:** about 200 × 104 µs per IGH read; rejected on cost (A20).
- **Summed k-mer orientation:** misorients heavily mutated reads (above).
- **Top 3 alleles instead of top 3 genes:** alleles of one gene often fill the top 3, so a better gene is never scored. Scoring up to 8 genes for every read raised coverage only from 1,843 to 1,849 of 1,871 and doubled the cost, so wider scoring is used only below the threshold.
- **Cutting the V-derived bases with the alignment score (+2/−3), earliest peak:** biased toward the projecting allele. On 15–35% mutated reads its primary accuracy was 93.9%, against 96.6% with the ±1 cut.
- **No cut:** N bases aligned to the germline's 3' end decide ties. Coverage on 5'-truncated reads was 97.5%.

## Consequences

Measured with the built-in `human` database (release build, one thread):

| Set (recoverable reads) | Primary in `v_call_valid` | `v_call_valid` ⊆ ties | No V hit |
|---|---|---|---|
| naive, seed 1 (2,000) | 100.00% | 99.95% | 0 |
| mutation 0–30%, seed 1 (1,871) | 99.20% | 99.57% | 0 |
| memory + 20% indels (2,913) | 99.83% | 100.00% | 0 |
| memory, 5' truncation 0–250 nt (2,886) | 99.34% | 99.79% | 2 |
| mutation 15–35% (2,628) | 96.65% | 96.77% | 46 |

- Ties are not complete by construction. The read alone cannot reproduce origin-based truth when:
  - an N base matches a sibling's 3' tail beyond the 3-base slack;
  - above about 20% mutation, the aligner's compensating insertion/deletion pairs misplace bases.
  The test `prefilter_never_drops_a_tie` therefore requires ≥ 99% coverage and prints the exact rate.
- Very heavily mutated reads (≥ 35%) can fall below the score threshold and give `no_v_hit` rather than a forced call.
- V step cost: about 0.28 ms per naive read, mostly the one traceback alignment; refinement and fallback add cost only on reads that need them.
- An allele whose IMGT gapping is shifted (#138), or that lacks gaps (macaque, #137), is compared at the wrong positions, as in the truth sets.
- Changing any constant here (k, scores, 90%, 3, 12, 50, ±1 cut, 3-base slack) changes calls, so it needs a superseding record.
