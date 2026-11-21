# 0029: Orientation, V search and tie sets by IMGT projection

## Status

Accepted (2026-11-20)

## Context

The annotator's first steps (spec *Pipeline* steps 2–3) orient the read, find the V gene and then the J gene, and report tie sets by the same fewest-mismatches rule the generator's truth uses (record 0025). Truth sets are computed from each read base's known origin; an annotator only sees the read.

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
- `place_gaps` then runs.
- **5′ start extension** (added 2026-11-21, ruling A28). The cut drops mutated first V bases: a prefix with a negative score, which is a single mutated first base 30% of the time at 30% mutation. So the V start then extends along the alignment's first diagonal to the read's start or the allele's start, whichever comes first. The extension is taken when it is at most 20 bases long and either at most 3 bases long or matching the allele on at least a quarter of its bases (as often as random sequence would).
  - A 5′ flank before the allele's start (leader, 5′ UTR, a primer upstream of IMGT 1) lies past the allele's start on the diagonal. The geometry protects it: the extension stops at the allele's start. The quarter-identity guard does not, since a leader matches the allele no worse than random sequence. When a deletion sits in the first ~5 V bases, though, the first diagonal is shifted by it and the extension can absorb 1–3 leader bases. (2026-11-21, final review: this bullet said the flank "lies past the allele's start on the diagonal, so the extension cannot absorb it".)
  - The extended bases are V-derived bases like the rest: they count in the tie set, as truth counts every V-origin base.
- Each aligned read base takes the allele's IMGT gapped index; inserted bases take none.

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

**J assignment** (`assign_j`, after a V hit; no V hit gives no J for now).
- Candidates are the J alleles of the V primary's locus, across all species in a mixed-species database. Each is scored with `score_semiglobal` against the read from 10 nt before the V alignment's end.
- The 3 best are aligned, gap pairs are collapsed, and the alignment's 5' end is cut to the J-derived bases. The cut keeps the maximum-scoring run of columns that ends at the alignment's 3' end, with the V cut's scores (+2/−2, gaps as above). Among equal scores the shortest run wins, so N bases that only break even are left out.
  - The 3' end is not cut: it is the J's own end or the read's end.
  - Ranking by the cut score rather than the DP score matters. The overlap alignment must start at the J's first base, so the DP charges each allele for its own 5'-trimmed bases.
- `place_gaps` runs with `protect` = the germline up to 5 nt past the first base of the Trp/Phe118 codon. This is the generator's J indel zone (A18).
- The best cut score becomes the projecting allele. Below 30 there is no J hit (`no_j_hit`).
- **Tie set (A21 rule, truth's sliding comparison).** The J-derived bases are taken as (projecting-allele index, read base). Inserted bases are skipped, and deleted positions have no read base. Every J allele of the locus is counted at its best shift, where read index `u` faces allele base `u + k`. A base the shifted allele does not cover, or a non-ACGT read base, is a mismatch. This is how short J siblings separate.
- The record 0008 pick among the fewest cuts the J-derived bases itself, and the count is repeated once over its cut. The tie set is exactly the fewest-mismatch alleles there, with no end slack. The primary is their record 0008 pick, and the reported alignment is the primary's own.

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

- **J: rank by DP score, gate the DP score at 30.** The forced 5' overlap gives N bases to alleles with long 5' ends (IGHJ1, IGHJ5 over IGHJ4). Naive tie coverage was 98.75%, and 11% of 0–30%-mutation reads had no J hit.
- **J: longest cut among equal scores** (the V rule): naive tie coverage 98.90%, against 99.15%.
- **J: among equal cut scores, the shortest projects:** within noise. Naive primary rose by 2 reads, bnab tie coverage fell by 15.
- **J: cut every locus allele instead of 3:** within noise (mutation 0–30%: 98.27% against 98.17% primary), at about 1.5 times the cost.
- **J floor on the cut score:** random 250-nt tails after a V (100,000 each, after IGHV3-23, IGKV3-20 and TRAV1-1) reach it at these rates:

  | Floor | Random tails reaching it (IGH / IGK / TRA) | No J hit on mutation 0–30% | No J hit on bnab |
  |---|---|---|---|
  | 20 | 194 / 630 / 1207 | 6 | 92 |
  | 30 | 3 / 14 / 18 | 11 | 118 |
  | 40 | 0 / 1 / 0 | 62 | 254 |

  The table was measured on an earlier build of the J step (before its final cut and refinement rules); the build these results describe gives no J hit on 7 reads at 0–30% and 107 on bnab at the floor of 30 (the J table under *Consequences*), not 11 and 118. (2026-11-21 label, final review.)

  Calls below 30 are mostly wrong: 4 of 10 at 0–30% mutation and 32 of 109 on bnab had the primary in `j_call_valid`. So 30 costs few correct calls and keeps random tails at about 10⁻⁴.

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

J assignment, same database and build. These rates are over reads whose `j_call_valid` is non-empty and holds the true allele. A read with no V hit counts as no J hit.

| Set | Primary in `j_call_valid` | `j_call_valid` ⊆ ties | No J hit |
|---|---|---|---|
| naive, seed 1 (2,000) | 99.85% | 99.15% | 0 |
| mutation 0–30%, seed 1 (1,968) | 98.17% | 95.83% | 7 |
| memory + 20% indels, seed 11 (2,995) | 99.40% | 98.20% | 1 |
| bnab preset, seed 1 (1,906) | 92.39% | 87.25% | 107 |

The J step costs about 71 µs per read.
- J tie misses have the same cause as V's: the read cannot show where J starts. IGHJ alleles share their 3' 30–35 nt and differ in the 5' bases next to the N region. An N base that happens to match one allele's 5' tail counts for it, and a cut that starts late drops the J bases that separate siblings. At 0–30% mutation the J start agrees with truth on 46% of reads and is off by up to 10 nt either way.
- On bnab reads, most of the loss is the floor: 107 reads (5.6%) score below 30 against a J mutated 15–35%. Among called reads, 97.9% have the primary in `j_call_valid`.

**5′ start extension** (ruling A28). The table compares the cut alone with the extension, over reads whose V and J calls are the true alleles (`tests/derive.rs`). `v_sequence_start` is read position agreement with truth.

| Set | `v_sequence_start` | FWR1 | `sequence_vdj` |
|---|---|---|---|
| naive, seed 1 | 100.00% → 100.00% | 100.00% → 100.00% | 100.00% → 100.00% |
| mutation 0–30%, seed 1 | 78.22% → 99.40% | 78.22% → 99.40% | 77.92% → 99.03% |
| memory + 20% indels, seed 11 | 90.21% → 99.68% | 89.92% → 99.28% | 90.03% → 99.50% |
| bnab, seed 1 | 62.41% → 97.00% | 62.41% → 96.77% | 62.34% → 96.77% |

- **Flanks** (`tests/assign_v.rs`, 0–60 random bases prepended, 0–30% mutation).
  - Before a V that starts at the allele's start: no read's V starts inside the flank (0 of 1,742; the test allows 0.5%).
  - Before a 5′-truncated V (0–60 nt cut), the worst case, the cut alone already keeps chance flank matches next to the V. It absorbs 1.4 flank bases per read (40% of reads); with the extension that is 3.4 (55%; the test allows a mean of 4).
  - An unguarded 20-base extension would absorb 6.6 bases per read here. A 10-base one would absorb 3.0, but V starts would agree on only 93% of bnab reads.
- **V calls**, over recoverable reads whose D and J are recoverable too (cut alone → extension):
  - primary in `v_call_valid` changes within noise: 0–30% 99.62 → 99.69%, memory with 5′ truncation 99.39 → 99.49%, bnab 99.50 → 99.59%;
  - `v_call_valid` ⊆ ties is also within noise: 0–30% 99.75 → 99.69%, memory with 5′ truncation 98.27 → 98.41%, bnab and memory + indels unchanged.
  - One 0–30% read lost coverage: a spurious 3-nt gap near the 5′ end left the first diagonal shifted, and the extension followed it.

Spurious gaps, measured as V alignments with a gap on indel-free reads, with `collapse_gap_pairs` off and on:

| Set | Off | On |
|---|---|---|
| mutation 0–30% | 1.87% | 0.43% |
| bnab | 8.14% | 1.11% |

On reads with a single true V indel where the primary is the true allele, the placement matches truth on 494 of 518 with collapse off and on. The remaining mismatches are indels within a few bases of a V end that the cut removes.

- **Ties are exact but not always truth's set.** The read alone cannot show where V ends. N bases that continue an allele's 3' tail count for that allele, while truth counts only V-origin bases. This dominates the misses on 5'-truncated reads, whose short V leaves few positions. Above about 25% mutation, the remaining compensating gaps and alignments to the wrong CDR length misplace bases.
- **Junk-read cost.** No-hit input costs about 1.3 ms per read, both strands included.
- **Typical cost:** 0.3 ms per naive read.
- **Gapping limits.** An allele whose IMGT gapping is shifted (#138), or that lacks gaps (macaque, #137), is compared at the wrong positions, as in the truth sets. (Amended 2026-11-21: a database without IMGT gaps is now refused at build and load, record 0036; the macaque case no longer arises.)
- **Constants.** Changing any constant here changes calls, so it needs a superseding record: k, scores, 90%, 3, 12, bin 8, 8 diagonal hits, J 30, near V 20, floor 50, the cut scores, the pair span of 40, the V start extension's 20 and 3 bases and quarter identity, or for J the 10-nt overlap, 3 cut alleles, the J floor of 30 and the J protect margin of 5.
- **2026-11-21 (fix wave, rulings B1 and B4; measured on the Task 12 benchmark sets, 5,000 reads each, seeds 20261201–09).**
  - **B1: short 3′ J fragments.** Below the J floor of 30, the best gapless alignment of the read's 3′ end to any J allele of the V's locus is called when all of these hold:
    - it has at most 2 mismatches;
    - it is at least 7 + ⌈log4(the locus's total J nt)⌉ bases long, plus 1 for TR loci (12 for IGH, IGK and IGL; 13 for TRB; 14 for TRA), plus 4 bases per mismatch;
    - the junction it implies (Cys104 on V's last diagonal through the Trp/Phe118 codon) is 6–180 nt; the bnab truth sets reach 162. Without a V or J anchor the junction cannot be bounded, so the fragment is not called.

    The tie set, the record 0008 pick and the re-cut then follow the usual steps, on that tail alignment. The joint junction step (record 0030) places the V end, D and J start of these reads as usual.
    - On igh_truncated (memory preset, 0–60 nt 5′ and 0–40 nt 3′ truncation), J call accuracy rose from 77.19% to 80.91% (`scripts/benchmark/short_j_null.py`), by J nt left:

      | J nt left | before | after | IgBLAST |
      |---|---|---|---|
      | < 10 (n = 415) | 0.0% | 0.0% | 39.3% |
      | 10–14 (n = 463) | 2.4% | 31.7% | 78.8% |
      | 15–19 (n = 510) | 61.2% | 70.4% | 91.6% |
      | ≥ 20 (n = 3,513) | 98.2% | 98.2% | 98.9% |

      D rose from 71.28% to 74.89% and the junction from 92.23% to 94.46%. IgBLAST reaches 91.1% on J by calling many more of the fragments under 15 nt; its J start reaches only 77.7%.
    - **Null.** `scripts/benchmark/short_j_null.py`: the full V allele plus 100,000 seeded random tails, per V.
      - With 250-nt tails (this record's null), J calls stay at today's rate for every locus. Per 100k: IGH 5, IGK 8, IGL 12, TRA 16, TRB 6. The junction cap removes every tail match from these reads; without it, an exact 12-nt rule raised TRA to 29.
      - With 1–150-nt tails, where a fragment can imply a plausible junction, the rule is held within +2 per 100k of the rate before B1 on every locus:

        | Locus | before | after |
        |---|---|---|
        | IGH | 2 | 2 |
        | IGK | 4 | 6 |
        | IGL | 0 | 2 |
        | TRA | 7 | 9 |
        | TRB | 1 | 3 |
    - **Rules swept**, measured with temporarily overridable constants: the five-locus short-tail null and the igh_truncated buckets (10–14 / 15–19 nt).

      | Rule | igh_truncated J | Short-tail null IGH / IGK / IGL / TRA / TRB |
      |---|---|---|
      | shortest exact 10 on every locus | 84.30% | 20 / 18 / – / 234 / – |
      | shortest exact 11 on every locus | 82.54% | 2 / 7 / – / 67 / – |
      | 12 / 12 (TRB) / 13 (TRA), the first B1 commit | 80.91% (31.7 / 70.4) | 2 / 6 / 2 / 15 / 8 |
      | **TR + 1: 13 (TRB) / 14 (TRA) (chosen)** | **80.91% (31.7 / 70.4)** | **2 / 6 / 2 / 9 / 3** |
      | TR + 2 | 80.91% | 2 / 6 / 2 / 8 / 3 |
      | 5 nt per mismatch | 79.86% (30.9 / 63.5) | 2 / 6 / 2 / 15 / 8 |

      TR + 1 is the least change that meets the +2 bound; the IGH truncated gain is untouched because IGH's minimum is unchanged.
  - **B4: later J start among equal cuts.** When the projecting J cut has at most 3 mismatched or gapped bases, the J alleles are counted only over read positions from the latest start among the top-3 cuts with the same score. Chance N bases matching a sibling's 5′ end, such as IGLJ3\*02's against IGLJ2\*01, therefore don't decide the allele.

    J call, before → after:

    | Set | J call | Tie coverage |
    |---|---|---|
    | igl | 98.70% → 99.82% | 98.70% → 99.82% |
    | igh_naive | 99.76% → 99.84% | 98.86% → 98.92% |
    | igh_memory | 99.56% → 99.62% | |
    | igh_bnab | 93.65% → 93.63% | |
    | igh_truncated | 80.89% → 80.91% | |
    | igk, tra, trb, igh_indels | unchanged | |

    It is a net gain, so it is kept. Variants that lost:
    - Counting from the latest start of all three top cuts: IGL 90.7%, igh_memory 93.5%, igh_bnab 74.4%.
    - Equal-score cuts on any J: igh_bnab −0.65 pp.
    - Cuts within 4 of the best score: IGK 94.3%, IGL 93.3%.

    The cut's own tie-break was already the shortest J among equal cut scores (`trim` keeps the latest minimal prefix), so that form of the ruling changes nothing.
  - **New constants:**
    - B1: the 7 + ⌈log4⌉ (+1 for TR) shortest fragment, 4 nt per mismatch, at most 2 mismatches, the 6–180-nt junction;
    - B4: the 3-difference limit.
- **2026-11-21 (Wave 2, user ruling: J may accept weaker and shorter evidence; settings are chosen by generated-truth accuracy).** J is searched only after a V call, so random non-antibody reads cannot get a J. The random-tail null (`scripts/benchmark/short_j_null.py`) is reported for information and **no longer gates the choice**; this supersedes the "+2 per 100k" bound of B1 above.
  - **Old rule (kept for reference).** Below the J floor of 30 the only fallback was B1: a gapless 3′ fragment with "at most 2 mismatches" and "at least 7 + ⌈log4(the locus's total J nt)⌉ bases long, plus 1 for TR loci (12 for IGH, IGK and IGL; 13 for TRB; 14 for TRA), plus 4 bases per mismatch", implying a 6–180-nt junction.
  - **New rule.** The floor of 30 on the cut score is unchanged. Below it, two fallbacks are tried in order:
    1. **Weak J.** Every J allele of the V's locus is realigned with mismatch −1 (match +2, gaps −10/−1 as before) and cut as usual. The best cut that holds its allele's Trp/Phe118 codon and implies a 6–180-nt junction projects when its cut score (the ±2 scale) is at least 16. This recovers full-length J segments mutated 30–40%, which the −3 mismatch alignment cuts short. The B4 later-start count does not apply to it.
    2. **Short 3′ fragment (B1, amended).** The best gapless alignment of the read's 3′ end is called when its matches minus mismatches reach 5 + ⌈log4(the locus's total J nt)⌉, plus 1 for TR loci: 10 for IGH, IGK and IGL, 11 for TRB, 12 for TRA. The mismatch cap and the per-mismatch length are gone. The 6–180-nt junction bound is kept.
  - **Selection.** Measured with temporarily overridable constants on the benchmark's 12 generated sets (5,000 reads each: the nine Task 12 sets plus `igh_memory_c`, `igh_bnab_c`, `igk_c`), graded with `abstar.grade` on `j_call`, `j_sequence_start`, `junction`, `junction_aa`, `cdr3`, `fwr4`, `productive`, `vj_in_frame` and `d_call`. Objective: the summed `j_call` + `junction` gain over all sets. Constraint: no set's `all` cell for those fields drops by more than 0.1 pp. Selected rows (Σ gain in pp; worst drop; J / junction in %; truncated J by J nt left < 10 / 10–14 / 15–19 / ≥ 20):

    | Setting | Σ gain | Worst drop | igh_bnab | igh_bnab_c | igh_truncated | truncated J by nt left |
    |---|---|---|---|---|---|---|
    | old rule | 0 | – | 93.63 / 94.50 | 93.23 / 94.73 | 80.91 / 94.46 | 0.0 / 31.7 / 70.4 / 98.2 |
    | floor 22 (no other change) | +11.44 | −0.20 bnab_c fwr4 | 94.46 / 95.33 | 94.08 / 95.49 | 86.20 / 97.27 | 0.2 / 53.3 / 92.9 / 99.5 |
    | floor 14 | +19.02 | −0.23 bnab_c fwr4 | 94.55 / 95.38 | 94.17 / 95.51 | 92.37 / 98.45 | 29.2 / 86.6 / 98.2 / 99.7 |
    | floor 1 | +22.71 | −0.69 truncated cdr3 | 94.84 / 95.38 | 94.48 / 95.51 | 95.44 / 98.47 | 71.1 / 92.7 / 98.8 / 99.7 |
    | weak ≥ 16, mismatch −2 (anchor, junction) | +23.01 | −0.13 bnab_c fwr4 | 96.98 / 98.18 | 97.00 / 98.48 | 85.25 / 98.53 | 5.3 / 49.9 / 87.8 / 98.8 |
    | weak ≥ 16, mismatch −1 (anchor, junction) | +23.28 | −0.097 bnab_c fwr4 | 97.04 / 98.21 | 97.10 / 98.58 | 85.25 / 98.51 | 5.3 / 49.9 / 87.6 / 98.8 |
    | weak ≥ 12, mismatch −1 | +24.01 | −0.22 truncated cdr3 | 97.06 / 98.23 | 97.17 / 98.60 | 85.85 / 98.51 | 10.8 / 51.6 / 88.0 / 98.8 |
    | tail score ≥ 10 only (no weak J) | +17.92 | −0.021 | 96.84 / 97.94 | 93.23 / 94.73 | 88.44 / 98.16 | 3.1 / 70.0 / 96.1 / 99.7 |
    | **weak ≥ 16, −1 (anchor, junction) + tail score ≥ 10 (chosen)** | **+27.15** | **−0.099 truncated fwr4** | **97.06 / 98.21** | **97.10 / 98.58** | **89.06 / 98.55** | **7.2 / 73.2 / 96.3 / 99.6** |
    | same, tail score ≥ 9 | +28.83 | −0.15 truncated fwr4 | 97.06 / 98.21 | 97.10 / 98.58 | 90.72 / 98.57 | 16.1 / 81.2 / 97.5 / 99.7 |
    | same, weak J without the anchor guard | +30.08 | −0.14 truncated fwr4 | 97.06 / 98.23 | 97.08 / 98.55 | 91.93 / 98.64 | 21.0 / 88.1 / 98.6 / 99.7 |
    | same, weak ≥ 20 | +25.03 | −0.070 bnab_c fwr4 | 96.88 / 98.01 | 96.77 / 98.27 | 88.42 / 98.10 | 3.1 / 70.0 / 95.9 / 99.7 |
    | weak ≥ 16 (−1) + old tail + floor 18 for read-end, gapless cuts inside J | +28.20 | −0.097 truncated fwr4 | 97.04 / 98.21 | 97.10 / 98.58 | 90.09 / 98.59 | 13.5 / 78.4 / 96.9 / 99.6 |

    - The sweep had 66 rows. Besides the rows above, it covered:
      - **Floors 26 to 1 with no other change:** Σ +6.3 to +22.7, worst drop −0.11 to −0.69.
      - **Floors 6–26 with junction, read-end, gapless and order guards:** Σ +16.1 to +32.3. Only "floor 18, gapless, read end inside J, after the weak J" stays within 0.1 pp.
      - **Weak thresholds 8–24:** Σ +15.5 to +24.3. Only 16 or more stays within 0.1 pp.
      - **Tail rules 9–12 nt with 2–4 nt per mismatch:** Σ +24.1 to +28.5. Rules with fewer nt than the chosen score-10 rule break the bound.
      - **Weak J with mismatch −2, and without its guards:** worse on the worst drop.
    - **Why lowering the floor alone fails.** On bnab the −3 mismatch alignment of a heavily mutated J stops early, so lower floors mostly add wrong alleles or 1–18-nt cuts that miss Trp118: at floor 1, 229 bnab reads gain a J but only 50 gain the right junction.
    - **The binding cells are `fwr4` and `cdr3`.** They are graded only where the J allele is truth's, so each newly right J adds a graded read. On igh_truncated `fwr4`, correct reads rise from 3,526 to 3,858 while wrong ones rise from 2 to 6. The accuracy cell falls 0.099 pp, although 332 more reads are right. Every setting with more truncated calls than the chosen one (tail ≥ 9, the unguarded weak J, the floor-18 path) breaks the 0.1-pp bound through these cells.
    - **Choice.** The floor-18 path meets the bound and gains 1.05 pp more than the chosen rule, nearly all of it truncated J (90.09 vs 89.06%). It needs a third acceptance path with four guards, so the simpler rule is kept (the brief: prefer the simpler of similar settings).
    - **The junction guard** changes no benchmark cell. It halves the 250-nt null (below), so it is kept.
    - The 54 igh_truncated reads with no recoverable J get one J call (none before).
  - **Null (for information; per 100,000 reads, 250-nt tails / 1–150-nt tails).**

    | Locus | old rule | chosen | chosen without the weak J's junction guard | tail score ≥ 10 only | weak J only (old tail) |
    |---|---|---|---|---|---|
    | IGH | 5 / 2 | 690 / 1,226 | 1,678 / 1,226 | 5 / 294 | 1,678 / 1,094 |
    | IGK | 8 / 6 | 995 / 1,521 | 1,829 / 1,521 | 8 / 192 | 1,829 / 1,451 |
    | IGL | 12 / 2 | 881 / 1,119 | 1,983 / 1,119 | 12 / 180 | 1,983 / 1,011 |
    | TRA | 16 / 9 | 6,073 / 8,038 | 12,782 / 8,038 | 16 / 364 | 12,782 / 7,886 |
    | TRB | 6 / 3 | 2,293 / 2,675 | 4,702 / 2,675 | 6 / 189 | 4,702 / 2,595 |

    Most false calls come from the weak J. A random tail after a V gets a J on 0.7–1.5% of IG reads and up to 8% of TRA reads (61 TRAJ alleles). The requirement that the weak cut reach the J's or the read's 3′ end changed neither accuracy nor the null, so it is not used.
  - **New constants:** the weak J's mismatch −1 and score 16; the fragment score 5 + ⌈log4⌉ (+1 for TR). Removed: B1's 2-mismatch cap and 4 nt per mismatch.
- **2026-11-21 (Wave 2, review fix round 1).** The review found that the weak J's guard was checked only on the projecting cut, and that the weak J pre-empted B1 even when B1 had a long exact match.
  - **Guard on the reported alignment.** Each fallback's guard (weak J: score ≥ 16, Trp/Phe118 codon held, junction 6–180; B1: junction 6–180) is re-applied to the alignment reported after the record 0008 pick and re-cut. A candidate that fails is dropped, and the other fallback is tried. Example: bnab_c g3609's weak cut held IGHJ4\*02's Trp118, but its re-cut on the picked IGHJ2\*01 started past it; the read now has no J.
  - **Order of the fallbacks.** Both candidates are computed. B1 goes first when its doubled score (the ±2 scale) is at least the weak cut's score; otherwise the weak J goes first.
    - Example: tra g0440's exact 47-nt TRAJ8\*01 tail (B1 score 47, 94 on the ±2 scale) is called as that fragment again. Before this fix, the weak J stretched J's 5′ end through a 4-nt insertion and 18 mismatches, and the joint step then started J past Phe118.
    - Preferring an anchor-holding weak J over an anchorless B1 fragment was tried. It lost 0.2 pp, so it is not used.
  - **Weak prefilter.** The locus alleles are ranked by score-only `score_semiglobal` with the weak scoring, and only the best `WEAK_TOP` = 8 are cut.
    - Top 8 matches cutting every allele on all 12 sets; top 3 cost 0.2 pp of summed J + junction.
    - On 30,000 V + random 250-nt tails the run takes 2.29 s, against 2.34 s for Wave 2's every-allele cut.
  - **Joint step (record 0030).** When the J alignment holds its Trp/Phe118 codon and the read has TGG, TTT or TTC there, J may not start after it.
  - **Effect (Σ gain over the old rule; the 12 sets):** +27.15 → +27.83 pp, worst `all` drop −0.099 → −0.047 (igh_truncated `fwr4`, still dilution).
    - igh_bnab J / junction: 97.06 / 98.21 → 97.09 / 98.25.
    - igh_bnab_c: 97.10 / 98.58 → 97.12 / 98.69.
    - igh_truncated: 89.06 / 98.55 → 89.20 / 98.70.
    - The igh_bnab_c `fwr4` dilution (−0.097) shrinks to −0.009.
  - **Null after the fix** (per 100k, 250-nt / 1–150-nt tails): IGH 673 / 1,214; IGK 987 / 1,515; IGL 871 / 1,105; TRA 5,430 / 6,281; TRB 2,243 / 2,579.
  - **Accuracy-gate cost.** One read, igh_bnab_c g1091: its weak J was the wrong allele but bounded the right junction. Its reported re-cut now fails the guard, so it has no J, C or junction. That lowered 15 floors (C call, J start, junction, productive and vj_in_frame cells that hold it).
  - **New constant:** `WEAK_TOP` = 8.
- **2026-11-21 (Wave 3; user rulings).**
  - **Rulings.**
    - False J calls on junk after a V are accepted: false positives are preferred over discarding real recombinations.
    - Settings lean permissive: where two are close, take the more permissive.
    - The 0.1-pp no-drop bound applies only to real losses. A cell graded only where the J allele is right (`fwr4`, `cdr3`, J mutation fields, …) that falls only because more reads are graded is dilution, and does not count. In a dilution cell, the number of correct reads does not fall.
  - **Old values (Wave 2 fix round 1):** weak J at "score 16"; B1 "10 for IGH, IGK and IGL, 11 for TRB, 12 for TRA" (5 + ⌈log4⌉, +1 TR).
  - **New values:** weak J at score **12**; B1 matches minus mismatches ≥ 3 + ⌈log4⌉ (+1 TR), i.e. **8** for IGH, IGK and IGL, **9** for TRB, **10** for TRA. Nothing else changed: order, guards, re-check after the record 0008 pick, `WEAK_TOP`, and the joint-step cap.
  - **The third acceptance path was swept and is not added.** The path is a strong-scoring (−3) cut below the floor of 30, tried after the weak J and B1. It must score at least N, imply a 6–180-nt junction, and optionally be gapless and end at the read's end inside J. On top of fix round 1 it gains only below the B1 threshold's doubled score.
    - The binding test is `missing_j_is_none`: a read ending in N/D bases gets no J at any truncation point. Floor 14 already calls IGHJ2\*01 on the full NP read, whose last 7 nt (`CTGGTAC`) are an exact IGHJ2 stretch.
    - Floors 16 and up keep that test green, but once B1 is at score 8 they add no call. A low path at 16 changed none of the 5,000 igh_truncated reads, none of the 12 sets and none of 6,000 synthetic mutated or indel-bearing J fragments.
    - Its gapless and read-end guards never changed a result in any setting tried.
    - Lowering B1 to 8 captures what the path gave: floor 16 with B1 at 10 reaches +30.87; B1 at 8 reaches +31.75.
  - **Sweep** (Σ = summed `j_call` + `junction` gain over the pre-Wave-2 build on the 12 benchmark sets, pp; "real loss" = a listed `all` cell with fewer correct reads, against the pre-Wave-2 build and against fix round 1; "N/D" = the N/D-ending test):

    | Setting | Σ | Worst real loss vs fix round 1 | Truncated J / junction | J on the 54 no-J truncated reads | N/D test |
    |---|---|---|---|---|---|
    | fix round 1 (weak 16, B1 10) | +27.83 | – | 89.20 / 98.70 | 1 | green |
    | + low floor 18 | +29.20 | none | 90.51 / 98.76 | 1 | green |
    | + low floor 16 | +30.35 | none | 91.66 / 98.76 | 1 | green |
    | + low floor 14 | +31.58 | none | 92.82 / 98.84 | 2 | fails (NP 33 nt) |
    | + low floor 10 | +33.08 | none | 94.29 / 98.86 | 8 | fails |
    | + low floor 2 | +33.99 | none | 95.20 / 98.86 | 29 | fails (8 NP prefixes) |
    | weak 12, B1 10, low 16 | +30.87 | none | 92.05 / 98.72 | 4 | green |
    | B1 8 | +31.13 | bnab_c `j_call` −0.083 (4,662 → 4,658 of 4,800) | 92.53 / 98.80 | 1 | green |
    | **weak 12, B1 8 (chosen)** | **+31.75** | **bnab_c `junction` −0.047 (4,157 → 4,155 of 4,212)** | **92.98 / 98.84** | **4** | **green** |
    | same + low floor 16 or 18 (any guards) | +31.75 | same | same | 4 | green |
    | weak 10, B1 9, low 2 | +34.21 | truncated `junction` −0.165 (4,774 → 4,766 of 4,837) | 95.40 / 98.53 | 32 | fails (low 2) |
    | B1 7, low 2 | +34.81 | bnab_c `j_call` −0.083 | 96.06 / 98.95 | 31 | fails |

    - No chosen-setting cell loses correct reads against the pre-Wave-2 build.
    - Against fix round 1 the only real losses are 2 bnab_c junctions and 4 bnab_c J calls, each under 0.1 pp. Both come from the lower B1 threshold taking a fragment before the weak J.
    - Dilution: igh_truncated `cdr3` 99.55% (3,512 of 3,528) → 99.38% (4,003 of 4,028).
    - Weak 10 loses 8 truncated junctions (−0.165 pp) against fix round 1, so 12 is the most permissive weak threshold within the bound.
    - **Marginal precision of the low path** (new igh_truncated calls right / wrong allele / on no-J reads, B1 at 10):

      | floor | right | wrong allele | on no-J reads |
      |---|---|---|---|
      | 18 | 66 | 1 | 0 |
      | 14 | 117 | 5 | 1 |
      | 10 | 80 | 13 | 6 |
      | 6 | 16 | 9 | 1 |
      | 2 | 17 | 15 | 8 |

      Below floor 6 the added calls are more often wrong than right, and the J is 1–3 nt.
  - **Accuracy gate:** 32 floors rose and 24 fell. All 24 are igh_bnab J mutation cells that hold one more correct read (dilution).
  - **Null** (per 100k, 250-nt / 1–150-nt tails): IGH 2,153 / 5,941; IGK 3,527 / 7,616; IGL 3,119 / 6,042; TRA 14,555 / 19,346; TRB 7,349 / 10,539. The ruling accepts it.

