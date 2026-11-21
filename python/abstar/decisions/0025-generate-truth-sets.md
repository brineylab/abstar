# 0025: Truth sets for generated reads: fewest mismatches, joint boundary ties

## Status

Accepted (2026-11-20)

## Context

`abstar generate` writes a truth record for every synthetic read so annotators, the rewritten abstar first, can be graded against it. The generative value is not always the right answer. Sibling alleles can be identical over the read, mutation can make another allele closer, and N bases can match germline. Grading against the generative value alone would penalise answers no annotator could avoid.

The generator records the origin of every read base (`Origin`: a gapped V index, an ungapped D or J index, N1, N2, P, a C index, or an inserted base). So the equally valid answers can be computed exactly, without aligning anything.

## Decision

**Fewest mismatches.** For each ambiguous field, the correct answer is the set of answers tied for the fewest mismatches against the observed read. The true generative values (`v_call`, `v_sequence_end`, …) are always recorded alongside these sets.

**Allele sets.** Each set holds the sorted allele names of every candidate in the database for the locus and species, anchored or not, that ties for the fewest mismatches over the read bases derived from that segment. Mutated bases count. N, P and inserted bases never count for any segment.
- **V:** each read base of origin `V(g)` is compared with the candidate's IMGT-gapped sequence at `g`.
  - Inserted bases have no V origin, so they are skipped.
  - Bases removed by a deletion are not in the read, so no read base maps into the deleted span.
- **D and J: sliding comparison.** The derived bases are compared with every candidate at every shift that covers at least one of them, and a candidate's score is its best shift.
  - J bases are indexed by their true J index. Read base `J(u)` is compared with `candidate[u + k]`, so an indel in J is skipped as in V.
  - D bases are indexed by their place in the D run.
  - A candidate's reverse complement is also compared, but only when the model's `d_inversion` is above 0 (`Candidates::d_inversion`). This follows the spec ("where inversion is enabled"). A model that never inverts D does not make a forward-only annotator compete with chance reverse-complement matches.
  - **Amended 2026-11-21 (ruling A24):** the reverse complement is now always compared, whatever the model's `d_inversion`, and `Candidates::d_inversion` is removed. `d_inversion` governs generation only. The valid set is what the read supports, and an annotator cannot know that a model forbids inversion. Record 0031 sets the annotator's side.
- **C:** a prefix comparison. Read base `C(u)` is compared with `candidate[u]`.
- **Coverage rule.** A derived read base that the candidate does not cover counts as a mismatch. That includes a gap `.`, a base past either end, or a shift that leaves it uncovered. A partial allele cannot explain bases it lacks. Candidate positions outside the read are ignored, which handles truncation.
- **Absent segments.** A segment with no derived bases in the read has an empty set and empty boundary lists, and is never unrecoverable. For example, a V lost to truncation, a D trimmed away, or no C appended.

**Unrecoverable.** `unrecoverable_v`, `unrecoverable_d` and `unrecoverable_j` are set when the true allele is not in its non-empty set: mutation made another allele strictly closer. No annotator can get these right, and they count toward the reported error floor.

**No D.** `d_null_valid` is true when the read's D-derived run is shorter than **5 nt**, and always true for VJ loci. The D set is still computed for short runs, so a D call from the set and a null call are both valid. For a short run, the V ends and J starts of the best no-D configurations are also added to the boundary lists (see below), so an annotator that legitimately omits D is graded fairly.

**Boundaries.** `v_sequence_end`, `d_sequence_start`, `d_sequence_end` and `j_sequence_start` are found by a joint maximum of a separable score.
- **Score.** In a configuration, each read base given to V, D or J scores +1 if it matches that segment's placed germline and −1 if it does not. Uncovered bases score −1. Bases left to N/P score 0, and inserted bases inside V or J score 0.
- **Placements.** Each placement fixes where each read position falls in the allele. Every allele in the segment's valid set is placed:
  - **V** by its true IMGT mapping. Past the last V-derived base, V continues with the allele's own ungapped bases after that gapped position, for as long as the allele has bases.
  - **D** at each tied best shift, in each tied orientation.
  - **J** at each tied best shift of the true J index (shift 0 is the true mapping). J can extend 5' with the allele's bases before the first J-derived base, for as long as the allele and the read have them.

  The spec says V and J are placed by the true mapping. For J, a sibling of a different length, such as one missing 5' bases, needs its own shift to line up. Using its best shift places the allele exactly as it was compared for the set.
- **Configurations.** V covers the read from its start to the V end. J covers the J start through the last J-derived base. Cuts are ordered: V end ≤ D start < D end ≤ J start. Each segment present in the read gets at least one base. Without a D (VJ loci, or a D run of length 0), the D terms are dropped.
- **Ties per allele combination.** Ties are computed separately for every combination of valid alleles: v ∈ `v_valid`, d ∈ `d_valid`, j ∈ `j_valid`. Each combination keeps the configurations tying its own best score over its placements (D and J shifts and orientations). No D is one more combination per (v, j) when the D run is under 5 nt, contributing V ends and J starts only. This follows the spec ("ties are computed for every allele in that segment's valid set, and their union is taken"; rulings G7 and G8, 2026-11-20).
  - For example, a V sibling missing the true V's last 3 bases cannot claim an N base that matches the true V's next germline base. So its own best end, the true end, is valid alongside the true allele's extended end.
  - A single global maximum across combinations would drop it, and would grade an annotator that called the short sibling, a valid V call, wrong on its V end.
- **Ties stored as lists.** Each boundary stores every value it takes in any configuration kept above, as a sorted list. It is a list, not a range, because tied cuts need not be contiguous. Values are 1-based positions in the final read, per AIRR. Each boundary is graded on its own.
- **Computation.** For each combination, a DP over prefix maxima of the V score and suffix maxima of the J score is run, with one linear pass per D placement. Tied configurations are recovered only for the placements that reach the combination's best.
- **How these cases come out under this score:**
  - An N base that matches V's next germline base moves the V end onto it: +1 beats 0. The true end stays valid only if the neighbouring segment can claim that base equally well. For example, it may also match J's 5' extension.
  - A mutated base at a segment edge is given to N. Its end moves inward by one: 0 beats −1.

**Indels.** Allele sets and boundaries use the true indel placement. Each indel records its range of equivalent placements (`equiv`, VCF-style normalisation on the unmutated sequence), and an annotator's placement anywhere in that range counts as correct. An annotator that picks a different but equally good placement could, rarely, rank alleles differently. This is accepted.

**Generator conventions the sets depend on** (Task 5):
- Indels are placed only in V upstream of Cys104 and in J downstream of Trp/Phe118, at least 3 nt from either codon, with germline on both sides. So they never touch the junction, the CDR3 or the boundary region.
- When `d_inverted`, `Origin::D(i)`, `d_trim_5` and `d_trim_3` index the reverse-complemented D.
- Trims are clamped to the segment's length, so a D can be trimmed to nothing (no D run).
- The C prefix is appended before indels and truncation, and is never mutated.

**Known limit (#138).** V candidates are compared at their IMGT-gapped positions. A V allele whose gapping is shifted, such as human IGLV3-2\*03, is compared at the wrong positions. It is a candidate in the set but may be under-counted. Fixing #138 fixes the sets. Macaque is excluded until its V genes are gapped (#137).

## Alternatives considered

- **Grade against the generative value only:** penalises answers no annotator could distinguish, such as identical siblings, alleles a mutation made closer, and N bases that match germline.
- **Align the read to every candidate:** slower and circular. The truth would depend on an aligner, the very thing being graded. The origin map makes alignment unnecessary.
- **Count mismatches alone for boundaries:** degenerate. N bases are free, so giving every base to N scores 0.
- **Count the true value as valid in addition to the tied maxima:** it would grade an unobservable end (N bases that match germline) as correct. The true value is recorded separately for anyone who wants it.
- **Always compare D reverse complements:** with short D runs, chance reverse-complement matches would mark forward-only annotators wrong on models that never invert D. (Adopted on 2026-11-21 by ruling A24; see the amendment above.)
- **Place J siblings by the true J index:** misplaces siblings of different lengths, and disagrees with the shift used to admit them to the set.
- **Store boundary sets as `[min, max]`:** tied cuts need not be contiguous.
- **One global maximum across all allele combinations:** drops the best cuts of valid alleles that score lower overall, such as shorter siblings (rejected in ruling G7).
- **Boundaries for a short D only with the D placed:** grades an annotator that legitimately calls no D against cuts that assume a D (rejected in ruling G8).

## Consequences

- Grading is per field: an answer is correct if it is in its `_valid` set or list.
- The unrecoverable rows give the error floor that a perfect annotator still hits.
- With `d_null_valid`, a null D call is graded correct, and the V-end and J-start lists include the cuts of the best no-D configurations as well as those with the D placed.
- Changing the 5 nt threshold, the score, or any comparison rule changes truth for the same seed, so it needs a superseding record.
- **Cost profile.** Candidate sequences are cached once per run (`Candidates::new`).
  - Per read, the V set compares the read's V bases (placed at their IMGT positions) with every V candidate in one vectorized pass.
  - D and J sets count matches at all shifts at once: candidates × run length × allele length.
  - Boundaries cost (distinct V score profiles inside the D windows) × (J placements) × (D placements) × D length, plus a linear pass per V allele.
    - V alleles that agree inside the D windows share one pass. This is common in large, ambiguous V sets from 5'-truncated reads.
    - Short D runs (1–3 nt) tie at many shifts across many D alleles in both orientations, so they are the expensive case.
  - Measured on the human IGH candidate set (198 V, 31 D, 7 J, 88 C, inversion on), single-threaded, release build:

    | Read | Time |
    |---|---|
    | Typical (15 nt D) | about 37 µs |
    | 3 nt D | 45 µs |
    | 1 nt D | 80 µs |
    | 1 nt D, 5'-truncated: IGHV1-69\*01 with 60 V nt kept (V set of 19) | 160 µs |
    | 1 nt D, 5'-truncated: IGHV1-69\*01 with 25 V nt kept (V set of 63) | 260 µs |
    | 1 nt D, 5'-truncated: IGHV3-30\*01 with 10 V nt kept (V set of 114, the worst case found) | 310 µs |

    The test is `truth_speed_human_igh`, ignored by default.
  - The V comparison is the largest fixed cost, so it grows linearly with the number of V alleles in the database.
