# 0030: Joint junction scoring for V end, D and J start, and C assignment

## Status

Accepted (2026-11-20)

## Context

After the V and J calls (record 0029) the annotator must place the V end, the D and the J start (spec *Pipeline* step 5). Truth (record 0025) defines these boundaries by a joint maximum of a separable score, computed from each read base's known origin. The annotator sees only the read, so it reproduces the same score with its own placements of the called alleles. The V and J tie sets are not changed by this step.

The generator's truth code is not a library dependency of the annotator (record 0014), so the rules are reimplemented here and mirror `crates/generate/src/truth.rs`.

## Decision

**Score** (the same as truth's).
- A configuration gives V the read up to the V end, optionally a D placement, and J the read from the J start.
- Cuts are ordered: V end ≤ D start < D end ≤ J start. V and J each get at least one base.
- Each read base given to V, D or J scores +1 if it is ACGT and equals that segment's placed germline base, and −1 otherwise. A germline position that does not exist is a −1.
- Bases left to N/P score 0. Bases inserted in V or J (alignment `Ins`) score 0.

**Placements.**
- **V:** the primary allele, through the V alignment's ops. Past the alignment's end, V continues on its last diagonal while the allele has bases. Truth continues V past its last V-derived base in the same way.
- **J:** the primary allele, through the J alignment's ops. Before the alignment's start, J continues on its first diagonal while the allele and the read have bases, as truth's J extends 5'.
- **D:** every D allele of the V primary's locus, forward and reverse-complemented (record 0031), at every offset.
- So the V end can move back from the alignment's end (a mutated edge base goes to N) or forward (an N base matching V's next germline base). The J start can move either way in the same way.

**Search range.**
- V ends cover every position V is placed at: from one base past the alignment's start to the allele's end. J starts cover every position J is placed at: from the read's or allele's start to the alignment's last base. This is the whole range truth allows for the placement, so nothing is approximated.
- D placements are restricted exactly, with no loss of the optimum or of its ties. Let `pa[c]` be the best V score with V end ≤ c, `sb[c]` the best J score with J start ≥ c, and `S0` the best no-D score.
  - A D of length L starting or ending at cut c gives a configuration of at most `pa[c] + L + sb[c]`.
  - So a D configuration reaching `S0` lies between the first and last cuts where that bound reaches `S0`. Only offsets overlapping that zone are scored, each clipped to it.

**Computation.**
- Prefix maxima of the V scores and suffix maxima of the J scores, with one linear pass per D placement (truth's DP).
- The tie-break below is then applied by enumerating the D starts and ends of the primary D's placements that reach the best score.

**D call.**
- The D orientation is chosen first, with the inverted-D margin of record 0031: inverted only when the best reverse-complement score beats the best forward score by at least 2.
- `Sd` is the best score in that orientation over every D allele and placement.
- The D tie set holds the alleles reaching `Sd` in that orientation. The primary is their record 0008 pick.
- No D (`d` is `None`) when:
  - the locus has no D (`locus_has_d` false: light chains, TRA, TRG);
  - `Sd ≤ S0`: no D configuration beats the best without a D;
  - or the chosen configuration's D run is shorter than **5 nt** (record 0025's threshold).
- With no D, the boundaries are the best no-D configuration's.

**Boundaries.**
- They come from the primary combination's best configuration: the primary V, the primary D (when called) and the primary J.
- Ties break toward the longest V, then the longest J, then the earliest D start, then the longest D.
- Truth keeps every tied configuration, so for the primary combination these choices are in truth's lists whenever the primaries and their placements are truth's.

**C.**
- Each C allele of the J primary's locus is compared base for base with the read from the J alignment's end: `allele[u]` against read base `end + u`.
- The tie set is the alleles with the fewest mismatches over all bases after J. As in truth, a base the allele does not cover, or a non-ACGT read base, is a mismatch. The primary is the record 0008 pick.
- No C call when fewer than **10 nt** follow J.
- `CHit` covers the primary's bases `0..g_end`, so `c_germline_start` is 1 (ruling A14).

**Interface.** `resolve_junction(read, v, j, alleles, d, locus_has_d)` and `assign_c(read, j, alleles, c)`.
- One database allele slice replaces the brief's separate V, D and J (and C) slices, since every hit's indices address the database slice.
- `d` and `c` are the D and C `SegmentIndex`, so uppercase sequences are prepared once per run, as for V and J.
- `assign_c` takes no `Junction`, because the J end is not moved by the junction step.

## Alternatives considered

- **A mismatch penalty of 2 or 3 inside D for choosing the D**, with the boundaries still scored ±1 (measured before ruling A24, both orientations, no margin):
  - On naive reads it helps. Correct D-and-boundary reads rose from 1,856 to 1,889 of 2,000 at 3.
  - On reads with 0–30% mutation it hurts: 1,272 fell to 1,002 of 1,621.
  - A penalty tied to the read's V mutation level could get both. It departs from truth's ±1 score and is left for a later record.
- **D search only between the V and J alignments' ends:** misses D bases that the alignments' trims gave to V or J. The exact zone above costs little more.
- **Null D on the longest tied D run instead of the chosen one:** no naive read was called null with a D expected, so the simpler rule is kept.
- **A single global maximum across allele combinations** (truth's G7 issue): the annotator scores one combination, the primaries, so it does not arise.

## Consequences

Measured with the built-in `human` database (release build), with truth comparing both D orientations (record 0025 as amended by ruling A24) and the inverted-D margin of 2 (record 0031). Rates are over reads where V, J and D are each absent or have the true allele in their valid set.
- A boundary is correct when it is in its `_valid` list. A null D's boundaries, like the null call itself, are correct iff `d_call_null_valid`.
- The D call is correct when its primary is in `d_call_valid`, or it is null and `d_call_null_valid`.
- "D placed right" narrows that: a null call with `d_call_null_valid`, or a called D in `d_call_valid`, in the true orientation, with a true D run of at least 5 nt. A shorter run puts every allele containing its few bases in `d_call_valid`, wherever the call is placed.

| Set | V end | D start | D end | J start (J alignment start) | All four | D call | All four, D placed right |
|---|---|---|---|---|---|---|---|
| naive, seed 1 (2,000) | 99.90% | 96.00% | 95.75% | 99.75% (98.55%) | 95.10% | 96.15% | 99.52% |
| mutation 0–30%, seed 1 (1,596) | 99.37% | 86.34% | 86.59% | 98.37% (96.93%) | 83.83% | 86.09% | 97.39% |
| memory + 20% indels, seed 11 (2,810) | 99.64% | 92.70% | 92.60% | 99.04% (97.79%) | 91.03% | 92.70% | 98.63% |
| bnab, seed 1 (1,208) | 95.53% | 71.94% | 71.69% | 93.87% (92.22%) | 67.80% | 72.35% | 94.56% |
| naive, d_inversion 0.5, seed 3 (2,000) | 99.90% | 94.50% | 94.65% | 99.70% (98.15%) | 93.30% | 94.40% | 99.51% |
| mutation 0–30%, d_inversion 0.5, seed 3 (1,587) | 98.87% | 82.36% | 82.73% | 98.36% (96.72%) | 80.09% | 82.42% | 98.48% |

- **The D call limits the D boundaries.** Truth's D set counts only the D-derived bases, which the read does not show. A short true D (5–10 nt) is often outscored by another allele placed partly over N bases. For example, 12 matches and 2 mismatches (+10) beat an exact 8-nt D (+8). Truth cannot accept such a D. This is a read-level limit under the origin-based truth (ruling A26).
- When the D is placed right, the boundaries agree with truth on 99.52% of naive reads. The rest are mostly V or J misplacements from the earlier steps, such as a spurious V insertion near Cys104 or a wrong short-J call.
- **J start.** The joint step moves the J start into `j_sequence_start_valid` for most reads whose J alignment start was outside it (0–30%: 96.93% → 98.37%). Most of the remainder follow a wrong D.
- **C:** on naive reads with 40 nt of C (seed 2, 2,000 reads) the primary is in `c_call_valid` on 100%.
- **Cost:** about 160 µs per naive read and 230 µs per bnab read, single thread. Without the D search it is about 2.5 µs. The D search scores 62 orientations × about 100 offsets × a window of up to the D's length, as one serial DP pass each. SIMD across offsets is the upgrade path.
- **Constants.** The ±1 score, the 5-nt D threshold, the 10-nt C threshold and the tie-break order are part of this decision. Changing any of them needs a superseding record.
- **2026-11-21 (fix wave, ruling B3): D margin of 6.** A D is called only when its best configuration beats the best no-D configuration by at least **6**; it was 1 (`Sd > S0`). The 5-nt D threshold stays.
  - **Effect.** A chance D placed over N bases rarely gains 6. An exact 5-nt true D gains only 5, but truth's 5-nt runs are where `d_call_valid` lists every allele containing those bases.
  - **Sweep.** D call accuracy on the Task 12 sets, 5,000 reads each; a null D is correct when `d_call_null_valid` (A11):

    | Setting | igh_naive | igh_memory | igh_bnab | igh_indels | igh_truncated | trb |
    |---|---|---|---|---|---|---|
    | margin 1 (before) | 95.84 | 93.20 | 68.32 | 92.54 | 74.69 | 85.72 |
    | margin 2–4 | 95.84 | 93.20 | 68.32 | 92.54 | 74.69 | 85.72–87.48 |
    | margin 5 | 95.90 | 93.22 | 68.38 | 92.58 | 74.71 | 90.28 |
    | **margin 6 (chosen)** | **96.22** | **93.33** | **68.57** | **92.76** | **74.89** | **91.94** |
    | margin 7 | 96.58 | 93.22 | 68.35 | 92.91 | 74.93 | 89.78 |
    | margin 8 | 96.28 | 92.13 | 65.45 | 92.02 | 74.54 | 83.56 |
    | min D run 7 (TRB only) | 95.84 | 93.20 | 68.32 | 92.54 | 74.69 | 87.58 |
    | min D run 7 (all loci) | 96.06 | 93.20 | 68.46 | 92.85 | 74.83 | 87.58 |
    | min D run 9 (TRB only) | 95.84 | 93.20 | 68.32 | 92.54 | 74.69 | 80.00 |

    Margin 6 maximizes the summed D accuracy with one locus-independent constant. Margin 7 for IGH alone is within 10 reads of it in total, trading naive for mutated sets. D boundaries move with the call: trb D start goes from 82.02% to 90.06%; IgBLAST reaches 92.8% on the TRB D call.
  - **Trade-off.**
    - Exact true D runs of 5 nt are no longer called. Two TRB reads in the accuracy gate are examples.
    - With no D, the no-D configuration can give the D's bases to V or J. On TRB, `v_sequence_end` fell from 99.56% to 99.44% and the V mutation fields from 99.80% to 99.51%. In one igh_memory gate read the J start moved 16 nt upstream onto IGHJ5\*02's 5′ end.
    - On the real BCR corpus, 829 reads lose a D call. Reads where only legacy calls a D rose from 182 to 980; in 972 of them legacy's D scores 6 or less against the read.
  - The margin of 6 is a constant of this decision.
