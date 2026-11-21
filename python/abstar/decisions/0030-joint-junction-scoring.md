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
- `Sd` is the best score over every D allele, orientation and placement.
- The D tie set holds the alleles reaching `Sd` in either orientation. The primary is their record 0008 pick.
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

- **A mismatch penalty of 2 or 3 inside D for choosing the D**, with the boundaries still scored ±1:
  - On naive reads it helps. Correct D-and-boundary reads rose from 1,856 to 1,889 of 2,000 at 3.
  - On reads with 0–30% mutation it hurts: 1,272 fell to 1,002 of 1,621.
  - A penalty tied to the read's V mutation level could get both. It departs from truth's ±1 score and is left for a later record.
- **D search only between the V and J alignments' ends:** misses D bases that the alignments' trims gave to V or J. The exact zone above costs little more.
- **Null D on the longest tied D run instead of the chosen one:** no naive read was called null with a D expected, so the simpler rule is kept.
- **A single global maximum across allele combinations** (truth's G7 issue): the annotator scores one combination, the primaries, so it does not arise.

## Consequences

Measured with the built-in `human` database (release build), over reads where V, J and D are each absent or have the true allele in their valid set. A boundary is correct when it is in its `_valid` list. A null D's boundaries, like the null call itself, are correct iff `d_call_null_valid`. The D call is correct when its primary is in `d_call_valid` in the true orientation.

| Set | V end | D start | D end | J start (J alignment start) | All four | D call | All four, D call correct |
|---|---|---|---|---|---|---|---|
| naive, seed 1 (2,000) | 99.85% | 94.60% | 94.65% | 99.80% (98.55%) | 93.50% | 93.35% | 99.41% |
| mutation 0–30%, seed 1 (1,621) | 99.38% | 83.10% | 83.41% | 98.40% (96.98%) | 79.64% | 80.75% | 97.17% |
| memory + 20% indels, seed 11 (2,817) | 99.65% | 90.27% | 90.38% | 98.86% (97.76%) | 88.18% | 88.57% | 98.56% |
| bnab, seed 1 (1,246) | 95.99% | 66.93% | 66.77% | 93.42% (92.54%) | 62.36% | 63.72% | 95.97% |

- **The D call limits the D boundaries.** Truth's D set counts only the D-derived bases, which the read does not show. A short true D (5–10 nt) is often outscored by another allele placed partly over N bases. On naive reads that is 48 of 2,000. For example, 12 matches and 2 mismatches (+10) beat an exact 8-nt D (+8). Truth cannot accept such a D.
- Searching the reverse complement (record 0031) adds 80 more naive misses, because the human-igh model never inverts D.
- When the D call is correct, the boundaries agree with truth on 99.41% of naive reads. The rest are mostly V or J misplacements from the earlier steps, such as a spurious V insertion near Cys104 or a wrong short-J call.
- **J start.** The joint step moves the J start into `j_sequence_start_valid` for most reads whose J alignment start was outside it (0–30%: 96.98% → 98.40%). Most of the remainder follow a wrong D.
- **C:** on naive reads with 40 nt of C (seed 2, 2,000 reads) the primary is in `c_call_valid` on 100%.
- **Cost:** about 165 µs per naive read and 230 µs per bnab read, single thread. Without the D search it is about 2.5 µs. The D search scores 62 orientations × about 100 offsets × a window of up to the D's length, as one serial DP pass each.
- **Constants.** The ±1 score, the 5-nt D threshold, the 10-nt C threshold and the tie-break order are part of this decision. Changing any of them needs a superseding record.
