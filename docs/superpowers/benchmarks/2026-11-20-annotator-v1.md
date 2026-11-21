# Annotator v1 benchmark: new abstar vs IgBLAST vs legacy abstar

Date run: 2026-11-21 (the file name keeps the plan's date); refreshed the same day after the fix wave and its review round (§0), after the final-review fix wave, after Wave 2, and after Wave 3 (both in the next sections). The report was accepted as v1 after Wave 2. The headline tables (summary, §1, §2) and §6 show the numbers after Wave 3; the dated sections keep the numbers of their time. The annotator is branch `rewrite` after the round-1 fixes (see §8). Machine: AMD Ryzen AI Max+ 395, 16 cores / 32 threads, Linux. All three tools ran with 32 threads.

**Summary.** On generator truth, new abstar ties or beats both other tools on most primary fields. It is far ahead on boundaries and on the junction wherever truth defines one. After Wave 2, three gaps remain, all with measured causes:

1. **3′-truncated reads.** The J call now leads IgBLAST: 93.0% against 91.1% (80.9% before Wave 2). Fragments under 10 nt are still called on only 30% of reads (IgBLAST 39%), and truncated D trails IgBLAST (86.6% vs 90.8%). The truncated junction is 98.8%, against IgBLAST's 82.0%.
2. **bnab J allele.** Under heavy mutation the J allele is often ambiguous (J 97.1%, against IgBLAST's 88.4% and legacy's 93.6%). The bnab junction now leads both tools: 98.2%, against legacy's 97.5% (94.5% before Wave 2).
3. **TRB D.** The D call is close to IgBLAST (91.9% vs 92.8%).

IGL J now leads both tools (99.8%).

**Constant region.** Since the final-review fix wave the benchmark grades C: three sets carry 60 nt of C, and the real-data comparison has a C row. New abstar's C call is 100% on IGH memory and IGK with C, and agrees with legacy on 99.6% of real reads (33.7% before the fix wave, when C was read one base off). On bnab with C it is 99.1% against legacy's 97.9%. It was 94.6% before Wave 2, because a C call needs a J call.

The 99.5% target holds for new abstar on IGK, IGL, TRA (except FWR3) and the non-D fields of IGH naive. It does not hold on D in any set with a D segment, or on the bnab, indel, truncated or TRB D-related fields.

On real data, new abstar agrees with the legacy baseline on 99.5% of V alleles (99.9% counting ties) and 99.1% of J alleles (99.6%); Wave 2 gave J calls to 3 of the 4 reads that had none, all agreeing with legacy. D agrees on 94.2% of reads, counting reads where both leave D null. The D figure is lower than before the fix wave because new abstar now abstains on chance-level D calls that legacy makes. Every disagreement class is explained below.

## Wave 2 (2026-11-21)

Three changes followed user rulings. Only the first moves benchmark numbers.

**1. J calling accepts weaker and shorter evidence** (record 0029, Wave 2 amendment and fix round 1). J is searched only after a V call, so a read without a V cannot get a J. The settings were chosen by generated-truth accuracy. The random-tail null is reported for information and no longer gates the choice; this supersedes B1's "+2 per 100k" bound. The J floor of 30 on the cut score is unchanged. Below it, two candidates are formed:

- **Weak J (new).** The locus J alleles are ranked by an alignment with mismatch −1 instead of −3, and the best 8 are cut with it. The best cut is a candidate if it holds its allele's Trp/Phe118 codon, implies a 6–180-nt junction and scores at least 16. Heavily mutated full-length J segments are called this way; the −3 alignment cut them short.
- **Short 3′ fragment (B1, amended).** The fragment is a candidate when its matches minus mismatches reach 10 (IG), 11 (TRB) or 12 (TRA). B1 needed an exact 12/13/14 nt plus 4 nt per mismatch, with at most 2 mismatches.
- **Order and guard.** The fragment goes first when its doubled score is at least the weak cut's score (both on the ±2 scale); otherwise the weak J goes first. A candidate whose reported alignment, after the record 0008 pick and re-cut, fails its guard is dropped, and the other is tried.
- **Joint step (record 0030).** When the J alignment holds its Trp/Phe118 codon and the read has TGG, TTT or TTC there, J may not start after it.

The order, the re-checked guard, the top-8 prefilter and the joint-step cap came from the review's fix round 1. Before them:
- a weak J's re-cut could lose its anchor;
- the joint step could start a mutated J past Trp118;
- tra g0440's exact 47-nt TRAJ8\*01 tail lost to a weak cut and its junction.

New abstar, subset `all`, before Wave 2 → after fix round 1 (5,000 reads per set; bold is a gain):

| set | `j_call` | `j_sequence_start` | `junction` | `d_call` | `v_sequence_end` | `c_call` |
|---|---|---|---|---|---|---|
| igh_bnab | 93.63 → **97.09** | 92.65 → **96.23** | 94.50 → **98.25** | 68.57 → **70.33** | 97.76 → 97.74 |  |
| igh_bnab_c | 93.23 → **97.12** | 92.37 → **96.48** | 94.73 → **98.69** | 69.31 → **70.81** | 98.17 → **98.29** | 94.62 → **99.12** |
| igh_truncated | 80.91 → **89.20** | 80.61 → **88.86** | 94.46 → **98.70** | 74.89 → **82.83** | 99.16 → **99.41** |  |
| igh_indels | 99.60 → **99.64** | 99.20 → **99.24** | 99.15 → **99.17** | 92.76 → **92.81** | 99.57 |  |
| igk | 100.00 | 99.98 | 99.94 → **99.98** | 100.00 | 100.00 |  |
| igk_c | 100.00 | 99.98 | 99.90 → **99.98** | 100.00 | 100.00 → 99.98 | 100.00 |
| igl | 99.82 | 99.80 | 99.96 → **100.00** | 100.00 | 99.98 |  |

igh_naive, igh_memory, igh_memory_c, tra and trb are unchanged on every field. The IGK and IGL junction gains come from the joint-step cap.

- **No `all` cell of `j_call`, `j_sequence_start`, `junction`, `junction_aa`, `cdr3`, `fwr4`, `productive`, `vj_in_frame` or `d_call` fell by more than 0.05 pp.**
  - The largest falls are on fields graded only where the J allele is truth's, so each newly right J adds a graded read:
    - igh_truncated `fwr4` 99.94 → 99.90: 340 more reads right, 2 more wrong;
    - igh_truncated `cdr3` 99.55 → 99.51;
    - igh_bnab_c `fwr4` 99.34 → 99.33.
  - Real one-read losses:
    - **igk_c, one V end.** The joint-step cap moved this read's V/J cut, and its V end left the valid list.
    - **One bnab read.** Its new J lets the joint V/D/J step (record 0030) pull the V end back.
- **J by J nt left on igh_truncated** (< 10 / 10–14 / 15–19 / ≥ 20): 0.0 / 31.7 / 70.4 / 98.2 → 7.5 / 73.7 / 96.9 / 99.7. IgBLAST: 39.3 / 78.8 / 91.6 / 98.9.
- **Reads truncated before J** (the watch case). Of the 54 igh_truncated reads whose truth holds no J base, 1 now gets a J (0 before). Junctions reported where truth has none went from 1 to 4.
- **Junction shortfalls.** No-J-hit junction misses fell from 169 to 11 on igh_bnab, from 175 to 10 on igh_bnab_c and from 254 to 40 on igh_truncated (shortfalls TSV).
- **Reads with a J call but no junction** (all 12 sets plus the corpus): 2,061 before Wave 2, 2,272 after it, 2,187 after fix round 1.
  - The rise is igh_truncated (424 → 622), whose added fragments mostly end before Trp118.
  - On bnab/bnab_c (52 → 64 → 54) the rest are strong or B1 cuts that start past Trp118, or reads with no Cys104.
  - The corpus fell from 805 to 739 (65 junctions now identical to legacy's).
- **Random-tail null** (`scripts/benchmark/short_j_null.py`; per 100k reads, 250-nt / 1–150-nt tails after the full V), for information:

  | locus | before | after |
  |---|---|---|
  | IGH | 5 / 2 | 673 / 1,214 |
  | IGK | 8 / 6 | 987 / 1,515 |
  | IGL | 12 / 2 | 871 / 1,105 |
  | TRA | 16 / 9 | 5,430 / 6,281 |
  | TRB | 6 / 3 | 2,243 / 2,579 |

  Most of these come from the weak J. Its junction bound halves the 250-nt rate (IGH 1,678 without it) at no cost on any benchmark cell. The sweep and the per-rule nulls are in record 0029.
- **Real corpus.** Three reads that had no J (c30155, c30367, c36768) are now annotated, and all three agree with legacy on J, junction and C.
  - J same allele 49,542 → 49,545; C same allele 49,800 → 49,803.
  - Junction identical 49,141 → 49,209. The joint-step cap gave 65 reads a junction identical to legacy's and one, ordinal 32471, a length-different one (49 → 50).
  - Otherwise no J call changed, and one D call moved: it agreed with legacy after tag stripping and now differs, with the read favoring legacy by its gapless score.
- **Speed.** Builds were timed alternately (3 runs each, 32 threads, same machine), median wall time:
  - Pre-Wave 2 → Wave 2, machine load average 6–9:
    - 12 generated sets (60,000 reads): 3.77 → 3.80 s;
    - corpus: 2.90 → 2.92 s.
  - Pre-Wave 2 → fix round 1, load average 34–45 from other jobs:
    - generated: 7.12 → 7.14 s;
    - corpus: 5.83 → 5.77 s;
    - 30,000 V + random 250-nt tails, where every read takes the fallbacks: 3.52 → 4.39 s (+25%; 4.46 s for Wave 2 before the fix).

  The speed TSV and §6 keep the Wave 2 run's times, because the fix-round benchmark ran under that load (its times roughly doubled for every tool run then).
- **Accuracy gate.** The floors were regenerated twice.
  - **Wave 2:** 507 rose and 178 fell, all on igh_bnab and igh_bnab_c.
    - 136 of the 178 hold more correct reads than before. Their n grew because more reads have the right J (and so the right D) allele and are now graded: J and D mutation fields, `fwr4`, `cdr3` and masks.
    - The other 42 are igh_bnab_c V-end and V mutation cells. They all trace to one read, g0785: its new J (junction now right) lets the joint step pull its V end from 292 back to 276.
  - **Fix round 1:** 65 rose and 15 fell. The 15 trace to one read, igh_bnab_c g1091. Its weak J was the wrong allele but bounded the right junction; its re-cut now fails the guard, so it has no J, C or junction.

**2. The D rule matches the amended spec** (commit `2ed76b5`). The spec's D call rule now states the B3 margin of 6 that the annotator already used (§0), not a 5-nt minimum. No code changed, so no D number moved.

**3. Macaque removed; IMGT-gapped databases required** (commits `b1afeaa`, `c960d29`; record 0036). The built-in macaque database is gone: `-g macaque` now fails with "database not found". Building or loading any database whose V alleles lack IMGT gaps (or where fewer than 95% of a species' V alleles carry TGT/TGC at IMGT 309–312) is refused. No benchmark set used macaque, so no number moved.

## Wave 3 (2026-11-21, after acceptance of v1)

**User rulings.**
- False J calls on junk after a V are accepted: false positives are preferred over discarding real recombinations.
- Settings lean permissive.
- The 0.1-pp bound applies only to real losses. A cell graded only on right-allele reads that falls because more reads are graded, with no fewer of them right, is dilution and does not count.

**Change** (record 0029, Wave 3 amendment). Two thresholds are lowered:
- weak J score 16 → **12**;
- B1 fragment score (matches minus mismatches) 10 / 11 / 12 → **8 / 9 / 10** (IG / TRB / TRA).

Order, guards, the re-check after the record 0008 pick and the joint-step cap are unchanged.

**The third acceptance path was not added.** It is a strong-scoring cut below 30, gapless and ending at the read's end inside J, tried after the weak J (+28.20 in record 0029's Wave 2 sweep).
- Re-swept on the fix-round-1 code, it helps only at floors of 14 or less. At 14 it calls IGHJ2\*01 on a read that ends in N/D bases whose last 7 nt match IGHJ2 exactly, which breaks the existing N/D test.
- At floors of 16 and up it adds no call once B1 is at 8: none on the 12 sets and none on 6,000 synthetic fragments.
- Lowering B1 to 8 gets more than the path did: +31.75 summed J + junction, against +30.87 for the path at 16 with B1 at 10.
- The sweep, with marginal precision by floor, is in record 0029.

New abstar, subset `all`, fix round 1 → Wave 3 (5,000 reads per set; bold is a gain). Correct reads are shown where a cell fell:

| set | `j_call` | `j_sequence_start` | `junction` | `d_call` | `cdr3` | `fwr4` |
|---|---|---|---|---|---|---|
| igh_truncated | 89.20 → **92.98** | 88.86 → **92.51** | 98.70 → **98.84** | 82.83 → **86.58** | 99.51 → 99.38 (3,851 → 4,003 right of 3,870 → 4,028; dilution) | 99.90 → 99.78 (3,866 → 4,019 of 3,870 → 4,028; dilution) |
| igh_bnab | 97.09 → **97.11** | 96.23 → **96.34** | 98.25 → **98.30** | 70.33 → **70.36** | 99.25 → 99.25 | 99.25 → 99.25 |
| igh_bnab_c | 97.12 → 97.10 (4,662 → 4,661) | 96.48 → 96.46 (4,631 → 4,630) | 98.69 → 98.65 (4,157 → 4,155; −0.047) | 70.81 → **70.90** | 99.59 → 99.56 (3,434 → 3,433) | 99.33 → 99.30 (3,425 → 3,424) |

- **The other nine sets are unchanged** on every listed field.
- **Real losses.** The only ones are the igh_bnab_c reads above: 1–2 reads per cell, at most −0.047 pp, plus C call 4,956 → 4,954. In each, a B1 fragment at the new threshold now goes before the weak J.
- **J by J nt left on igh_truncated** (< 10 / 10–14 / 15–19 / ≥ 20): 7.5 / 73.7 / 96.9 / 99.7 → 30.4 / 90.9 / 99.2 / 99.7. IgBLAST: 39.3 / 78.8 / 91.6 / 98.9.
- **Reads truncated before J** (the watch case). 4 of the 54 igh_truncated reads whose truth holds no J base now get a J (1 after Wave 2). Junctions reported where truth has none went from 4 to 6. A read ending in N/D bases (the NP test read) still gets no J at any truncation point.
- **Reads with a J call but no junction** (12 sets plus corpus): 2,187 → 2,360. The rise is all in igh_truncated (622 → 791, short fragments ending before Trp118), plus igh_bnab 31 → 34 and igh_bnab_c 23 → 24.
- **Random-tail null** (per 100k, 250-nt / 1–150-nt tails after the full V), accepted by the ruling:

  | locus | fix round 1 | Wave 3 |
  |---|---|---|
  | IGH | 673 / 1,214 | 2,153 / 5,941 |
  | IGK | 987 / 1,515 | 3,527 / 7,616 |
  | IGL | 871 / 1,105 | 3,119 / 6,042 |
  | TRA | 5,430 / 6,281 | 14,555 / 19,346 |
  | TRB | 2,243 / 2,579 | 7,349 / 10,539 |

- **Real corpus.** No call or junction changed (the real-data TSV is identical but for row order).
- **Speed** (benchmark run, load average 2–10): 9 generated sets 3.56 s (12,600 reads/s); corpus 2.84 s (17,600 reads/s).
- **Accuracy gate.** Floors regenerated: 32 rose and 24 fell. All 24 are igh_bnab J mutation cells that hold one more correct read (dilution), e.g. `j_mutations/all` 1,721 of 1,778 → 1,722 of 1,780.

## Final-review fix wave (2026-11-21)

The final branch review found that C was assigned one base off on real reads and was called on any tail. The fixes are in records 0030 (C note), 0025/0026 (generator) and 0032 (productivity).

- **C at the J/C boundary.** Built-in C alleles start with the J's last base (the codon the exons share), and real reads continue the J with C base 2. `assign_c` now tries the C allele's first base within 3 nt of J's last base, and C starts after the base J supplies (`c_germline_start` 2). If J's alignment reached into C, J is cut back. The generator now appends C from base 2, as real transcripts read; before, it duplicated the shared base, so truth and annotator agreed on a layout no real read has. Both sides take the C codon frame from that shared base.
- **Evidence gate.** No C call below 0.8 identity over the C bases (or fewer than 10 C bases). On random 10–100-nt tails after J (12,000 reads, six sets) 5 reach 0.8, all with 10–12-nt tails. On the real corpus, 28 C-bearing reads fall below it.
- **Productivity.** A non-ACGT base in the junction makes a read nonproductive. Generated reads hold only ACGT, so no generator cell moved.

C fields, before → after. Generator sets of 2,000 reads (accuracy-gate seed 20261121, truth from the fixed generator) were annotated with the annotator before the fix wave (7f27996) and after it. The C mutation fields are graded where the C call is truth's allele (n after the fix):

| set | `c_call` | `c_mutations` | `c_mutations_aa` | `c_mutation_count` | `fwr4` |
|---|---|---|---|---|---|
| igh_memory + 60 nt C | 23.50 → 100.00 | 0.00 → 100.00 (260) | 0.00 → 100.00 | 0.00 → 100.00 | 99.68 → 99.73 |
| igh_bnab + 60 nt C | 22.65 → 94.95 | 0.00 → 99.25 (265) | 0.00 → 99.25 | 0.00 → 99.25 | 97.48 → 98.63 |
| igk + 60 nt C | 100.00 → 100.00 | 0.00 → 100.00 (386) | 0.00 → 100.00 | 0.00 → 100.00 | 99.95 → 99.95 |

The IGK C allele set shares one sequence, so IGK's call was right before too, but every C mutation was wrong. Of the 101 bnab C misses after the fix, 93 have no J (86) or no V (7) call; C needs J. The other 8 have a J end more than 3 nt from truth's.

Benchmark sets (5,000 reads, §8 seeds), all three tools; IgBLAST runs without a C database:

| set | `c_call` new / legacy | C mutation fields new / legacy (n new / legacy) |
|---|---|---|
| igh_memory_c | **100.0** / 99.96 | **100.0** / 84.8 (664 / 664) |
| igh_bnab_c | 94.6 / **97.9** | **99.4** / 42.6 (656 / 679) |
| igk_c | 100.0 / 100.0 | 100.0 / 100.0 (1,019 / 1,019) |

Legacy's C mutation misses on IGH are not traced. On the nine original sets, which have no C, only `c_call` moved: the gate removed spurious C calls on a few junk tails (igh_bnab 99.94 → 100.0, igh_indels 99.76 → 100.0).

Real corpus (49,999 reads plus legacy's expected failure), new abstar vs legacy's `c_call`:

| class | before | after |
|---|---|---|
| same allele (new's primary is legacy's first allele) | 16,874 | **49,800** |
| tie overlap | 38 | 86 |
| same gene, other allele | 155 | 0 |
| different gene | 32,841 | 0 |
| only legacy calls C | 4 | 26 |
| only new calls C | 24 | 18 |
| both null | 64 | 70 |
| median `c_identity` / `c_mutation_count` | 0.29 / 106 | 1.00 / 0 |

The 86 tie overlaps are short C prefixes that IGHG1–4 or IGHA1/2 share. The 26 reads only legacy calls are gated out: an indel inside C, or a J end several bases off, defeats the gapless comparison. `c_germline_start` is 2 on 49,898 of new's calls; legacy's is 2 on 46,379 of its 49,912.

The accuracy gate gained `igh_memory_c`, `igh_bnab_c` and `igk_c`. Cells graded on fewer than 5 reads carry no floor (32 trb `ambiguous_j` cells lost theirs). No floor was lowered; four igh_bnab `c_call` floors rose.

## 0. Fix wave B1–B5 (2026-11-21)

Commit `7b4d18d`, tightened in review round 1 (B1's TR minimums, B2's distance bound, pinning tests). The amended records are 0029 (B1, B4), 0030 (B3) and 0032 (B2); their sweeps are recorded there.

| ruling | change | measured effect (new abstar, subset `all`, before → after) |
|---|---|---|
| **B1** short 3′ J fragments | Below the J floor of 30, a gapless 3′ J fragment is called when all of these hold: at most 2 mismatches; length ≥ 7 + ⌈log4(locus J nt)⌉, plus 1 for TR loci (12 for IG, 13 for TRB, 14 for TRA), plus 4 nt per mismatch; and the junction it implies is 6–180 nt, so both anchors must be known. Such reads then get a jointly placed V end, D and J start. | igh_truncated: J 77.19 → 80.91, J start 76.91 → 80.61, D 71.28 → 74.89, junction 92.23 → 94.46, V end 98.98 → 99.16. J by J nt left is in the next table. The 250-nt random-tail null (record 0029) is unchanged: per 100k, IGH 5, IGK 8, IGL 12, TRA 16, TRB 6. The 1–150-nt tail null stays within +2 per 100k of the pre-B1 rate: IGH 2 → 2, IGK 4 → 6, IGL 0 → 2, TRA 7 → 9, TRB 1 → 3. IGK, IGL, TRA and TRB all sit exactly at that +2 bound, on counts this small. The first B1 commit allowed TRA 15 and TRB 8; round 1 raised the TR minimums by 1, with no change on any benchmark set (`scripts/benchmark/short_j_null.py` reproduces both). |
| **B2** Cys104 projection | When V ends before Cys104, a TGT/TGC codon on V's last diagonal, at most 50 nt past V's last base and ending before J, is the anchor. All 107 bnab projections lie 1–49 nt out and match truth. | igh_bnab: junction 92.30 → 94.50, FWR3 96.87 → 99.03, CDR3 97.10 → 99.22, productive 92.30 → 94.50. igh_memory and igh_truncated move by about one read each; naive, light-chain and TCR sets are unchanged. Without the TGT/TGC requirement, about 126 naive and 70 TRA reads would wrongly get a junction. |
| **B3** D margin 6 | A D must beat the no-D configuration by at least 6 (was 1). | D call: trb 85.72 → 91.94, igh_naive 95.84 → 96.22, memory 93.20 → 93.33, bnab 68.32 → 68.57, indels 92.54 → 92.76, truncated 74.69 → 74.89. trb D start goes from 82.02 to 90.06, and gene_segment_mask from 88.3 to 95.2. Cost on TRB: V end 99.56 → 99.44, V mutation fields 99.80 → 99.51 (the no-D V end can take D bases); exact 5-nt true D runs are no longer called. |
| **B4** later J start among equal cuts | For J cuts with at most 3 differences, alleles are counted from the latest start among equal-score top cuts. | igl J 98.70 → 99.82, J start 98.96 → 99.80; igh_naive J 99.76 → 99.84; memory 99.56 → 99.62; bnab −1 read. Broader variants lost up to 19 pp, so this narrow form is the one kept (record 0029). |
| **B5** legacy tag | The benchmark's legacy converter now strips `__homo_sapiens` from legacy's calls. This is a formatting quirk, not an annotation error. | Legacy: truncated V 74.8 → 98.4, D +2.7 to +32 pp (trb 44.4 → 76.4). Its region and mask fields are now graded on the tagged reads too, which lowers some of them: truncated FWR3 98.1 → 74.4 at n 3,352 → 4,420; TRB gene_segment_mask 99.2 → 84.8 at n 1,872 → 3,293. §4 keeps the as-shipped numbers. |

J call on igh_truncated by J nt left:

| tool | < 10 (n=415) | 10–14 (n=463) | 15–19 (n=510) | ≥ 20 (n=3,513) |
|---|---|---|---|---|
| new abstar, before | 0.0 | 2.4 | 61.2 | 98.2 |
| new abstar, after B1 | 0.0 | 31.7 | 70.4 | 98.2 |
| new abstar, after Wave 2 | 7.2 | 73.2 | 96.3 | 99.6 |
| IgBLAST | 39.3 | 78.8 | 91.6 | 98.9 |
| legacy | 0.7 | 11.9 | 68.6 | 92.4 |

**Accuracy gate.** Cells below the 99.5% target fell from 630 to 579 of 2,276. The floors were regenerated: 41 were lowered, 293 were raised. Seven of the lowered cells had fallen below their old floors (breaches); the other 34 stayed above their old floors.

All 41 trace to B3's D margin except two bnab `nongermline_mask_aa` cells (one read each; not traced):
- **igh_memory/igh_naive `j_sequence_start` and the igh_naive J mutation fields on `ambiguous_d`.** In g1089 and g1677, B3 drops a D. The no-D configuration then extends the J start onto IGHJ5\*02's long 5′ end (316 vs valid [318, 319]; 298 vs [313, 314]), so the J span and its mutation fields fail. The two memory `ambiguous_v`/`ambiguous_j` breaches are g1677.
- **D cells on bnab and trb `ambiguous_j`/`ambiguous_d`.** These are true D runs of 5 nt (gain 5 < 6) that are now left null. g1372 and g1646 are the trb breaches.
- **TRB cap-effect cells.** `v_sequence_end` and the ten V mutation fields on `all` and `ambiguous_d`, plus `cdr_mask`. These cells measured 1.000 (floor capped at 0.995) and now measure 0.999 and 0.998. B3's no-D V end takes one or two D bases on one or two reads. The floor rule (measured − 0.5 pp, capped at 0.995) lowers the floor although no floor was breached.

The round-1 B1 tightening moved no gate cell.

## 1. At a glance (subset `all`, recoverable reads, % correct)

Each cell is **new / IgBLAST / legacy**; bold is the best when the tools differ. 5,000 reads per set, fixed seeds. Light chains and TRA have no D, so their D cells are "null D correct". Legacy is graded with its species tag stripped (B5).

| set | v_call | d_call | j_call | v_sequence_end | d_sequence_start | d_sequence_end | j_sequence_start | junction | cdr3 | productive |
|---|---|---|---|---|---|---|---|---|---|---|
| igh_naive | 100.0 / 100.0 / 100.0 | 96.2 / **96.8** / 95.6 | 99.8 / **99.9** / **99.9** | **99.9** / 99.8 / 99.7 | **96.1** / 92.9 / 94.6 | **96.3** / 93.3 / 94.4 | **99.6** / 96.6 / 98.6 | **99.8** / 97.3 / 97.3 | **99.8** / 97.3 / 97.3 | **100.0** / 99.1 / **100.0** |
| igh_memory | 99.9 / 99.9 / **100.0** | **93.3** / 90.6 / 91.0 | 99.6 / 99.7 / **99.8** | **99.9** / **99.9** / 98.7 | **93.5** / 80.2 / 91.1 | **93.2** / 78.9 / 90.2 | **99.2** / 84.0 / 96.5 | **99.4** / 97.6 / 97.6 | **99.5** / 97.7 / 97.7 | **100.0** / 98.8 / 99.8 |
| igh_bnab | 98.6 / **99.2** / 96.5 | **70.4** / 54.4 / 62.7 | **97.1** / 88.4 / 93.6 | 97.7 / **99.1** / 80.5 | **71.1** / 43.1 / 61.9 | **71.3** / 41.7 / 61.7 | **96.3** / 56.4 / 80.0 | **98.3** / 91.7 / 97.5 | 99.2 / 97.4 / **99.6** | **98.3** / 95.7 / 97.5 |
| igh_indels | **99.7** / 99.6 / 98.9 | **92.8** / 90.0 / 90.9 | 99.6 / 99.6 / **99.8** | **99.6** / 98.9 / 97.5 | **93.0** / 79.4 / 91.0 | **92.9** / 78.8 / 90.1 | **99.2** / 84.0 / 95.9 | **99.2** / 96.1 / 96.8 | **99.4** / 96.4 / 96.9 | **100.0** / 98.6 / 99.7 |
| igh_truncated | **99.9** / **99.9** / 98.4 | 86.6 / **90.8** / 69.0 | **93.0** / 91.1 / 74.8 | 99.5 / **99.7** / 73.4 | **86.8** / 79.5 / 68.6 | **86.7** / 77.9 / 68.0 | **92.5** / 77.7 / 72.1 | **98.8** / 82.0 / 90.5 | **99.4** / 95.3 / 95.7 | **100.0** / 96.0 / 75.1 |
| igk | 100.0 / 100.0 / 100.0 | 100.0 / 100.0 / 100.0 | 100.0 / 100.0 / 100.0 | **100.0** / 99.5 / 99.5 | 100.0 / 100.0 / 100.0 | 100.0 / 100.0 / 100.0 | **100.0** / 97.9 / 99.4 | 100.0 / 100.0 / 100.0 | 100.0 / 100.0 / 100.0 | 100.0 / 100.0 / 100.0 |
| igl | 99.9 / 99.8 / **100.0** | 100.0 / 100.0 / 100.0 | **99.8** / 99.7 / 99.5 | **100.0** / 99.4 / 99.6 | 100.0 / 100.0 / 100.0 | 100.0 / 100.0 / 100.0 | **99.8** / 97.6 / 99.1 | 100.0 / 100.0 / 100.0 | 100.0 / 100.0 / 100.0 | 100.0 / 100.0 / 100.0 |
| tra | **100.0** / 99.7 / 99.9 | 100.0 / 100.0 / 100.0 | 100.0 / 100.0 / 100.0 | **100.0** / 99.2 / 99.2 | 100.0 / 100.0 / 100.0 | 100.0 / 100.0 / 100.0 | **99.9** / 96.1 / 99.0 | **99.9** / 97.5 / 97.5 | **99.9** / 97.5 / 97.5 | **100.0** / 98.4 / 99.0 |
| trb | **100.0** / 99.6 / **100.0** | 91.9 / **92.8** / 76.4 | 100.0 / 100.0 / 100.0 | **99.4** / 99.0 / 98.8 | **90.1** / 88.9 / 71.7 | **90.2** / 90.0 / 72.7 | **99.7** / 97.2 / 99.6 | **99.9** / 98.6 / 98.5 | **99.9** / 98.6 / 98.6 | **100.0** / 99.2 / 99.9 |

How to read the table:

- **Grading.** `abstar.grade(truth, calls)` scores each tool only on the fields it emits; the rest are "n/a" in the TSVs. The primary call is the first allele of a comma-separated AIRR call.
  - New abstar is graded on 78 fields (86 on the C sets, which add the eight C mutation fields).
  - Legacy is graded on 54 (60 on the C sets). It has no `_length` fields, no D/J mutation fields and no mutation frequencies; its C mutation fields are graded on the C sets.
  - IgBLAST is graded on 40: calls, boundaries, junction, regions and productivity. It has no mutation strings, no masks and no `_length`.
- **Regions (`cdr3`, `fwr*`) and masks** are graded only on reads where that tool's alleles are right, so their n differs by tool. The masks need the V, D and J alleles all right, which on igh_bnab leaves 1,515 of 5,000 reads for new abstar (1,288 for legacy): its bnab mask accuracies describe that easier ~30% of the set, not the whole, and are not comparable across tools.
- **IgBLAST and legacy junctions.** Both report a junction on reads where the anchor codon was trimmed away. Truth and new abstar leave the junction null there (record 0032). This costs both tools 136 reads on igh_naive and 125 on tra (all of their junction misses there).

### Pooled by difficulty subset (all sets, % correct, new / IgBLAST / legacy)

| subset | n (v_call) | v_call | d_call | j_call | junction |
|---|---|---|---|---|---|
| high_mutation | 3183 | 98.1 / 98.9 / 95.4 | 66.1 / 48.4 / 57.6 | 96.2 / 85.5 / 91.8 | 97.9 / 89.0 / 96.8 |
| long_cdr3 | 6689 | 99.6 / 99.8 / 99.1 | 93.5 / 89.1 / 89.6 | 99.4 / 97.7 / 98.1 | 99.5 / 98.5 / 98.9 |
| has_indel | 2877 | 99.2 / 99.0 / 97.0 | 90.7 / 85.9 / 87.5 | 99.2 / 98.0 / 98.7 | 98.2 / 94.0 / 95.9 |
| truncated_5p | 4838 | 99.9 / 99.9 / 98.4 | 86.5 / 90.8 / 69.0 | 93.0 / 91.1 / 74.7 | 98.9 / 82.0 / 90.4 |
| truncated_3p | 4816 | 99.9 / 99.9 / 98.4 | 86.4 / 90.7 / 68.7 | 92.9 / 91.0 / 74.4 | 98.8 / 81.8 / 90.4 |
| ambiguous_v | 2394 | 99.3 / 99.4 / 98.1 | 82.7 / 79.4 / 72.6 | 96.7 / 92.4 / 89.5 | 98.7 / 88.6 / 93.9 |
| ambiguous_d | 8015 | 99.6 / 99.7 / 98.7 | 83.1 / 80.4 / 78.0 | 98.1 / 95.5 / 94.6 | 99.3 / 94.0 / 96.6 |
| ambiguous_j | 1786 | 99.0 / 99.5 / 98.5 | 71.2 / 76.5 / 64.2 | 89.0 / 85.9 / 85.0 | 96.8 / 88.5 / 91.4 |

## 2. Full tables (committed alongside this report)

- [`2026-11-20-annotator-v1-accuracy.tsv`](2026-11-20-annotator-v1-accuracy.tsv): every (set, tool, field, subset) cell with n, correct and accuracy (10,144 rows, including the three C sets).

The committed TSVs include the three C sets (final-review fix wave). The counts quoted in §1–§4 (tables, convention evidence, tallies) are over the nine original sets unless a line says otherwise; there, only `c_call` moved.
- [`2026-11-20-annotator-v1-shortfalls.tsv`](2026-11-20-annotator-v1-shortfalls.tsv): every new-abstar miss on a primary field, classified, with counts and example read ids.
- [`2026-11-20-annotator-v1-legacy-naming.tsv`](2026-11-20-annotator-v1-legacy-naming.tsv): legacy calls absent from abstar's DB, and legacy accuracy as shipped and with the species tag stripped.
- [`2026-11-20-annotator-v1-realdata.tsv`](2026-11-20-annotator-v1-realdata.tsv): real-data agreement classes, with read evidence.
- [`2026-11-20-annotator-v1-speed.tsv`](2026-11-20-annotator-v1-speed.tsv): wall time per set and tool.
- [`2026-11-20-annotator-v1-conventions.tsv`](2026-11-20-annotator-v1-conventions.tsv): convention-only candidates per tool and field, with the evidence counts.
- [`2026-11-20-annotator-v1-tallies.tsv`](2026-11-20-annotator-v1-tallies.tsv): at-target and head-to-head tallies, raw and without convention-only fields.

**Convention-only fields.** Some fields follow a tool's output convention that differs from truth's while the underlying nt call agrees. They are left out of the tallies below; their raw numbers stay in the accuracy TSV.

A (tool, field) is convention-only when at least 90% of the tool's misses on it, over reads with the true V and J alleles and at least 50 misses, fit the rule's evidence test. `run_benchmark.py::CONV_RULES` holds the rules, and `2026-11-20-annotator-v1-conventions.tsv` holds every candidate with its counts.

| tool | fields | evidence | misses fitting / misses |
|---|---|---|---|
| IgBLAST | `fwr4`, `fwr4_end`, `fwr4_aa` | FWR4 ends before truth's J end with the same start: one base short on nearly every read (igh_naive 4,956 of 4,986 −1) | 37,460 / 37,640 (`fwr4_end`); 99.3% (`fwr4`), 90.7% (`_aa`) |
| IgBLAST | `junction_aa`, `cdr3_aa` | nt junction or CDR3 identical, aa one trailing residue longer: IgBLAST translates the partial codon of out-of-frame junctions (`…CLTS` vs `…CLTSG`). Counted over aa misses whose nt field is right. | 9,599 / 9,600; 9,627 / 9,628 |
| IgBLAST, legacy, new | `fwr1`, `fwr1_start`, `fwr1_aa` | FWR1 starts after truth's V start with the same end: the V 5′ alignment drops mutated first bases. Truth has no valid list for the V start. 100% on igh_naive for all three. For new abstar this is not an output convention but its own V 5′ limit (the A28 extension stops short, record 0029); excluding these cells takes them out of its denominator. | IgBLAST 3,716 / 3,737; legacy 4,091 / 4,137; new 131 / 131 |
| legacy | `fwr4`, `fwr4_end`, `fwr4_aa` | FWR4 ends before truth's J end with the same start: mutated last J bases are trimmed (memory −1/−2 nt; 100% on naive) | 3,592 / 3,636 (`fwr4_end`); 96.4% (`fwr4`), 94.7% (`_aa`) |
| new (Wave 2 only) | `fwr4_end` | Same rule. After Wave 2 it had 51 misses (47 before, under the 50-miss minimum), 47 of them a mutated J 3′ end that the alignment trims. As with FWR1 this is new abstar's own J 3′ limit, not an output convention. Since Wave 3 its 57 misses fit only 82% (47 / 57), under the 90% bar, so it is graded again. | 47 / 51 (Wave 2), 47 / 57 (Wave 3) |

Candidates that did not pass stay in the tallies:
- legacy `gene_segment_mask_aa`: in 81% of its misses the nt mask is graded right, so it is mostly a codon-labeling convention, but below the bar;
- legacy `cdr_mask`/`_aa`: 43% / 35% are 3′-trimmed or have no truth junction;
- legacy `nongermline_mask`/`_aa` and `gene_segment_mask`;
- legacy `junction_aa`/`cdr3_aa`: no aa-only misses beyond 26;
- every other new-abstar candidate.

**Against the 99.5% target** (raw → without convention-only fields):

| tool | cells at target, all subsets | cells at target, subset `all` |
|---|---|---|
| new abstar | 2,102 of 3,311 → 2,012 of 3,182 | 504 of 678 → 483 of 651 |
| IgBLAST | 825 of 1,732 → 767 of 1,387 | 183 of 360 → 168 of 288 |
| legacy | 935 of 2,335 → 836 of 2,077 | 223 of 486 → 195 of 432 |

Legacy's raw count fell from 1,078 after B5, because stripping the tag grades its region fields on the tagged reads too. New abstar's raw count was 2,132 (515 on `all`) before Wave 2, 2,090 after it, 2,134 after its fix round 1 and 2,102 after Wave 3. The moves are almost all igh_truncated region, mask and J-mutation cells near the target. They grade more reads as more J alleles are called right, with more of those reads right (Wave 3: `fwr4` 3,866 of 3,870 → 4,019 of 4,028).

**Head to head.** These counts use cells with n ≥ 30 that both tools grade.

A field that is convention-only for either tool in a pair is left out of that pair's adjusted count.

- **vs IgBLAST:** new abstar is ≥ IgBLAST in 1,145 of 1,276 shared cells (89.7%) without convention-only fields. Raw it is 1,460 of 1,592 (91.7%). Earlier raw counts: 1,453 after Wave 2, 1,418 before Wave 2, 1,413 before the final-review fix wave, 1,369 before the B fix wave. It is more than 0.5 pp behind in 43 cells (44 raw; 50 after Wave 2, 79 before it, 124 before the B fix wave).
- **vs legacy:** new abstar is ≥ legacy in 1,790 of 1,908 (93.8%) without convention-only fields. Raw it is 2,019 of 2,145 (94.1%). Earlier raw counts: 1,974 before Wave 2, 1,960 before the final-review fix wave, 1,910 before the B fix wave. It is more than 0.5 pp behind in 27 (55 before Wave 2, 57 before the final-review fix wave, 115 before the B fix wave).
- **All 12 sets (with C)**, as in the tallies TSV: vs IgBLAST 1,508 of 1,724 adjusted (1,935 of 2,152 raw), 58 behind (64 after Wave 2, 106 before it); vs legacy 2,465 of 2,658 adjusted (2,776 of 2,979 raw), 37 behind (36 after Wave 2, 107 before it). The C sets repeat their base sets' gaps (igh_bnab_c behaves like igh_bnab) and add legacy's C mutation misses.

Raw and adjusted tallies are in `2026-11-20-annotator-v1-tallies.tsv`.

The `all`-subset fields where new abstar still trails by more than 0.5 pp are listed below. Each is explained in §3.

| set | field | new | IgBLAST | legacy |
|---|---|---|---|---|
| igh_naive | d_call | 96.2 | 96.8 | |
| igh_bnab | v_call / v_sequence_end | 98.6 / 97.7 | 99.2 / 99.1 | |
| igh_truncated | d_call | 86.6 | 90.8 | |
| trb | d_call | 91.9 | 92.8 | |

The fix wave closed several of the earlier gaps:

- bnab FWR3 (99.0 vs IgBLAST 99.3, legacy 97.8) and bnab CDR3 (99.2 vs legacy 99.6) are now within 0.5 pp.
- IGL J leads both tools.
- TRB D boundaries now lead IgBLAST.

Wave 2 closed three more: bnab junction, productive and vj_in_frame (now 98.2 / 98.2 / 98.3, ahead of legacy); truncated D boundaries (83.1 / 83.1, ahead of IgBLAST); and truncated V end (99.4, within 0.5 pp).

## 3. New-abstar shortfalls, classified

### Method

Each failing (read, field) is placed in one class by `run_benchmark.py::shortfalls`. Where the class depends on what the read shows, both alleles are scored against the read with one rule: the best gapless local score, +1 per match and −1 per mismatch. V and J are scored over the whole read. D is scored in both orientations, over the span between the V end and J start of truth and of the call. Ties and calls that score higher than truth mean the read cannot favor truth: **genuine ambiguity the truth set does not capture**. A call that scores lower than truth is an **annotator error**.

| class | sets (n reads) | verdict | examples |
|---|---|---|---|
| **No J hit, < 20 J nt left: below the J floor, the weak J and the short-fragment rule.** After Wave 2 these are mostly fragments under 10 nt (matches minus mismatches below 10), plus 3 reads with 20+ J nt (was 55). J accuracy by J nt left is in the Wave 2 section. These reads cause most remaining truncated-set D misses (560, was 964) and junction misses (40, was 254). | igh_truncated (510 + 3; was 881 + 55) | **Annotator limit.** Since Wave 2 the rule is chosen by truth accuracy, not the null. A looser rule (score 9, or an unguarded weak J) gains J calls here but adds wrong-offset fragments, whose `fwr4`/`cdr3` errors break the 0.1-pp cell bound (record 0029). Below about 10 nt the call is close to a guess, so part of this is genuine ambiguity. | g0003, g0011, g0020 |
| **No J hit, heavy mutation** | igh_bnab (14; was 203), igh_bnab_c (11; was 215) | **Known annotator limit**, mostly closed by the weak J (Wave 2). 11 bnab junction misses remain (was 169). | g0100, g0191, g0271 |
| **D called on a read that holds no D base.** Truth's `d_call_valid` is empty, so only a null call is correct (A11). The call is a chance placement over N. | trb 162 (was 567), igh_bnab 60, igh_indels 48, igh_naive 43, igh_memory 33, igh_truncated 32 | **Annotator limit**, reduced by B3 | trb g0054, g0099, g0118 |
| **D not called, true D present** | trb 144 + 8 (true D best score ≤ 6 / 7–9), bnab 28, indels 12, memory 9, truncated 8, naive 5 | Mostly **genuine ambiguity**: true D runs of ≤ 6 matching nt that B3's margin no longer calls. This is the trade-off for the class above. | trb g0007, g0123, g0127 |
| **D call scores above the true D in the read.** A chance D in N/P bases beats a short or mutated true D. | igh_bnab 714 (576 high_mutation), igh_indels 189, igh_memory 163, igh_truncated 142, igh_naive 104, trb 73 | **Genuine ambiguity** (A26: truth is origin-based) | igh_naive g0011, g0032, g0146 |
| **D call ties the true D** | igh_bnab 201, igh_memory 113, igh_indels 93, igh_truncated 70, igh_naive 33, trb 13 | **Genuine ambiguity.** Tie-breaking differs from the generator's draw. | igh_memory g0032, g0134, g0187 |
| **D call: true D scores higher** | 3–6 per set | **Annotator D error** (small) | trb g1785, g3418 |
| **J call ties or outscores the true allele; truth's `_valid` misses it** | bnab 92 (69 by the same count before Wave 2: the weak J calls more heavily mutated J segments, whose allele the read often cannot decide), memory 16, indels 14, igl 9, truncated 20 (was 8), naive 4 | **Truth-definition limit** (the J tie-coverage miss in the ledger). The IGLJ2/IGLJ3 ties (56 reads before) were fixed by B4. | igl g0535, g1387 |
| **V/J call: true allele scores higher** | bnab V 31, J 12; memory V 3, J 3; indels V 3, J 4; truncated V 5 | **Annotator error.** Almost all are heavy mutation, the known bnab V/J limit. | bnab g0184, g0593, g0739 |
| **V/J cut outside `_valid` with the right allele** | J start: naive 14, memory 22, indels 25, truncated 17, bnab 71, trb 14. V end: trb 26, indels 18, truncated 38, bnab 64 | **Truth-definition limit** (A27 valid-cut) in naive, memory and TCR; heavy-mutation cuts in bnab. On TRB, B3's no-D configuration adds a few V-end/J-start moves into what truth calls D. | naive g0201, g0441; trb g0173 |
| **Junction null: Cys104 beyond every valid V end, and the projected codon is not TGT/TGC** | memory 5, truncated 5, indels 4, bnab 1 | **Convention** (record 0032 as amended by B2): a mutated anchor codon is not projected, because without the Cys check about 126 naive reads would wrongly get a junction | memory g0800, g2357 |
| **Junction null: anchor codon mutated, V/J cut stops short** | memory 10, indels 4, igk 2, igl 2, bnab 1, truncated 1, tra 1, trb 1 | **Convention** (record 0032) | memory g1198, g1588 |
| **Junction follows a V-end/J-start miss** | naive 11, memory 10, indels 22, bnab 31, truncated 7 | Inherits the cut class above | naive g0201 |
| **Region fields on indel reads.** igh_indels has_indel: v_mutations 95.3%, v_insertions 95.6%, v_deletions 96.7%, fwr1 98.3%, fwr3 98.2%. | igh_indels | **Known annotator limits:** split indels (record 0033) and V 5′ extension after a 5′ indel (A28) | see the shortfalls TSV |
| **FWR3 end** (naive 99.4, memory 99.1, tra 99.3, truncated 98.8) | | **Known residual:** V trimmed just short of Cys104 (Task 11 deferred) | Task 11 follow-up report |

### IgBLAST and legacy failure modes (for context)

- **IgBLAST.**
  - D/J starts sit +3 to +5 nt inside truth's valid cuts. This alignment-trimming convention costs it J start (84.0% on memory) and D boundaries.
  - Junctions are reported where truth's anchor was trimmed away (above).
  - Of its 45 productivity misses on igh_naive, 25 overlap its junction misses.
- **Legacy.**
  - The species tag on calls (§4, stripped for grading since B5).
  - V ends often fall ≥ 5 nt short of truth's valid cuts: 80.5% on bnab, 98.7% on memory.
  - Junctions are projected (as IgBLAST).
  - Its regions on 5′-truncated reads are often wrong (truncated FWR3 74.4%).
  - It returns no row for failed records: 16 bnab and 79 truncated reads, all graded wrong.

## 4. Legacy naming effect

Legacy's bundled human DBs carry the same alleles as abstar's built-in `human` DB, with identical sequences. The IG set is the same OGRDB release; the TR set is the same IMGT download, but holds only functional alleles (338 vs abstar's F+ORF+P).

Every legacy FASTA header carries a `__homo_sapiens` suffix. Legacy usually strips it from calls but not always. It is a formatting quirk, so since B5 the benchmark strips it before grading (§1 and §2). As-shipped accuracy, for reference:

| set | segment | legacy calls not in abstar's DB (all are tagged) | accuracy as shipped | with the tag stripped (graded) |
|---|---|---|---|---|
| igh_truncated | v | 1,151 of 4,804 | 74.8 | 98.4 |
| igh_bnab | d | 844 of 3,570 | 49.0 | 62.7 |
| trb | d | 2,739 of 4,832 | 44.4 | 76.4 |
| igh_memory | d | 372 of 4,836 | 86.5 | 91.0 |
| igh_naive | d | 277 of 4,996 | 92.9 | 95.6 |

Every legacy call absent from abstar's DB was a tagged name, and every one was graded wrong before B5. Stripping the tag explains nearly all of legacy's truncated-set V deficit and roughly a third to a half of its D deficit. The rest of its D gap, and its J and junction gaps, are not naming.

## 5. Real data: legacy BCR corpus

**Corpus.** `legacy/test_data/bcr_corpus/` holds 50,000 real BCR reads. `baseline.parquet` is legacy abstar 0.8.0's committed output: 49,999 rows plus one expected failure, ordinal 11457.

**Setup.** New abstar (after Wave 2), IgBLAST and a fresh legacy run were joined on `corpus_ordinal`. The fresh legacy run reproduces the baseline's V, D and J calls on all 49,999 rows, so the baseline is current legacy behavior.

**Read evidence.** Same gapless scoring as §3, comparing new's allele with the baseline's allele. The species tag is stripped on both sides; 2,111 D calls agree only after stripping.

**Fix-wave changes on these reads** (new abstar, before → after):

- 197 J calls changed, mostly B4 on IGLJ2/IGLJ3.
- 829 reads lost a D call (B3).
- 17 junctions became defined (14 identical to legacy's; B1 and B4 J changes).
- One junction changed: ordinal 26844, IGHJ6\*03 → IGHJ6\*02, whose anchor sits one codon later.

### V (50,000 reads): 99.48% same allele, 99.94% counting tie overlap (unchanged)

| class | count | read evidence | explanation | example ordinals |
|---|---|---|---|---|
| same allele | 49,742 | | | |
| tie overlap: new's call is in legacy's list, or legacy's primary is in new's ties | 226 | 205 tie, 12 favor legacy, 9 favor new | Genuine ties broken differently, e.g. IGKV1D-13\*02 vs IGKV1-13\*01 or IGHV3-30\*03 vs \*18 | tie 28, 68, 122; legacy better 1850, 4372; new better 1218, 1233 |
| different gene | 20 | 8 new, 10 legacy, 2 tie | 28 of the 31 conflicting V calls (this row and the next) have a V indel in at least one tool's alignment. The tools explain the read with different allele + indel combinations, so the gapless score cannot decide. Examples: IGLV1-44 vs IGLV1-36 with an IMGT196 CAG deletion; IGKV2-28 vs IGKV2D-40 where legacy places a TGA deletion. | 883, 1598, 36494, 43978 |
| same gene, other allele | 11 | 2 new, 3 legacy, 6 tie | Same indel mechanism (IGLV1-44\*01 vs \*02, IGHV3-21\*01 vs \*07) | 2758, 6934 |
| only new called | 1 | | Ordinal 11457, legacy's expected failure (JunctionAnchorError). New abstar annotates it: IGHV4-39\*08 / IGHJ4\*02, junction CAKHRAPLVTGPADFDSW. | 11457 |

### J (50,000 reads): 99.09% same allele (98.93% before the fix wave, 99.08% before Wave 2), 99.60% counting tie overlap

| class | count | read evidence | explanation | example ordinals |
|---|---|---|---|---|
| same allele | 49,545 (was 49,542) | | | |
| tie overlap | 255 | 149 favor new, 96 tie, 10 favor legacy | Legacy lists both alleles (mostly IGLJ2\*01,IGLJ3\*02); new picks the better-matching one in 149 | 132, 530, 426 |
| different gene | 188 (was 268) | 114 favor new, 67 tie, 7 favor legacy | **New better:** IGHJ4\*02/IGHJ5\*02 (69) and others, where the read matches new's allele by 1+ more bases. **Ties:** 48 of 67 are IGLJ2\*01/IGLJ3\*02, identical over the aligned J and not in each other's tie list (J tie-coverage limit). B4 removed 85 of these. | new better 581, 728, 1339; ties 887, 1206; legacy better 4538, 7417 |
| same gene, other allele | 10 | 4 tie, 5 legacy, 1 new | IGHJ6\*02 vs \*03, IGKJ2\*01 vs \*04 | 28437, 3358 |
| one null | 2 (was 5) | | 1 is new no_j_hit (36217); 1 is ordinal 11457, legacy's expected failure. Wave 2 called 30155, 30367 and 36768, all as legacy does. | 36217, 11457 |

### D (50,000 reads): 94.2% agree (13,666 same allele + 33,411 both null; was 95.0%)

| class | count | read evidence | explanation | example ordinals |
|---|---|---|---|---|
| only legacy called | 980 (was 182) | 972 with legacy's D best score ≤ 6, 8 with 7–9 | Legacy calls a D from ≤ 6 matching nt; new abstains (B3 margin). These are chance-level D calls, the class that cost new abstar TRB D accuracy against truth. | 86, 168, 200 |
| tie overlap | 893 | 804 tie, 57 new, 32 legacy | Short D segments that several D genes match equally | 6, 62, 159 |
| different gene | 1,038 | **760 favor new**, 145 tie, 133 favor legacy | Where they differ, new's D usually matches the read better | new better 20, 124, 149; tie 355, 1132; legacy better 597, 665 |
| only new called | 5 (was 36) | 3 with score ≥ 10, 2 with 7–9 | Includes 11457 | 12555, 22247, 11457 |
| same gene, other allele | 7 | 7 favor new | IGHD2-21\*02 vs \*01, IGHD3-10\*03 vs \*01 | 786, 1538 |

### C (50,000 reads): 99.60% same allele (was 33.75%), 99.77% counting tie overlap

| class | count | read evidence | explanation |
|---|---|---|---|
| same allele | 49,803 (was 16,874; 49,800 before Wave 2) | | |
| tie overlap | 86 | 86 tie | Short C prefixes that IGHG1–4 or IGHA1/2 share; the read cannot decide |
| only legacy called | 23 (was 26) | | New's evidence gate (identity < 0.8): an indel inside C, or a J end more than 3 nt off, defeats the gapless comparison. Legacy aligns C with gaps. |
| only new called | 18 | | Reads legacy leaves without C |
| both null | 70 | | |

The fresh legacy run reproduces the baseline's C call on all 49,912 rows that have one. The before → after classes are in the final-review section above.

### Junction (amino acid): 98.42% identical (49,209; 49,141 before Wave 2, 49,144 before its fix round 1)

| class | count | explanation | example ordinals |
|---|---|---|---|
| only legacy has a junction | 740 (826 before the fix wave, 809 before Wave 2, 806 before its fix round 1) | 739 are new `junction_undefined`; 1 is no_j_hit. In 736 the J alignment starts after the Trp/Phe118 codon's first base. | |
| …legacy flags it nonproductive, junction does not end in W/F | 653 | The anchor codon is not in the J alignment. Legacy projects the J frame and reports a junction, then flags it nonproductive; new abstar leaves it null (record 0032 convention). | 47 (CQQYGSSP), 134 |
| …junction ends in W/F | 83 (149 before Wave 2's fix round 1, whose joint-step cap now gives 65 of them legacy's junction) | Legacy places the J start where W/F is included, e.g. ordinal 956: legacy J from germline 10, new from germline 39 on a repetitive IGHJ6 5′ end. **Annotator J-start limit** (the J start outside `_valid` class). | 380, 531, 956, 1264 |
| …V alignment stops before Cys104 | 3 | Projected codon not TGT/TGC (B2's rule); 4 of the 7 such reads before the fix wave now have a junction | 5074 |
| length differs | 50 (one more since fix round 1: ordinal 32471, where new's junction stops at an F codon before legacy's W) | Different anchor choices, mostly where one tool calls an indel near an anchor. In 9234 new's junction runs into FWR4 (a later W). In 8150 and 9369 the tools take the Cys one codon apart. In 26844, B4's IGHJ6\*02 call puts the anchor one codon later. | 9234, 8150, 26844 |
| only new has a junction | 1 | Ordinal 11457 (above) | 11457 |

**Locus** agrees on all 49,999 rows; the 50,000th is ordinal 11457.

**Cross-check with IgBLAST** (same germline DB as new abstar). New vs IgBLAST:

- V same allele 49,769, plus 216 tie overlap.
- J same allele 49,403 (49,268 before the fix wave, 49,400 before Wave 2), plus 276 tie overlap.
- D same allele 13,274 and both-null 33,775. IgBLAST alone calls a D on 616 reads, new alone on 126.

IgBLAST vs the legacy baseline:

- V same allele 49,897, plus 15 tie overlap.
- J same allele 49,611, plus 117 tie overlap.

## 6. Speed (same inputs, 32 threads, wall time including process start-up)

| input | new abstar | IgBLAST | legacy abstar |
|---|---|---|---|
| generated, 9 sets × 5,000 = 45,000 reads | 3.56 s (12,600 reads/s) | 21.9 s (2,050 reads/s) | 49.2 s (910 reads/s) |
| real corpus, 50,000 BCR reads | 2.84 s (17,600 reads/s) | 25.0 s (2,000 reads/s) | 28.7 s (1,740 reads/s) |

**Against the spec's targets** (informational: the spec names a 10-core M-series laptop, this is a 16-core AMD machine): ≥ 10k sequences/s is met here on every set but bnab (7,200–8,300 reads/s) and on the real corpus; a single sequence takes 107 ms end to end (`abstar run` on a one-read FASTA, median of 7, including Python start-up), within the ≤ 300 ms target.

Per-set times are in the speed TSV. Notes on the measurements:

- New abstar was re-timed after Wave 3 (load average 2–10); IgBLAST and legacy times on the nine original sets are from the first run (same inputs and machine). The C sets' times for all three tools are in the speed TSV. Before Wave 2 these rows read 3.55 s and 2.88 s, and after Wave 2 3.64 s and 2.94 s (under more load).
- New abstar runs 0.3–0.7 s per 5,000-read set. bnab is the slowest at 8,300 reads/s (7,200 with C).
- Legacy's per-set time includes about 1.3 s of Python import.
- New abstar ran with `--threads 32`, IgBLAST with `-num_threads 32`, and legacy with `--mmseqs_threads 32 --n_processes 32`.

## 7. Mouse

Mouse: annotated; no generator validation; no mouse real-data corpus in the repo.

## 8. Setup and versions

| component | version / source |
|---|---|
| new abstar | 1.0.0.dev0, `rewrite`: Wave 3 (J rule, record 0029) on Wave 2 and its review fix round 1 (base c960d29), release build in `.venv` (Python 3.12.14, polars 2.1.0) |
| IgBLAST | igblastn 1.22.0 (bioconda), makeblastdb 2.17.0+ |
| legacy abstar | 0.9.0, installed from a temporary copy of `legacy/` (`legacy/` is untouched); abutils 0.6.1 with its bundled MMseqs2 (commit 7e28409); Python 3.12.16 |
| MMseqs2 (env) | 19.ab57c (bioconda; installed per the brief, though legacy uses abutils' bundled binary) |
| env | micromamba 2.9.0, env `abstar-bench` (`scripts/benchmark/setup_env.sh`) |

### Truth sets

Built with `abstar generate --n 5000 --format parquet`:

| set | locus | seed | options |
|---|---|---|---|
| igh_naive | IGH | 20261201 | `--preset naive` |
| igh_memory | IGH | 20261202 | `--preset memory` |
| igh_bnab | IGH | 20261203 | `--preset bnab` (productive only) |
| igh_indels | IGH | 20261204 | `--preset memory --indels 0.5` |
| igh_truncated | IGH | 20261205 | `--preset memory --truncate-5p uniform:0-60 --truncate-3p uniform:0-40` |
| igk | IGK | 20261206 | `--preset naive` |
| igl | IGL | 20261207 | `--preset naive` |
| tra | TRA | 20261208 | `--preset tcr` |
| trb | TRB | 20261209 | `--preset tcr` |
| igh_memory_c | IGH | 20261210 | `--preset memory --constant 60` |
| igh_bnab_c | IGH | 20261211 | `--preset bnab --constant 60` |
| igk_c | IGK | 20261212 | `--preset naive --constant 60` |

The IGK/IGL/TCR presets match the accuracy gate's.

### How each tool was run

- **IgBLAST.**
  - Databases come from abstar's built-in `human` alleles via `scripts/benchmark/make_igblast_db.py`: `makeblastdb -parse_seqids -dbtype nucl`, IG and TR V/D/J kept separate.
  - FWR/CDR boundaries come from a custom internal-data file (`-custom_internal_data abstar.ndm`), derived from each V allele's IMGT gaps. FWR3 includes Cys104 (IMGT 196–312).
  - J frame and CDR3 end come from an auxiliary file derived from abstar's J anchors: frame = anchor mod 3, CDR3 end = anchor − 1.
  - Spot checks match NCBI's `human.ndm.imgt` and `human_gl.aux` exactly (IGHV1-18\*01, TRAV1-1\*01, IGKV1-39\*01, IGHJ1\*01, IGKJ1\*01, TRBJ2-7\*01).
  - Options: `-organism human -domain_system imgt -outfmt 19`, `-ig_seqtype Ig` or `TCR` per set. IgBLAST is told the receptor type; new abstar is not.
- **Legacy.** `abstar run FASTA DIR -o parquet --receptor bcr|tcr`.
  - Calls have the `__homo_sapiens` tag stripped (B5).
  - Its list fields are converted to abstar's string format:
    - `70:A>C` → `IMGT70:A>C`
    - `77-85:9>SEQ` → `IMGT77:del:SEQ`
    - `118:9>SEQ` → `IMGT118:ins:SEQ`
  - `|`-joined lists become lists. Column names are already AIRR.

Commands:

```
scripts/benchmark/setup_env.sh
MAKEBLASTDB=~/micromamba/envs/abstar-bench/bin/makeblastdb .venv/bin/python scripts/benchmark/make_igblast_db.py WORK/igblast_db
.venv/bin/python scripts/benchmark/run_benchmark.py --work WORK --igblast-db WORK/igblast_db [--out docs/superpowers/benchmarks --overwrite]
.venv/bin/python scripts/benchmark/short_j_null.py   # random-tail J null (information only) and truncated J buckets
```

`run_benchmark.py` keeps the raw outputs and timings in WORK. A rerun regrades without rerunning the tools; deleting a tool's output directory reruns that tool. Tables go to WORK/tables unless `--out` is given, and existing tables are only replaced with `--overwrite`. It also asserts that its per-read failure extraction matches `abstar.grade`'s counts for every primary field.

### Limits of this benchmark

- The read-evidence score is gapless, so it is approximate on reads with indels. This matters mainly for the real-data V conflicts, 28 of 31 of which involve an indel.
- Only IGH has truncation, indel and bnab subsets; IGK, IGL, TRA and TRB are single preset sets.
- C is graded on three sets with 60 nt of C (IGH memory, IGH bnab, IGK). Until the final-review fix wave the benchmark had no C in its truth sets, so every `c_call` row only measured that no C was called, and no C mutation field was graded. IgBLAST runs without a C database and is not graded on C. Generated C is never mutated, so C mutation accuracy on real reads (sequencing errors, indels in C) is measured only through the real-data agreement.
- The random-tail false-J rate is not part of any score. Since Wave 2, J settings are chosen by truth accuracy. Since Wave 3 the user accepts false J calls on junk after a V, and a V followed by random bases gets a J on 2–8% of IG reads and 15–19% of TRA reads (Wave 3 section). The generated sets hold no reads with a V and junk where J should be, so the benchmark does not charge for this.
- The `fwr4`, `cdr3` and mask cells are graded only on reads whose alleles are right, so a change that calls more J alleles right can lower these cells while more reads are right (Wave 2: igh_truncated `fwr4` −0.04 pp with 340 more reads right).
- The built-in macaque database was removed in Wave 2, and databases without IMGT gaps are refused (record 0036). Its re-addition from IMGT-gapped sources is #137. TRA/TRD pairing is out of scope.
