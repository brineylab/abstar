# Annotator v1 benchmark: new abstar vs IgBLAST vs legacy abstar

Date run: 2026-11-21 (the file name keeps the plan's date); refreshed the same day after the fix wave (§0). The annotator is branch `rewrite` at `7b4d18d`. Machine: AMD Ryzen AI Max+ 395, 16 cores / 32 threads, Linux. All three tools ran with 32 threads.

**Summary.** On generator truth, new abstar ties or beats both other tools on most primary fields. It is far ahead on boundaries and on the junction wherever truth defines one. After the fix wave, three gaps remain, all with measured causes:

1. **3′-truncated reads.** J fragments shorter than about 15 nt, or mutated ones, still get no call (J 80.9% vs IgBLAST 91.1%). The short-fragment rule is held back by the random-tail null.
2. **bnab junction.** New abstar finds no J under heavy mutation (junction 94.5% vs legacy 97.5%).
3. **TRB D.** The D call is now close to IgBLAST (91.9% vs 92.8%).

IGL J now leads both tools (99.8%).

The 99.5% target holds for new abstar on IGK, IGL, TRA (except FWR3) and the non-D fields of IGH naive. It does not hold on D in any set with a D segment, or on the bnab, indel, truncated or TRB D-related fields.

On real data, new abstar agrees with the legacy baseline on 99.5% of V alleles (99.9% counting ties) and 99.1% of J alleles (99.6%). D agrees on 94.2% of reads, counting reads where both leave D null. The D figure is lower than before the fix wave because new abstar now abstains on chance-level D calls that legacy makes. Every disagreement class is explained below.

## 0. Fix wave B1–B5 (2026-11-21)

Commit `7b4d18d`. The amended records are 0029 (B1, B4), 0030 (B3) and 0032 (B2); their sweeps are recorded there.

| ruling | change | measured effect (new abstar, subset `all`, before → after) |
|---|---|---|
| **B1** short 3′ J fragments | Below the J floor of 30, a gapless 3′ J fragment is called when all of these hold: at most 2 mismatches; length ≥ 7 + ⌈log4(locus J nt)⌉ (12 for IG and TRB, 13 for TRA), plus 4 nt per mismatch; and the junction it implies is 6–180 nt. Such reads then get a jointly placed V end, D and J start. | igh_truncated: J 77.19 → 80.91, J start 76.91 → 80.61, D 71.28 → 74.89, junction 92.23 → 94.46, V end 98.98 → 99.16. J by J nt left is in the next table. The 250-nt random-tail null (record 0029) is unchanged: per 100k, IGH 5, IGK 8, IGL 12, TRA 16, TRB 6. The 1–150-nt tail null rises: IGH 2 → 2, IGK 4 → 6, IGL 0 → 2, TRA 7 → 15, TRB 1 → 8. |
| **B2** Cys104 projection | When V ends before Cys104, a TGT/TGC codon on V's last diagonal, ending before J, is the anchor. | igh_bnab: junction 92.30 → 94.50, FWR3 96.87 → 99.03, CDR3 97.10 → 99.22, productive 92.30 → 94.50. igh_memory and igh_truncated move by about one read each; naive, light-chain and TCR sets are unchanged. Without the TGT/TGC requirement, about 126 naive and 70 TRA reads would wrongly get a junction. |
| **B3** D margin 6 | A D must beat the no-D configuration by at least 6 (was 1). | D call: trb 85.72 → 91.94, igh_naive 95.84 → 96.22, memory 93.20 → 93.33, bnab 68.32 → 68.57, indels 92.54 → 92.76, truncated 74.69 → 74.89. trb D start goes from 82.02 to 90.06, and gene_segment_mask from 88.3 to 95.2. Cost on TRB: V end 99.56 → 99.44, V mutation fields 99.80 → 99.51 (the no-D V end can take D bases); exact 5-nt true D runs are no longer called. |
| **B4** later J start among equal cuts | For J cuts with at most 3 differences, alleles are counted from the latest start among equal-score top cuts. | igl J 98.70 → 99.82, J start 98.96 → 99.80; igh_naive J 99.76 → 99.84; memory 99.56 → 99.62; bnab −1 read. Broader variants lost up to 19 pp, so this narrow form is the one kept (record 0029). |
| **B5** legacy tag | The benchmark's legacy converter now strips `__homo_sapiens` from legacy's calls. This is a formatting quirk, not an annotation error. | Legacy: truncated V 74.8 → 98.4, D +2.7 to +32 pp (trb 44.4 → 76.4). Its region and mask fields are now graded on the tagged reads too, which lowers some of them: truncated FWR3 98.1 → 74.4 at n 3,352 → 4,420; TRB gene_segment_mask 99.2 → 84.8 at n 1,872 → 3,293. §4 keeps the as-shipped numbers. |

J call on igh_truncated by J nt left:

| tool | < 10 (n=415) | 10–14 (n=463) | 15–19 (n=510) | ≥ 20 (n=3,513) |
|---|---|---|---|---|
| new abstar, before | 0.0 | 2.4 | 61.2 | 98.2 |
| new abstar, after B1 | 0.0 | 31.7 | 70.4 | 98.2 |
| IgBLAST | 39.3 | 78.8 | 91.6 | 98.9 |
| legacy | 0.7 | 11.9 | 68.6 | 92.4 |

**Accuracy gate.** Cells below the 99.5% target fell from 630 to 579 of 2,276. Seven cells fell below their old floors, all by one read and all from B3 (see the commit body); the floors were regenerated.

## 1. At a glance (subset `all`, recoverable reads, % correct)

Each cell is **new / IgBLAST / legacy**; bold is the best when the tools differ. 5,000 reads per set, fixed seeds. Light chains and TRA have no D, so their D cells are "null D correct". Legacy is graded with its species tag stripped (B5).

| set | v_call | d_call | j_call | v_sequence_end | d_sequence_start | d_sequence_end | j_sequence_start | junction | cdr3 | productive |
|---|---|---|---|---|---|---|---|---|---|---|
| igh_naive | 100.0 / 100.0 / 100.0 | 96.2 / **96.8** / 95.6 | 99.8 / **99.9** / **99.9** | **99.9** / 99.8 / 99.7 | **96.1** / 92.9 / 94.6 | **96.3** / 93.3 / 94.4 | **99.6** / 96.6 / 98.6 | **99.8** / 97.3 / 97.3 | **99.8** / 97.3 / 97.3 | **100.0** / 99.1 / **100.0** |
| igh_memory | 99.9 / 99.9 / **100.0** | **93.3** / 90.6 / 91.0 | 99.6 / 99.7 / **99.8** | **99.9** / **99.9** / 98.7 | **93.5** / 80.2 / 91.1 | **93.2** / 78.9 / 90.2 | **99.2** / 84.0 / 96.5 | **99.4** / 97.6 / 97.6 | **99.5** / 97.7 / 97.7 | **100.0** / 98.8 / 99.8 |
| igh_bnab | 98.6 / **99.2** / 96.5 | **68.6** / 54.4 / 62.7 | **93.6** / 88.4 / **93.6** | 97.8 / **99.1** / 80.5 | **69.1** / 43.1 / 61.9 | **69.4** / 41.7 / 61.7 | **92.7** / 56.4 / 80.0 | 94.5 / 91.7 / **97.5** | 99.2 / 97.4 / **99.6** | 94.5 / 95.7 / **97.5** |
| igh_indels | **99.7** / 99.6 / 98.9 | **92.8** / 90.0 / 90.9 | 99.6 / 99.6 / **99.8** | **99.6** / 98.9 / 97.5 | **92.9** / 79.4 / 91.0 | **92.8** / 78.8 / 90.1 | **99.2** / 84.0 / 95.9 | **99.2** / 96.1 / 96.8 | **99.4** / 96.4 / 96.9 | **100.0** / 98.6 / 99.7 |
| igh_truncated | **99.9** / **99.9** / 98.4 | 74.9 / **90.8** / 69.0 | 80.9 / **91.1** / 74.8 | 99.2 / **99.7** / 73.4 | 75.1 / **79.5** / 68.6 | 75.2 / **77.9** / 68.0 | **80.6** / 77.7 / 72.1 | **94.5** / 82.0 / 90.5 | **99.6** / 95.3 / 95.7 | **99.8** / 96.0 / 75.1 |
| igk | 100.0 / 100.0 / 100.0 | 100.0 / 100.0 / 100.0 | 100.0 / 100.0 / 100.0 | **100.0** / 99.5 / 99.5 | 100.0 / 100.0 / 100.0 | 100.0 / 100.0 / 100.0 | **100.0** / 97.9 / 99.4 | 99.9 / **100.0** / **100.0** | 99.9 / **100.0** / **100.0** | 100.0 / 100.0 / 100.0 |
| igl | 99.9 / 99.8 / **100.0** | 100.0 / 100.0 / 100.0 | **99.8** / 99.7 / 99.5 | **100.0** / 99.4 / 99.6 | 100.0 / 100.0 / 100.0 | 100.0 / 100.0 / 100.0 | **99.8** / 97.6 / 99.1 | 100.0 / 100.0 / 100.0 | 100.0 / 100.0 / 100.0 | 100.0 / 100.0 / 100.0 |
| tra | **100.0** / 99.7 / 99.9 | 100.0 / 100.0 / 100.0 | 100.0 / 100.0 / 100.0 | **100.0** / 99.2 / 99.2 | 100.0 / 100.0 / 100.0 | 100.0 / 100.0 / 100.0 | **99.9** / 96.1 / 99.0 | **99.9** / 97.5 / 97.5 | **99.9** / 97.5 / 97.5 | **100.0** / 98.4 / 99.0 |
| trb | **100.0** / 99.6 / **100.0** | 91.9 / **92.8** / 76.4 | 100.0 / 100.0 / 100.0 | **99.4** / 99.0 / 98.8 | **90.1** / 88.9 / 71.7 | **90.2** / 90.0 / 72.7 | **99.7** / 97.2 / 99.6 | **99.9** / 98.6 / 98.5 | **99.9** / 98.6 / 98.6 | **100.0** / 99.2 / 99.9 |

How to read the table:

- **Grading.** `abstar.grade(truth, calls)` scores each tool only on the fields it emits; the rest are "n/a" in the TSVs. The primary call is the first allele of a comma-separated AIRR call.
  - New abstar is graded on 78 fields.
  - Legacy is graded on 54. It has no `_length` fields, no D/J/C mutation fields and no mutation frequencies.
  - IgBLAST is graded on 40: calls, boundaries, junction, regions and productivity. It has no mutation strings, no masks and no `_length`.
- **Regions (`cdr3`, `fwr*`) and masks** are graded only on reads where that tool's alleles are right, so their n differs by tool.
- **IgBLAST and legacy junctions.** Both report a junction on reads where the anchor codon was trimmed away. Truth and new abstar leave the junction null there (record 0032). This costs both tools 136 reads on igh_naive and 125 on tra (all of their junction misses there).

### Pooled by difficulty subset (all sets, % correct, new / IgBLAST / legacy)

| subset | n (v_call) | v_call | d_call | j_call | junction |
|---|---|---|---|---|---|
| high_mutation | 3183 | 98.1 / 98.9 / 95.4 | 63.6 / 48.4 / 57.6 | 91.5 / 85.5 / 91.8 | 92.6 / 89.0 / 96.8 |
| long_cdr3 | 6689 | 99.6 / 99.8 / 99.1 | 92.7 / 89.1 / 89.6 | 98.3 / 97.7 / 98.1 | 98.4 / 98.5 / 98.9 |
| has_indel | 2877 | 99.2 / 99.0 / 97.0 | 90.3 / 85.9 / 87.5 | 98.4 / 98.0 / 98.7 | 97.5 / 94.0 / 95.9 |
| truncated_5p | 4838 | 99.9 / 99.9 / 98.4 | 74.8 / 90.8 / 69.0 | 80.9 / 91.1 / 74.7 | 94.5 / 82.0 / 90.4 |
| truncated_3p | 4816 | 99.9 / 99.9 / 98.4 | 74.6 / 90.7 / 68.7 | 80.7 / 91.0 / 74.4 | 94.4 / 81.8 / 90.4 |
| ambiguous_v | 2394 | 99.3 / 99.4 / 98.1 | 77.9 / 79.4 / 72.6 | 90.8 / 92.4 / 89.5 | 95.1 / 88.6 / 93.9 |
| ambiguous_d | 8015 | 99.6 / 99.7 / 98.7 | 81.3 / 80.4 / 78.0 | 95.4 / 95.5 / 94.6 | 97.7 / 94.0 / 96.6 |
| ambiguous_j | 1786 | 99.0 / 99.5 / 98.5 | 65.0 / 76.5 / 64.2 | 83.3 / 85.9 / 85.0 | 93.2 / 88.5 / 91.4 |

## 2. Full tables (committed alongside this report)

- [`2026-11-20-annotator-v1-accuracy.tsv`](2026-11-20-annotator-v1-accuracy.tsv): every (set, tool, field, subset) cell with n, correct and accuracy (7,378 rows).
- [`2026-11-20-annotator-v1-shortfalls.tsv`](2026-11-20-annotator-v1-shortfalls.tsv): every new-abstar miss on a primary field, classified, with counts and example read ids.
- [`2026-11-20-annotator-v1-legacy-naming.tsv`](2026-11-20-annotator-v1-legacy-naming.tsv): legacy calls absent from abstar's DB, and legacy accuracy as shipped and with the species tag stripped.
- [`2026-11-20-annotator-v1-realdata.tsv`](2026-11-20-annotator-v1-realdata.tsv): real-data agreement classes, with read evidence.
- [`2026-11-20-annotator-v1-speed.tsv`](2026-11-20-annotator-v1-speed.tsv): wall time per set and tool.

**Against the 99.5% target.** New abstar has 2,130 of 3,311 graded cells (all subsets) at or above 99.5%, up from 2,114; 515 of its 678 cells in subset `all` meet it. IgBLAST has 825 of 1,732 (183 of 360 in `all`); legacy has 935 of 2,335 (223 of 486). Legacy's count fell from 1,078 because stripping the tag grades its region fields on the tagged reads too.

**Head to head.** These counts use cells with n ≥ 30 that both tools grade.

- **vs IgBLAST (1,592 shared cells):** new abstar is ≥ IgBLAST in 1,418 (was 1,369) and more than 0.5 pp behind in 79 (was 124).
- **vs legacy (2,145 shared cells):** new abstar is ≥ legacy in 1,967 (was 1,910) and more than 0.5 pp behind in 56 (was 115).

The `all`-subset fields where new abstar still trails by more than 0.5 pp are listed below. Each is explained in §3.

| set | field | new | IgBLAST | legacy |
|---|---|---|---|---|
| igh_naive | d_call | 96.2 | 96.8 | |
| igh_bnab | v_call / v_sequence_end | 98.6 / 97.8 | 99.2 / 99.1 | |
| igh_bnab | junction / productive / vj_in_frame | 94.5 / 94.5 / 94.6 | – / 95.7 / 95.7 | 97.5 / 97.5 / 97.9 |
| igh_truncated | d_call, d boundaries, j_call, v_sequence_end | 74.9, 75.1/75.2, 80.9, 99.2 | 90.8, 79.5/77.9, 91.1, 99.7 | |
| trb | d_call | 91.9 | 92.8 | |

The fix wave closed several of the earlier gaps:

- bnab FWR3 (99.0 vs IgBLAST 99.3, legacy 97.8) and bnab CDR3 (99.2 vs legacy 99.6) are now within 0.5 pp.
- IGL J leads both tools.
- TRB D boundaries now lead IgBLAST.

## 3. New-abstar shortfalls, classified

### Method

Each failing (read, field) is placed in one class by `run_benchmark.py::shortfalls`. Where the class depends on what the read shows, both alleles are scored against the read with one rule: the best gapless local score, +1 per match and −1 per mismatch. V and J are scored over the whole read. D is scored in both orientations, over the span between the V end and J start of truth and of the call. Ties and calls that score higher than truth mean the read cannot favor truth: **genuine ambiguity the truth set does not capture**. A call that scores lower than truth is an **annotator error**.

| class | sets (n reads) | verdict | examples |
|---|---|---|---|
| **No J hit, < 20 J nt left: below both the J floor and B1's short-fragment rule.** These fragments are under 12 nt (13 for TRA), or are mutated. Another 55 reads have 20–25 J nt with 2–9 mutations. J accuracy by J nt left is in §0's table. These reads cause most remaining truncated-set D misses (910) and junction misses (232 + 22). | igh_truncated (881 + 55) | **Annotator limit, bounded by the null.** A shorter or looser rule raises the random-tail false-J rate (record 0029: 10 nt would give TRA 234 per 100k on 1–150-nt tails). Below about 10 nt the call is close to a guess, so part of this is genuine ambiguity. | g0003, g0006, g0011 |
| **No J hit, heavy mutation** | igh_bnab (203) | **Known annotator limit** (bnab heavy-mutation J). With the junction projection fixed, this is the main reason new trails legacy on bnab junction (169 reads). | g0032, g0074, g0100 |
| **D called on a read that holds no D base.** Truth's `d_call_valid` is empty, so only a null call is correct (A11). The call is a chance placement over N. | trb 162 (was 567), igh_bnab 60, igh_indels 48, igh_naive 43, igh_memory 33, igh_truncated 32 | **Annotator limit**, reduced by B3 | trb g0054, g0099, g0118 |
| **D not called, true D present** | trb 144 + 8 (true D best score ≤ 6 / 7–9), bnab 28, indels 12, memory 9, truncated 8, naive 5 | Mostly **genuine ambiguity**: true D runs of ≤ 6 matching nt that B3's margin no longer calls. This is the trade-off for the class above. | trb g0007, g0123, g0127 |
| **D call scores above the true D in the read.** A chance D in N/P bases beats a short or mutated true D. | igh_bnab 714 (576 high_mutation), igh_indels 189, igh_memory 163, igh_truncated 142, igh_naive 104, trb 73 | **Genuine ambiguity** (A26: truth is origin-based) | igh_naive g0011, g0032, g0146 |
| **D call ties the true D** | igh_bnab 201, igh_memory 113, igh_indels 93, igh_truncated 70, igh_naive 33, trb 13 | **Genuine ambiguity.** Tie-breaking differs from the generator's draw. | igh_memory g0032, g0134, g0187 |
| **D call: true D scores higher** | 3–6 per set | **Annotator D error** (small) | trb g1785, g3418 |
| **J call ties or outscores the true allele; truth's `_valid` misses it** | bnab 67, memory 16, indels 14, igl 9, truncated 8, naive 4 | **Truth-definition limit** (the J tie-coverage miss in the ledger). The IGLJ2/IGLJ3 ties (56 reads before) were fixed by B4. | igl g0535, g1387 |
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

**Setup.** New abstar (after the fix wave), IgBLAST and a fresh legacy run were joined on `corpus_ordinal`. The fresh legacy run reproduces the baseline's V, D and J calls on all 49,999 rows, so the baseline is current legacy behavior.

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

### J (50,000 reads): 99.08% same allele (was 98.93%), 99.59% counting tie overlap

| class | count | read evidence | explanation | example ordinals |
|---|---|---|---|---|
| same allele | 49,542 | | | |
| tie overlap | 255 | 149 favor new, 96 tie, 10 favor legacy | Legacy lists both alleles (mostly IGLJ2\*01,IGLJ3\*02); new picks the better-matching one in 149 | 132, 530, 426 |
| different gene | 188 (was 268) | 114 favor new, 67 tie, 7 favor legacy | **New better:** IGHJ4\*02/IGHJ5\*02 (69) and others, where the read matches new's allele by 1+ more bases. **Ties:** 48 of 67 are IGLJ2\*01/IGLJ3\*02, identical over the aligned J and not in each other's tie list (J tie-coverage limit). B4 removed 85 of these. | new better 581, 728, 1339; ties 887, 1206; legacy better 4538, 7417 |
| same gene, other allele | 10 | 4 tie, 5 legacy, 1 new | IGHJ6\*02 vs \*03, IGKJ2\*01 vs \*04 | 28437, 3358 |
| one null | 5 | | 4 are new no_j_hit (J fragments too short or mutated even for B1); 1 is ordinal 11457, legacy's expected failure | 30155, 30367, 36217 |

### D (50,000 reads): 94.2% agree (13,667 same allele + 33,411 both null; was 95.0%)

| class | count | read evidence | explanation | example ordinals |
|---|---|---|---|---|
| only legacy called | 980 (was 182) | 972 with legacy's D best score ≤ 6, 8 with 7–9 | Legacy calls a D from ≤ 6 matching nt; new abstains (B3 margin). These are chance-level D calls, the class that cost new abstar TRB D accuracy against truth. | 86, 168, 200 |
| tie overlap | 893 | 804 tie, 57 new, 32 legacy | Short D segments that several D genes match equally | 6, 62, 159 |
| different gene | 1,037 | **760 favor new**, 145 tie, 132 favor legacy | Where they differ, new's D usually matches the read better | new better 20, 124, 149; tie 355, 1132; legacy better 597, 665 |
| only new called | 5 (was 36) | 3 with score ≥ 10, 2 with 7–9 | Includes 11457 | 12555, 22247, 11457 |
| same gene, other allele | 7 | 7 favor new | IGHD2-21\*02 vs \*01, IGHD3-10\*03 vs \*01 | 786, 1538 |

### Junction (amino acid): 98.28% identical (49,141)

| class | count | explanation | example ordinals |
|---|---|---|---|
| only legacy has a junction | 809 (was 826) | 805 are new `junction_undefined`; 4 are no_j_hit. In 802 the J alignment starts after the Trp/Phe118 codon's first base. | |
| …legacy flags it nonproductive, junction does not end in W/F | 653 | The anchor codon is not in the J alignment. Legacy projects the J frame and reports a junction, then flags it nonproductive; new abstar leaves it null (record 0032 convention). | 47 (CQQYGSSP), 134 |
| …junction ends in W/F | 149 | Legacy places the J start where W/F is included, e.g. ordinal 956: legacy J from germline 10, new from germline 39 on a repetitive IGHJ6 5′ end. **Annotator J-start limit** (the J start outside `_valid` class). | 380, 531, 956, 1264 |
| …V alignment stops before Cys104 | 3 | Projected codon not TGT/TGC (B2's rule); 4 of the 7 such reads before the fix wave now have a junction | 5074 |
| length differs | 49 | Different anchor choices, mostly where one tool calls an indel near an anchor. In 9234 new's junction runs into FWR4 (a later W). In 8150 and 9369 the tools take the Cys one codon apart. In 26844, B4's IGHJ6\*02 call puts the anchor one codon later. | 9234, 8150, 26844 |
| only new has a junction | 1 | Ordinal 11457 (above) | 11457 |

**Locus** agrees on all 49,999 rows; the 50,000th is ordinal 11457.

**Cross-check with IgBLAST** (same germline DB as new abstar). New vs IgBLAST:

- V same allele 49,769, plus 216 tie overlap.
- J same allele 49,400 (was 49,268), plus 276 tie overlap.
- D same allele 13,274 and both-null 33,775. IgBLAST alone calls a D on 616 reads, new alone on 126.

IgBLAST vs the legacy baseline:

- V same allele 49,897, plus 15 tie overlap.
- J same allele 49,611, plus 117 tie overlap.

## 6. Speed (same inputs, 32 threads, wall time including process start-up)

| input | new abstar | IgBLAST | legacy abstar |
|---|---|---|---|
| generated, 9 sets × 5,000 = 45,000 reads | 3.57 s (12,600 reads/s) | 21.9 s (2,050 reads/s) | 49.2 s (910 reads/s) |
| real corpus, 50,000 BCR reads | 3.24 s (15,400 reads/s) | 25.0 s (2,000 reads/s) | 28.7 s (1,740 reads/s) |

Per-set times are in the speed TSV. Notes on the measurements:

- New abstar was re-timed after the fix wave; IgBLAST and legacy times are from the first run (same inputs and machine).
- New abstar runs 0.3–0.6 s per 5,000-read set. bnab is the slowest at 8,100 reads/s.
- Legacy's per-set time includes about 1.3 s of Python import.
- New abstar ran with `--threads 32`, IgBLAST with `-num_threads 32`, and legacy with `--mmseqs_threads 32 --n_processes 32`.

## 7. Mouse

Mouse: annotated; no generator validation; no mouse real-data corpus in the repo.

## 8. Setup and versions

| component | version / source |
|---|---|
| new abstar | 1.0.0.dev0, `rewrite` @ 7b4d18d (fix wave on 10c3612), release build in `.venv` (Python 3.12.14, polars 2.1.0) |
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
.venv/bin/python scripts/benchmark/run_benchmark.py --work WORK --igblast-db WORK/igblast_db
```

`run_benchmark.py` keeps the raw outputs and timings in WORK; a rerun regrades without rerunning the tools. It also asserts that its per-read failure extraction matches `abstar.grade`'s counts for every primary field.

### Limits of this benchmark

- The read-evidence score is gapless, so it is approximate on reads with indels. This matters mainly for the real-data V conflicts, 28 of 31 of which involve an indel.
- Only IGH has truncation, indel and bnab subsets; IGK, IGL, TRA and TRB are single preset sets.
- Macaque is deferred (#137), and TRA/TRD pairing is out of scope.
