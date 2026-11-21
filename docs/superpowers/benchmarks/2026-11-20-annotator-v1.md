# Annotator v1 benchmark: new abstar vs IgBLAST vs legacy abstar

Date run: 2026-11-21 (the file name keeps the plan's date). Branch `rewrite` at `10c3612`. Machine: AMD Ryzen AI Max+ 395, 16 cores / 32 threads, Linux. All three tools ran with 32 threads.

**Summary.** On generator truth, new abstar ties or beats both other tools on most primary fields. It is far ahead on boundaries and on the junction wherever truth defines one. There are four exceptions, all with measured causes:

1. **3′-truncated reads.** J and D calls fail when fewer than about 20 J nt remain. The cause is the J score floor of 30. IgBLAST is better here (J 91.1% vs 77.2%).
2. **TRB D.** New abstar calls a D on reads that hold no D base. IgBLAST abstains on those reads, and a null D is correct there (D 92.8% vs 85.7%).
3. **IGL J.** New abstar breaks IGLJ2\*01/IGLJ3\*02 ties the wrong way (98.7% vs IgBLAST 99.7% and legacy 99.5%).
4. **bnab junction.** No J hit under heavy mutation, plus a Cys104 that lies beyond every read-supported V end (92.3% vs legacy 97.5%).

The 99.5% target holds for new abstar on IGK, TRA (except FWR3), IGL (except the J call and J start) and the non-D fields of IGH naive. It does not hold on D in any set with a D segment, or on bnab, indel, truncated or TRB D-related fields. On real data, new abstar agrees with the legacy baseline on 99.5% of V alleles (99.9% counting ties) and 98.9% of J alleles (99.4%). D agrees on 95.0% of reads, counting reads where both leave D null. Every disagreement class is explained below.

## 1. At a glance (subset `all`, recoverable reads, % correct)

Each cell is **new / IgBLAST / legacy**; bold is the best when the tools differ. 5,000 reads per set, fixed seeds. Light chains and TRA have no D, so their D cells are "null D correct".

| set | v_call | d_call | j_call | v_sequence_end | d_sequence_start | d_sequence_end | j_sequence_start | junction | cdr3 | productive |
|---|---|---|---|---|---|---|---|---|---|---|
| igh_naive | 100.0 / 100.0 / 100.0 | 95.8 / **96.8** / 92.9 | 99.8 / **99.9** / **99.9** | **99.9** / 99.8 / 99.7 | **95.6** / 92.9 / 94.6 | **95.8** / 93.3 / 94.4 | **99.6** / 96.6 / 98.6 | **99.8** / 97.3 / 97.3 | **99.8** / 97.3 / 97.3 | **100.0** / 99.1 / **100.0** |
| igh_memory | 99.9 / 99.9 / **100.0** | **93.2** / 90.6 / 86.5 | 99.6 / 99.7 / **99.8** | **99.9** / **99.9** / 98.7 | **93.3** / 80.2 / 91.1 | **93.1** / 78.9 / 90.2 | **99.2** / 84.0 / 96.5 | **99.4** / 97.6 / 97.6 | **99.5** / 97.7 / 97.7 | **100.0** / 98.8 / 99.8 |
| igh_bnab | 98.6 / **99.2** / 94.9 | **68.3** / 54.4 / 49.0 | **93.7** / 88.4 / 93.6 | 97.7 / **99.1** / 80.5 | **68.8** / 43.1 / 61.9 | **69.1** / 41.7 / 61.7 | **92.6** / 56.4 / 80.0 | 92.3 / 91.7 / **97.5** | 97.1 / 97.4 / **99.6** | 92.3 / 95.7 / **97.5** |
| igh_indels | **99.7** / 99.6 / 98.9 | **92.5** / 90.0 / 86.6 | 99.6 / 99.6 / **99.8** | **99.6** / 98.9 / 97.5 | **92.7** / 79.4 / 91.0 | **92.6** / 78.8 / 90.1 | **99.2** / 84.0 / 95.9 | **99.2** / 96.1 / 96.8 | **99.4** / 96.4 / 96.9 | **100.0** / 98.6 / 99.7 |
| igh_truncated | **99.9** / **99.9** / 74.8 | 71.3 / **90.8** / 65.6 | 77.2 / **91.1** / 74.8 | 99.0 / **99.7** / 73.4 | 71.5 / **79.5** / 68.6 | 71.6 / **77.9** / 68.0 | 76.9 / **77.7** / 72.1 | **92.2** / 82.0 / 90.5 | **99.6** / 95.3 / 95.7 | **99.7** / 96.0 / 75.1 |
| igk | 100.0 / 100.0 / 100.0 | 100.0 / 100.0 / 100.0 | 100.0 / 100.0 / 100.0 | **100.0** / 99.5 / 99.5 | 100.0 / 100.0 / 100.0 | 100.0 / 100.0 / 100.0 | **100.0** / 97.9 / 99.4 | 99.9 / **100.0** / **100.0** | 99.9 / **100.0** / **100.0** | 100.0 / 100.0 / 100.0 |
| igl | 99.9 / 99.8 / **100.0** | 100.0 / 100.0 / 100.0 | 98.7 / **99.7** / 99.5 | **100.0** / 99.4 / 99.6 | 100.0 / 100.0 / 100.0 | 100.0 / 100.0 / 100.0 | 99.0 / 97.6 / **99.1** | 100.0 / 100.0 / 100.0 | 100.0 / 100.0 / 100.0 | 100.0 / 100.0 / 100.0 |
| tra | **100.0** / 99.7 / 99.9 | 100.0 / 100.0 / 100.0 | 100.0 / 100.0 / 100.0 | **100.0** / 99.2 / 99.2 | 100.0 / 100.0 / 100.0 | 100.0 / 100.0 / 100.0 | **99.9** / 96.1 / 99.0 | **99.9** / 97.5 / 97.5 | **99.9** / 97.5 / 97.5 | **100.0** / 98.4 / 99.0 |
| trb | **100.0** / 99.6 / **100.0** | 85.7 / **92.8** / 44.4 | 100.0 / 100.0 / 100.0 | **99.6** / 99.0 / 98.8 | 82.0 / **88.9** / 71.7 | 82.4 / **90.0** / 72.7 | **99.8** / 97.2 / 99.6 | **99.9** / 98.6 / 98.5 | **99.9** / 98.6 / 98.6 | **100.0** / 99.2 / 99.9 |

How to read the table:

- **Grading.** `abstar.grade(truth, calls)` scores each tool only on the fields it emits; the rest are "n/a" in the TSVs. The primary call is the first allele of a comma-separated AIRR call.
  - New abstar is graded on 78 fields.
  - Legacy is graded on 54. It has no `_length` fields, no D/J/C mutation fields and no mutation frequencies.
  - IgBLAST is graded on 40: calls, boundaries, junction, regions and productivity. It has no mutation strings, no masks and no `_length`.
- **Regions (`cdr3`, `fwr*`) and masks** are graded only on reads where that tool's alleles are right, so their n differs by tool.
- **Legacy is graded as shipped.** Its allele names carry a `__homo_sapiens` tag on some calls (see §4), which costs it most of its truncated-set V and TRB D score.
- **IgBLAST and legacy junctions.** Both report a junction on reads where the anchor codon was trimmed away. Truth and new abstar leave the junction null there (record 0032). This costs both tools 136 reads on igh_naive and 125 on tra (all of their junction misses there).

### Pooled by difficulty subset (all sets, % correct, new / IgBLAST / legacy)

| subset | n (v_call) | v_call | d_call | j_call | junction |
|---|---|---|---|---|---|
| high_mutation | 3183 | 98.1 / 98.9 / 93.3 | 63.3 / 48.4 / 42.1 | 91.6 / 85.5 / 91.8 | 89.7 / 89.0 / 96.8 |
| long_cdr3 | 6689 | 99.6 / 99.8 / 98.5 | 92.2 / 89.1 / 85.5 | 98.0 / 97.7 / 98.1 | 97.9 / 98.5 / 98.9 |
| has_indel | 2877 | 99.2 / 99.0 / 96.8 | 90.0 / 85.9 / 82.4 | 98.4 / 98.0 / 98.7 | 97.2 / 94.0 / 95.9 |
| truncated_5p | 4838 | 99.9 / 99.9 / 74.7 | 71.2 / 90.8 / 65.5 | 77.1 / 91.1 / 74.7 | 92.2 / 82.0 / 90.4 |
| truncated_3p | 4816 | 99.9 / 99.9 / 74.5 | 71.0 / 90.7 / 65.2 | 76.9 / 91.0 / 74.4 | 92.1 / 81.8 / 90.4 |
| ambiguous_v | 2394 | 99.3 / 99.4 / 89.8 | 76.4 / 79.4 / 63.9 | 89.9 / 92.4 / 89.5 | 93.9 / 88.6 / 93.9 |
| ambiguous_d | 8015 | 99.6 / 99.7 / 94.8 | 80.7 / 80.4 / 53.4 | 94.9 / 95.5 / 94.6 | 96.8 / 94.0 / 96.6 |
| ambiguous_j | 1786 | 99.0 / 99.5 / 84.7 | 64.3 / 76.5 / 55.7 | 82.6 / 85.9 / 85.0 | 91.7 / 88.5 / 91.4 |

## 2. Full tables (committed alongside this report)

- [`2026-11-20-annotator-v1-accuracy.tsv`](2026-11-20-annotator-v1-accuracy.tsv): every (set, tool, field, subset) cell with n, correct and accuracy (7,378 rows).
- [`2026-11-20-annotator-v1-shortfalls.tsv`](2026-11-20-annotator-v1-shortfalls.tsv): every new-abstar miss on a primary field, classified, with counts and example read ids.
- [`2026-11-20-annotator-v1-legacy-naming.tsv`](2026-11-20-annotator-v1-legacy-naming.tsv): legacy calls absent from abstar's DB, and legacy accuracy with the species tag stripped.
- [`2026-11-20-annotator-v1-realdata.tsv`](2026-11-20-annotator-v1-realdata.tsv): real-data agreement classes, with read evidence.
- [`2026-11-20-annotator-v1-speed.tsv`](2026-11-20-annotator-v1-speed.tsv): wall time per set and tool.

**Against the 99.5% target.** New abstar has 2,114 of 3,311 graded cells (all subsets) at or above 99.5%, and 515 of 678 in subset `all`. IgBLAST has 825 of 1,732 (183 of 360 in `all`); legacy has 1,078 of 2,335 (243 of 486).

**Head to head.** These counts use cells with n ≥ 30 that both tools grade.

- **vs IgBLAST (1,592 shared cells):** new abstar is ≥ IgBLAST in 1,369 and more than 0.5 pp behind in 124.
- **vs legacy (2,145 shared cells):** new abstar is ≥ legacy in 1,910 and more than 0.5 pp behind in 115.

The `all`-subset fields where new abstar trails by more than 0.5 pp are listed below. Each is explained in §3.

| set | field | new | IgBLAST | legacy |
|---|---|---|---|---|
| igh_naive | d_call | 95.8 | 96.8 | |
| igh_bnab | v_call / v_sequence_end | 98.6 / 97.7 | 99.2 / 99.1 | |
| igh_bnab | fwr3 / fwr3_end | 96.9 / 97.1 | 99.3 / 99.5 | 99.3 / 99.7 |
| igh_bnab | junction / productive / vj_in_frame | 92.3 / 92.3 / 92.4 | – / 95.7 / 95.7 | 97.5 / 97.5 / 97.9 |
| igh_bnab | cdr3 (+_aa, _start, _end, _length) | 97.1 | | 99.6 |
| igh_truncated | d_call, d boundaries, j_call, j_sequence_start, v_sequence_end | 71.3, 71.5/71.6, 77.2, 76.9, 99.0 | 90.8, 79.5/77.9, 91.1, 77.7, 99.7 | |
| igl | j_call | 98.7 | 99.7 | 99.5 |
| trb | d_call / d_sequence_start / d_sequence_end | 85.7 / 82.0 / 82.4 | 92.8 / 88.9 / 90.0 | |
| trb | gene_segment_mask / nongermline_mask | 88.3 / 98.2 | n/a | 99.2 / 99.9 (n = 1,872 vs 3,163; graded only where legacy's untagged D is right) |

## 3. New-abstar shortfalls, classified

### Method

Each failing (read, field) is placed in one class by `run_benchmark.py::shortfalls`. Where the class depends on what the read shows, both alleles are scored against the read with one rule: the best gapless local score, +1 per match and −1 per mismatch. V and J are scored over the whole read. D is scored in both orientations, over the span between the V end and J start of truth and of the call. Ties and calls that score higher than truth mean the read cannot favor truth: **genuine ambiguity the truth set does not capture**. A call that scores lower than truth is an **annotator error**.

| class | sets (n reads) | verdict | examples |
|---|---|---|---|
| **No J hit: J fragment below the J score floor (30).** 3′-truncated reads with < 20 J nt (1,063). Another 58 have 20–25 J nt carrying 2–9 mutations. J accuracy by J nt left: < 10: 0% (IgBLAST 39%); 10–14: 2% (79%); 15–19: 61% (92%); ≥ 20: 98.2% (98.6%). This causes most truncated-set D misses (1,080) and junction misses (339 + 24). It also causes 44 of the 47 V-end misses there: without J, the V 3′ end is not jointly placed (the Task 8 deferred minor). | igh_truncated (1,121) | **Annotator limit (new finding at this scale).** A design floor from Task 6. Truth counts any J base as recoverable. Below ~10 nt the call is close to a guess, so part of it is genuine ambiguity, but IgBLAST shows calls at 10–19 nt are feasible. | g0003, g0006, g0011 |
| **No J hit, heavy mutation** | igh_bnab (203) | **Known annotator limit** (bnab heavy-mutation J) | g0032, g0074, g0100 |
| **D called on a read that holds no D base.** Truth's `d_call_valid` is empty, so only a null call is correct (A11). The call is a chance placement over N. | trb 567, igh_naive 57, igh_indels 57, igh_memory 40, igh_truncated 42, igh_bnab 71 | **Annotator limit:** a permissive D threshold (A26 family). IgBLAST abstains on most of these, which is why it leads on TRB D (399 of new abstar's 714 TRB D misses are reads IgBLAST gets right, almost all with a null D). | trb g0012, g0014, g0030 |
| **D call scores above the true D in the read.** A chance D in N/P bases beats a short or mutated true D. | igh_bnab 720 (581 high_mutation), igh_indels 192, igh_memory 167, igh_truncated 139, igh_naive 109, trb 98 | **Genuine ambiguity** (A26: truth is origin-based) | igh_naive g0011, g0032, g0146 |
| **D call ties the true D** | igh_bnab 210, igh_memory 116, igh_indels 104, igh_truncated 72, igh_naive 39, trb 34 | **Genuine ambiguity.** Tie-breaking differs from the generator's draw. | igh_memory g0032, g0134, g0187 |
| **D call: true D scores higher** | 3–15 per set (trb 15, bnab 6, naive 3) | **Annotator D error** (small) | trb g1785, g2163, g3418 |
| **D not called, true D present** | trb 8, bnab 9 (best true score ≤ 9) | Mostly **genuine ambiguity** (true D ≤ 6 nt matching) | trb g0683, g1708 |
| **J call ties the true allele; truth's `_valid` misses it.** IGL: IGLJ2\*01 called IGLJ3\*02 (48) or the reverse (17). The annotator's J start lies 1–5 nt upstream, where chance bases match IGLJ3\*02's 5′ end; truth's origin-based valid set holds only the true allele. | igl 56 + 9 above, igh 6–34 per set | **Truth-definition limit** (the J tie-coverage miss in the ledger). IgBLAST and legacy start J later (IgBLAST J starts are +3/+5 nt), so they happen to agree with origin truth. A more conservative J-start tie-break would recover most of these. | igl g0231, g0407, g0450 |
| **V/J call: true allele scores higher** | bnab V 31, J 12; memory V 3, J 4; indels V 3, J 4; truncated V 5 | **Annotator error.** Almost all are heavy mutation, the known bnab V/J limit. | bnab g0184, g0593, g0739 |
| **V/J cut outside `_valid` with the right allele** | J start: naive 12, memory 21, indels 26, truncated 16, bnab 73, trb 9. V end: trb 20, indels 18, bnab 66 | **Truth-definition limit** (A27 valid-cut) in naive/memory/TCR; heavy-mutation cuts in bnab | naive g0201, g0441; trb g0375 |
| **Junction null: Cys104 beyond every valid V end.** The read's V end (all of `v_sequence_end_valid`) stops 10–30 nt short of an intact Cys codon because the V 3′ bases are heavily mutated. Truth places the anchor from the origin; record 0032 requires the anchor inside the V alignment. Legacy and IgBLAST project it. | bnab 53, memory 5, truncated 5, indels 4 | **Annotator limit / convention (new finding).** A germline-projection fallback would recover these. This class plus no-J-hit explains new abstar trailing legacy on bnab junction and CDR3. | bnab g0062 (V end 262, valid [262], Cys at 289), g0164, g0222 |
| **Junction null: anchor codon mutated, V/J cut stops short** | memory 10, indels 4, trb 3, igk 2, igl 2 | **Convention** (record 0032): the anchor codon is not inside the alignment | memory g1198, g1588 |
| **Junction follows a V-end/J-start miss** | naive 11, memory 10, indels 22, bnab 50 | Inherits the cut class above | naive g0201 |
| **Region fields on indel reads.** igh_indels has_indel: v_mutations 95.3%, v_insertions 95.6%, v_deletions 96.7%, fwr1 98.3%, fwr3 98.2%. | igh_indels | **Known annotator limits:** split indels (record 0033) and V 5′ extension after a 5′ indel (A28) | see the shortfalls TSV |
| **FWR3 end** (naive 99.1, memory 99.1, tra 99.4, truncated 98.7) | | **Known residual:** V trimmed just short of Cys104 (Task 11 deferred) | Task 11 follow-up report |

Nothing outside the ledger's known classes appears except two findings at this scale:

1. The J-floor effect on 3′-truncated reads.
2. The Cys104-beyond-valid-V-end junction loss on bnab.

The TRB/IGH "D called with no D base" class was already known under A26, but this benchmark shows IgBLAST avoids most of it.

### IgBLAST and legacy failure modes (for context)

- **IgBLAST.**
  - D/J starts sit +3 to +5 nt inside truth's valid cuts. This alignment-trimming convention costs it J start (84.0% on memory) and D boundaries.
  - Junctions are reported where truth's anchor was trimmed away (above).
  - Of its 45 productivity misses on igh_naive, 25 overlap its junction misses.
- **Legacy.**
  - The species tag on calls (§4).
  - V ends often fall ≥ 5 nt short of truth's valid cuts: 80.5% on bnab, 98.7% on memory.
  - Junctions are projected (as IgBLAST).
  - It returns no row for failed records: 16 bnab and 79 truncated reads, all graded wrong.

## 4. Legacy naming effect

Legacy's bundled human DBs carry the same alleles as abstar's built-in `human` DB, with identical sequences. The IG set is the same OGRDB release; the TR set is the same IMGT download, but holds only functional alleles (338 vs abstar's F+ORF+P).

Every legacy FASTA header carries a `__homo_sapiens` suffix. Legacy usually strips it from calls but not always, so we graded it as shipped and measured the effect:

| set | segment | legacy calls not in abstar's DB (all are tagged) | accuracy as is | with the tag stripped |
|---|---|---|---|---|
| igh_truncated | v | 1,151 of 4,804 | 74.8 | 98.4 |
| igh_bnab | d | 844 of 3,570 | 49.0 | 62.7 |
| trb | d | 2,739 of 4,832 | 44.4 | 76.4 |
| igh_memory | d | 372 of 4,836 | 86.5 | 91.0 |
| igh_naive | d | 277 of 4,996 | 92.9 | 95.6 |

Every legacy call absent from abstar's DB was wrong, and every one was a tagged name. Naming therefore explains nearly all of legacy's truncated-set V deficit (74.8 → 98.4, against new abstar's 99.9) and roughly a third to a half of its D deficit. The rest of its D gap and its J and junction gaps are not naming. Region and mask fields are also graded only where the untagged allele matches, which shrinks legacy's n on TRB masks (1,872 vs 3,163).

## 5. Real data: legacy BCR corpus

**Corpus.** `legacy/test_data/bcr_corpus/` holds 50,000 real BCR reads. `baseline.parquet` is legacy abstar 0.8.0's committed output: 49,999 rows plus one expected failure, ordinal 11457.

**Setup.** New abstar, IgBLAST and a fresh legacy run were joined on `corpus_ordinal`. The fresh legacy run reproduces the baseline's V, D and J calls on all 49,999 rows, so the baseline is current legacy behavior.

**Read evidence.** Same gapless scoring as §3, comparing new's allele with the baseline's allele. The species tag is stripped on both sides; 2,571 D calls agree only after stripping.

### V (50,000 reads): 99.48% same allele, 99.94% counting tie overlap

| class | count | read evidence | explanation | example ordinals |
|---|---|---|---|---|
| same allele | 49,742 | | | |
| tie overlap: new's call is in legacy's list, or legacy's primary is in new's ties | 226 | 205 tie, 12 favor legacy, 9 favor new | Genuine ties broken differently, e.g. IGKV1D-13\*02 vs IGKV1-13\*01 or IGHV3-30\*03 vs \*18 | tie 28, 68, 122; legacy better 1850, 4372; new better 1218, 1233 |
| different gene | 20 | 8 new, 10 legacy, 2 tie | 28 of the 31 conflicting V calls (this row and the next) have a V indel in at least one tool's alignment. The tools explain the read with different allele + indel combinations, so the gapless score cannot decide. Examples: IGLV1-44 vs IGLV1-36 with an IMGT196 CAG deletion; IGKV2-28 vs IGKV2D-40 where legacy places a TGA deletion. | 883, 1598, 36494, 43978 |
| same gene, other allele | 11 | 2 new, 3 legacy, 6 tie | Same indel mechanism (IGLV1-44\*01 vs \*02, IGHV3-21\*01 vs \*07) | 2758, 6934 |
| only new called | 1 | | Ordinal 11457, legacy's expected failure (JunctionAnchorError). New abstar annotates it: IGHV4-39\*08 / IGHJ4\*02, junction CAKHRAPLVTGPADFDSW. | 11457 |

### J (50,000 reads): 98.93% same allele, 99.43% counting tie overlap

| class | count | read evidence | explanation | example ordinals |
|---|---|---|---|---|
| same allele | 49,466 | | | |
| different gene | 268 | 152 tie, 114 favor new, 2 favor legacy | **Ties:** IGLJ2\*01 vs IGLJ3\*02 (116 of 152), identical over the aligned J; the rest are IGHJ1/IGHJ4 and IGHJ4/IGHJ5 pairs. New's tie list misses the alternative; this is the J tie-coverage limit, the same mechanism as the IGL generated-set shortfall. **New better:** mostly IGHJ4\*02/IGHJ5\*02 (69) and IGLJ2\*01/IGLJ3\*02 (25), where the read matches new's allele by 1+ more bases. | ties 623, 662, 690; new better 581, 728, 1339; legacy better 7417, 39483 |
| tie overlap | 249 | 150 favor new, 89 tie, 10 favor legacy | Legacy lists both alleles; new picks the better-matching one in 150 | 132, 530, 426 |
| same gene, other allele | 12 | 6 tie, 5 legacy, 1 new | IGHJ6\*02 vs \*03, which differ at a single base | 12353, 3358 |
| one null | 5 | | 4 new no_j_hit (short or mutated J fragments, the J floor); 1 is ordinal 11457, legacy's expected failure | 30155, 30367, 36217, 36768 |

### D (50,000 reads): 95.0% agree (14,127 same allele + 33,380 both null)

| class | count | read evidence | explanation | example ordinals |
|---|---|---|---|---|
| tie overlap | 1,154 | 1,034 tie, 78 new, 42 legacy | Short D segments that several D genes match equally | 6, 62, 159 |
| different gene | 1,114 | **825 favor new**, 154 tie, 135 favor legacy | Where they differ, new's D usually matches the read better | new better 20, 124, 149; tie 355, 1132; legacy better 597, 665 |
| only legacy called | 182 | 181 with legacy's D best score ≤ 6 | Legacy calls a D from ≤ 6 matching nt; new abstains (chance level) | 168, 229, 285 |
| only new called | 36 | 23 with score ≤ 6, 10 with 7–9, 3 with ≥ 10 | New's chance-level D calls (the same class as the TRB "no D base" shortfall); includes 11457 | 1752, 2376, 3585 |
| same gene, other allele | 7 | 7 favor new | IGHD2-21\*02 vs \*01 | 786, 1538 |

### Junction (amino acid): 98.26% identical (49,128)

| class | count | explanation | example ordinals |
|---|---|---|---|
| only legacy has a junction | 826 | 822 are new `junction_undefined`; 4 are no_j_hit. In 815 of the 822 the J alignment starts after the Trp/Phe118 codon's first base. | |
| …legacy flags it nonproductive, junction does not end in W/F | 665 | The anchor codon is not in the J alignment. Legacy projects the J frame and reports a junction, then flags it nonproductive; new abstar leaves it null (record 0032 convention). | 47 (CQQYGSSP), 134 |
| …junction ends in W/F | 157 | 110 productive per legacy. Legacy places the J start where W/F is included, e.g. ordinal 956: legacy J from germline 10, new from germline 39 on a repetitive IGHJ6 5′ end. **Annotator J-start limit** (the J start outside `_valid` class). | 380, 531, 956, 1264 |
| …V alignment stops before Cys104 | 7 | Same as the bnab class above | 5074, 11567 |
| length differs | 45 | Different anchor choices, mostly where one tool calls an indel near an anchor. In 9234 new's junction runs into FWR4 (a later W). In 8150 and 9369 the tools take the Cys one codon apart. | 9234, 8150, 9369 |
| only new has a junction | 1 | Ordinal 11457 (above) | 11457 |

**Locus** agrees on all 49,999 rows; the 50,000th is ordinal 11457.

**Cross-check with IgBLAST** (same germline DB as new abstar). New vs IgBLAST:

- V same allele 49,769, plus 216 tie overlap.
- J same allele 49,268, plus 263 tie overlap.
- D same allele 13,674 and both-null 33,553.

IgBLAST vs the legacy baseline:

- V same allele 49,897, plus 15 tie overlap.
- J same allele 49,611, plus 117 tie overlap.

## 6. Speed (same inputs, 32 threads, wall time including process start-up)

| input | new abstar | IgBLAST | legacy abstar |
|---|---|---|---|
| generated, 9 sets × 5,000 = 45,000 reads | 3.65 s (12,300 reads/s) | 21.9 s (2,050 reads/s) | 49.2 s (910 reads/s) |
| real corpus, 50,000 BCR reads | 3.30 s (15,100 reads/s) | 25.0 s (2,000 reads/s) | 28.7 s (1,740 reads/s) |

Per-set times are in the speed TSV. Notes on the measurements:

- New abstar runs 0.3–0.6 s per 5,000-read set. bnab is the slowest at 8,250 reads/s.
- Legacy's per-set time includes about 1.3 s of Python import.
- New abstar ran with `--threads 32`, IgBLAST with `-num_threads 32`, and legacy with `--mmseqs_threads 32 --n_processes 32`.

## 7. Mouse

Mouse: annotated; no generator validation; no mouse real-data corpus in the repo.

## 8. Setup and versions

| component | version / source |
|---|---|
| new abstar | 1.0.0.dev0, `rewrite` @ 10c3612, release build in `.venv` (Python 3.12.14, polars 2.1.0) |
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
- **Legacy.** `abstar run FASTA DIR -o parquet --receptor bcr|tcr`. Its list fields are converted to abstar's string format:
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
