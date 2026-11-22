# IMGT gapping of V germline sequences: design

Date: 2026-11-21. Status: draft for review. Board cards: #152 (this tool), #137 (macaque, its first use).

## Intent

abstar refuses a germline database unless its V alleles are IMGT-gapped (record 0036), because regions, numbering, mutations and masks all read IMGT positions. Users who have only ungapped V sequences, or V sequences gapped by another scheme, have no way to make them valid. This tool adds IMGT gaps to V sequences so that such a database can be built.

The first use is macaque. The shipped macaque source (`germline-sources/macaque/v.fasta`, IGH from KIMDB 1.1, IGK/IGL from IMGT GENE-DB) has `.` gaps, but not IMGT's: only 8 of 1,071 V alleles have Cys104 at IMGT nt 310–312. Macaque is therefore no longer shipped (#137).

## Goals

- Gap V sequences from primates and rodents, the near-term targets, accurately enough that a gapped database passes validation and its regions match IMGT's.
- Work for any species when the user supplies a close enough IMGT-gapped template database. Bovine and camelid are possible later targets.
- Never guess. An allele the tool cannot gap with confidence is reported and left out, never written with doubtful gaps.

## Non-goals

- Gapping D, J or C sequences. Only V is gapped. D/J/C pass through unchanged, as they do today.
- Bundling IMGT sequences for species abstar does not already ship. Templates for new species are the user's to supply, which keeps IMGT's licence out of the wheel.
- A profile-model (ANARCI-style) aligner. It can come later behind the same interface if single templates prove too weak for distant species.
- IMGT insertion codes (e.g. 32A). A fixed-width gapped string cannot hold them.

## Anti-goals

- Silently wrong gaps. A sequence with plausible-looking but wrong gaps is worse than a refusal, because every downstream region and mutation inherits the error.
- Copying another scheme's gaps. Any `.` in the input is stripped first; the input's gapping is never trusted.

## Approaches considered

- **A. Best single template, aligned as protein (chosen by the user).** Align each translated V to every template V allele, transfer framework positions from the best template, and place CDR1/CDR2 gaps by IMGT's length rule. Reuses existing ideas, simple to test, and strong for close species.
- **B. Profile model from all templates (ANARCI-style).** More robust for distant species, but several times the code for a case not yet needed. Deferred.
- **C. Call ANARCI.** Adds Python and HMMER dependencies, which the rewrite dropped. It also numbers protein only, so the result would still have to be mapped back to nucleotides.

## IMGT facts this design relies on

Verified against the 1,033 anchored V alleles in the shipped human, c57bl6 and balbc databases.

- IMGT positions are amino-acid positions. Position p occupies gapped nt `3(p−1) .. 3p`. Regions: FR1 1–26, CDR1 27–38, FR2 39–55, CDR2 56–65, FR3 66–104 (Cys104 at nt 310–312), then CDR3 from 105.
- Conserved residues: Cys23, Trp41 and Cys104 (and a hydrophobic residue at 89). Cys104 is what the annotator requires (`find_v_anchor`). Cys23 and Trp41 are checks only, since some functional alleles lack them.
- **CDR gaps depend only on CDR length.** In the shipped data each length has exactly one pattern:

  | CDR1 length | Gapped positions | CDR2 length | Gapped positions |
  |---|---|---|---|
  | 12 | none | 10 | none |
  | 11 | 33 | 9 | 61 |
  | 10 | 32–33 | 8 | 60–61 |
  | 9 | 32–34 | 7 | 60–62 |
  | 8 | 31–34 | 6 | 59–62 |
  | 7 | 31–35 | 5 | 59–63 |
  | 6 | 30–35 | 4 | 58–63 |
  | 5 | 30–36 | 3 | 58–64 |
  | 4 | 29–36 | | |

  Gaps open from the loop's centre outwards: CDR1 in the order 33, 32, 34, 31, 35, 30, 36, 29, then 37, 28, 38, 27; CDR2 in the order 61, 60, 62, 59, 63, 58, 64, then 57, 65, 56. The orders beyond the observed lengths follow IMGT's published rule (Lefranc et al. 2003). A CDR longer than 12 (CDR1) or 10 (CDR2) cannot be represented.
- **Framework gaps vary by locus and family** (e.g. IGH has gaps at 10 and 73; IGK at 73, 81 and 82; IGL at 10 and 73 in some families and 81–82 in others; TRA/TRB vary more). These come from the template.
- **5′-partial alleles** start with as many `.` as needed to put their first base at its IMGT nt position, even when that is not a whole codon (e.g. 16 or 44 leading dots). The 3′ end after Cys104 is contiguous, with no gaps.

## Design

### Where it lives

- A new module `crates/germline/src/gap.rs` holds the algorithm and a small scalar protein aligner. `abstar_germline` cannot depend on `abstar_annotate`, whose aligner is nucleotide-only anyway. The scale is small: about 1,000 query alleles against a few hundred same-locus templates, each about 100 aa. Parallelism is optional.
- It is exposed to Python through `crates/py`, the same way the other germline functions are.

### Interfaces

- **Command:** `abstar germline gap V_FASTA [--templates DB ...] [--output PATH] [--report PATH] [--json]`
  - Reads V sequences, gapped or not. Writes IMGT-gapped FASTA to `--output`, or to stdout.
  - `--templates` names template databases by name or path, and may be repeated. It defaults to a fixed list: `human`, `c57bl6` and `balbc`. These are the databases whose gaps come from IMGT or OGRDB. The default never includes a database this tool gapped (macaque, once it ships), so a rebuild can't inherit its own earlier mistakes. `human+c57bl6` is left out because it only duplicates the other two. When `--templates` is given, it replaces the default, so `--templates c57bl6` uses mouse alone (needed for the round-trip tests and for bovine/camelid with a user database). *(Builder check, my call.)*
  - The databases used as templates, with their manifest checksums, are recorded in the output database's provenance (recipe builds) and in the report, so a gapping can be reproduced. *(Builder check, my call.)*
  - The report (`--report` TSV, otherwise a summary on stderr; `--json` for machine-readable output) has one row per input allele: allele name, status (`gapped` or `failed`), template used, its species, the protein identity to it, CDR1/CDR2 lengths, warnings, and the failure reason.
  - Exit status follows the existing error conventions (`abstar manual errors`). Having some failed alleles is not an error. Having no alleles gapped is.
- **Recipe option:** a per-source key `gap_v = true` (and optionally `gap_templates = [...]`, with the same meaning as `--templates`). The source's V sequences are stripped and gapped during `abstar germline build`. Failed alleles are left out of the database and listed in the manifest's build warnings, the existing mechanism. Database validation (record 0036) then runs as usual.
- **Python:** `abstar.germline.gap(sequences, templates=None) -> (gapped, report)`, mirroring the command.

### Algorithm (per input allele)

1. **Clean.** Strip `.` and whitespace and uppercase the sequence. Infer its locus the way the germline build does today. Template candidates are the template databases' V alleles of the same locus that have a valid Cys104 anchor. If the locus is unknown or has no templates, use every anchored V allele.
2. **Translate.** Translate the sequence in all three frames. Each template's protein is its gapped sequence translated by IMGT position, so each template residue carries its IMGT position.
3. **Align.** Do a semi-global protein alignment (free end gaps on both sequences, since 5′-partial alleles exist) of each frame against each candidate template, with BLOSUM62 and affine gaps. Rank (frame, template) pairs by score.
4. **Transfer framework positions.** Take the best-ranked pair. Every query residue aligned to a template residue in FR1, FR2 or FR3 takes that residue's IMGT position. A template framework residue with no query partner becomes a gap at that position. That's correct, since IMGT gaps are absent positions. *(Builder check, my calls, below.)*
   - **Extra query residues:** framework positions differ within a locus (IGL: 10/73 versus 10/73/81/82; TRA: extra gaps at 7, 50 or 83–84). So query residues inserted between two template framework residues fill the empty IMGT positions between those residues, when there are enough of them. They go into the positions that templates of that locus most often fill; ties go to the leftmost position. The pair is rejected (next pair tried) only when there are more inserted residues than empty positions.
   - **Missing query residues:** in a repetitive stretch the aligner may put a deletion anywhere within it. Among the placements that score the same, choose the position where templates of that locus most often have a gap (e.g. 81–82 rather than 80), and the leftmost if none does.
   - **Near-tied templates:** if a passing pair scoring within a small margin of the chosen one would give a different gapped string, the chosen (top) pair stands and the report flags the allele with a `template_disagreement` warning. It does not fail: the round-trip tests measure how often the top pair is wrong, and failing would drop alleles the tool most likely gets right.
5. **Place CDR gaps by rule.** The query residues between its last FR1 residue and first FR2 residue are CDR1; those between FR2 and FR3 are CDR2. Lay each CDR out on its IMGT positions by the length rule above, whatever the template's CDR looked like. A CDR longer than its positions allow rejects the pair.
6. **Build the gapped nucleotide string.** Each IMGT position takes its codon, or `...` when absent. Add leading dots for a 5′-partial start, including a partial first codon, and append the contiguous 3′ tail after position 104.
7. **Check.**
   - **Required:** the gapped string, with dots removed, equals the cleaned input; every internal gap is a whole `...` codon. If either fails, try the next-ranked pair.
   - **Cys104 preferred** *(builder check, my call)*: a pair that puts TGT/TGC at nt 310–312 beats any pair that doesn't, even a higher-scoring one. Pairs are tried in rank order until one passes the required checks and places Cys104. If none places Cys104, the best pair that passes the required checks is written with a `no_cys104` warning, as IMGT does for alleles with a mutated Cys. The database-level 95% rule (record 0036) still applies.
   - **Warnings only:** a stop codon in the IMGT-positioned translation (`stop_codon`: pseudogenes are gapped as IMGT gaps them, and kept, as the legacy macaque database kept them; 32 of the 1,071 macaque V alleles have a stop in every frame) *(builder check, my call)*; Cys23 and Trp41 missing; protein identity to the template below a threshold (suggested 60%) as a signal of a distant template.
8. **Give up honestly.** If no pair passes within the top N (suggested 10), the allele fails, with the reason from its best pair.

### Macaque (#137)

- Change the macaque recipe to `gap_v = true` with the default templates. Restore it to `scripts/build_builtin_germlines.py` and the built-in list.
- Records 0019 and 0036 get dated amendments, and the built-in database tests get macaque counts again.
- Record 0019 requires built-in sequences to match legacy exactly (`test_sequences_match_legacy`). For macaque that comparison uses ungapped sequences, since the gaps change by design. Any allele the tool fails is listed in the test as an allowed dropout, with its reason. *(Builder check, my call.)*
- Done means at least 99% of the macaque V alleles with an open reading frame are anchored (the #137 bar). Pseudogenes with a stop in every frame (32 of 1,071) are gapped and kept but don't count toward the bar. The database must also pass validation. *(Builder check, my call.)*

## Validation

- **Rule table:** CDR gap placement reproduces every anchored shipped allele's CDR1/CDR2 gaps exactly.
- **Self round trip:** each shipped database's V alleles, stripped and regapped with that same database as templates (excluding the allele itself), match IMGT's gapping exactly in at least 99% of alleles.
- **Cross-species round trip (the real test):**
  - c57bl6 and balbc V alleles regapped with human templates only.
  - Human IG V alleles regapped with mouse templates only.
  - Report exact gapped-string agreement and Cys104 placement per locus.
  - Floors are set from the first measurement, with Cys104 placement at least 99% wherever IMGT's own gapping has it. Human TR has no mouse counterpart, so it gets only the self round trip.
- **Macaque:** at least 99% anchored. A spot check of IGK/IGL against IMGT GENE-DB's own macaque gapping, if those sequences are at hand.
- **Failure paths:** a framework insertion, an over-long CDR, a non-V sequence and an empty input each fail with their stated reason.

## Left to the builder

- Protein gap penalties and whether to add a cheap k-mer prefilter, if speed needs it.
- The identity warning threshold and the top-N retry limit (suggested values above), set from the round-trip measurements.
- How report rows map into the manifest's warnings.
- Whether the gapping runs multithreaded.
