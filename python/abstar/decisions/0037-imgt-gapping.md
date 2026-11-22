# 0037: IMGT gapping of V germline sequences

## Status

Accepted (2026-11-21)

## Context

A germline database must have IMGT-gapped V alleles (record 0036): regions, numbering, mutations and masks read IMGT positions. Users with ungapped V sequences, or V sequences gapped by another scheme, had no way to make them valid. The built-in macaque V genes (KIMDB) are the first case: only 8 of 1,071 have Cys104 at IMGT nt 310–312 (#137). This record fixes how `abstar germline gap` (#152, `crates/germline/src/gap`) adds IMGT gaps to V sequences, and the constants it uses. The design is `docs/superpowers/specs/2026-11-21-imgt-gapping-design.md`.

IMGT position p occupies gapped nt `3(p−1)..3p`. Regions: FR1 1–26, CDR1 27–38, FR2 39–55, CDR2 56–65, FR3 66–104 (Cys104 at nt 310–312), then CDR3 from 105. In the 1,033 anchored V alleles of the shipped human, c57bl6 and balbc databases, CDR gaps depend only on CDR length. Framework gaps depend on locus and family.

## Decision

**Method: best single template, aligned as protein** (`gap_v`). Per input allele:

1. **Clean.** Strip `.` and whitespace and uppercase. The input's own gaps are never trusted. The locus comes from the first three characters of the name. Candidates are the template V alleles of that locus with a Cys104 anchor. If the locus is unknown or has no templates, every anchored template V allele is a candidate.
2. **Translate** the query in its three frames. A codon with a non-ACGT base is `X`, which scores −1 against every residue. Each template residue carries its IMGT position, read from the template's gapped string.
3. **Align** every (frame, template) pair. The aligner is semi-global (both ends free on both sequences) with BLOSUM62 and affine gaps. A gap opening in a template CDR column (`TemplateOnly`), or next to one (`QueryOnly`), costs less than one in a framework column. CDR lengths vary between alleles while framework lengths rarely do, and the CDR's residues are laid out by rule in step 6 anyway. Pairs are ranked by score, then template order, then frame.
4. **Floor.** The best pair must score above 0 and align at least 50 query residues; otherwise the allele fails `no_alignment` (it is not a V). Lower pairs are held to the same 50-residue floor before they are tried, so a short wrong-frame hit is never written in place of a refusal.
5. **Transfer framework positions** from the template, for each of the first `TOP_N` pairs in rank order:
   - A query residue paired with a template residue takes its position. A template framework residue with no query partner leaves its position as a gap.
   - **Framework occupancy rules.** Occupancy is, per IMGT position, how many candidate templates have a residue there.
     - Missing query residues: among the placements of a framework deletion with the same alignment score, the one whose emptied positions have the lowest summed occupancy wins (ties leftmost). The deletion stays within one framework region.
     - Extra query residues between two template residues of the same framework fill the empty positions between them that have the highest occupancy, in order (ties leftmost). More extra residues than empty positions fails the pair with `framework_insertion`.
   - **Overhangs.** Query residues before the first pair run contiguously up to it (`unplaceable_overhang` if they would pass position 1).
   - **Tail.** Positions after 104 are contiguous whatever the template's tail does. Pairs beyond 104 are dropped, and the residues after the last placed residue at ≤ 104 run on from 105. A 3′ overhang is allowed only when the alignment reaches 104 or beyond (`unplaceable_overhang` otherwise).
6. **CDR rule.** The query residues between the last FR1 residue and the first FR2 residue are CDR1; those between FR2 and FR3 are CDR2. Each is laid out by its length, whatever the template's CDR looks like. Gaps open from the loop's centre outwards:
   - CDR1 order: 33, 32, 34, 31, 35, 30, 36, 29, 37, 28, 38, 27.
   - CDR2 order: 61, 60, 62, 59, 63, 58, 64, 57, 65, 56.
   - A CDR of length k leaves the first (12 − k) or (10 − k) positions of its order empty. This reproduces every shipped allele. The orders beyond the observed lengths follow Lefranc et al. (2003).
   - A CDR1 longer than 12 or a CDR2 longer than 10 fails the pair with `cdr_too_long`.
   - A query that starts inside a CDR (5′-partial) has no known CDR length. Its residues keep the template's positions. They are right-aligned on the CDR's last positions only if one of them has no position.
7. **Build** the gapped nucleotide string: a codon per position, `...` for an absent one. Leading dots put the first base at its IMGT nt position, including a partial first codon. Trailing bases that don't fill a codon are kept.
8. **Check.** The string with dots removed must equal the cleaned input, and every internal dot run must be whole codons. This holds by construction, so a pair that fails it has no named reason. A `debug_assert!` fires if every tried pair fails that way.
9. **Choose.** The first passing pair (in rank order) with TGT/TGC at nt 310–312 wins. If no passing pair has one, the first passing pair is written with `no_cys104`, as IMGT writes alleles with a mutated Cys104. If no pair passes, the allele fails with the first named reason in rank order.

**Default templates.** The defaults are exactly `human`, `c57bl6` and `balbc` (`DEFAULT_TEMPLATES`), whose gaps come from IMGT or OGRDB. `human+c57bl6` would only duplicate two of them. A database this tool gapped (macaque, once it ships) is never a default template, so a rebuild cannot inherit its own earlier mistakes. `--templates` / `gap_templates` / `templates=` replace the default list.

**Warning codes** (the allele is written):

| Code | Meaning |
|---|---|
| `no_cys104` | no passing pair puts TGT/TGC at nt 310–312 |
| `stop_codon` | a stop in the translated frame (pseudogenes are gapped and kept, as IMGT gaps them) |
| `no_cys23` | no Cys at position 23 |
| `no_trp41` | no Trp at position 41 |
| `low_identity` | protein identity to the template, over aligned pairs, below `LOW_IDENTITY` |
| `template_disagreement` | another passing pair scoring within `TIE_MARGIN` of the winner gives a different gapped string |

**Failure reasons** (the allele is not written): `empty_sequence` (nothing left after cleaning), `no_alignment` (best pair at or below score 0, or under 50 aligned residues), `framework_insertion`, `cdr_too_long` and `unplaceable_overhang` (the first named reason among the tried pairs).

**Constants** (`gap/mod.rs`, `gap/align.rs`):

| Constant | Value | Why |
|---|---|---|
| gap open (framework) | 13 | sweep, below |
| gap extend | 1 | sweep |
| CDR gap open (`cdr_gap_open`) | 5 | sweep |
| `TOP_N` | 10 | 1 loses 2 self and 2 cross-species exact matches, 3 loses 1 and 2, 20 gains none |
| `LOW_IDENTITY` | 0.60 | flags 5 of the 10 wrongly written alleles and 84 of 1,804 exact ones. 0.70 would catch 3 more wrong ones and flag 408 exact ones |
| `TIE_MARGIN` | 10 score units | flags 2 of the 10 wrongly written alleles and 12 of 1,804 exact ones. At 5: 1 and 4. At 0: none |
| `no_alignment` floor | 50 aligned residues | plan ruling: below 50 aligned residues is not a V. It also applies to each tried pair |

**Measurement** (`crates/germline/tests/gap_roundtrip.rs`; `cargo test -p abstar_germline --release --test gap_roundtrip -- --nocapture`, about 5 s):

- Each anchored V allele of a shipped database is stripped and regapped, then compared with IMGT's string.
  - Self round trip: leave-one-out, with the allele's own database minus the allele.
  - Cross-species: mouse from human templates only; human IG from c57bl6 + balbc only. Human TR has no mouse counterpart.
- *Exact* counts all alleles, so a failure counts as not exact. *Cys104* counts the written alleles whose Cys104 codon is the one IMGT puts at nt 310–312. A refusal places no Cys104, right or wrong.

**Sweep.** Exact matches per set; Σ is the summed cross-species exact rate in percent. The rule: keep the highest Σ whose self round trip holds ≥ 99% in every database. All 70 combinations of open 9–13, extend 1–2 and CDR open {1, 2, 3, 4, 5, 7, 9} were run; excerpt:

| open | extend | CDR open | human | c57bl6 | balbc | c57bl6 ← human | balbc ← human | human IG ← mouse | self ≥ 99% | Σ |
|---|---|---|---|---|---|---|---|---|---|---|
| 11 | 1 | none (= 11) | 572/577 | 191/193 | 261/263 | 189/193 | 257/263 | 306/340 | no | 285.6 |
| 11 | 1 | 3 | 572/577 | 192/193 | 262/263 | 190/193 | 261/263 | 321/340 | yes | 292.1 |
| 12 | 1 | 4 | 572/577 | 192/193 | 262/263 | 192/193 | 262/263 | 323/340 | yes | 294.1 |
| 13 | 1 | 4 | 572/577 | 192/193 | 262/263 | 192/193 | 262/263 | 323/340 | yes | 294.1 |
| **13** | **1** | **5** | **573/577** | **192/193** | **262/263** | **192/193** | **262/263** | **323/340** | **yes** | **294.1** |
| 12 | 2 | 4 | 573/577 | 191/193 | 261/263 | 192/193 | 262/263 | 327/340 | no | 295.3 |
| 13 | 2 | 5 | 573/577 | 191/193 | 261/263 | 192/193 | 262/263 | 326/340 | no | 295.0 |

- Three settings tie on Σ. 13/1/5 has the most self matches.
- Every extend-2 setting holds c57bl6 to 191/193 (98.96%), so none qualifies, though some reach a higher Σ.
- Without the CDR open cost (re-measured with the fixes below), c57bl6 stays at or below 191/193 for every open 9–13 / extend 1–2, so no setting holds the self round trip.
- That first row is the starting point (open 11, extend 1, no CDR cost). It already includes three fixes made during tuning: the contiguous tail, the per-pair 50-residue floor and the 5′-partial CDR rule. Before them, human was 571/577.
- `TIE_MARGIN` changes only the warning, not the gapped strings, so it was chosen from warning counts (constants table).

**Final round trips** (open 13, extend 1, CDR open 5):

| Set | Locus | n | exact | failed | Cys104 (of written) |
|---|---|---|---|---|---|
| human (self) | IGH | 198 | 198 | 0 | 198 |
| | IGK | 63 | 63 | 0 | 63 |
| | IGL | 79 | 79 | 0 | 79 |
| | TRA | 103 | 99 | 0 | 103 |
| | TRB | 117 | 117 | 0 | 117 |
| | TRG | 11 | 11 | 0 | 11 |
| | TRD | 6 | 6 | 0 | 6 |
| c57bl6 (self) | IGH | 99 | 99 | 0 | 99 |
| | IGK | 91 | 91 | 0 | 91 |
| | IGL | 3 | 2 | 1 | 2 |
| balbc (self) | IGH | 162 | 162 | 0 | 162 |
| | IGK | 98 | 98 | 0 | 98 |
| | IGL | 3 | 2 | 1 | 2 |
| c57bl6 from human | IGH | 99 | 99 | 0 | 99 |
| | IGK | 91 | 90 | 0 | 91 |
| | IGL | 3 | 3 | 0 | 3 |
| balbc from human | IGH | 162 | 162 | 0 | 162 |
| | IGK | 98 | 97 | 0 | 98 |
| | IGL | 3 | 3 | 0 | 3 |
| human IG from mouse | IGH | 198 | 198 | 0 | 198 |
| | IGK | 63 | 63 | 0 | 63 |
| | IGL | 79 | 62 | 13 | 66 |

Self totals: human 573/577 (99.3%), c57bl6 192/193 (99.5%), balbc 262/263 (99.6%). Floors in the test:

- self ≥ 99% exact per database;
- cross-species exact per locus at the measured rate rounded down to a whole percent: mouse from human IGH 100, IGK 98, IGL 100; human from mouse IGH 100, IGK 100, IGL 78;
- Cys104 ≥ 99% of the written alleles everywhere.

**What is still wrong** (10 alleles written with gaps that differ from IMGT's, 16 refused):

- Human self (4 wrong):
  - TRAV5\*01: FR1 gap at 9 instead of 8.
  - TRAV40\*01: FR1 gap at 7 instead of 8, and its CDR2/FR3 gaps shifted. Its nearest templates (TRAV9-2) have neither.
  - TRAV8-4\*06 and \*07 are 5′-partial alleles that start in or just before CDR1, and IMGT places their first residues apart from the rest. TRAV8-4\*07's four CDR1 residues follow the length rule (27, 28, 37, 38). TRBV7-9\*07's follow its full allele's positions, which the template rule reproduces. No one rule gives both.
- Mouse self (2 refused): the one IGLV3-like allele per strain (IGLV0-TXCK, IGLV0-FBPO). Its strain's only other lambda templates are IGLV1/2, so it fails `framework_insertion`. Without the per-pair 50-residue floor, both were written from a short wrong-frame alignment.
- Mouse IGK from human (1 wrong each): IGKV0-MAQO and IGKV0-JCV5 take IMGT's FR1 gap at 4–10 as a 5′-partial start (positions 1–7 empty).
- Human IGL from mouse:
  - 13 refused: human IGLV5 and IGLV11 have no mouse counterpart (`framework_insertion`).
  - 4 wrong: IGLV3-25\*02, \*03 and \*03_g74c put a CDR1/FR2 boundary residue in CDR1. IGLV3-21\*04 has the CDR2/FR3 boundary one residue off.

## Alternatives considered

- **A profile model from all templates (ANARCI-style).** More robust for distant species, but several times the code for a case not yet needed. It can come later behind the same interface.
- **Calling ANARCI.** It adds Python and HMMER dependencies, which the rewrite dropped, and it numbers protein only.
- **One gap open everywhere** (the starting point). It puts the indel that a CDR length difference needs into the flanking framework about as often as into the CDR. That gives a framework gap IMGT doesn't have, or a spurious `framework_insertion`. No setting of open 9–13 / extend 1–2 holds the self round trip at 99% without the CDR cost.
- **Right-aligning every partial CDR start.** It got neither of the two shipped alleles that start inside CDR1 right. Template positions get one.
- **Failing an allele with `template_disagreement`.** That would drop alleles the tool mostly gets right: of the 14 alleles flagged at margin 10, 12 are exact.
- **Counting failures against the Cys104 rate.** A refusal never writes a misplaced anchor, which is what the rate guards. Refusals are counted in the exact rate and listed in the report instead.

## Consequences

- An IMGT-gapped database can be built from ungapped or differently gapped V sequences of a species close to human or mouse. An allele the tool cannot place is reported and left out, never written with doubtful gaps.
- Accuracy depends on template coverage of the query's families. A family with no counterpart among the templates (human IGLV5/IGLV11 against mouse) fails rather than guessing.
- `template_disagreement` and `low_identity` are weak signals: together they flag 5 of the 10 wrongly written alleles above. A clean report does not prove the gaps are right.
- **Constants.** Changing any constant here changes the gapped output, so it needs a superseding record: the framework gap open 13, extend 1, the CDR gap open 5, BLOSUM62 (with `X` at −1), `TOP_N` 10, `LOW_IDENTITY` 0.60, `TIE_MARGIN` 10, the 50-residue floor, the CDR gap orders, and the default template list.
