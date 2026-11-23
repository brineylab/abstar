# 0032: IMGT regions, truncation, alignment strings and productivity

## Status

Accepted (2026-11-21)

## Context

After the V, J, junction and C calls (records 0029–0031) the annotator derives the IMGT regions, the junction, the amino-acid fields, `sequence_alignment` and `germline_alignment`, productivity and the assembled sequences (spec *Pipeline* step 7). The generator derives the same fields from each base's origin (`crates/generate/src/derive.rs`, record 0026; productivity in `crates/generate/src/recombine.rs`). The grader compares regions exactly. The annotator reimplements those rules on its called alleles (record 0014): it never uses generator code.

Coordinates here are 0-based half-open on the oriented read. They become AIRR 1-based closed positions when the record is built (record 0016).

## Decision

**Region mapping** (`derive::regions`).
- Each read base given to V by the junction step (record 0030) maps to an ungapped position of the V primary. This uses the V alignment's ops, and past the alignment's end the V's last diagonal. The ungapped position then maps to the primary's IMGT gapped position. These are `VHit.imgt`'s positions, extended to the junction's V end.
- The IMGT gapped positions (1-based) are FWR1 1–78, CDR1 79–114, FWR2 115–165, CDR2 166–195, and FWR3 196–312, which ends with the Cys104 codon (310–312). A V base past 312 is CDR3.
- An inserted base takes the region of the base before it. So an inserted base at a region boundary goes to the upstream region (ruling A5).
- **The Cys104 codon** is the 3 read bases from the V anchor's first base. Its position is counted from the read's last V base, so it holds only when the V keeps that base (ruling A7). Those 3 bases are FWR3, whatever segment they came from.
- N and D bases after the Cys104 codon are CDR3. J bases before the J anchor (Trp/Phe118) are CDR3, and J bases from the anchor on are FWR4. Bases after J's end (C) have no region.
- A region is the run of bases with its label. It is null when no base of it is in the read. It is empty when its neighbours on both sides are in the read but every base of it is deleted (start and end are null, record 0026).
- **CDR3** is null when the V keeps no anchor base, or the V allele has no anchor. It is also null when the J allele has no anchor, or the J's first base lies past the anchor and no Trp/Phe118 codon is placed there (below). Otherwise it is the covered part: it is partial when the read stops inside it (ruling A4).
- **Without a J call** (the read stops before J): V ends at its alignment's end, CDR3 runs to the read's end, and FWR4 and the junction are null.
- **The junction** runs from the V anchor's first base through the Trp/Phe118 codon. It is null unless both anchor bases are in their segments and the Trp/Phe118 codon is in the read (`junction_undefined` in Task 10). The CDR3 can be partial while the junction is null.
- *Amended 2026-11-22 (#162): a J that starts past its anchor.* A mutated or inserted J 5′ end, or an aligner indel next to the codon, can move the J alignment's start 1–15 nt past the Trp/Phe118 codon. The codon is then placed on the read:
  - first, at the read position J's first diagonal projects the anchor to, if the 3 bases there are TGG, TTT or TTC;
  - otherwise, J from its anchor through its alignment's end is realigned (record 0029's scores, gaps anywhere) to the read from after the Cys104 codon through J's end. The codon is where J's anchor lands, if J's first base aligns, the 3 bases there are TGG, TTT or TTC in the V frame, they are at most 15 nt from the projection, and at least one more J base aligns before any gap. *(Amended again 2026-11-22, before 1.0.0: this replaces a rule that took the nearest V-frame codon starting the FWR4 motif `[WF]G.G`.)*

  Either way the codon starts after the Cys104 codon and ends within J. It ends the junction and CDR3, and FWR4 starts there, even where V's alignment reached into it by chance. J's start does not move. Without a placed codon the junction stays null.
- **Amino acids.** The reading frame is the V's: codons start in phase with the first V base's IMGT gapped position (`Regions::frame`).
  - `translate` drops a partial terminal codon. A codon with a base other than ACGT is `X`.
  - A region's translation is the codons whose third base lies in the region, so a codon crossing a boundary belongs to the downstream region. `R_length` is that translation's length in amino acids.
  - `junction_aa` translates the junction from its first base.

**Assembled sequences.** `sequence_vdj` runs from the first V base through J's end, or through the read's end when there is no J call (ruling A8). `sequence_vdjc` runs through C's end and is null without C.

**Alignment strings** (`derive::alignments`). Truth does not define these, so they are ungraded.
- They cover `sequence_vdj`. V carries the primary's IMGT gaps: `.` in both strings at every gapped position between two aligned V positions.
- N bases align to `-` in `germline_alignment` (the legacy convention). The D aligns to the called allele in its called orientation (record 0031), and J to the J primary, as the junction step placed them.
- An inserted base is `-` in `germline_alignment`. A deleted base is `-` in `sequence_alignment`.
- The `_aa` strings group the same columns in threes from the first V-frame codon, in both strings. Three `.` translate to `.`, and three `-` to `-`. A codon mixing gaps and bases is `X`. A partial terminal codon is dropped.

**Productivity** (`derive::productivity`), as the generator defines it.
- `vj_in_frame`: a junction exists, its length is a multiple of 3, and its Trp/Phe118 codon is in the V frame.
- `stop_codon`: TAA, TAG or TGA in the V frame over `sequence_vdj` (no C).
- `productive`: `vj_in_frame`, no stop codon, a Cys104 codon of TGT or TGC, and a Trp/Phe118 codon of TGG, TTT or TTC. (2026-11-21: also every junction base ACGT; see the final-review note at the end.)
- Without a junction, `vj_in_frame` and `productive` are F.

## Alternatives considered

- **Null CDR3 whenever the junction is null:** truth keeps the covered part (ruling A4), and the grader compares exactly.
- **Regions from the V primary's anchor alone, ignoring the junction step's V end:** a V trimmed into or past the Cys104 codon would still get a full FWR3 and a CDR3 that truth does not have.
- **Stop codons from the read's start**, as the generator scans: generated reads start with V, so this is the same on truth sets. On real reads it would count a leader before V.
- **Per-region frames:** one frame per read is what truth uses, so an out-of-frame indel shows in the downstream translations.
- **Truth's fallback for a V allele without IMGT gaps** (FWR4 only): every V allele in the built-in databases is IMGT-gapped, and record 0005 requires gapped sequences, so this path is not mirrored.

## Consequences

Agreement with truth (`tests/derive.rs`, seeds as listed) is measured on reads whose V and J calls are the true alleles, the grader's condition for region fields. With truth's V start, V end and J start, the regions match exactly: 100% on naive reads (918/918). The misses on the other sets come from upstream calls (below).

| Set | Graded reads | All regions | With truth's cuts | Junction | `productive` | `sequence_vdj` |
|---|---|---|---|---|---|---|
| naive, seed 1 | 1,986 | 98.29% | 100.00% (918/918) | 98.79% | 99.90% | 100.00% |
| mutation 0–30%, seed 1 | 1,654 | 94.68% | 99.60% (492/494) | 96.74% | 100.00% | 99.03% |
| memory + 20% indels, seed 11 | 2,788 | 96.59% | 98.77% (1,043/1,056) | 98.89% | 100.00% | 99.50% |
| bnab, seed 1 | 1,301 | 93.01% | 98.66% (221/224) | 97.54% | 97.46% | 96.77% |
| naive + 40 nt C, seed 2 | 1,985 | 98.24% | 100.00% (946/946) | 98.74% | 100.00% | 100.00% |
| naive, truncated 0–120 nt 5′ and 0–60 nt 3′, seed 4 | 1,369 | 98.17% | 99.46% (732/736) | 99.27% | 99.93% | 100.00% |

- **Valid cuts that change a region (about 1.7% of naive reads).** When the V is trimmed near Cys104, or the J near Trp/Phe118, a chance N match lets the read support another V end or J start. These cuts are valid under record 0025. When that cut decides whether the anchor base is in its segment, or where FWR3 ends, the annotator's regions and junction differ from truth's origin-based ones. Like ruling A26, this is a read-level limit. `region_accuracy_naive` asserts exact regions with truth's cuts (≥ 99.5%) and at least 98% overall.
- **The V start (FWR1, `sequence_vdj`).** Record 0029's 5′ start extension (ruling A28) gives V back the mutated first bases its cut drops. Without it, FWR1 agreed on 78% of 0–30% reads and 62% of bnab reads; with it, 99.4% and 96.8%. The remaining V-start misses are reads whose extension is too long (over 20 bases) or too unlike the allele, mostly bnab. Truth has no valid list for the V start.
- **Indel placement and the J end.** The remaining misses with truth's cuts are indels placed differently at a region boundary (memory + indels, bnab). There are also mutated last J bases cut from the J alignment's 3′ end, and J anchors trimmed away in reads that stop before J, which the read cannot show.
- The `R_length` fields cite this record. `sequence_alignment` and `germline_alignment` are ungraded until truth defines them.
- **2026-11-21 (rulings A31–A33, truth side):** truth now does three things the annotator already did. `vj_in_frame` is null without a junction. A trimmed anchor codon is placed from a valid V end or J start that the true allele supports, so the junction and CDR3 are defined. Regions follow each segment's normalized indel placement. On the accuracy gate (seed 20261121, 2,000 reads), the annotator rule is unchanged and agreement rose: `junction` on igh_naive from 98.75% to 99.75% and on tra from 99.15% to 99.95%; `fwr3` on igh_naive from 98.40% to 99.15%; `cdr2` on igh_bnab from 99.14% to 99.60%. The remaining `fwr3` misses are reads whose V is trimmed just short of Cys104 and has a valid V end one or two bases later. The annotator takes that end, so the extra base is FWR3 in its call, while truth counts it as an N base outside FWR3.
- **2026-11-21 (fix wave, ruling B2): Cys104 projection.** When the V's last base lies before the V anchor (every V end the read supports stops short of Cys104, typically a heavily mutated V 3′ end), the Cys104 codon is projected along V's last diagonal, at most 50 nt past V's last base. The projection is kept when the 3 read bases there are TGT or TGC and end before the J start. The 50-nt bound covers every projection measured on igh_bnab: 107 reads, distances 1–49, all at truth's anchor; 104 of them are within 30.
  - The projected codon is FWR3, as an anchor kept by V would be, and the junction and CDR3 start there. The V end itself does not move.
  - **Why the Cys requirement.** Without it, every V trimmed past Cys104 (truth: junction null) would get a junction, the IgBLAST and legacy convention. On the Task 12 sets that would wrongly add about 126 junctions on igh_naive and 70 on tra. With it, 110 of 111 bnab reads whose projected codon is TGT/TGC have truth's junction at the projected position.
  - **Before → after:**

    | Set | junction | FWR3 | CDR3 |
    |---|---|---|---|
    | igh_bnab | 92.30% → 94.50% | 96.87% → 99.03% | 97.10% → 99.22% |
    | igh_memory | 99.44% → 99.46% | | |
    | igh_truncated | 94.46% → 94.44% (one read) | | |
    | igh_naive, igk, igl, tra, trb | unchanged | | |
- **2026-11-21 (final review, ruling I5): ambiguity and motifs in productivity.**
  - **Ambiguous junction bases.** `productive` is F when the junction holds a base other than A, C, G or T (an `N` or IUPAC code from the sequencer): the read cannot show that the junction encodes a functional chain. AGENTS.md lists ambiguity among the productivity rules, and legacy abstar calls such reads nonproductive. Before, an `N` in CDR3 left `productive` T. `vj_in_frame` and `stop_codon` are unchanged. The generator emits only ACGT bases, so truth's productivity needs no change and the accuracy gate is unaffected.
  - **Motifs.** The Trp/Phe118 codon may be TGG, TTT or TTC on every locus, heavy or light, BCR or TCR. A locus-specific motif (Trp for IGH, Phe for light chains and TCRs) would call real rearrangements whose J allele carries the other residue nonproductive (2 human TRAJ alleles carry Trp118), and the generator's selection uses the same universal set. So the rule stays universal.
- **2026-11-22 (#162): the Trp/Phe118 codon when J starts past it.** On real data (the rc1 run over all 2,372,992 `jaffe-paired` chains), 1,356 chains (0.057%; 1,333 IGH) had V and J called but `junction_undefined`. Their J alignment started 1–12 nt past the codon, while legacy abstar gave each a junction.
  - **What fixes them.** With the placement rule (above), 1,334 now have legacy's junction, 1 differs, and 21 stay null. Seven of the 21 are IGKV2-29\*01, which has no Cys104 anchor in the database. Of the other 9,612 test-split chains, none changed.
  - **Why the motif.** Taking the nearest V-frame TGG, TTT or TTC without the motif matched legacy on 1,349. On the accuracy gate, though, it added wrong junctions where the gate's truth (the generator's records) has another one or none: a TGG in the D where truth's J is trimmed past its anchor, and V-frame codons in out-of-frame reads. Requiring the exact projection or the motif removed all of them.
  - **The gate.** No cell went down and 215 went up. The fallback fires on 15 gate reads, and all 15 match truth. `junction`: igh_bnab 98.57% → 98.81%, igh_bnab_c 98.43% → 98.73%, igk, igk_c and tra 99.95% → 100%, trb 99.85% → 99.90%.
- **2026-11-22 (#162, rc2 rerun): realigning J replaces the motif.** The rc2 run over all `jaffe-paired` chains left 21 of the 1,356 null: the 7 IGKV2-29\*01 chains and 14 others.
  - **Which 14.** Ten had J's own Trp/Phe with the next residues changed, which breaks the motif: a Gly119 substitution whose J alignment carried a spurious 1-nt insertion (4), a 3-nt insertion after the Trp (4), Gly→Asp with a spurious 3-nt insertion (1), and a 15-nt tandem duplication in a kappa FWR4 (1). In all ten the read stays in frame through C, so legacy's junction and productivity are right. In the other four, J is deleted through its Trp, and legacy took a Trp encoded by D or N 6–9 nt upstream.
  - **What realigning finds.** Realigning finds legacy's codon in 1,342 of the 1,356 (the motif found 1,334), including every chain the motif fixed. It leaves 6 null and the 1 difference is unchanged. Two of the six keep J's own codon, but another alignment scores as well or better: a TTG with a frameshift scores the same as legacy's TGG with three substitutions, and a mutated Trp beats legacy's 3-nt insertion after an N-region TGG. The other four are the D- or N-Trp chains, kept null on purpose (below).
  - **Why the conditions.** Without the V-frame condition the realignment took out-of-frame TGGs in two D/N-Trp chains. Without the aligned-base condition it took the D's TGG in two more, followed by a 6-nt insertion. A TGG with an insertion right after it is as likely an N or D codon as J's, and the gate's truth (the generator) has no junction when J is trimmed past its anchor. The projection window is 15 nt so it reaches the kappa duplication.
  - **The gate.** No cell went down and 48 went up, all in igh_bnab_c: `junction` 98.73% → 98.85%, `cdr3` 99.33% → 99.48%.
