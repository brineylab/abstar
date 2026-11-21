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
- **CDR3** is null when the V keeps no anchor base, or the V allele has no anchor. It is also null when the J allele has no anchor, or the J's first base lies past the anchor. Otherwise it is the covered part: it is partial when the read stops inside it (ruling A4).
- **Without a J call** (the read stops before J): V ends at its alignment's end, CDR3 runs to the read's end, and FWR4 and the junction are null.
- **The junction** runs from the V anchor's first base through the Trp/Phe118 codon. It is null unless both anchor bases are in their segments and the Trp/Phe118 codon is in the read (`junction_undefined` in Task 10). The CDR3 can be partial while the junction is null.
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
- `productive`: `vj_in_frame`, no stop codon, a Cys104 codon of TGT or TGC, and a Trp/Phe118 codon of TGG, TTT or TTC.
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
| mutation 0–30%, seed 1 | 1,653 | 74.53% | 99.52% (413/415) | 96.73% | 99.88% | 77.92% |
| memory + 20% indels, seed 11 | 2,788 | 87.41% | 98.95% (944/954) | 98.89% | 99.86% | 90.03% |
| bnab, seed 1 | 1,301 | 60.11% | 99.35% (153/154) | 97.54% | 97.46% | 62.34% |
| naive + 40 nt C, seed 2 | 1,985 | 98.24% | 100.00% (946/946) | 98.74% | 100.00% | 100.00% |
| naive, truncated 0–120 nt 5′ and 0–60 nt 3′, seed 4 | 1,369 | 98.17% | 99.46% (732/736) | 99.27% | 99.93% | 100.00% |

- **Valid cuts that change a region (about 1.7% of naive reads).** When the V is trimmed near Cys104, or the J near Trp/Phe118, a chance N match lets the read support another V end or J start. These cuts are valid under record 0025. When that cut decides whether the anchor base is in its segment, or where FWR3 ends, the annotator's regions and junction differ from truth's origin-based ones. Like ruling A26, this is a read-level limit. `region_accuracy_naive` asserts exact regions with truth's cuts (≥ 99.5%) and at least 98% overall.
- **The V start (FWR1, `sequence_vdj`).** Mutated first bases are cut from the V alignment's 5′ end (record 0029's trim). FWR1 and `sequence_vdj` then start late, which accounts for most misses at 0–30% mutation and on bnab reads. Truth has no valid list for the V start.
- **Indel placement and the J end.** The remaining misses with truth's cuts are indels placed differently at a region boundary (memory + indels, bnab). There are also mutated last J bases cut from the J alignment's 3′ end, and J anchors trimmed away in reads that stop before J, which the read cannot show.
- The `R_length` fields cite this record. `sequence_alignment` and `germline_alignment` are ungraded until truth defines them.
