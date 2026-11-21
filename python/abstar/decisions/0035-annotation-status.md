# 0035: Annotation status and record assembly

## Status

Accepted (2026-11-21)

## Context

The annotator's steps (records 0029–0034) each stop when they cannot continue: no V, no J, a missing anchor. Records 0004 and 0016 fix one row per input record, the statuses, and AIRR coordinates converted once. The output needs a rule for which status each stopping point gives, which fields each status still fills, and what the score and identity fields measure. The V alignment's score (`VHit.aln.score`) mixes the aligner's DP score with the 5′ extension's columns (record 0029), so it is not a usable `v_score`.

## Decision

*Amended 2026-11-21 (record 0036): a database whose V genes are not IMGT-gapped is refused when it is built or loaded, so the status table has no `macaque_not_imgt_gapped` row and there is no V-only path. The removed row read: "V allele not IMGT-gapped (see below) / `unassigned` / `macaque_not_imgt_gapped` / `rev_comp`, `locus`, the V call, ties, species, V coordinates, CIGAR, score, identity".*

**Statuses and failure reasons** (`Engine::annotate`, `record.rs`):

| Stopping point | `annotation_status` / `failure_reason` | Fields filled |
|---|---|---|
| Empty sequence; a character outside `ACGTN` + IUPAC | `invalid_input` / `empty_sequence`, `invalid_characters` | identity only |
| No V on either strand | `unassigned` / `no_v_hit` | identity only |
| V found, no J | `unassigned` / `no_j_hit` | everything V supports: the V call and alignment, FWR1–FWR3 (and the partial CDR3, ruling A4), V mutations, `sequence_vdj` (to the read's end, ruling A8), `gene_segment_mask` and `nongermline_mask` (and their `_aa` forms), and `cdr_mask` only when V reaches the Cys104 codon. With no J to end it, the V(D)J region runs to the read's end, so the partial CDR3, `cdr_mask` and `sequence_vdj` can include J or C bases the J step missed. (2026-11-21, final review: this cell ended "and its masks".) |
| V and J found, junction undefined (an anchor missing or truncated) | `unassigned` / `junction_undefined` | every call and every region the read covers; `junction` null |
| Everything | `annotated` / null | all |
| A panic while annotating | `error` / `internal_error: <msg>` | identity only (record 0004) |

So a record never invents what a step could not establish, and never drops what an earlier step did establish. There is no V-less J path: without V there is no frame or locus to constrain J (record 0029).

**IMGT gapping.** *Amended 2026-11-21 (record 0036): every V allele of a loaded database is IMGT-gapped, so the annotator checks nothing per allele; the text below is the superseded rule.* ~~A V allele is IMGT-gapped when its gapped sequence has `.` gaps, except a macaque V allele without an anchor. The built-in macaque IGH V genes (KIMDB) carry gaps that do not follow IMGT numbering: the database builder finds no Cys104 at IMGT 310–312 for 1,063 of 1,071 of them. Regions, mutations and masks from those gaps would be plausible and wrong, so those reads stop with the V call. The species clause goes when #137 regaps them. The 8 macaque V alleles with anchors (IMGT GENE-DB) annotate normally.~~

**Coordinates.** Internal positions are 0-based half-open on the oriented read. `record.rs` converts each once: `start = a + 1`, `end = b`. `*_germline_start/end` are 1-based on the ungapped allele as used (the reverse complement for an inverted D). An empty region (every base deleted) has `""`, length 0 and null start and end.

**Per-segment alignment fields.** Each segment's `_sequence_start/end`, `_germline_start/end`, `_cigar`, `_score`, `_identity` and its mutations come from one set of columns: V and J as the junction step cut them, the D as it placed it, C over its covered prefix (`derive::segments`).
- `_score` is the score of those columns under record 0029's scoring: +2 per match, −3 per mismatch (a non-ACGT read base mismatches), −10 to open and −1 to extend a gap. It is not the aligner's DP score, which covers a different stretch than the reported one.
- `_identity` is matches over aligned (non-gap) columns, so it equals 1 − the nt mutation frequency over the same bases.

**Calls.** `<seg>_call` and `<seg>_call_ties` hold plain allele names; `<seg>_species` and `<seg>_species_ties` come from the alleles, entry by entry (record 0007).

**Productivity.** `productive` and `stop_codon` are filled whenever there is a V(D)J region; `productive` is false without a junction, as in truth. `vj_in_frame` is null without a junction.

**Loci.** The `loci` run option (`--locus`) restricts the indexes the engine builds, so other loci are never searched. J, D and C candidates still come from the V's locus.

**Engine.** `Engine` borrows the database and builds the V, D, J and C indexes once per run. It is `Sync` and allocates nothing per allele per read.

## Alternatives considered

- **`unassigned` records with only the status:** throws away a V call (or V and J calls) that the steps established and that users filter on.
- **`annotated` with null fields for truncated reads:** a downstream filter on `annotation_status == "annotated"` would then admit records without a junction.
- **The aligner's DP score as `v_score`:** it covers the overlap alignment before the cut, so it disagrees with the reported coordinates and identity.
- **Flagging every V allele without an anchor that reaches IMGT 312:** would also stop reads on four human alleles whose gaps are IMGT's but whose Cys104 is mutated (IGKV2-29*01, IGLV3-2*03, TRBV7-3*02/03), under a macaque-named reason. (2026-11-21: moot, record 0036 refuses a non-IMGT-gapped database; an anchorless V allele of a valid database takes the normal path.)

## Consequences

- Smoke check (Python `abstar.grade`, 2,000 generated human IGH reads, seed 1): naive 1,959 annotated and 41 `junction_undefined`, all with a null truth junction. At 0–30% mutation: 1,885 annotated, 107 `junction_undefined` and 8 `no_j_hit`. V, J and D call accuracy match the step-level measurements in records 0029–0031.
- ~~A `macaque` run annotates only the reads on its 8 anchored V alleles until #137.~~ (2026-11-21: macaque is no longer shipped; record 0036.)
- The stub annotator remains for the engine seam tests only.
