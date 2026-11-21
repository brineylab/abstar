# 0026: Generate outputs: files, truth columns, event encodings, per-read seeding

## Status

Accepted (2026-11-20)

## Context

`abstar generate` writes synthetic reads and, for each one, a truth row that annotators are graded against by joining on `sequence_id`. The truth row mixes AIRR fields (the true generative values), the valid-answer sets of 0025, the recorded events, and difficulty flags. Each needs a fixed meaning and encoding. The output must also be reproducible: the same seed must give byte-identical files at any thread count or batch size.

## Decision

**Files.** `generate_to_files` writes into the output directory:
- `<name>.fasta`: one record per read, header `>sequence_id`, the sequence on one line (no wrapping).
- `<name>.truth.tsv` and/or `<name>.truth.parquet`: one row per read, in read order, with the columns of `crates/generate/truth.toml` in its order. The shared writers (0022) encode them: `T`/`F`, empty cells for null, lists joined with `,`.
- `generate.json`, written last. Its keys are `abstar_version`, `started` and `finished` (RFC 3339 UTC), `seed`, `parameters` (every `GenerateParams` field, with distributions in their `fixed:`/`uniform:`/`beta:` string form, plus `preset`, `formats`, `germline_db`, `db_path` and `model_path`), `model` (`name`, `sha256` of the model file), `database` (`name` and `alleles_sha256` from the database's `manifest.toml`), `counts` (`n`, `productive`, `attempts`) and `warnings`.

As with `run.json` (0018), an existing `generate.json` is deleted before any output file is created or truncated, so a failed rerun cannot look finished. Invalid parameters fail before that point, leaving the previous run's files intact. `generate_fasta_to_writer` writes FASTA only and treats a closed pipe as a normal end. It skips the truth sets, which FASTA does not need. `generate_in_memory` returns the truth rows as Arrow batches, at least one batch even for `n = 0`.

**Warnings.** Model-binding warnings are listed as they are. A per-read warning is listed once, as the text before its first `:`, with the number of reads that raised it, for example `SHM safety cap reached (12 sequences)`. This keeps `generate.json` small for million-read runs.

**`sequence_id`.** `g` followed by the 1-based read number, zero-padded to the number of digits in `n` (`g001` … `g500`). Ids sort in read order and never collide with each other.

**Seeding.** Read `i` (0-based) uses its own `ChaCha8Rng`, seeded with `seed_for(seed, i)`, where `seed_for(g, i) = splitmix64(g ^ splitmix64(i))`. Reads are simulated in parallel within each batch and written in index order, so output does not depend on threads or batch size. With `productive_only`, attempt `j` (0-based) of read `i` uses `seed_for(seed, i * 1000 + j)`. The first productive attempt is kept. If none of 1000 attempts is productive, the run fails with `INVALID_ARGUMENT` ("productive rate too low"), because the model and parameters cannot produce the requested reads. `counts.attempts` is the number of reads simulated, so `1 − n / attempts` is the rejection rate. A drawn seed is recorded in `generate.json`; when FASTA streams to stdout (no `generate.json`) the CLI prints it to stderr as `seed: N` before generating, even with `-q`.

**Python layout** (final review, ruling G15). `abstar.generate` is a plain function and `abstar.GenerateSummary` its return type, both exported lazily from `abstar/__init__.py` (implementation in the private `abstar._generate`), so `inspect.signature` and `help()` show the function. The recombination models have their own public module, `abstar.model`, mirroring `abstar.germline`: `list()`, `info(name_or_path)`, `resolve(name_or_path, species, locus)` and `ModelInfo`. `list()` lists a model that fails to load with its reason in `error`, as `germline.list()` does. An earlier callable-module trick (`abstar.generate` as a module with `__call__`, so `abstar.generate.list_models()` also worked) was dropped: type checkers and `help()` saw a module, not a function.

**Coordinates.**
- `*_sequence_start`/`*_sequence_end` are 1-based closed positions in `sequence`: the first and last read base derived from that segment. They are null when no base of the segment is in the read (truncated away, or a D trimmed to nothing).
- `*_germline_start`/`*_germline_end` are 1-based positions of those same bases in the ungapped allele. For an inverted D they are positions in the D's reverse complement, the sequence actually used, so start ≤ end and `d_trim_5` keeps its meaning.
- The `_valid` boundary lists of 0025 use the same read coordinates.

**True values and AIRR columns.** `v_call`, `d_call`, `j_call` and `c_call` are the alleles actually used, and `v_gene`… their genes. These have the same names and meanings as `abstar run`'s columns (`v_gene` is defined by 0008), so grading is a join. A call is the generative value even when truncation or trimming removed every base of the segment. `rev_comp` is always F.
- `junction` runs from the first base of the true V's Cys104 codon to the last base of the true J's Trp/Phe118 codon, as read bases (after mutation). It is null when either codon is not wholly in the read; the row is still written.
- `junction_aa` translates the junction from its first base, dropping a trailing partial codon.
- `productive`, `vj_in_frame` and `stop_codon` follow AIRR, from the read. Without a junction, `productive` and `vj_in_frame` are F.
- **`np1`/`np2` are the read bases strictly between the true segments' read spans**, after mutation: `np1` between the last V base and the first D base (or the first J base when the read has no D base), `np2` between the last D base and the first J base. They therefore include P-nucleotides and any mutation of N bases, which is what an annotator sees and reports. The N bases as generated are recoverable from the events and the database (`replay`). A value is null when a flanking segment is not in the read, and an empty string when the segments abut.
- `v_identity` is the fraction of the read's V-derived bases that match the true V allele, which is `1 − v_mutation_frequency`, so it is counted in the same V alignment. Inserted and deleted bases are not counted. `high_mutation` likewise follows that `v_mutation_frequency`.

**`_valid` columns.** Each ambiguous field `F` has a list column `F_valid` holding every equally valid answer (0025): allele and gene sets as sorted names, boundaries as sorted positions. A null answer's validity is the separate boolean `d_call_null_valid`. Lists are never null; a segment absent from the read has an empty list.

**Events.**
- Trims: `v_trim_3`, `d_trim_5`, `d_trim_3`, `j_trim_5`. A negative value is that many P-nucleotides added. D trims refer to the D as used (reverse-complemented when `d_inverted`). D trims and `d_inverted` are null without a D.
- `mutations`: `"POS:FROM>TO"`, sorted by read position, one entry per position with the net change. FROM is the unmutated base and TO the base in the read. A position hit twice appears once, and a position mutated back to its original base does not appear. The ordered hit list is internal; graders compare observed substitutions.
- `indels`: `"ins|del:POS:LEN:BASES:MIN-MAX"`.
  - POS is the 1-based read position of the first inserted base, or of the first base after a deletion.
  - BASES are the inserted bases, or the deleted germline bases.
  - MIN-MAX is the range of equivalent POS values (0025).
  - An indel is listed only when the read has a base on both sides of it. An indel lost to, or cut by, truncation is unobservable: it is dropped, and the equivalence range is clipped to positions that stay in the read.
- `truncate_5p`, `truncate_3p`: bases cut from each end. `v_mutation_frequency`: the realized divergence of the read's V-derived bases (null when there are none).

**Annotator truth fields** (regions, amino acids, masks, per-segment mutations). These carry the names and encodings of `abstar run`'s v1 fields, so each is graded by a join. They are valued for the true alleles and derived from the origin map and the alleles only (`crates/generate/src/derive.rs`), independently of the annotator. Where the annotator plan's rules left a case open, this record fixes it (controller rulings A3–A10 of the annotator plan's Task 1 review), and the annotator follows it.
- **Reading frame.** Every translation uses the V frame: codons start at read positions in phase with the IMGT-gapped position of the read's first V base. With a frameshifting indel the frame does not follow the germline codons after it; this is the continuous V(D)J query frame.
- **Regions.** On the true V's IMGT-gapped nt positions: FWR1 1–78, CDR1 79–114, FWR2 115–165, CDR2 166–195, FWR3 196–312 (ending with the Cys104 codon); CDR3 the bases after the Cys104 codon and before the Trp/Phe118 codon (V tail, P, N, D and J bases); FWR4 the true J's bases from the Trp/Phe118 codon to the J end.
  - An inserted base belongs to the region of the germline base it follows (the base its `ins` position names).
  - A region cut by truncation is partial. A region inside the read whose every base is deleted is empty: `R` and `R_aa` are `""`, `R_length` 0, `R_start`/`R_end` null. A region with no base in the read is null.
  - The Cys104 codon is the 3 read bases from the V anchor's first base, whatever their origin: when `v_trim_3` cuts into the codon but keeps its first base, the N or P bases filling IMGT 311–312 are FWR3. CDR3 runs from there to the J's Trp/Phe118 base, so `cdr3 = junction[3..len−3]` whenever `junction` is set (productivity needs only the anchors' first bases).
  - CDR3's ends are known when both anchors' first bases survive trimming. A partially covered region holds the covered part, so when truncation removes one anchor codon CDR3 holds the covered bases between the known boundary and the read end (start, end and `cdr3_length` of what is present, partial terminal codons dropped), while `junction` stays null. CDR3 is null when trimming removed either anchor's first base or an allele has no anchor.
  - FWR4 is null when the J has no anchor. With no V base in the read it is still valued, translated in the J frame (phased by the Trp/Phe118 codon); every other region field is then null.
  - `R_aa` holds the codons (first V base to the last read base that is not a dangling insertion) whose 3rd base is in R, so a codon crossing a boundary belongs to the downstream region and partial codons are not translated. `R_length` is the length of `R_aa`, so `cdr3_length` counts the CDR3 in the V frame, not `junction_aa` length − 2 on an out-of-frame read.
  - Every region except FWR4, `sequence_aa`, the `sequence_vdj*` fields and the masks are null when the read has no V base or the V allele is not IMGT-gapped.
- **Sequences.** `sequence_aa` translates the whole read in the V frame. `sequence_vdj` runs from the first V base to the last J base, or to the last V/N/D-derived base when the read has no J base (null without a V base); `sequence_vdjc` runs to the last constant-region base (null without one). Each `_aa` translates its sequence in the V frame.
- **Masks** cover only the V(D)J region (user decision A12): nt masks index `sequence_vdj` (one character per base) and aa masks `sequence_vdj_aa` (one per codon); they are null when that sequence is null. No mask has a constant-region base or a `C` label.
  - `cdr_mask`: `1`/`2`/`3` for CDR1/CDR2/CDR3 (a partial CDR3 included), `0` for framework bases. Null when CDR3's ends are undefined (as for a null CDR3 above, apart from truncation).
  - `gene_segment_mask`: `V`, `D`, `J`; `N` for N and P nucleotides. An inserted base takes the segment it sits in (always V or J).
  - `nongermline_mask`: `1` for N, P and inserted bases and for bases differing from the germline base behind them.
  - **Amended 2026-11-21 (ruling A29):** `nongermline_mask` and `nongermline_mask_aa` now follow each segment's alignment, the one the mutation lists use (indels at the indel strings' placement, rulings A16, A17). A base is `0` only when that alignment pairs it with an equal germline base. The bases it sets aside as inserted, and N and P bases, are `1`. So a segment's `1`s are its substitutions plus its inserted bases, on both truth's side and an annotator's. Before, the masks followed the generated placement: on 36 of 577 indel reads (`memory`, 20% indels, seed 11) truth's V nongermline count differed from its own `v_mutation_count` plus inserted bases. `gene_segment_mask`, `cdr_mask` and the regions still follow the origin map. An inserted base takes its segment either way, so `gene_segment_mask` does not depend on the placement.
  - The `_aa` masks have one character per codon of the matching translation: the 3rd base's character for `cdr_mask_aa` and `gene_segment_mask_aa`; for `nongermline_mask_aa`, `1` when a base is N, P or inserted or the codon's amino acid differs from the translation of the germline bases behind it.
- **Per-segment mutations** (`v`, `d`, `j`, `c`).
  - All of a segment's fields are null when no base of it is in the read (truncated away, trimmed to nothing, or no such segment).
  - nt: `IMGT{pos}:REF>ALT` for V (IMGT-gapped position), `{pos}:REF>ALT` for D, J and C (1-based in the allele; for an inverted D, in the reverse complement), in germline order. REF is the germline base, ALT the read base.
  - aa: over germline codons whose three bases are all in the read; REF translates the germline codon, ALT the read bases aligned to it (an inserted base inside the codon is ignored). Codons: V by IMGT codon (`IMGT{n}`; gap codons skipped); J in phase with the Trp/Phe118 codon; D in the V frame at the read's last V base, the frame D occupies in the junction (null without a V base in the read); C continuing the J's Trp/Phe118 frame across the J end, as the spliced transcript reads (null when the J has no anchor). Codons are numbered from 1 at the allele's first complete codon in that frame.
  - Counts are list lengths. Every per-segment mutation field (nt and aa lists, counts, frequencies and their denominators) and the indel strings come from one alignment per segment, the one with the indel strings' placement (rulings A13, A15, A16, A17): the segment's read bases, from its first to its last and minus the inserted bases at the insertion strings' placement, keep their order and fill the germline positions from its first aligned one on, skipping those deleted at the deletion strings' placement. That is what an aligner-based annotator sees. Between an indel's generated placement and that one, bases are therefore numbered (and compared) at their aligned germline positions, not their origins, so a list can differ there from the recorded events (`mutations`); every aa change still has an nt change in its codon. `S_mutation_frequency` is the count ÷ the segment's germline positions in the read (deleted positions and inserted bases excluded), null when there are none; `S_mutation_frequency_aa` is the aa count ÷ codons counted, null when none is. `c_germline_start`/`c_germline_end` are the 1-based first and last constant-allele positions in the read (null without a C base), so the C fields are gradeable. The existing `v_mutation_frequency` has this definition and serves as the V member; it too is null without a V base or when the V allele is not IMGT-gapped, like every V field (`high_mutation` stays F).
  - Indels: `{pos}:ins:BASES` (pos the germline position the inserted bases follow) and `{pos}:del:BASES` (pos the first deleted germline position, BASES the deleted germline bases), with V positions as `IMGT{pos}`. Each is placed where an optimal read-vs-germline aligner puts it (ruling A17, superseding A3's VCF-style left normalization against the unmutated germline): among every placement of that indel (same kind and length) in the segment, with the segment's other indels fixed and its indels processed left to right, the one with the fewest mismatches between the read and the germline over the segment, the leftmost on ties. A placement must have an aligned read base of the segment on both sides and must lie in the zone where the generator places indels (ruling A18; never in CDR3, record 0025): in V at least 3 nt upstream of the Cys104 codon, in J at least 3 nt downstream of the Trp/Phe118 codon. BASES of an insertion are the read bases that placement sets aside. The `indels` column keeps the generated placement and its MIN-MAX equivalence range. An indel is listed only with a read base on both sides, as in `indels`. Indels are generated only in V and J, so D and C lists are empty when present. Masks and regions still follow the bases actually inserted or deleted (except the nongermline masks since ruling A29, above).

**Difficulty flags** (extending the labels of 0025):
- `high_mutation`: `v_mutation_frequency` ≥ 0.20.
- `long_cdr3`: CDR3 length (`junction_aa` length − 2) ≥ 24 aa for IGH, ≥ 16 aa for other loci. Null without a junction.
- `has_indel`: `indels` is not empty.
- `truncated_5p` / `truncated_3p`: that end was cut.
- `ambiguous_v` / `ambiguous_j`: the allele set has more than one member.
- `ambiguous_d`: more than one D answer is valid, counting null when `d_call_null_valid`. A short D with one tied allele is ambiguous, but a VJ read is not.

## Alternatives considered

- **`np1`/`np2` as generated (before mutation, without P):** these are not what an annotator can observe, and they would be unrelated to the boundaries graded. They remain recoverable from the events.
- **Raw mutation hits in order:** a position hit twice would appear twice, and a back-mutation would list a change the read does not have. Grading wants observed substitutions; `replay` keeps the exact hit list internal.
- **Keep indels outside the read, with read positions out of range:** these would be unobservable events that no annotator could report, and they would set `has_indel` on reads with no indel.
- **D germline coordinates in the forward allele for an inverted D:** start would exceed end, and they would disagree with `d_trim_5`/`d_trim_3`.
- **One RNG stream per batch or per thread:** output would depend on the batch size or the thread count.
- **Retry productive-only reads forever:** a model that cannot produce productive reads would hang the run.
- **Delete `generate.json` before validating parameters:** a typo would destroy a finished run's provenance without producing anything.
- **The generated indel placement in the per-segment lists:** an aligner places an indel at one end of its equivalence range, so exact grading would fail on equivalent placements. The masks and regions keep the generated placement (the nongermline masks moved to the lists' alignment on 2026-11-21, ruling A29; see the amendment above).
- **VCF-style left normalization against the unmutated germline (ruling A3, superseded by A17):** it ignores the mutations, so where shifting the indel moves a mutated base onto a germline base it matches, an optimal aligner places the indel elsewhere (on `abstar.generate(4000, "IGH", seed=11, preset="bnab", indels=0.5)`, 188 of 957 V deletions and 142 of 836 V insertions had a placement with fewer mismatches), and a correct annotator would be graded wrong. Ties within the mutated read also extend past the unmutated equivalence range; the leftmost tie is the aligner convention.
- **A null CDR3 when truncation removes one anchor codon:** it would break the rule that a partially covered region holds the covered part. A CDR3 whose anchor codon was trimmed stays null, because its boundary is not in the germline as used.
- **Constant-region bases in the masks:** they are outside the annotated V(D)J region, so no mask covers them (the earlier draft indexed `sequence_vdjc` and labelled them `C`/`0`).
- **Germline-frame codons for region and mask amino acids:** on frameshifted reads they would not match the read's own translation (`sequence_vdj_aa`), which an annotator reports.

## Consequences

- Grading is a join on `sequence_id`. A true value is compared with the same-named column, and an ambiguous field with its `_valid` list.
- The same seed and parameters give byte-identical FASTA and TSV at any thread count and batch size. Changing the seeding scheme, the id form, or the order in which a read consumes random numbers changes the output for a given seed and needs a superseding record.
- `productive_only` costs about 1/(productive rate) simulations per read, and is refused when that rate is below roughly 1 in 1000.
- The truth schema lives in `truth.toml`, so the manual's field reference and the Arrow schema cannot drift from the rows.
- Adding the annotator truth fields changed the truth TSV bytes for a given seed (not the FASTA): the generator's pinned TSV hash was updated with them.
- Ruling A29 (2026-11-21) changed the truth TSV bytes for a given seed (not the FASTA) on reads with an indel: the pinned TSV hash was updated.
