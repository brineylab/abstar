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
- `v_identity` is the fraction of the read's V-derived bases that match the true V allele, which is `1 − v_mutation_frequency`. Inserted and deleted bases are not counted, so an aligner's identity on an indel read can differ slightly.

**`_valid` columns.** Each ambiguous field `F` has a list column `F_valid` holding every equally valid answer (0025): allele and gene sets as sorted names, boundaries as sorted positions. A null answer's validity is the separate boolean `d_call_null_valid`. Lists are never null; a segment absent from the read has an empty list.

**Events.**
- Trims: `v_trim_3`, `d_trim_5`, `d_trim_3`, `j_trim_5`. A negative value is that many P-nucleotides added. D trims refer to the D as used (reverse-complemented when `d_inverted`). D trims and `d_inverted` are null without a D.
- `mutations`: `"POS:FROM>TO"`, sorted by read position, one entry per position with the net change. FROM is the unmutated base and TO the base in the read. A position hit twice appears once, and a position mutated back to its original base does not appear. The ordered hit list is internal; graders compare observed substitutions.
- `indels`: `"ins|del:POS:LEN:BASES:MIN-MAX"`.
  - POS is the 1-based read position of the first inserted base, or of the first base after a deletion.
  - BASES are the inserted bases, or the deleted germline bases.
  - MIN-MAX is the range of equivalent POS values (0025).
  - An indel is listed only when the read has a base on both sides of it. An indel lost to, or cut by, truncation is unobservable: it is dropped, and the equivalence range is clipped to positions that stay in the read.
- `truncate_5p`, `truncate_3p`: bases cut from each end. `v_mutation_frequency`: the realized divergence of the read's V-derived bases (0 when there are none).

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

## Consequences

- Grading is a join on `sequence_id`. A true value is compared with the same-named column, and an ambiguous field with its `_valid` list.
- The same seed and parameters give byte-identical FASTA and TSV at any thread count and batch size. Changing the seeding scheme, the id form, or the order in which a read consumes random numbers changes the output for a given seed and needs a superseding record.
- `productive_only` costs about 1/(productive rate) simulations per read, and is refused when that rate is below roughly 1 in 1000.
- The truth schema lives in `truth.toml`, so the manual's field reference and the Arrow schema cannot drift from the rows.
