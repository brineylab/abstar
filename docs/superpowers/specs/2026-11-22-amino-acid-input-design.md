# Amino-acid input: design

Status: draft (2026-11-22). Decision record to write: `0039-amino-acid-input.md`.

## Intent

Let `abstar run` annotate antibody and TCR sequences given as protein, not nucleotides. Many
sequences exist only as protein: structures, patents, therapeutics, published mAbs. Today the
only option is to back-translate them, which invents codons, mutations and D calls. The
primary users are structure and patent datasets.

The output stays an AIRR Rearrangement table with the same columns. The one hard break with
AIRR is that `sequence` (required, non-nullable) is null, because there is no nucleotide
sequence to report. Every field that is about nucleotides, or whose coordinates are defined on
`sequence`, is null too. Every amino-acid field is filled the way nucleotide mode fills it, as
far as the protein supports it.

## Goals

- Gene and allele calls (V, J), locus, IMGT regions, CDR3/junction, amino-acid mutations,
  amino-acid masks and productivity for protein input.
- Same entry points (`abstar run`, `abstar.run`), outputs (TSV, Parquet, DataFrame), statuses
  and `run.json` as nucleotide mode.
- An accuracy gate for amino-acid mode built from the existing `abstar generate` truth sets.

## Non-goals (v1)

- **D calls.** A D contributes a few residues; a call from them is noise. `d_*` fields are null.
- **C calls.** `c_*` fields are null, even though structure data often includes CH1 (deferred
  by decision, 2026-11-22). See Follow-ups.
- **Nucleotide indel fields in amino-acid mode.** `*_insertions`/`*_deletions` use IMGT
  nucleotide positions and bases; null. Amino-acid indels get their own fields (below).
- **Amino-acid coordinate fields** (`*_start_aa`-style). Not in v1.
- **Auto-detecting** protein input. The mode is always explicit.
- **SIMD alignment and a k-mer prefilter** for protein. Scalar alignment first; see Performance.
- **`abstar generate` protein output.** The gate translates nucleotide truth instead.

## Anti-goals

- **Any change to nucleotide-mode results.** Nucleotide values must stay identical. The existing
  accuracy gate and engine tests must pass with no floor rewrites. The one allowed change is
  the schema: `sequence` becomes an optional (nullable) column in every output, nucleotide
  included. `schema_matches_toml` in `crates/annotate/tests/writers.rs` is updated to match.
- **Overloading a field's meaning by mode.** A field whose AIRR meaning is nucleotide
  coordinates is null in amino-acid mode, never re-used for residue positions.
- **Fabricated nucleotide data** (back-translation in any form).

## Constraints

- AGENTS.md rewrite rules: new code in `crates/` and `python/abstar/`; `legacy/` untouched.
- Fields are single-sourced in `crates/annotate/fields.toml` (record 0012); descriptions there
  are the docs.
- Lower-bound dependency versions only. No new dependency is needed: `block_aligner` (already a
  dependency) exports `scores::BLOSUM62` with `Matrix::get(a, b)`.
- Germline database format unchanged; built-in databases are not rebuilt.

## Approaches

1. **A second `Annotate` implementation, `AaEngine` (chosen).** The pipeline calls annotators
   through `trait Annotate` (`crates/annotate/src/engine.rs`). `AaEngine` sits beside `Engine`;
   the pipeline, writers, schema, `run.json`, Python API and CLI are shared. The nucleotide path
   is not touched, so it cannot change. Cost: amino-acid versions of assignment and derivation
   (regions, mutations, masks, record assembly), each much smaller than its nucleotide
   counterpart because there is no frame, strand, D or codon logic.
2. **Thread an alphabet parameter through `Engine`, `assign` and `derive`.** Rejected: nearly
   every nucleotide function carries codon or strand logic (`orient`, `KmerIndex`'s 2-bit
   encoding, the SIMD sentinels, frame and Cys projection in `regions`, D placement), so every
   one would grow a branch, and every branch is a chance to change nucleotide results.
3. **Back-translate protein and run the nucleotide engine.** Rejected: invents codons, so
   nucleotide mutations, D calls, NP regions and allele calls would all be artifacts.

Within approach 1, the aligner: a new scalar protein aligner (chosen) rather than generalizing
the nucleotide aligner (its scalar backend is held to tie-break parity with the SIMD backend,
and `is_match`/`place_gaps` are ACGT-specific) or using block-aligner's `AAMatrix` SIMD path
(its padding and sentinel bytes `E` and `H` are amino acids; the trick needs rework).

## Design

### Mode selection

- Rust: `RunOptions.sequence_type: SequenceType { Nt, Aa }`, default `Nt`. The pipeline builds
  `Engine` or `AaEngine` from it. `run.json` records it under `parameters.sequence_type`.
- Python: `abstar.run(..., sequence_type="nt")`, values `"nt"` / `"aa"`, anything else
  `INVALID_ARGUMENT`. In `"aa"` mode a bare string that is not an existing path is a sequence
  when it consists only of amino-acid characters (below); `_classify` takes the mode.
- CLI: `abstar run --sequence-type {nt,aa}`, default `nt`.
- One run is one type. No per-record detection.

### Input validation (`invalid_input` in amino-acid mode)

Uppercased. Allowed: the 20 standard residues, `X`, `B`, `Z`, `J`, `U`, `O`, and `*` (stop).
Anything else (including `-` and `.`) is `invalid_characters`; empty is `empty_sequence`.
`invalid_encoding` is unchanged.

### Germline translation (once, in `AaEngine::new`)

From the existing nucleotide alleles, honoring `loci` and `include_nonfunctional` exactly as
`SegmentIndex::for_loci` does:

- **V:** translate `sequence_gapped` codon by codon; a `...` codon becomes `.`, keeping IMGT
  numbering (gapped residue index `k` is IMGT position `k + 1`). A trailing partial codon is
  dropped. A codon mixing `.` and bases, or a V allele with no `.` at all (not IMGT-gapped),
  excludes the allele from amino-acid mode. Cys104 is gapped residue index 103.
- **J:** translate the ungapped sequence from `anchor % 3`; partial codons at either end
  dropped. The Trp/Phe118 residue is index `anchor / 3`. Residue `k` is labeled `k + 1`, the
  same labels nucleotide mode's J codons use. A J allele without an anchor is excluded.
- **D, C:** not used.
- Non-ACGT bases in a codon translate to `X`.

Suggested: align against the distinct V amino-acid sequences (many alleles are synonymous) and
expand each hit back to its alleles. This cuts alignment work and makes synonymous alleles tie
by construction.

### Alignment

New module `crates/annotate/src/align/protein.rs`: scalar semi-global affine-gap (Gotoh)
alignment of query vs germline residues, query ends and germline ends free, scored with
BLOSUM62, gap open −11 and extend −1 (BLAST protein defaults; calibrate if the gate says so).
`.` (IMGT gap) characters are removed before aligning and mapped back afterwards. Returns the
existing `Alignment`/`Op` types (coordinates in residues). `Match` means identical residues,
`Mismatch` any other aligned pair. Gap placement is deterministic: among equal-scoring
alignments, gaps go leftmost (record 0029's convention).

### Assignment

- **V.** Align the query against every candidate V. Best score wins if it reaches
  `MIN_V_SCORE_AA` (calibrated); else `unassigned` / `no_v_hit`. Trim the alignment to its
  maximum-scoring contiguous stretch (the residues that came from V), as `trim` does for
  nucleotides. Ties: same-locus V alleles with the fewest mismatches over those residues
  (record 0029); primary by record 0008's tiebreak (`primary_by_tiebreak`). Locus is the
  primary's.
- **J.** Align the query after V's last residue (or after Cys104, whichever is earlier, so V
  can be cut back) against the locus's J alleles; best score must reach `MIN_J_SCORE_AA`
  (calibrated), else `no_j_hit`. Trim the start only (`keep_end`). Ties and primary as for V.
- **V/J boundary.** If V's trimmed end overlaps J's trimmed start, cut at the position that
  maximizes V-prefix plus J-suffix score; equal scores take the earliest cut. Residues between
  V's end and J's start are non-germline (`N` in `gene_segment_mask_aa`).

### Regions and junction

V residues take their region from their IMGT position: FWR1 1–26, CDR1 27–38, FWR2 39–55,
CDR2 56–65, FWR3 66–104 (ends with Cys104, as in nucleotide mode), then CDR3 to the residue
before Trp/Phe118, then FWR4 from Trp/Phe118 to J's end. An inserted residue takes the region
of the residue before it.

- **Cys104:** the query residue aligned to the V allele's IMGT 104. If V's alignment ends before
  104 and within 16 residues of it, project along V's last diagonal and accept the residue only
  if it is `C` and before J's start (the residue-level form of ruling B2).
- **Trp/Phe118:** the query residue aligned to J's anchor residue. If J's alignment starts past
  its anchor, project along J's first diagonal and accept the residue if it is the J allele's
  own (`W` or `F`); otherwise look back up to 5 residues from the last-diagonal projection
  (the residue-level form of #162 / record 0032).
- **CDR3** is null if either anchor is missing; **junction_aa** is Cys104 through Trp/Phe118
  inclusive, null unless both are present.

### Fields

Filled in amino-acid mode:

| Field(s) | Value |
|---|---|
| `sequence_id`, `annotation_status`, `failure_reason`, `source_file` | as in nucleotide mode |
| `sequence_aa` | the input, uppercased |
| `locus`, `v_*`/`j_*` `_call`, `_gene`, `_species`, `_call_ties`, `_species_ties` | as in nucleotide mode |
| `v_score`, `v_identity`, `j_score`, `j_identity` | BLOSUM62 score and residue identity over the reported residues |
| `fwr1_aa` … `fwr4_aa`, `cdr1_aa` … `cdr3_aa` | region residues |
| `fwr1_length` … `cdr3_length` | residue counts (already amino-acid lengths in nucleotide mode) |
| `junction_aa` | Cys104 through Trp/Phe118 |
| `sequence_vdj_aa` | first V residue through J's last residue (the read's end without J) |
| `sequence_alignment_aa`, `germline_alignment_aa` | V with the primary's IMGT gaps (`.` in both), J as aligned, an inserted residue `-` in germline, a deleted one `-` in sequence, non-germline residues against `-` |
| `v_mutations_aa`, `j_mutations_aa` | `IMGT{pos}:REF>ALT` (V), `{pos}:REF>ALT` (J), as in nucleotide mode |
| `v_/j_mutation_count_aa`, `_mutation_frequency_aa` | count; count over aligned germline residues |
| `v_insertions_aa`, `v_deletions_aa`, `j_insertions_aa`, `j_deletions_aa` (new) | `IMGT{pos}:ins:RESIDUES` / `IMGT{pos}:del:RESIDUES` (V), `{pos}:ins:…` / `{pos}:del:…` (J); POS is the germline residue the insertion follows, or the first deleted germline residue. Same run rule as nucleotide mode: a gap run with an aligned residue on each side. Empty list when the segment has none |
| `cdr_mask_aa`, `gene_segment_mask_aa`, `nongermline_mask_aa` | one character per residue of `sequence_vdj_aa`; segments are `V`, `N`, `J` only |
| `stop_codon` | `*` in `sequence_vdj_aa` |
| `productive` | junction present, no stop, `junction_aa` starts with `C` and ends with `W` or `F`, and has no `X` |

Null in amino-acid mode: `sequence`, `rev_comp`, `vj_in_frame`, every cigar, every
`*_sequence_start/end` and `*_germline_start/end`, every region `_start/_end`, `junction`,
`junction_length`, `np1`, `np1_length`, `np2`, `np2_length`, `sequence_vdj`, `sequence_vdjc`,
`sequence_vdjc_aa`, `sequence_alignment`, `germline_alignment`, every nucleotide
region/mask/mutation field, every nucleotide `*_insertions`/`*_deletions`, and all `d_*` and
`c_*` fields.

`germline_alignment_aa` and `gene_segment_mask_aa` differ from nucleotide mode for heavy chains
by design: there are no D residues, so the CDR3 middle is non-germline.

### Schema and records

- `fields.toml`: `sequence` becomes `nullable = true`, description noting it is null for
  amino-acid input (record 0039). Descriptions of fields whose meaning changes in amino-acid
  mode say so: `sequence_aa`, `productive`, `stop_codon`, `vj_in_frame`,
  `germline_alignment_aa`, `gene_segment_mask_aa`.
- `Record.sequence` becomes `Option<String>`; constructors update mechanically.
- Four new `abstar`-origin fields, `v_insertions_aa`, `v_deletions_aa`, `j_insertions_aa`,
  `j_deletions_aa` (`list<string>`, nullable, `decision = "0039"`), each placed after its
  nucleotide counterpart. Null in nucleotide mode in v1 (see Follow-ups). Adding columns is
  the only other nucleotide-output change.
- This makes `sequence` nullable in nucleotide-mode Parquet/Arrow schemas too: nucleotide
  values are unchanged, but the column is no longer marked required. A per-mode schema would
  break the single-sourced `fields.toml` (record 0012). The decision record notes this for
  Parquet readers that check schemas strictly.
- Failure records (`invalid_input`, `invalid_encoding`, panics in `annotate_guarded`) put the
  raw input in `sequence_aa` in amino-acid mode. Suggested: give `Annotate` a `failed(id, seq,
  status, reason)` method that `Engine` and `AaEngine` each implement, so the pipeline stays
  mode-blind.
- Amino-acid output does not pass `airr` validation (null `sequence`). Documented, not tested.

### Accuracy gate

`tests/python/test_accuracy.py` gains an amino-acid section. It reuses the truth frames the
nucleotide gate already generates (every set in `SETS`, IG and TR; no extra generation),
filtered to truth `vj_in_frame == True`. Out-of-frame reads are dropped because truth
`sequence_aa` past their junction is not real protein. In-frame reads with a stop codon stay,
so `productive` is graded in both directions. Each kept row's truth `sequence_aa` is annotated
with `sequence_type="aa"`.

Grading is a small helper in the test file, not `abstar.grade`. `abstar.grade` and `_grade.py`
are unchanged: their span, same-allele and nucleotide tie-set rules do not fit amino-acid
output. Per read, the helper checks:

- `v_call` / `j_call`: truth's allele is in the amino-acid `*_call_ties`.
- `v_gene` / `j_gene`: truth's gene is among the genes of the amino-acid tie set.
- `locus`, `productive`, `stop_codon`: equality.
- `cdr3_aa`, `junction_aa`, every region `_aa` and `_length`: equality, graded only on reads
  whose V call (for V regions) or J call (for FWR4) is correct by the rule above, so that
  synonymous alleles still count.
- `v_mutations_aa`: the entries at IMGT 1–100 equal truth's entries at IMGT 1–100. This is
  away from the V/N boundary, where amino-acid mode places the cut itself.

Floors go in a separate file, `tests/python/data/accuracy_floors_aa.toml`, with the same
measured − 0.5 pp rule and the same `ABSTAR_WRITE_FLOORS` flow.

Thresholds (`MIN_V_SCORE_AA`, `MIN_J_SCORE_AA`, gap costs) are chosen by true accuracy on these
sets, leaning permissive; null rates are informational only. A negative control (a handful of
non-immunoglobulin proteins) should come back `no_v_hit`.

### Performance

Scalar alignment against every distinct V amino-acid sequence is a few million cells per query;
fine for antibody panels. Mark it with a `ponytail:` comment naming the ceiling; the upgrade is
an amino-acid k-mer prefilter (5-bit residues, k ≈ 3) and then block-aligner's `AAMatrix`. Report
measured throughput on a 2,000-sequence gate set in the PR.

### Tests

- Rust: protein aligner unit tests (scoring, leftmost gaps, free ends); `AaEngine` end-to-end
  on a translated IGH read checking calls, regions, `junction_aa`, mutations, masks; failure
  statuses (`no_v_hit`, `no_j_hit`, `junction_undefined`, `invalid_input`); null-field
  contract (every field in the null list is null); amino-acid indels on a V with a known
  inserted and a known deleted residue (strings and positions). Indel fields are not on the
  accuracy gate: truth records nucleotide indels only.
- Python: `sequence_type` validation and string classification; CLI flag; `run.json` records
  the mode.
- Nucleotide mode unchanged: existing suites pass with no floor rewrites.

## Left to the builder

Internal module layout; how ties are computed mechanically; the overlap tiebreak details beyond
"earliest cut"; exact threshold values (by the calibration above); whether the gate covers TR
sets in v1.

## Follow-ups (not v1)

- C calls: translate C in the frame continuing J (nucleotide mode's `c_frame`), giving isotype
  for Fab and full-length heavy-chain protein.
- Amino-acid indel fields in nucleotide mode (needs a rule for indels that aren't a multiple
  of 3 or that split a codon).
- Amino-acid coordinate fields, if users ask.
- Prefilter and SIMD for large protein repertoires.
