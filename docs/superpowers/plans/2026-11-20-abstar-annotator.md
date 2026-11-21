# abstar Annotator Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace `abstar run`'s stub annotator with real V(D)J assignment and annotation. It is built test-first against `abstar generate` answer keys and benchmarked against IgBLAST, old abstar and real data.

**Architecture:**
- **Search:** staged search finds V and J. Orientation and a k-mer prefilter come first, then SIMD affine semi-global alignment of every allele of the candidate genes.
- **Junction:** a joint junction step settles V end, D (both orientations) and J start, using matches-minus-mismatches scoring. Then C is found.
- **Derivations:** regions, alignments, junction, productivity, mutations and masks are derived from those placements.
- **Engine:** a run-level `Engine` (index built once) replaces `annotate_stub` behind the existing streaming pipeline.
- **Ground truth:** the generator grows to provide truth for every new field. A Python grader scores output against truth.

**Tech Stack:** Rust (`block-aligner` 0.5 with AVX2/NEON features as in `spikes/vsearch`, rayon, existing arrow/parquet via `abstar_table`), Python ≥ 3.11 (Polars, Typer), micromamba/bioconda for the external benchmark only.

**Spec:** `docs/superpowers/specs/2026-11-20-abstar-annotator-design.md`. Read it, `AGENTS.md`'s scientific-correctness rules, and records 0004, 0007, 0008, 0016, 0018, 0025 and 0026 before any task.

## Global Constraints

Foundation and generate constraints still hold:
- no subprocesses in `crates/` or `python/abstar/` (scripts outside may use them);
- lower-bound-only dependencies;
- tests set `ABSTAR_HOME` to a tmp path;
- `legacy/` is read-only;
- every consequential decision gets a decision record (next free id: 0028) in the same commit;
- `abstar_generate` never depends on `abstar_annotate` (enforced by `tests/python/test_crate_graph.py`). The annotator **library** never uses generator code, and truth rules are implemented independently on each side. `abstar_annotate`'s **tests** may list `abstar_generate` under `[dev-dependencies]` to generate test reads and truth in memory. The library code must not use it; `test_crate_graph.py` should assert `abstar_generate` appears in annotate's `dev-dependencies` only.

Additions:

**Coordinates**
- Output `*_start`/`*_end` are AIRR 1-based closed positions in `sequence`. When `rev_comp` is T they address the reverse complement.
- `*_germline_start`/`*_germline_end` are 1-based positions in the ungapped called allele. For an inverted D they address the reverse-complemented allele.
- Internally everything is 0-based half-open on the oriented read, converted once when the `Record` is built (record 0016).

**IMGT regions (nt positions on the V's gapped coordinates)**

| Region | IMGT nt positions |
|---|---|
| FWR1 | 1–78 |
| CDR1 | 79–114 |
| FWR2 | 115–165 |
| CDR2 | 166–195 |
| FWR3 | 196–312, including the Cys104 codon (310–312) |
| CDR3 | after the Cys104 codon to before the Trp/Phe118 codon |
| FWR4 | the Trp/Phe118 codon to the J end |

- The junction is CDR3 plus both anchor codons. It starts at the Cys104 codon, so FWR3 and the junction share that codon (AGENTS.md: FWR3 and CDR3 share the junction-start-codon endpoint).
- A codon crossing a region boundary belongs to the downstream region. Partial terminal codons don't translate (AGENTS.md).

**Mutation strings and counts**
- V nt: `IMGT{pos}:{ref}>{alt}`; D/J/C nt: `{pos}:{ref}>{alt}`. `pos` is the IMGT-gapped position (V) or the 1-based position in the called allele (D/J/C).
- aa: V uses `IMGT{aa_pos}:{ref}>{alt}` (IMGT codon number); D/J/C use `{codon_pos}:{ref}>{alt}`.
- J aa frame is set by the Trp/Phe118 codon. D aa frame is the frame D occupies in the junction (V frame).
- Insertions are `{pos}:ins:{bases}`, where `pos` is the germline position the inserted bases follow (V: `IMGT{pos}`). Deletions are `{pos}:del:{bases}`, with `pos` the first deleted germline position.
- Frequencies are mutations ÷ aligned germline positions of that segment (codons for aa), excluding deleted positions. Insertions and deletions are not substitutions and are not counted in them.

**Masks**
- `cdr_mask`: `0` FWR, `1`/`2`/`3` CDR1/2/3.
- `gene_segment_mask`: `V`/`N`/`D`/`N`/`J`. P-nucleotides count as `N`.
- `nongermline_mask`: `1` for N/P bases, mismatches and inserted bases, else `0`.
- Masks cover only the V(D)J region (user decision): nt masks index `sequence_vdj`, one character per base, never constant-region bases. aa masks index `sequence_vdj_aa`, one per codon. A codon takes the segment or region of its 3rd base for `gene_segment_mask_aa`, and the downstream-region rule for `cdr_mask_aa`. `nongermline_mask_aa` is `1` if the codon's aa differs from germline or any of its bases is N/P/inserted.

**Status and failure reasons** (`annotation_status` / `failure_reason`)

| Status | `failure_reason` |
|---|---|
| `annotated` | null |
| `unassigned` | `no_v_hit`, `no_j_hit` (V-derived fields filled), `junction_undefined` (an anchor missing or truncated; V/J calls filled), `macaque_not_imgt_gapped` (V allele lacks IMGT gaps, #137) |
| `invalid_input` | `empty_sequence` / `invalid_characters` |
| `error` | `internal_error: <msg>` |

**Tie-break for primary calls (record 0008):** score (fewest mismatches), then F > ORF > P > NF > unknown, then usage, then natural sort of names.

**Quality targets** (spec *Goals*, user-agreed as targets):
- ≥ 99.5% per field and difficulty subset on recoverable reads, on generator truth;
- match or beat IgBLAST and old abstar per subset;
- ≥ 10k sequences/s on a 10-core M-series laptop; one sequence ≤ 300 ms end to end; memory flat in input size.

Shortfalls are reported, not hidden.

**Loci held to the targets:** human IGH, IGK, IGL, TRA, TRB.
- Mouse databases annotate, with no generator check.
- Macaque returns `unassigned`/`macaque_not_imgt_gapped` for V alleles without gaps.
- No TRA/TRD cross-locus pairing: J candidates come from the V's locus.

## Review Focus

1. **Reverse-complemented input.** A read and its reverse complement annotate identically apart from `rev_comp`, with coordinates addressing the oriented sequence. Test in Task 5.
2. **N and IUPAC bases in the read.** No crash. A read base that isn't A/C/G/T counts as a mismatch for every candidate, so ties are unaffected. Test in Task 5.
3. **Non-antibody or very short input.** Random DNA, adapters, or a TCR read against an IG-only database must give `unassigned`/`no_v_hit`, never a forced call. Test in Task 5.
4. **Truncated reads.** A read starting inside CDR1/FWR2, or ending inside CDR3, must give partial regions and null (never invented) region values and junction, with the right status. Test in Task 8.
5. **A large in-frame indel near a region boundary.** For example 12 codons deleted in FWR3, or 9 inserted in CDR1. The IMGT regions, junction and mutation positions must stay correct. Test in Task 8.

---

### Task 1: Generator truth for the new fields

**Files:**
- Modify: `crates/generate/truth.toml`, `crates/generate/src/pipeline.rs` (or a new `crates/generate/src/derive.rs`), `crates/generate/tests/pipeline.rs`
- Modify: `python/abstar/decisions/0026-generate-outputs.md` (new truth fields)

**Interfaces:**
- Consumes: `Simulated` (origin map, true alleles, mutations, indels, truncation), `Truth`.
- Produces truth columns, using the **same names and encodings as the annotator output** (Task 3 field list). Each is valued for the true allele and true boundaries:
  - for each region R in fwr1, cdr1, fwr2, cdr2, fwr3, cdr3, fwr4: `R`, `R_aa`, `R_start`, `R_end`, `R_length`;
  - `sequence_aa`, `sequence_vdj`, `sequence_vdj_aa`, `sequence_vdjc`, `sequence_vdjc_aa`;
  - `cdr_mask`, `cdr_mask_aa`, `gene_segment_mask`, `gene_segment_mask_aa`, `nongermline_mask`, `nongermline_mask_aa`;
  - for each segment S in v, d, j, c: `S_mutations`, `S_mutations_aa`, `S_mutation_count`, `S_mutation_count_aa`, `S_mutation_frequency`, `S_mutation_frequency_aa`, `S_insertions`, `S_deletions`.

  All encodings follow Global Constraints. Fields that need an absent segment or anchor are null.

- [ ] **Step 1: Write the failing tests** in `crates/generate/tests/pipeline.rs`:
  - `region_truth_hand_traced`: extend the existing hand-traced IGK row.
    - FWR1–FWR3 starts and ends follow IMGT gapped positions mapped to read positions.
    - CDR3 = `junction[3..len-3]`, FWR4 starts at the Trp/Phe codon.
    - `cdr3_length == len(cdr3_aa)`.
  - `masks_consistent`: over 2,000 reads with every option on:
    - each nt mask's length equals `len(sequence_vdj)`;
    - `gene_segment_mask` matches the origin map;
    - `nongermline_mask[i] == 1` iff origin is N/P/inserted or the base differs from the true germline base;
    - each aa mask has `len(sequence_vdj_aa)` characters.
  - `segment_mutations_from_events`: over 2,000 mutated reads (with inversion enabled on a test model copy), every entry in `v/d/j_mutations` corresponds to a recorded mutation at that germline position. Their union over segments, plus N/P hits, equals the `mutations` event list. Frequencies use the aligned-germline denominator.
- [ ] **Step 2: Run and confirm they fail.** Run `cargo test -p abstar_generate --test pipeline`.
- [ ] **Step 3: Implement** the derivations in the generate crate, independently of `crates/annotate`, and add the fields to `truth.toml` (abstar fields cite decision 0026).
- [ ] **Step 4: Run** `cargo test --workspace`. The existing hash pin `neutral_output_pinned` must be updated, since the truth schema grew. Update it in this commit, with a note in the commit message.
- [ ] **Step 5: Update record 0026 and commit** with `feat(generate): truth for regions, masks, aa fields and per-segment mutations`.

---

### Task 2: Python grader

**Files:**
- Create: `python/abstar/grade.py`, `tests/python/test_grade.py`
- Modify: `python/abstar/__init__.py` (lazy export `grade`), `python/abstar/manual/topics/generate.md` (grading section)

**Interfaces:**
- Produces:

```python
def grade(truth: "pl.DataFrame", calls: "pl.DataFrame") -> "pl.DataFrame"
# One row per (field, subset). Columns: field, subset, n, correct, accuracy.
# subsets: "all", plus each truth flag that is True for a row
#   (high_mutation, long_cdr3, has_indel, truncated_5p, truncated_3p, ambiguous_v/d/j).
```

**Rules** (spec *Evaluation harness*):
- **Join and scope:** join on `sequence_id`. A row counts for a segment's fields only when that segment is recoverable (`not unrecoverable_<seg>`).
- **Calls:** a primary call is correct if it's in `<seg>_call_valid`. A null `d_call` is correct if `d_call_null_valid`.
- **Boundaries:** `v_sequence_end`, `d_sequence_start`, `d_sequence_end`, `j_sequence_start` are correct if in their `_valid` lists.
- **Exact match:** `junction`, `junction_aa`, `productive`, `vj_in_frame`, `stop_codon`.
- **Allele-dependent fields** (regions, masks, `<seg>_mutations*`, indels): graded only where the primary call equals the true allele.
- **Boundary-dependent fields:** use the shared-span rule.
  - Mutation and indel lists are compared after restricting both to germline positions inside the intersection of the truth's and the annotator's `<seg>_germline_start..end`. Counts and frequencies are recomputed from the restricted lists on both sides.
  - Masks are compared position by position in read coordinates, only where both `gene_segment_mask`s agree.
- **Import:** polars is imported inside `grade`.

- [ ] **Step 1: Write the failing tests** in `tests/python/test_grade.py`. Use hand-built tiny truth and call DataFrames:
  - `test_call_in_valid_set`: a call in `v_call_valid` → correct; a call outside → wrong; an unrecoverable row is excluded from V fields.
  - `test_null_d_when_allowed`.
  - `test_boundary_lists`.
  - `test_shared_span_mutations`: an annotator that put a mutated edge base in N matches truth on the shared span.
  - `test_subsets`: per-flag rows exist with correct `n`.
  - `test_grade_round_trip`: `abstar.grade(truth, truth_as_calls)` gives 1.0 on every row. Use real `abstar.generate(500, "IGH", seed=1)` output, renaming the truth columns that match the output schema.
- [ ] **Step 2: Run and confirm they fail.** Run `.venv/bin/python -m pytest tests/python/test_grade.py -q`.
- [ ] **Step 3: Implement** `grade.py`.
- [ ] **Step 4: Run and confirm they pass.** Run `.venv/bin/python -m pytest tests/python -q`.
- [ ] **Step 5: Commit** with `feat: abstar.grade() scores annotations against generator truth`.

---

### Task 3: Output schema v1 and the run-level engine seam

**Files:**
- Modify:
  - `crates/annotate/fields.toml`
  - `crates/annotate/src/{pipeline.rs,stub.rs,lib.rs}`
  - tests in `crates/annotate/tests/` and `tests/python/` that assert the field list
- Create: `crates/annotate/src/engine.rs`, `python/abstar/decisions/0028-annotator-field-naming.md`

**Interfaces:**
- Produces:

```rust
pub trait Annotate: Sync { fn annotate(&self, id: &str, seq: &str) -> Record; }
pub struct EngineOptions { pub threads: Option<usize> }        // extended by later tasks if needed
pub struct Engine { /* indexes and germline data built once per run */ }
impl Engine { pub fn new(db: &Database, opts: &EngineOptions) -> Result<Engine, RunError> }
impl Annotate for Engine { /* until Task 10: delegates to annotate_stub */ }
pub struct Stub<'a>(pub &'a Database);  impl Annotate for Stub<'_> { /* annotate_stub */ }
```

- `run_to_files`, `run_to_writer` and `run_in_memory` build an `Engine` once per run. The `_with` variants take `&dyn Annotate` instead of `&dyn Fn(&str, &str, &Database) -> Record`.

**`fields.toml` additions**, appended in this order after the existing fields and before `annotation_status`. Origins: AIRR fields `airr-optional`; others `abstar`, citing 0028 unless noted.
1. `locus`, `vj_in_frame`, `stop_codon`, `junction_length`.
2. For each of v, d, j, c: `<s>_sequence_start`, `<s>_sequence_end`, `<s>_germline_start`, `<s>_germline_end`, `<s>_score`, `<s>_identity`.
3. `np1`, `np1_length`, `np2`, `np2_length`, `d_inverted` (abstar).
4. `sequence_aa`, `sequence_alignment_aa`, `germline_alignment_aa`.
5. For each region R (fwr1, cdr1, fwr2, cdr2, fwr3, cdr3, fwr4): `R`, `R_aa`, `R_start`, `R_end` (AIRR) and `R_length` (abstar, aa length).
6. `sequence_vdj`, `sequence_vdj_aa`, `sequence_vdjc`, `sequence_vdjc_aa`.
7. `cdr_mask`, `cdr_mask_aa`, `gene_segment_mask`, `gene_segment_mask_aa`, `nongermline_mask`, `nongermline_mask_aa` (abstar).
8. For each of v, d, j, c: `<s>_mutations`, `<s>_mutations_aa`, `<s>_insertions`, `<s>_deletions` (`list<string>`); `<s>_mutation_count`, `<s>_mutation_count_aa` (`int64`); `<s>_mutation_frequency`, `<s>_mutation_frequency_aa` (`float64`) (abstar).

Each description states its coordinate space (fill `coordinates`) and encoding.

**Decision ids.** Records must stay contiguous (`tests/python/test_decisions.py`), so in this task **every** new abstar field cites 0028. Later tasks re-point fields when they add their own records:
- `d_inverted` → 0031;
- `R_length` → 0032;
- mutation fields → 0033;
- masks → 0034.

- [ ] **Step 1: Write the failing tests:**
  - `crates/annotate/tests/writers.rs::schema_matches_toml` still passes with the new list.
  - A new `engine_seam`: `run_in_memory` produces identical output whether called through `Engine` (still the stub) or `run_in_memory_with(&Stub(&db))`.
  - Python `test_run_api.py`: columns equal the `fields.toml` names.
- [ ] **Step 2: Run and confirm they fail.**
- [ ] **Step 3: Implement** the schema additions and the seam. Existing stub behaviour is unchanged.
- [ ] **Step 4: Verify.** Run `cargo test --workspace` and `.venv/bin/python -m pytest tests/python -q` (rebuild the extension first). `abstar manual fields` lists every new field.
- [ ] **Step 5: Add record 0028** (naming policy: AIRR names where AIRR defines them, legacy names reused, the new-field pattern). Commit with `feat(annotate): v1 output schema and run-level engine seam`.

---

### Task 4: Alignment primitives

**Files:**
- Create: `crates/annotate/src/align/{mod.rs,kmer.rs,aligner.rs}`, `crates/annotate/tests/align.rs`
- Modify: `crates/annotate/Cargo.toml` (block-aligner with target-specific SIMD features, as `spikes/vsearch/Cargo.toml` does)

**Interfaces:**
- Produces:

```rust
pub struct KmerIndex { /* k-mer → (allele index, offset) postings */ }
impl KmerIndex {
    pub fn new(seqs: &[&[u8]], k: usize) -> KmerIndex;
    pub fn hits(&self, query: &[u8]) -> Vec<(usize /*allele*/, u32 /*shared k-mers*/)>;   // sorted by allele index
}
pub enum Op { Match, Mismatch, Ins, Del }                     // Ins = base in query only; Del = base in germline only
pub struct Alignment { pub score: i32, pub q_start: usize, pub q_end: usize, pub g_start: usize, pub g_end: usize, pub ops: Vec<(Op, u32)> }
pub struct Scoring { pub matched: i8, pub mismatch: i8, pub gap_open: i8, pub gap_extend: i8 }
pub fn align_semiglobal(query: &[u8], germline: &[u8], s: &Scoring) -> Alignment   // free query ends; free germline 5' and 3' ends
pub fn cigar(a: &Alignment, query_len: usize) -> String                          // AIRR CIGAR incl. soft clips
```

- **Defaults:** k = 7, match 2, mismatch −3, gap open −10, gap extend −1. The builder may tune these within the indel requirement. Record the final values in record 0029.
- **Non-ACGT bases** in the query never score as a match.

- [ ] **Step 1: Write the failing tests:**
  - `kmer_hits_rank_true_allele`: a read sampled from allele X ranks X first.
  - `align_exact`: a substring of a germline aligns with all matches and correct coordinates.
  - `align_mismatches`.
  - `align_deletion_36nt` and `align_insertion_27nt`: a 12-codon deletion and a 9-codon insertion inside a 300-nt V are recovered as single gap ops at the right positions.
  - `align_truncated_5p`.
  - `n_never_matches`.
  - `cigar_format`: e.g. `5S290M3D10M`.
  - `deterministic`.
- [ ] **Step 2: Run and confirm they fail.** Run `cargo test -p abstar_annotate --test align`.
- [ ] **Step 3: Implement.** block-aligner, as in the spike, with traceback. A scalar fallback is acceptable on targets without SIMD.
- [ ] **Step 4: Run and confirm they pass.**
- [ ] **Step 5: Commit** with `feat(annotate): k-mer index and SIMD semi-global alignment`.

---

### Task 5: Orientation and V assignment

**Files:**
- Create: `crates/annotate/src/assign/{mod.rs,orient.rs,v.rs,ties.rs}`, `crates/annotate/tests/assign_v.rs`, `python/abstar/decisions/0029-tie-sets-and-prefilter.md`

**Interfaces:**
- Consumes: Task 4; `abstar_germline::{Database, Allele, Locus, Segment}`.
- Produces:

```rust
pub struct Oriented { pub seq: Vec<u8> /* uppercased oriented read */, pub rev_comp: bool }
pub fn orient(read: &[u8], v_index: &KmerIndex, j_index: &KmerIndex) -> Option<Oriented>   // None ⇒ no_v_hit
pub struct VHit { pub primary: usize /*allele idx*/, pub ties: Vec<usize>, pub aln: Alignment,
                  pub imgt: Vec<Option<usize>> /* read pos → IMGT gapped index of primary (None = inserted) */, pub mismatches: u32 }
pub fn assign_v(read: &Oriented, alleles: &[Allele], index: &KmerIndex, scoring: &Scoring) -> Option<VHit>
pub fn primary_by_tiebreak(tied: &[usize], alleles: &[Allele]) -> usize          // record 0008 order
```

**Rules** (spec *Pipeline* steps 2–3):
- **Orientation:** sum the k-mer hits for V plus J on each strand. The larger total wins, with ties going forward. Below a minimum evidence (no V k-mer hits on either strand) → `None`.
- **Prefilter completeness:** rank V genes by their best allele's k-mer count. Align every allele of every gene whose best count ≥ (best gene count − margin), plus the top gene in all cases. Derive the margin so that no allele outside the set can tie the best on mismatches. Bounding mismatches by shared k-mers is acceptable. State and justify the rule in record 0029, and test it.
- **Tie set:** project the best-scoring alignment onto IMGT gapped positions. Then count mismatches for each aligned candidate at those positions over the read's V-derived bases (indel positions skipped; uncovered bases count as mismatches). The tie set is the candidates with the minimum count. The primary comes from `primary_by_tiebreak`.

- [ ] **Step 1: Write the failing tests** (built-in human database, reads built in the test from germline sequences):
  - `exact_v_read_called`.
  - `sibling_tie_reported`: a read truncated so two alleles are identical over it gives both in ties, with the primary per tie-break.
  - `reverse_complement_identical`: Review Focus 1.
  - `iupac_bases_no_crash`: Review Focus 2.
  - `random_dna_no_hit` and `tcr_read_against_ig_only_db_no_hit`: Review Focus 3. Build an IG-only database in the test from the built-in human alleles.
  - `prefilter_never_drops_a_tie`: over 2,000 reads generated in-test with `abstar_generate` (dev-dependency; `generate_in_memory` with the built-in `human-igh` model, seed 1, mutation uniform 0–0.3), the V tie set contains every allele of `v_call_valid` whenever the true allele is recoverable.
  - `v_call_accuracy_naive`: on 2,000 naive IGH reads (seed 1, `naive` preset), the primary V is in `v_call_valid` for ≥ 99% of recoverable reads. Print the exact rate.
  - Put the shared helper that builds generated reads and truth in `crates/annotate/tests/common/mod.rs`. Tasks 6–9 reuse it.
- [ ] **Step 2: Run and confirm they fail.**
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run and confirm they pass.**
- [ ] **Step 5: Add record 0029 and commit** with `feat(annotate): orientation and V assignment with complete tie sets`.

---

### Task 6: J assignment

**Files:**
- Create: `crates/annotate/src/assign/j.rs`, `crates/annotate/tests/assign_j.rs`

**Interfaces:**
- Produces: `pub struct JHit { pub primary: usize, pub ties: Vec<usize>, pub aln: Alignment, pub mismatches: u32 }` and `pub fn assign_j(read: &Oriented, v: &VHit, alleles: &[Allele], scoring: &Scoring) -> Option<JHit>`.
- J candidates are the J alleles of the V's locus only.

**Rules:**
- Align every J allele of the locus to the read from (V alignment end − 10) onward.
- Tie set: fewest mismatches over the J-derived bases by sliding comparison, with uncovered bases counting as mismatches. Primary by tie-break.
- `None` → `no_j_hit`.

- [ ] **Step 1: Write the failing tests:**
  - `exact_j_called`.
  - `short_j_sibling_tie`.
  - `j_from_v_locus_only`: an IGK V with an IGH J-like tail never gets an IGHJ.
  - `missing_j_is_none`.
  - `j_call_accuracy_naive` on the Task 5 generated set: ≥ 99% in `j_call_valid` on recoverable reads, exact rate printed.
- [ ] **Step 2: Run and confirm they fail.**
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run and confirm they pass.**
- [ ] **Step 5: Commit** with `feat(annotate): J assignment`.

---

### Task 7: Joint junction step (D, both orientations) and C assignment

**Files:**
- Create: `crates/annotate/src/assign/{junction.rs,c.rs}`, `crates/annotate/tests/junction.rs`, `python/abstar/decisions/0030-joint-junction-scoring.md`, `python/abstar/decisions/0031-d-both-orientations.md`

**Interfaces:**
- Produces:

```rust
pub struct Junction { pub v_end: usize, pub d: Option<DHit>, pub j_start: usize }   // 0-based half-open on the oriented read
pub struct DHit { pub primary: usize, pub ties: Vec<usize>, pub inverted: bool, pub start: usize, pub end: usize, pub g_start: usize, pub g_end: usize }
pub fn resolve_junction(read: &Oriented, v: &VHit, j: &JHit, d_alleles: &[Allele], v_alleles: &[Allele], j_alleles: &[Allele], locus_has_d: bool) -> Junction
pub struct CHit { pub primary: usize, pub ties: Vec<usize>, pub start: usize, pub end: usize, pub g_end: usize }
pub fn assign_c(read: &Oriented, j: &JHit, junction: &Junction, c_alleles: &[Allele]) -> Option<CHit>
```

**Rules** (spec step 5):
- **Scoring:** cut points are scored by matches minus mismatches over germline-attributed bases (N/P score 0). V's germline continues past its alignment end; J's extends upstream of its alignment start.
- **D candidates:** every D allele of the locus, in forward and reverse-complement orientation, at every offset within (v_end_min, j_start_max).
- **D call:** the D tie set holds alleles (with orientation) reaching the best score for the primary V × primary J combination. `inverted` is set from the primary D. D is null when the best D run is < 5 nt.
- **Boundaries:** chosen from the primary combination's best configuration. Ties break toward the longest V, then the longest J, then the earliest D start.
- **D germline coordinates** address the reverse-complemented allele when inverted.
- **C:** prefix alignment of each C allele of the locus starting at J end. Tie set by fewest mismatches; `None` when fewer than 10 nt follow J.

- [ ] **Step 1: Write the failing tests:**
  - `n_matching_v_extends_v_end`: G5 semantics, hand-built.
  - `mutated_edge_base_goes_to_n`.
  - `short_d_is_null`: a 3-nt D run gives `d == None`.
  - `inverted_d_found`: a read built with the reverse complement of IGHD3-3*01 gives `inverted == true` and the right germline coordinates.
  - `vj_locus_has_no_d`.
  - `c_call_igg1`: a read with 30 nt of IGHG1 after J gives IGHG1, and the IGHG alleles identical over 30 nt are reported as ties.
  - `boundary_accuracy_naive` on the Task 5 generated set: the four boundaries are in their `_valid` lists for ≥ 99% of recoverable reads, exact rates printed.
- [ ] **Step 2: Run and confirm they fail.**
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run and confirm they pass.**
- [ ] **Step 5: Add records 0030 and 0031 (update `d_inverted`'s `decision` to 0031), and commit** with `feat(annotate): joint junction resolution with both D orientations; C assignment`.

---

### Task 8: Regions, alignment strings, junction, productivity, assembled sequences

**Files:**
- Create: `crates/annotate/src/derive/{mod.rs,regions.rs,alignment.rs,translate.rs,productivity.rs}`, `crates/annotate/tests/derive.rs`, `python/abstar/decisions/0032-regions-and-truncation.md`

**Interfaces:**
- Consumes: `Oriented`, `VHit`, `JHit`, `Junction`, `CHit`, alleles.
- Produces:

```rust
pub struct Regions { pub fwr1: Option<(usize, usize)>, /* … cdr1, fwr2, cdr2, fwr3, cdr3, fwr4 */ pub junction: Option<(usize, usize)> }
pub fn regions(read: &Oriented, v: &VHit, j: &JHit, junction: &Junction, v_alleles: &[Allele], j_alleles: &[Allele]) -> Regions
pub fn translate(nt: &[u8], frame: usize) -> String             // partial terminal codons dropped; codons with N → 'X'
pub struct Alignments { pub sequence: String, pub germline: String, pub sequence_aa: String, pub germline_aa: String }
pub fn alignments(read: &Oriented, v: &VHit, junction: &Junction, j: &JHit, alleles: …) -> Alignments
pub struct Productivity { pub productive: bool, pub vj_in_frame: bool, pub stop_codon: bool }
pub fn productivity(read: &Oriented, regions: &Regions, v: &VHit) -> Productivity
```

**Rules** (Global Constraints IMGT table; spec step 7):
- **Region boundaries:** map from IMGT gapped positions through `VHit.imgt`. FWR3 ends at the end of the Cys104 codon, CDR3 runs between the anchor codons, and FWR4 starts at the Trp/Phe118 codon.
- **Partial regions:** a region with no read bases is null. A partially covered region holds the covered part, with start/end of what's present. The junction is null if either anchor codon is missing.
- **Alignment strings:**
  - `sequence_alignment` and `germline_alignment` carry IMGT gaps (`.`) for V.
  - N bases align to `-`/gap in `germline_alignment`, the legacy convention: an N in the query aligns to a gap in germline.
  - Insertions in the query are gaps in germline; deletions are gaps in sequence.
- **Productivity:** same definition as the generator. `productive` = `vj_in_frame` and no V-frame stop (V start to J end) and conserved TGT/TGC and TGG/TTT/TTC.
- **Assembled sequences:** `sequence_vdj` = oriented read from V start to J end. `sequence_vdjc` extends it to C end.

- [ ] **Step 1: Write the failing tests:**
  - `regions_hand_built_igh`: exact starts, ends and sequences for a hand-built IGHV1-2*02 + IGHD3-3*01 + IGHJ4*02 read.
  - `regions_igk`.
  - `truncated_into_cdr1`: Review Focus 4. FWR1 null, CDR1 partial, the rest present.
  - `truncated_into_cdr3`: junction null, status later `junction_undefined`.
  - `large_indel_regions`: Review Focus 5. A 36-nt FWR3 deletion and a 27-nt CDR1 insertion leave the IMGT regions and junction correct.
  - `productivity_cases`: in frame, out of frame, stop codon, mutated anchor.
  - `translate_rules`.
  - `region_accuracy_naive` on the Task 5 generated set: exact match with truth regions where the V call is the true allele, ≥ 99%.
- [ ] **Step 2: Run and confirm they fail.**
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run and confirm they pass.**
- [ ] **Step 5: Add record 0032 (update the `R_length` fields' `decision` to 0032), and commit** with `feat(annotate): IMGT regions, alignments, translation and productivity`.

---

### Task 9: Mutations, indels and masks

**Files:**
- Create: `crates/annotate/src/derive/{mutations.rs,masks.rs}`, `crates/annotate/tests/mutations.rs`, `python/abstar/decisions/0033-mutation-numbering.md`, `python/abstar/decisions/0034-masks.md`

**Interfaces:**
- Produces:

```rust
pub struct SegmentMutations { pub nt: Vec<String>, pub aa: Vec<String>, pub insertions: Vec<String>, pub deletions: Vec<String>,
                              pub count: i64, pub count_aa: i64, pub frequency: f64, pub frequency_aa: f64 }
pub fn mutations(read: &Oriented, v: &VHit, junction: &Junction, j: &JHit, c: Option<&CHit>, alleles: …) -> [Option<SegmentMutations>; 4]   // V, D, J, C
pub struct Masks { pub cdr: String, pub cdr_aa: String, pub segment: String, pub segment_aa: String, pub nongermline: String, pub nongermline_aa: String }
pub fn masks(…) -> Masks
```

All encodings follow Global Constraints.

- [ ] **Step 1: Write the failing tests:**
  - `v_mutation_imgt_numbering`: a hand-built read with A>G at IMGT 153 gives `["IMGT153:A>G"]` and the matching aa string.
  - `d_and_j_mutations`: new for non-V segments. Inverted-D positions are in reverse-complement coordinates.
  - `frequency_denominator_truncated`: a 5'-truncated read's V frequency uses only aligned positions.
  - `indel_strings`.
  - `masks_hand_built`: exact strings for the Task 8 hand-built read.
  - `mutation_accuracy_memory`: on a `memory`-preset set generated in-test like Task 5's, the per-segment mutation lists match truth by the shared-span rule for ≥ 99% of reads with a correct primary call. Implement the rule in the Rust test helper as specified in Task 2.
- [ ] **Step 2: Run and confirm they fail.**
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run and confirm they pass.**
- [ ] **Step 5: Add records 0033 and 0034 (update the mutation fields' `decision` to 0033 and the mask fields' to 0034), and commit** with `feat(annotate): per-segment mutations, indels and masks`.

---

### Task 10: Assemble records; the engine replaces the stub

**Files:**
- Modify: `crates/annotate/src/{engine.rs,stub.rs,pipeline.rs}`, tests that assumed `unassigned`/`engine_not_implemented` (`crates/annotate/tests/pipeline.rs`, `tests/python/test_run_api.py`, `tests/python/test_cli.py`)
- Create: `crates/annotate/src/record.rs`, `crates/annotate/tests/engine.rs`, `python/abstar/decisions/0035-annotation-status.md`
- Modify: `python/abstar/manual/topics/conventions.md` (remove the "stub annotator" note; describe statuses)

**Interfaces:**
- Consumes: Tasks 4–9.
- Produces: `Engine::annotate` running the full pipeline. A `Record` builder converts 0-based internal coordinates to AIRR once. It fills calls, ties, `<seg>_gene`, `<seg>_species` and `<seg>_species_ties` (species from the alleles), scores and identities, regions, sequences, masks and mutations. It sets `annotation_status`/`failure_reason` per Global Constraints.

**Rules:**
- `no_j_hit` records still carry V-derived fields: V call/ties, V coordinates, FWR1–FWR3, V mutations.
- `junction_undefined` records carry calls and whatever regions exist.
- Macaque: a V primary allele whose gapped sequence lacks `.` gives `unassigned`/`macaque_not_imgt_gapped`, with V call fields filled.
- Mixed-species databases: `v_call` holds the plain name; `v_species` and `v_species_ties` are filled (record 0007).
- The `Engine` builds the k-mer indexes once per run, and is `Sync`.

- [ ] **Step 1: Write the failing tests:**
  - `engine_end_to_end_igh`: the hand-built read from Task 8 gives exact values for a chosen set of fields: calls, junction, regions, the first mutation entries, masks and status `annotated`.
  - `no_j_hit_keeps_v_fields`.
  - `macaque_reports_status`: built-in `macaque` database.
  - `mixed_species_fields`: `human+c57bl6`.
  - The existing pipeline tests are updated so valid reads expect `annotated` or a specific `unassigned` reason. Determinism across threads and batches still holds.
  - Python: `abstar.run("<a real IGHV1-2 + IGHJ4 read>")` gives `annotation_status == "annotated"` and `v_call` starting with `IGHV1-2`. `abstar run --seq …` likewise.
- [ ] **Step 2: Run and confirm they fail.**
- [ ] **Step 3: Implement.** Keep `annotate_stub` as `Stub` for the seam tests only.
- [ ] **Step 4: Verify.** Run `cargo test --workspace`, rebuild, `.venv/bin/python -m pytest tests/python -q`, and the sphinx build.
- [ ] **Step 5: Add record 0035 and commit** with `feat(annotate): real annotation replaces the stub`.

---

### Task 11: Accuracy gate, speed, and the AIRR check

**Files:**
- Create: `tests/python/test_accuracy.py`, `tests/python/data/accuracy_floors.toml`
- Modify: `crates/annotate/benches/pipeline.rs`, `scripts/bench.py`, `.github/workflows/rewrite.yml`

**Rules:**
- **Python accuracy test.** `test_accuracy.py` generates fixed-seed sets with `abstar.generate` for each held locus: IGH naive, IGH memory, IGH bnab (productive-only), IGK, IGL, TRA, TRB (2,000 reads each). It annotates them with `abstar.run` and grades with `abstar.grade`. It asserts each (field, subset) accuracy ≥ a committed floor table (`tests/python/data/accuracy_floors.toml`).
  - **Floors:** set to the measured value minus 0.5 percentage points, capped at the 99.5% target.
  - **Target table:** it prints a full table of measured vs the 99.5% target, so shortfalls are visible.
- **Speed.** The criterion bench `pipeline_10k` now runs the real engine on 10k corpus reads. `scripts/bench.py` reports `throughput_rps` for real annotation. Record both in the commit message against the targets: ≥ 10k sequences/s (M-series target; informational here), single sequence ≤ 300 ms.
- **AIRR check.** `airr.validate_rearrangement` must pass on a 2,000-read annotated TSV.

- [ ] **Step 1: Write the tests and the floors file**, initially with a placeholder floor of 0.0. Run the suite to measure.
- [ ] **Step 2: Set the floors** (measured − 0.5 pp, ≤ 99.5) and commit them together with the measured table in the commit message.
- [ ] **Step 3: Add the CI step**, in release, on ubuntu-latest/py3.13, like the generate validation step.
- [ ] **Step 4: Run the benches.** Run `cargo bench -p abstar_annotate --bench pipeline` and `.venv/bin/python scripts/bench.py`, and record the numbers.
- [ ] **Step 5: Commit** with `test(annotate): accuracy gate, AIRR validation and benchmarks`.

---

### Task 12: External benchmark and report

**Files:**
- Create:
  - `scripts/benchmark/setup_env.sh` (micromamba env `abstar-bench` with `igblast`, `mmseqs2`, and the legacy abstar deps from `legacy/requirements.txt`, installing legacy abstar from `legacy/` into that env only)
  - `scripts/benchmark/make_igblast_db.py` (builds IgBLAST V/D/J databases and an auxiliary file from abstar's built-in `human` germlines with `makeblastdb`, so names match)
  - `scripts/benchmark/run_benchmark.py` (generates sets, runs IgBLAST, old abstar and new abstar, converts each to AIRR columns, grades with `abstar.grade`, and compares new-abstar gene calls with `legacy/test_data/bcr_corpus/baseline.parquet`)
  - `docs/superpowers/benchmarks/2026-11-20-annotator-v1.md` (the dated report)

**Rules:**
- **Sets:** fixed seeds; 5,000 reads per held locus and subset: IGH naive, memory, bnab, indels-heavy (`--indels 0.5`), truncated (`--truncate-5p uniform:0-60 --truncate-3p uniform:0-40`), plus IGK, IGL, TRA, TRB.
- **The report has:**
  - per-field and per-subset accuracy for all three tools vs the 99.5% target;
  - each shortfall classified as genuine ambiguity the truth sets don't capture vs annotator error, with examples;
  - real-data gene-call agreement with legacy, every disagreement class explained with counts and examples;
  - speed;
  - mouse: "annotated; no generator validation; no mouse real-data corpus in the repo".
- **Placement:** scripts live outside the package and may use subprocesses. Nothing here is a CI dependency.

- [ ] **Step 1: Write and run** `setup_env.sh` and `make_igblast_db.py`. Record the tool versions in the report.
- [ ] **Step 2: Write and run** `run_benchmark.py`.
- [ ] **Step 3: Write the report.**
- [ ] **Step 4: Commit** with `docs: annotator v1 benchmark against IgBLAST, legacy abstar and real data`. Do not push. The controller hands the report to the user, whose acceptance defines "done" (spec).
