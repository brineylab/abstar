# abstar generate Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build `abstar generate`, a synthetic V(D)J generator. Every sequence comes with exact true values and fewest-mismatch valid-answer sets, for test-driven development and benchmarking of annotators.

**Architecture:**
- A new Rust crate `abstar_generate` depends only on `abstar_germline` and a new shared crate `abstar_table`.
- `abstar_table` provides schema code generation and TSV/Parquet writers. They are extracted from `abstar_annotate`, so the generator never sees annotation code.
- Models are TOML files with built-ins shipped in the package.
- Python exposes `abstar.generate()` and the Typer command `abstar generate` through `abstar._core`.

**Tech Stack:**
- Rust: `rand` + `rand_chacha`, `rayon`, `serde`/`toml`, and the existing arrow/parquet/needletail versions.
- Python ≥ 3.11 with Typer and Polars, as in the foundation.

**Spec:** `docs/superpowers/specs/2026-11-20-abstar-generate-design.md`. Read it, and the foundation spec it builds on, before any task. Section names below refer to the generate spec.

## Global Constraints

All foundation Global Constraints still hold (`docs/superpowers/plans/2026-11-20-abstar-foundation.md`). In particular:
- no subprocesses in `crates/` or `python/abstar/`;
- lower-bound-only dependencies;
- tests set `ABSTAR_HOME` to a tmp path;
- `legacy/` is read-only;
- every consequential decision gets a decision record in the same commit.

Additions for this plan:

**Crate graph**
- `crates/generate` (`abstar_generate`) may depend on `abstar_germline` and `abstar_table` only. It must never depend on `abstar_annotate`, directly or transitively.
- `crates/table` (`abstar_table`) depends on no abstar crate.

**New error codes** (both exit 2; add them to `errors.ERROR_CODES` and the manual):
- `MODEL_NOT_FOUND`: the recombination model does not exist.
- `MODEL_INVALID`: the recombination model is malformed or does not fit the germline database.

**Coordinates in the truth table**
- All `*_start`/`*_end` fields use AIRR convention: 1-based, closed, positions in `sequence`.
- `*_germline_start`/`*_germline_end` are 1-based positions in the ungapped germline allele.
- Records 0016 and 0018 still govern TSV encoding.

**Reproducibility**
- Sequence *i* is generated from a `ChaCha8Rng` seeded with `seed_for(global_seed, i)`.
- Output is byte-identical for a given seed regardless of thread count.
- When no seed is given, one is drawn from OS entropy and recorded.

**Thresholds** (each recorded in the truth decision record)
- `high_mutation`: realized V mutation ≥ 0.20.
- `long_cdr3`: junction_aa length − 2 ≥ 24 (IGH) or ≥ 16 (other loci).
- No-D threshold: best D match < 5 nt.
- Model/database coverage threshold: ≥ 0.90 of model V, D and J probability mass each.
- Productive-only attempt limit: 1000 attempts per sequence.

**Default names**
- Built-in models: `human-igh`, `human-igk`, `human-igl`, `human-tra`, `human-trb`, at `python/abstar/models/<name>/model.toml`.
- User models: `$ABSTAR_HOME/models/<name>/model.toml`.
- Default model: `f"{species}-{locus.lower()}"`, where `species` is the database's single species. Mixed-species databases require `--model`.

## Review Focus

1. **Model/database mismatch.** A model naming genes the database lacks warns and renormalizes. Coverage under 0.90 fails with `MODEL_INVALID` naming the missing genes; it never silently skews usage. Mixed-species databases without `--model` fail with `INVALID_ARGUMENT`. Tests are in Task 3.
2. **Productive-only with a low productive rate.** It must stop after the attempt limit with `INVALID_ARGUMENT` ("productive rate too low"), never loop forever. Test in Task 7.
3. **Same seed, different `--threads` or batch size.** Output must be byte-identical. Test in Task 7.
4. **Extreme parameters:**
   - mutation 0.0 and 0.5;
   - truncation longer than the read;
   - `--n 0`;
   - a junction removed by truncation (`junction` null, no crash).

   All must produce valid output. Tests are in Tasks 5 and 7.
5. **Closed stdout pipe.** `abstar generate --n 100000 --locus IGH | head -1` must exit 0 with no traceback. Rerunning into the same OUTPUT_DIR overwrites the outputs. Tests are in Task 10.

---

### Task 1: Shared table crate (`abstar_table`): schema codegen and writers

Extract the writers and the `fields.toml` code generator from `abstar_annotate` into a crate both tools can use. Behavior for `abstar run` must not change.

**Files:**
- Create: `crates/table/Cargo.toml`, `crates/table/src/{lib.rs,writers.rs,codegen.rs}`, `crates/table/tests/writers.rs`
- Modify:
  - `Cargo.toml` (workspace member)
  - `crates/annotate/{Cargo.toml,build.rs,src/schema.rs,src/writers.rs}`: `writers.rs` becomes thin wrappers bound to `arrow_schema()`, keeping today's public API
- Create: `python/abstar/decisions/0022-shared-table-crate.md`
- Test: `tests/python/test_crate_graph.py`

**Interfaces:**
- Produces:

```rust
// abstar_table::codegen, used from build.rs (abstar_table is both a dependency and a build-dependency)
pub fn generate(toml_path: &Path, out_path: &Path, record_type: &str, const_name: &str);
// Generates: struct <record_type>, fn arrow_schema() -> SchemaRef, fn to_record_batch(&[<record_type>]) -> RecordBatch,
// and pub const <const_name>: &str (include_str of the toml). Types: string | bool | int64 | float64 | list<string> | list<int64>.
// The validation rules of today's build.rs carry over unchanged (unknown type, duplicate name, abstar field without 4-digit decision).

// abstar_table::writers
pub enum Format { Airr, Parquet }                                   // FromStr "airr" | "parquet"
pub struct TsvWriter<W: Write>;  impl<W: Write> TsvWriter<W> { pub fn new(w: W, schema: &SchemaRef) -> io::Result<Self>; pub fn write(&mut self, b: &RecordBatch) -> io::Result<()>; pub fn finish(self) -> io::Result<W> }
pub struct ParquetSink;  impl ParquetSink { pub fn create(path: &Path, schema: SchemaRef) -> Result<Self, ParquetError>; pub fn write(&mut self, b: &RecordBatch) -> Result<(), ParquetError>; pub fn finish(self) -> Result<(), ParquetError> }
```

- `abstar_annotate::writers::{Format, TsvWriter::new(w), ParquetSink::create(path)}` keep their current signatures.
- TSV encodes `list<int64>` as integers joined by `,`.

- [ ] **Step 1: Write the failing tests**
  - `crates/table/tests/writers.rs`:
    - `list_int64_round_trip`: a schema with a `list<int64>` column holding `[3, 4]` and null gives TSV cell `3,4` and an empty cell. Parquet reads back the identical batch.
    - `header_follows_given_schema`: the header equals that schema's names.
  - `tests/python/test_crate_graph.py::test_generate_never_depends_on_annotate`:
    - parse `crates/generate/Cargo.toml` (skip if absent) and `crates/table/Cargo.toml` with `tomllib`;
    - assert neither has `abstar_annotate` in `dependencies` or `build-dependencies`;
    - assert `crates/table` has no `abstar_` dependency.
- [ ] **Step 2: Run and confirm failure.** Run `cargo test -p abstar_table`. Expected: the package does not exist.
- [ ] **Step 3: Implement.** Move the codegen logic out of `crates/annotate/build.rs` into `abstar_table::codegen`, and the writer bodies into `abstar_table::writers`. `crates/annotate/build.rs` becomes one call to `abstar_table::codegen::generate(...)`.
- [ ] **Step 4: Verify.** Run `cargo test --workspace` and `.venv/bin/python -m pytest tests/python -q`; both must pass unchanged. Rebuild with `maturin develop` first.
- [ ] **Step 5: Add record 0022** (why the writers moved: the generator must reuse them without seeing annotation code). Commit with `refactor: shared abstar_table crate for schema codegen and writers`.

---

### Task 2: Generate crate skeleton, model format and loading

**Files:**
- Create: `crates/generate/Cargo.toml`, `crates/generate/src/{lib.rs,error.rs,model.rs}`, `crates/generate/tests/model.rs`, `crates/generate/tests/data/{tiny_igh.model.toml,bad_*.model.toml}`
- Modify: `Cargo.toml` (workspace member)

**Interfaces:**
- Produces:

```rust
pub enum GenerateError { /* code(): "MODEL_NOT_FOUND" | "MODEL_INVALID" | "INVALID_ARGUMENT" | "IO_ERROR" | "INTERNAL" */ }
impl GenerateError { pub fn code(&self) -> &'static str }

pub struct Dist { pub values: Vec<i64>, pub probs: Vec<f64> }          // discrete; probs sum to 1 ± 1e-6
pub struct Source { pub citation: String, pub used_for: String }
pub struct Model {
    pub format_version: u32, pub name: String, pub species: String, pub locus: Locus, pub description: String,
    pub sources: Vec<Source>,
    pub v_usage: BTreeMap<String, f64>, pub d_given_j: BTreeMap<String, BTreeMap<String, f64>>, pub j_usage: BTreeMap<String, f64>,
    pub allele_usage: BTreeMap<String, f64>,                            // optional per-allele split; absent ⇒ equal within gene
    pub trim: Trims,                                                    // v_3, d_5, d_3, j_5: default Dist + per-gene BTreeMap<String, Dist>; negative values = P-nucleotides
    pub insertions: BTreeMap<String, Dist>,                             // keys "vd","dj" (loci with D) or "vj"
    pub insertions_long_tail: BTreeMap<String, Dist>,                   // same keys; mixed in with weight from GenerateParams
    pub n_markov: BTreeMap<String, [f64; 4]>,                           // keys "start","A","C","G","T" → P(A,C,G,T)
    pub d_inversion: f64,
    pub shm_model: String,                                              // "HH_S5F" | "HKL_S5F" | "none"
}
pub fn load_model(path: &Path) -> Result<Model, GenerateError>          // missing file ⇒ MODEL_NOT_FOUND; parse/validation ⇒ MODEL_INVALID naming the key
```

- Model file format (TOML, mirroring the struct): one `[[source]]` table per source; `[v_usage]`, `[j_usage]`, `[d_given_j.IGHJ4]`, `[trim.v_3.default]`, `[trim.v_3.genes."IGHV1-2"]`, `[insertions.vd]`, `[n_markov]`; each `Dist` is `values = [...]`, `probs = [...]`.
- Validation:
  - every probability table sums to 1 ± 1e-6;
  - `values` and `probs` have equal lengths;
  - `locus` parses;
  - D keys are present iff the locus has D (IGH, TRB, TRD);
  - `d_inversion` is in [0, 1].

- [ ] **Step 1: Write the failing tests** in `crates/generate/tests/model.rs`:
  - `tiny_model_loads`: the fixture loads with 3 V genes, 2 D, 2 J.
  - `missing_file_is_model_not_found`.
  - `probabilities_must_sum_to_one`: a `bad_sum` fixture gives `MODEL_INVALID`, and the message names the table key.
  - `d_tables_required_for_igh_forbidden_for_igk`.
- [ ] **Step 2: Run and confirm failure.** Run `cargo test -p abstar_generate`.
- [ ] **Step 3: Implement** `error.rs` and `model.rs` (serde, toml). The crate depends on `abstar_germline` and `abstar_table` only.
- [ ] **Step 4: Run and confirm the tests pass**, along with `tests/python/test_crate_graph.py`.
- [ ] **Step 5: Commit** with `feat(generate): crate skeleton and recombination model format`.

---

### Task 3: Bind a model to a germline database

**Files:**
- Create: `crates/generate/src/bind.rs`, `crates/generate/tests/bind.rs`, `python/abstar/decisions/0023-generate-model-format.md`

**Interfaces:**
- Consumes: `abstar_germline::{Database, Allele, Locus, Segment}`, Task 2 `Model`.
- Produces:

```rust
pub struct Bound { pub model: Model, pub v: Vec<(Allele, f64)>, pub d_given_j: BTreeMap<String, Vec<(Allele, f64)>>, pub j: Vec<(Allele, f64)>,
                   pub candidates: Candidates, pub warnings: Vec<String> }
pub struct Candidates { pub v: Vec<Allele>, pub d: Vec<Allele>, pub j: Vec<Allele>, pub c: Vec<Allele> }   // every allele of the locus (and species) in the database
pub fn bind(model: Model, db: &Database) -> Result<Bound, GenerateError>
```

**Rules:**
- **Matching.** Model genes match database alleles by `gene`, restricted to `allele.species == model.species` and `allele.locus == model.locus`. Mixed-species databases work because `name`/`gene` are plain and `species` filters.
- **Generatable alleles.** Only V and J alleles with `anchor.is_some()` can be generated. An anchorless allele's share of its gene goes to that gene's anchored alleles. A gene with none is dropped with a warning.
- **Missing genes.** A model gene with no generatable allele is dropped with a warning naming it, and the remaining mass is renormalized.
- **Coverage.** If retained mass is below 0.90 for V, D or J, the result is `MODEL_INVALID`, listing the dropped genes.
- **Candidates.** All alleles of the locus and species are candidates, anchored or not. V candidates must have a gapped sequence containing `.` (else `MODEL_INVALID`, which is why macaque fails); D and J candidates need no anchor.

- [ ] **Step 1: Write the failing tests** (a tiny database built in code with `abstar_germline::write`/`load`):
  - `drops_missing_gene_and_renormalizes`: model V {A: 0.5, B: 0.45, Z: 0.05}, database lacks Z. Result: one warning naming Z, and probabilities A ≈ 0.5263, B ≈ 0.4737.
  - `low_coverage_fails`: Z has 0.2 of the mass, so the result is `MODEL_INVALID` with "Z" in the message.
  - `anchorless_allele_not_generated_but_is_candidate`.
  - `species_filter_in_mixed_db`: hs/mm alleles share a plain name, and only the model's species is bound.
  - `ungapped_v_candidates_rejected`.
- [ ] **Step 2: Run and confirm failure.**
- [ ] **Step 3: Implement** `bind.rs`.
- [ ] **Step 4: Run and confirm the tests pass.**
- [ ] **Step 5: Add record 0023**, covering the model format, name mapping, the coverage threshold of 0.90, and the anchor rules. Commit with `feat(generate): bind models to germline databases`.

---

### Task 4: Somatic hypermutation (S5F)

**Files:**
- Create:
  - `crates/generate/src/shm.rs`
  - `crates/generate/data/s5f/{HH_S5F_mutability.tsv,HH_S5F_substitution.tsv}` and, if openly licensed, `HKL_S5F_*`
  - `crates/generate/data/s5f/SOURCES.md`
  - `crates/generate/tests/shm.rs`
  - `python/abstar/decisions/0024-s5f-mutation-model.md`

**Interfaces:**
- Produces:

```rust
pub struct S5f { /* 5-mer (ACGT only) → mutability: f64; 5-mer → [f64; 4] substitution probs (self = 0) */ }
pub fn s5f(name: &str) -> Result<&'static S5f, GenerateError>          // "HH_S5F" | "HKL_S5F"; tables embedded via include_str!
pub struct Mutation { pub pos: usize, pub from: u8, pub to: u8 }       // 0-based read position at time of mutation
pub fn mutate(seq: &mut Vec<u8>, mutable: &[bool], v_positions: &[usize], target_v_freq: f64, model: &S5f, rng: &mut impl Rng) -> Vec<Mutation>
```

**Rules:**
- **Picking a position.** Pick among `mutable` positions with probability ∝ the mutability of the 5-mer centred on that position. Pad the ends with N, and treat any 5-mer containing N or an IUPAC code as the mean mutability. Then pick the substitution from the substitution table.
- **Stopping.** Repeat until the realized V divergence (positions in `v_positions` differing from their original base) ≥ `target_v_freq × v_positions.len()` (rounded). Contexts update after each hit, and repeat hits and back-mutations are allowed.
- **Safety cap.** Stop with a warning after 10 × target hits.
- **Sources.** Obtain the tables from the openly licensed supplementary data of Yaari et al. 2013 (Front. Immunol.). For the light-chain table, use the published source only if its license allows redistribution. Otherwise use HH_S5F for light chains, and state that in record 0024 and `SOURCES.md`. Never copy from GPL/AGPL software packages.

- [ ] **Step 1: Write the failing tests**:
  - `tables_complete`: all 1024 5-mers have a mutability, and substitution rows sum to 1.
  - `realized_v_divergence_hits_target`: for targets {0.0, 0.05, 0.20, 0.35, 0.50} on a 300 nt V, the realized fraction ≥ target and < target + 2/300.
  - `hotspots_mutate_more`: over 2,000 runs at target 0.05 on a fixed sequence, the mean hits per WRC hotspot centre exceed the mean at SYC coldspot centres by at least 2×.
  - `non_mutable_positions_untouched`.
  - `deterministic_for_seed`.
- [ ] **Step 2: Run and confirm failure.**
- [ ] **Step 3: Implement** `shm.rs`, with the tables and `SOURCES.md` (URL, license, retrieval date and sha256 per file).
- [ ] **Step 4: Run and confirm the tests pass.**
- [ ] **Step 5: Add record 0024 and commit** with `feat(generate): S5F somatic hypermutation`.

---

### Task 5: Recombination, indels, truncation, constant region, productivity

**Files:**
- Create: `crates/generate/src/{recombine.rs,params.rs}`, `crates/generate/tests/recombine.rs`

**Interfaces:**
- Consumes: Task 3 `Bound`, Task 4 `mutate`.
- Produces:

```rust
pub enum Spec { Fixed(f64), Uniform(f64, f64), Beta(f64, f64) }        // FromStr: "0.1" | "fixed:0.1" | "uniform:0.15-0.35" | "beta:2,10"
pub struct GenerateParams {
    pub n: usize, pub seed: u64, pub mutation: Spec, pub long_cdr3_weight: f64, pub indel_rate: f64, pub frameshift_indels: bool,
    pub indel_duplication: bool, pub truncate_5p: Spec, pub truncate_3p: Spec /* nt */, pub constant_len: Option<usize>,
    pub isotype_weights: BTreeMap<String, f64>, pub productive_only: bool, pub shm_model: Option<String>, pub threads: Option<usize>, pub batch_size: usize /* 4096 */,
}
pub fn preset(name: &str) -> Result<GenerateParams, GenerateError>    // "naive" | "memory" | "bnab" | "tcr"; unknown ⇒ INVALID_ARGUMENT listing names
pub enum Origin { V(usize /* gapped index in true allele */), D(usize), J(usize), N1, N2, P, C(usize), Inserted }
pub struct Indel { pub kind: IndelKind /* Ins|Del */, pub read_pos: usize, pub germline_pos: usize, pub len: usize, pub bases: String, pub equiv: (usize, usize) }
pub struct Simulated {
    pub read: Vec<u8>, pub origin: Vec<Origin>,                        // origin[i] for read[i]
    pub v: Allele, pub d: Option<Allele>, pub d_inverted: bool, pub j: Allele, pub c: Option<Allele>,
    pub v_trim_3: i64, pub d_trim_5: i64, pub d_trim_3: i64, pub j_trim_5: i64,   // negative = P-nucleotides added
    pub np1: String, pub np2: String, pub mutations: Vec<Mutation>, pub indels: Vec<Indel>,
    pub truncate_5: usize, pub truncate_3: usize, pub v_mutation_frequency: f64,
}
pub fn simulate(bound: &Bound, params: &GenerateParams, rng: &mut ChaCha8Rng) -> Simulated
pub fn replay(s: &Simulated) -> Vec<u8>                                // rebuilds read from events + alleles (test oracle and validation item 3)
```

**Rules** (spec: *The recombination model*, *Somatic hypermutation*, *Indels*, *Truncation and constant region*, *Productivity*, *Presets*):
- **Order:** sample alleles → trim/P → N (Markov) → assemble (D reverse-complemented if inverted) → mutate (S5F over V(D)J, never C) → indels → append C prefix → truncate.
- **Indels:**
  - Allowed only in V or J, at least 3 nt away from the Cys104 and Trp/Phe118 codons.
  - Lengths 1–12 codons, in multiples of 3 unless `frameshift_indels`.
  - Position weight 3× inside CDR1/CDR2 (IMGT 27–38 and 56–65 aa, mapped through gapped positions) and their 3-codon flanks.
  - `equiv` is the leftmost-to-rightmost range of equivalent placements.
- **Presets:**

  | Preset | Mutation | Other settings |
  |---|---|---|
  | `naive` | fixed 0 | — |
  | `memory` | uniform 0.02–0.12 | — |
  | `bnab` | uniform 0.15–0.35 | `long_cdr3_weight` tuned (Task 11), `indel_rate` 0.10 |
  | `tcr` | fixed 0 | — |

- **Productivity** follows AIRR definitions computed on the final read:
  - `vj_in_frame`: the junction length is a multiple of 3 and the J anchor is in the V frame;
  - `stop_codon`: a stop codon in the V-frame translation;
  - `productive` = `vj_in_frame` and not `stop_codon`, with both anchor codons conserved.

- [ ] **Step 1: Write the failing tests** (`tests/recombine.rs`, using the tiny database and model):
  - `replay_rebuilds_read`: over 2,000 seeds with mutation, indels, truncation and constant enabled, `replay(s) == s.read`.
  - `p_nucleotides_are_palindromic`: `v_trim_3 = -2` yields the reverse complement of the V's last 2 nt.
  - `anchors_never_touched_by_indels`: over 2,000 seeds at indel rate 1.0, the Cys/Trp codons equal germline (barring mutation of other positions).
  - `truncation_longer_than_read`: truncation of 1,000 nt gives an empty read and no panic.
  - `spec_parsing`: each `Spec` form parses, and bad strings give `INVALID_ARGUMENT`.
  - `productive_flags_match_translation`: on hand-built junctions of 1 nt and 3 nt, and one with a stop codon.
- [ ] **Step 2: Run and confirm failure.**
- [ ] **Step 3: Implement** `params.rs` and `recombine.rs`.
- [ ] **Step 4: Run and confirm the tests pass.**
- [ ] **Step 5: Commit** with `feat(generate): recombination, indels, truncation and productivity`.

---

### Task 6: Truth sets

**Files:**
- Create: `crates/generate/src/truth.rs`, `crates/generate/tests/truth.rs` (including a brute-force oracle module), `python/abstar/decisions/0025-generate-truth-sets.md`

**Interfaces:**
- Consumes: Task 5 `Simulated`, Task 3 `Candidates`.
- Produces:

```rust
pub struct Truth {
    pub v_valid: Vec<String>, pub d_valid: Vec<String>, pub d_null_valid: bool, pub j_valid: Vec<String>, pub c_valid: Vec<String>, // allele names, sorted
    pub v_end_valid: Vec<i64>, pub d_start_valid: Vec<i64>, pub d_end_valid: Vec<i64>, pub j_start_valid: Vec<i64>,           // 1-based, sorted, deduped
    pub unrecoverable_v: bool, pub unrecoverable_d: bool, pub unrecoverable_j: bool,
}
pub fn truth(s: &Simulated, c: &Candidates) -> Truth
```

**Rules** (spec: *Truth* and *Computing the sets without alignment*):
- **V allele set.**
  - Compare the read bases with origin `V(g)`, skipping bases with `Inserted` origin and positions inside any deletion's span, against each candidate's gapped sequence at index `g`.
  - A read base where the candidate has `.` or no base counts as a mismatch (the coverage rule).
  - The set is all candidates with the minimum mismatches.
- **D and J allele sets.**
  - Take the contiguous run of read bases with origin `D(_)` (or `J(_)`).
  - Compare it against every candidate, and for D also its reverse complement, at every offset, including offsets where the candidate covers only part of the run (uncovered bases count as mismatches).
  - Each candidate's score is its best offset; the set is the candidates with the minimum score.
  - `d_null_valid` = the D run is shorter than 5 nt.
- **C allele set:** prefix comparison.
- **Unrecoverable:** the true allele is not in its set.
- **Boundaries.**
  - Score: +1 per matching germline-attributed base, −1 per mismatch, 0 for N/P.
  - For each allele in the segment's valid set, cut points range over positions where that allele's germline continuation is defined (V and J positioned by the true mapping, D by each tied best offset).
  - Maximize over ordered, non-overlapping (V end, D start, D end, J start) with DP over the separable score.
  - Collect every configuration tying the global maximum, and take each boundary's value set across all of them. Without D, skip the D terms.

- [ ] **Step 1: Write the failing tests:**
  - `oracle_agreement`: over 3,000 seeds on the tiny database (which has near-identical alleles by design) at mutation uniform 0–0.3 with indels, every `Truth` field equals the brute-force oracle. The oracle enumerates every candidate × offset × every boundary 4-tuple directly.
  - `n_matching_germline_widens_boundary`: hand-built case, V end at 100 and N1 starting with the next V germline base gives `v_end_valid == [100, 101]`.
  - `mutated_edge_base_goes_to_n`: a mutated last V base gives `v_end_valid` containing true end − 1.
  - `closer_allele_makes_unrecoverable`: mutate toward a sibling allele until strictly closer; the true allele is absent and `unrecoverable_v` is set.
  - `short_d_null_valid`.
  - `truncated_read_coverage`: a 5'-truncated read gives the same V set as computing only over the covered positions.
- [ ] **Step 2: Run and confirm failure.**
- [ ] **Step 3: Implement** `truth.rs`.
- [ ] **Step 4: Run and confirm the tests pass.**
- [ ] **Step 5: Add record 0025**, covering:
  - the fewest-mismatches rule and true values;
  - the unrecoverable label;
  - boundary scoring and joint ties stored as lists;
  - the coverage rule;
  - D/J sliding comparison;
  - the no-D threshold of 5 nt;
  - indel equivalence ranges;
  - the known limit for mis-gapped V alleles (#138).

  Commit with `feat(generate): fewest-mismatch truth sets`.

---

### Task 7: Truth schema and the generate pipeline

**Files:**
- Create: `crates/generate/truth.toml`, `crates/generate/build.rs`, `crates/generate/src/{schema.rs,pipeline.rs,meta.rs}`, `crates/generate/tests/pipeline.rs`, `python/abstar/decisions/0026-generate-outputs.md`

**Interfaces:**
- Consumes: Tasks 3–6; `abstar_table::{codegen, writers}`.
- Produces:

```rust
pub const TRUTH_TOML: &str;  pub struct TruthRecord;  pub fn arrow_schema() -> SchemaRef;  pub fn to_record_batch(&[TruthRecord]) -> RecordBatch   // generated
pub struct GenerateMeta { pub abstar_version: String, pub germline_db: String, pub db_path: PathBuf, pub model_name: String, pub model_path: PathBuf }
pub struct GenerateSummary { pub n: u64, pub productive: u64, pub attempts: u64, pub fasta: Option<PathBuf>, pub truth: BTreeMap<String, PathBuf>, pub generate_json: Option<PathBuf>, pub seed: u64, pub warnings: Vec<String> }
pub fn seed_for(global: u64, i: u64) -> u64                                       // splitmix64(global ^ splitmix64(i))
pub fn generate_to_files(bound: &Bound, params: &GenerateParams, out_dir: &Path, name: &str, formats: &[Format], meta: &GenerateMeta) -> Result<GenerateSummary, GenerateError>
pub fn generate_fasta_to_writer<W: Write>(bound: &Bound, params: &GenerateParams, out: W) -> Result<GenerateSummary, GenerateError>   // BrokenPipe ⇒ Ok
pub fn generate_in_memory(bound: &Bound, params: &GenerateParams) -> Result<Vec<RecordBatch>, GenerateError>   // ≥ 1 batch
```

**`truth.toml` fields, in this order.** Origin is `airr-required`/`airr-optional` where AIRR defines the field; otherwise `abstar`, with decision `0025` (sets, flags) or `0026` (events):
1. AIRR fields:
   - `sequence_id` (`"g{i}"`, 1-based, zero-padded to the width of `n`);
   - `sequence`, `rev_comp` (always F), `productive`, `vj_in_frame`, `stop_codon`, `locus`;
   - `v_call`, `d_call`, `j_call`, `c_call`;
   - `junction`, `junction_aa`, `junction_length`;
   - `v_sequence_start`, `v_sequence_end`, `d_sequence_start`, `d_sequence_end`, `j_sequence_start`, `j_sequence_end`, `c_sequence_start`, `c_sequence_end`;
   - `v_germline_start`, `v_germline_end`, `d_germline_start`, `d_germline_end`, `j_germline_start`, `j_germline_end`;
   - `np1`, `np1_length`, `np2`, `np2_length`, `v_identity`.
2. Gene fields: `v_gene`, `d_gene`, `j_gene`, `c_gene`.
3. Valid-answer sets:
   - `v_call_valid`, `d_call_valid`, `j_call_valid`, `c_call_valid`, `v_gene_valid`, `d_gene_valid`, `j_gene_valid` (`list<string>`);
   - `d_call_null_valid` (`bool`);
   - `v_sequence_end_valid`, `d_sequence_start_valid`, `d_sequence_end_valid`, `j_sequence_start_valid` (`list<int64>`).
4. Events:
   - `v_trim_3`, `d_trim_5`, `d_trim_3`, `j_trim_5` (`int64`, negative = P);
   - `d_inverted`;
   - `mutations` (`list<string>`, each `"{1-based read pos}:{from}>{to}"`);
   - `indels` (`list<string>`, each `"{ins|del}:{read_pos}:{len}:{bases}:{equiv_min}-{equiv_max}"`);
   - `truncate_5p`, `truncate_3p`;
   - `v_mutation_frequency` (`float64`).
5. Flags: `high_mutation`, `long_cdr3`, `has_indel`, `truncated_5p`, `truncated_3p`, `unrecoverable_v`, `unrecoverable_d`, `unrecoverable_j`, `ambiguous_v`, `ambiguous_d`, `ambiguous_j`.

**Rules:**
- **Files written:** `<name>.fasta` (header = `sequence_id`, 60-col wrap not required), `<name>.truth.tsv` and/or `.parquet`, then `generate.json`.
  - `generate.json` keys: `abstar_version`, `started`, `finished`, `seed`, `parameters` (every `GenerateParams` field plus preset name), `model` {name, sha256}, `database` {name, alleles_sha256}, `counts` {n, productive, attempts}, `warnings`.
  - Delete an existing `generate.json` before writing outputs (as in record 0018).
- **Productive-only:** attempt *j* of sequence *i* uses `seed_for(seed, i * 1000 + j)`. After 1000 failed attempts the run fails with `INVALID_ARGUMENT` ("productive rate too low: …").
- **Batching and order:** rayon `par_iter` within batches; output order = index order.

- [ ] **Step 1: Write the failing tests** (`tests/pipeline.rs`):
  - `byte_identical_across_threads_and_batches`: n=5000, seed 7. FASTA and TSV from (threads 1, batch 3) and (threads 4, batch 4096) are byte-identical.
  - `files_and_meta`: the files exist, `generate.json` parses, its seed is recorded, and TSV columns equal the `truth.toml` names.
  - `zero_n`: header-only TSV, empty FASTA, one 0-row batch.
  - `productive_only_limit`: a model with forced frameshift (insertions always length 1, no D) gives `INVALID_ARGUMENT`.
  - `fasta_broken_pipe`: a writer failing after 100 bytes gives `Ok`.
  - `junction_removed_by_truncation`: `junction` null, row present.
- [ ] **Step 2: Run and confirm failure.**
- [ ] **Step 3: Implement** `truth.toml`, `build.rs` (calling `abstar_table::codegen::generate`), `schema.rs`, `meta.rs` and `pipeline.rs`.
- [ ] **Step 4: Run and confirm the tests pass.**
- [ ] **Step 5: Add record 0026** (output files, `_valid` convention, event encodings, `sequence_id` form, per-sequence seeding). Commit with `feat(generate): truth schema and deterministic generate pipeline`.

---

### Task 8: Default human models

**Files:**
- Create:
  - `python/abstar/models/{human-igh,human-igk,human-igl,human-tra,human-trb}/model.toml`
  - `germline-sources/models/SOURCES.md`, with the per-table citations and the extraction method
  - `crates/generate/tests/builtin_models.rs`
  - `python/abstar/decisions/0027-default-recombination-models.md`

**Interfaces:**
- Consumes: Task 2 format, Task 3 `bind`, the built-in `human` database.

**Rules** (spec: *Default human models (source C)*):
- **Sources.** Draw each distribution from published, peer-reviewed summary data: gene usage, trim, insertion length and N composition per locus. Name the paper and table or figure in the model's `[[source]]` and `SOURCES.md`. Do not use IGoR/OLGA model files.
- **Fallbacks.** Where no usable published distribution exists, use a documented simple fallback and label it as such in `sources`:
  - gene usage uniform over functional (F) genes;
  - trims from a truncated normal with a published mean/SD;
  - N length geometric with a published mean.
- **SHM model:** `shm_model` = `HH_S5F` for IGH, the Task 4 light-chain choice for IGK/IGL, and `none` for TRA/TRB.
- **Long tail.** Give `insertions_long_tail` for IGH `vd`/`dj` enough mass at long lengths that Task 11 can tune `bnab` to the spec's CDRH3 targets.

- [ ] **Step 1: Write the failing test.** `each_builtin_model_binds`: for each of the five models, `load_model` + `bind` against the built-in `human` database succeed, warnings are listed, and coverage is ≥ 0.90 for each segment.
- [ ] **Step 2: Run and confirm failure.**
- [ ] **Step 3: Build the five model files and `SOURCES.md`.**
- [ ] **Step 4: Run and confirm the test passes.**
- [ ] **Step 5: Add record 0027** (source C, each locus's sources, fallbacks used, no GPL files). Commit with `feat(generate): default human recombination models`.

---

### Task 9: IGoR model importer

**Files:**
- Create: `crates/generate/src/igor.rs`, `crates/generate/tests/igor.rs`, `crates/generate/tests/data/igor_tiny/{model_params.txt,model_marginals.txt}` (a hand-written fixture in IGoR's file format, not copied from any IGoR model)

**Interfaces:**
- Produces:

```rust
pub struct ImportReport { pub kept: Vec<String>, pub marginalized: Vec<String>, pub dropped_genes: Vec<String> }
pub fn import_igor(dir: &Path, name: &str, species: &str, db: &Database) -> Result<(Model, ImportReport), GenerateError>   // MODEL_INVALID on parse errors
pub fn write_model(model: &Model, path: &Path) -> Result<(), GenerateError>      // TOML, atomic (temp sibling + rename)
```

**Rules:**
- Parse IGoR's `model_params.txt` (event definitions) and `model_marginals.txt` (event probabilities) as documented in the IGoR manual.
- Map each IGoR event onto our model: gene choice; deletions, with negative = palindromic; insertion lengths; and the dinucleotide Markov model to `n_markov`.
- Conditional distributions our format can't hold, such as V-deletion given V gene beyond per-gene tables, are marginalized and listed in `marginalized`.
- IGoR gene names are mapped to germline genes by the part before `*`.

- [ ] **Step 1: Write the failing tests:**
  - `imports_tiny_fixture`: the fixture imports, the expected probabilities appear in the `Model`, and `bind` succeeds against the tiny database.
  - `report_lists_marginalized_events`.
  - `round_trip_write_load`: `write_model` then `load_model` gives an equal `Model`.
- [ ] **Step 2: Run and confirm failure.**
- [ ] **Step 3: Implement** `igor.rs`.
- [ ] **Step 4: Run and confirm the tests pass.**
- [ ] **Step 5: Commit** with `feat(generate): IGoR model importer`.

---

### Task 10: Python bindings, `abstar.generate()`, CLI and manual

**Files:**
- Modify: `crates/py/{Cargo.toml,src/lib.rs}`, `python/abstar/{_core.pyi,__init__.py,errors.py,cli.py}`, `python/abstar/manual/__init__.py`
- Create:
  - `python/abstar/generate.py` (API plus model resolve/list/info)
  - `python/abstar/manual/topics/generate.md` (hand-written: presets, truth grading, the known limits)
- Test: `tests/python/test_generate_api.py`, `tests/python/test_generate_cli.py`

**Interfaces:**
- Consumes: Tasks 2, 3, 7, 9.
- Produces:

```python
# _core (GIL released during Rust work)
def generate_to_files(db_path: str, model_path: str, params: dict, output_dir: str, name: str, formats: list[str]) -> dict
def generate_fasta_to_stdout(db_path: str, model_path: str, params: dict) -> dict
def generate_in_memory(db_path: str, model_path: str, params: dict) -> object      # __arrow_c_stream__
def truth_toml() -> str
def model_summary(path: str) -> dict
def model_import_igor(igor_dir: str, name: str, species: str, db_path: str, out_path: str) -> dict
# python/abstar/generate.py
@dataclass(frozen=True)
class GenerateSummary: n: int; productive: int; attempts: int; seed: int; fasta: Path | None; truth: dict[str, Path]; generate_json: Path | None; warnings: list[str]
@dataclass(frozen=True)
class ModelInfo: name: str; kind: Literal["builtin", "user", "path"]; path: Path; species: str; locus: str; description: str; sources: list[dict]
def generate(n: int, locus: str, *, germline_db: str = "human", model: str | None = None, preset: str | None = None, seed: int | None = None,
             output_dir: str | Path | None = None, name: str = "generated", formats=("airr",), threads: int | None = None, **options) -> "pl.DataFrame | GenerateSummary"
def resolve_model(name_or_path, species: str, locus: str) -> tuple[Path, str]   # explicit model, else f"{species}-{locus.lower()}"
def list_models() -> list[ModelInfo]; def model_info(name_or_path) -> ModelInfo; def import_igor(igor_dir, *, name, species, germline_db="human") -> dict
```

- **`**options` keys:** `mutation`, `long_cdr3_weight`, `indels`, `frameshift_indels`, `indel_duplication`, `truncate_5p`, `truncate_3p`, `constant` (int nt or None), `isotype_weights`, `productive_only`, `shm_model` (overrides the model's, e.g. to mutate TCRs for testing). An unknown key raises `INVALID_ARGUMENT`. A user model whose name equals a built-in model name raises `BUILTIN_NAME_CONFLICT` (record 0010 applies to models).
- **Presets** are applied in Rust, and explicit options override them.

**CLI** (spec: *CLI and Python API*):
- `abstar generate [OUTPUT_DIR]` with `--n`, `--locus`, `-g/--germline-db`, `--model`, `--preset`, `--seed`, `--mutation`, `--long-cdr3`, `--indels`, `--frameshift-indels`, `--truncate-5p`, `--truncate-3p`, `--constant INT`, `--productive-only`, `--name`, `--format` (repeatable), `--threads`, `--json`, `-q`.
- `abstar model list|info NAME|import-igor DIR --name N --species S [-g DB]` (a top-level group like `abstar germline`).
- `abstar generate` is a plain command (not a group), so the positional OUTPUT_DIR never collides with a subcommand name.
- With no OUTPUT_DIR, FASTA goes to stdout and the summary to stderr, with the same flush and closed-pipe handling as `run` (`cli.main`).

**Manual:**
- New topics `generate` (the prose file plus the command tree) and `truth-fields` (rendered from `_core.truth_toml()` exactly as `fields` is).
- Mention both in `examples.md`.

- [ ] **Step 1: Write the failing tests:**
  - **API:**
    - `generate(100, "IGH", seed=1)` returns 100 rows with columns == `truth.toml` names.
    - The same seed gives an identical DataFrame.
    - `preset="bnab"` gives mean `v_mutation_frequency` within 0.15–0.35.
    - An unknown option or preset gives `INVALID_ARGUMENT`.
    - Mixed database `human+c57bl6` without `model` gives `INVALID_ARGUMENT`.
    - `model="nope"` gives `MODEL_NOT_FOUND`.
    - `output_dir` writes `generated.fasta`, `generated.truth.tsv` and `generate.json`.
    - `import_igor` on the Task 9 fixture lands in `ABSTAR_HOME/models/<name>`; importing as `human-igh` → `BUILTIN_NAME_CONFLICT`, nothing written.
  - **CLI:**
    - `abstar generate --n 3 --locus IGH --seed 1` prints a 6-line FASTA to stdout.
    - The closed-pipe test from Review Focus 5.
    - Re-running into an OUTPUT_DIR overwrites the files.
    - `abstar model list --json` includes `human-igh`.
    - `abstar manual truth-fields` lists every `truth.toml` name.
    - `ERROR_CODES` includes both new codes.
- [ ] **Step 2: Run and confirm failure.**
- [ ] **Step 3: Implement** the bindings, API, CLI and manual.
- [ ] **Step 4: Verify.** Run `cargo test --workspace`, `.venv/bin/python -m pytest tests/python -q` and `.venv/bin/sphinx-build -W -b html docs docs/_build`. `abstar --help` must stay light (the existing test).
- [ ] **Step 5: Commit** with `feat: abstar generate API, CLI and manual`.

---

### Task 11: Generator validation suite and bnab tuning

**Files:**
- Create: `crates/generate/tests/validation.rs` (marked `#[ignore]` where it needs over 10 s in debug; CI runs it with `--release -- --include-ignored`), `crates/generate/benches/generate.rs`
- Modify: `crates/generate/src/params.rs` (tuned `bnab` `long_cdr3_weight`), `.github/workflows/rewrite.yml` (add the release validation run)

**Interfaces:**
- Consumes: everything above, and the built-in `human-igh` model and `human` database.

- [ ] **Step 1: Write the validation tests** (spec: *Validation of the generator itself*):
  1. `distribution_recovery`: with n = 200,000 from `human-igh` (mutation 0), the empirical V-gene usage, J usage, each trim distribution, N lengths and N first-base composition match the bound model.
     - Statistic: chi-square with p > 0.001 per table, or total variation < 0.01 where expected counts are small.
     - The D-inversion rate is within 3 SE.
  2. `bnab_targets`: with n = 50,000, preset `bnab`, locus IGH:
     - realized `v_mutation_frequency` lies in [0.15, 0.35 + 2/300];
     - ≥ 25% of sequences have CDRH3 (junction_aa length − 2) ≥ 24 aa;
     - ≥ 2% have it ≥ 35 aa;
     - 10% ± 1% have `has_indel`.
  3. `event_replay_builtin`: `replay` equals the read for 20,000 sequences of each built-in model with all options on.
  4. `oracle_agreement_builtin`: the Task 6 oracle agrees on 2,000 sequences from `human-igh` and `human-igl`.
- [ ] **Step 2: Run them.** Run `cargo test -p abstar_generate --release --test validation -- --include-ignored`. Tune `bnab`'s `long_cdr3_weight` (and only that) until `bnab_targets` passes. If distribution recovery fails, fix the sampler, never the tolerance.
- [ ] **Step 3: Add the bench.** It times `generate_to_files` for 1M IGH sequences with truth TSV. Record the time in the commit message against the spec target (≤ 10 s on the M-series laptop target; informational on other hardware).
- [ ] **Step 4: Add the CI step and commit** with `test(generate): validation suite, bnab tuning and benchmark`.
