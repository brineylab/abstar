# abstar Foundation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the shared foundation of the rewritten abstar: a Rust workspace built into one Python wheel, the germline database format and `abstar germline` commands, the single-source output schema with TSV/Parquet writers, a streaming `run` pipeline (with a stub annotator), the Python API, the Typer CLI, `abstar manual`, decision records, and a benchmark harness. A throwaway alignment speed test comes first.

**Architecture:** Rust owns every hot path: FASTA parsing, multithreaded batches, writers, and germline building/validation. Python (Typer CLI and API) resolves names and paths, calls a few coarse-grained Rust functions, and turns Arrow data into Polars. Output fields are defined once in `crates/annotate/fields.toml`. Rust code is generated from it, and the manual reads it.

**Tech Stack:**
- Rust (stable), PyO3 (abi3-py311), maturin, pyo3-arrow, arrow and parquet, needletail, rayon, serde/toml/serde_json, thiserror, sha2
- Python ≥ 3.11, Typer, Polars
- Tests: pytest, `airr` (official AIRR validator), criterion

**Spec:** `docs/superpowers/specs/2026-11-20-abstar-foundation-design.md`. Read it before any task. Section names below refer to it.

## Prerequisites

- **Rust toolchain:** this machine has no Rust toolchain installed. Install rustup with the stable toolchain before Task 1. That needs the user's approval, because it installs software.
- **Python environment:** an isolated Python ≥ 3.11 environment with `maturin` and `uv` available. Invoke tools through that interpreter.
- **Development install:** `maturin develop --release` (or `pip install -e .`) after any Rust change, before running pytest.

## Global Constraints

**Versions and dependencies**
- Python `requires-python = ">=3.11"`. PyO3 feature `abi3-py311`, so one wheel per platform covers every supported Python.
- Python dependency versions are **lower bounds only**:
  - runtime: `typer>=0.12`, `polars>=1.44`
  - test: `pytest>=8`, `airr>=2.0`
- Package version `1.0.0-dev.0` in the Cargo workspace; maturin derives the Python version `1.0.0.dev0`. No releases or pre-releases from this plan.

**Process and I/O**
- No external binaries and no subprocesses anywhere in `crates/` or `python/abstar/`.
- No temp files on the run path. The germline build uses only a staging directory next to its destination.
- Output row order equals input order, whatever the thread count or batch size.
- Every input record produces exactly one output row.
- `annotation_status` ∈ {`annotated`, `unassigned`, `invalid_input`, `error`}.
- Foundation `failure_reason` values:
  - `engine_not_implemented` for `unassigned`
  - `empty_sequence` / `invalid_characters` for `invalid_input`
  - `internal_error: <message>` for `error`
- A record's `sequence_id` is the FASTA/FASTQ header's first whitespace-delimited token, byte-for-byte. It is always a string.
- Valid sequence characters, case-insensitive: `ACGTN` plus IUPAC `RYSWKMBDHV`.
- `sequence` is written exactly as read; it is not uppercased.

**Exit codes and error codes**
- Exit codes:

  | Exit code | Meaning |
  |---|---|
  | `0` | Success |
  | `1` | Fatal run error |
  | `2` | Invalid usage or input |
  | `3` | Strict-mode record failure |

- Error codes, each with its exit code:

  | Error code | Exit code |
  |---|---|
  | `INPUT_NOT_FOUND` | 2 |
  | `INPUT_UNREADABLE` | 2 |
  | `INVALID_ARGUMENT` | 2 |
  | `OUTPUT_NAME_CONFLICT` | 2 |
  | `DATABASE_NOT_FOUND` | 2 |
  | `DATABASE_INVALID` | 2 |
  | `RECIPE_INVALID` | 2 |
  | `DATABASE_EXISTS` | 2 |
  | `BUILTIN_NAME_CONFLICT` | 2 |
  | `IO_ERROR` | 1 |
  | `INTERNAL` | 1 |
  | `STRICT_RECORD_FAILURE` | 3 |

**CLI behavior**
- Results go to stdout only when stdout is the result destination. Logs and progress always go to stderr.
- Progress is shown only when stderr is a terminal.
- No interactive prompts, ever.

**Germline databases**
- Germline anchor positions are **0-based offsets into the ungapped allele sequence**: the first base of the V Cys104 codon, or of the J Trp118/Phe118 codon.
- User databases live at `$ABSTAR_HOME/germline/<name>/`, where `ABSTAR_HOME` defaults to `~/.abstar`. Tests must set `ABSTAR_HOME` to a pytest `tmp_path` and never touch the real home directory.

**AIRR**
- Target schema: AIRR 2.0, as reported by the `airr` package.
- TSV encoding: `T`/`F` booleans, empty cells for nulls, lists joined with `,`.
- Parquet: native types and native list columns.

**Repository rules**
- `legacy/` is a read-only reference. Never edit it, import from it, or depend on it, except to read data files named in this plan.
- Every consequential decision a task makes gets a decision record (Task 3 format) in the same commit.
- The future `generate` crate must never depend on `annotate`. This plan does not create it.

## Review Focus

1. **Messy FASTA input.** Real files arrive with:
   - multi-line records
   - CRLF line endings
   - lowercase bases
   - blank lines
   - gzip compression
   - a header with an empty sequence
   - duplicate IDs
   - IDs such as `10E8` and `0012`

   Every one must yield one row per record with the ID preserved exactly. The tests are in Task 7 and Task 10.
2. **Output name collisions.** Two inputs whose output stems collide (`a/s1.fasta` and `b/s1.fasta`, or `s1.fa` and `s1.fastq.gz`) must fail up front with `OUTPUT_NAME_CONFLICT` before writing anything, not silently overwrite each other. The test is in Task 7.
3. **Broken pipe on stdout.** `abstar run --seq … | head -1` and stdout closed early must exit 0 without a traceback, and stdout must contain only TSV. The test is in Task 11.
4. **Interrupted or conflicting germline builds.** A failed build must leave no partial database directory. Building a user database named `human` must fail with `BUILTIN_NAME_CONFLICT`. Tests must never write under the package directory or the real `~/.abstar`. The tests are in Task 5 and Task 8.
5. **A missing path that looks like a sequence.** `abstar.run("sample1.fasta")` must raise `INPUT_NOT_FOUND`, not annotate the string `sample1.fasta` as an invalid sequence. Zero inputs must return a full-schema, zero-row DataFrame. The tests are in Task 10.

---

### Task 1: Throwaway alignment speed test

This task confirms the throughput targets before anything is built on them (spec: *Milestones and order*, item 1). It is throwaway: no TDD, and it never joins the workspace.

**Files:**
- Create: `spikes/vsearch/` (a standalone Cargo project, outside any workspace)
- Create: `docs/superpowers/spikes/2026-11-20-alignment-speed.md` (the report)

**Interfaces:**
- Consumes:
  - `abstar/germline_dbs/bcr/human/ungapped/v.fasta` (342 human IGHV alleles)
  - `test_data/bcr_corpus/sequences.fasta.gz` (50,000 real BCR sequences, median 576 nt)
  - `test_data/bcr_corpus/baseline.parquet` (legacy `v_call` per record, for a sanity check)
- Produces: the report only.

- [ ] **Step 1: Build a V-assignment prototype in `spikes/vsearch/`.**
  - A k-mer prefilter over the 342 V alleles keeps the top few candidates per query (choose k and the number of candidates by experiment).
  - Semi-global affine-gap alignment of the query against each candidate, using SIMD. Any crate is allowed here, e.g. `block-aligner`, or hand-written NEON/AVX2. The point is to measure what is achievable.
  - Try both orientations; the best score wins.
  - Use rayon across records.
- [ ] **Step 2: Measure throughput.** Replicate the corpus ×20 (1M records). Measure V-step records/sec single-threaded and at all cores. Run in release mode with `RUSTFLAGS="-C target-cpu=native"`.
- [ ] **Step 3: Sanity-check accuracy.** Report the gene-level agreement rate with the legacy `v_call` (allele suffix stripped, ties counted as agreement). This is informational, not a gate.
- [ ] **Step 4: Write the report.** Include:
  - the machine, core count, and chip
  - k and candidate count
  - single-thread and all-core rates
  - scaling efficiency
  - the agreement rate

  **Gate:** the V step must reach ≥ 3,000 records/sec/core single-threaded, with ≥ 70% of linear scaling at all cores. The V step is the dominant cost, so this leaves room for J, D and annotation within the spec target of ≥ 10k records/sec on a 10-core M-series laptop.

  State pass or fail plainly. **On a fail, stop and bring the report to the user:** the spec says the design is revisited.
- [ ] **Step 5: Commit**

```bash
git add spikes/vsearch docs/superpowers/spikes/2026-11-20-alignment-speed.md
git commit -m "spike: measure V search and alignment throughput"
```

---

### Task 2: Move legacy code aside; workspace, wheel and CI skeleton

**Files:**
- Move (git mv, tracked files only) into `legacy/`:
  - `abstar/`, `scripts/`, `test_data/`
  - `pyproject.toml`, `requirements.txt`, `requirements-test.txt`, `requirements-corpus.txt`
  - `MANIFEST.in`, `.coveragerc`, `coverage-floors.json`
  - `Dockerfile`, `.travis.yml`, `.readthedocs.yaml`
  - `docs/source`, `docs/Makefile`, `docs/make.bat`, `docs/doc_requirements.txt`, `docs/EDGE-CASES.md`

  `docs/superpowers/` stays where it is.
- Delete: `.github/workflows/{corpus,pytest,python-publish,scheduled-linkcheck}.yml`. They keep running on `master`.
- Create:
  - `Cargo.toml` (workspace: members `crates/germline`, `crates/annotate`, `crates/py`; `exclude = ["spikes"]`; `workspace.package.version = "1.0.0-dev.0"`)
  - `crates/germline/Cargo.toml` + `src/lib.rs` (crate `abstar_germline`)
  - `crates/annotate/Cargo.toml` + `src/lib.rs` (crate `abstar_annotate`, which depends on `abstar_germline`)
  - `crates/py/Cargo.toml` + `src/lib.rs` (cdylib, PyO3 module `abstar._core`)
  - `pyproject.toml`:
    - maturin backend
    - `python-source = "python"`, `module-name = "abstar._core"`, `manifest-path = "crates/py/Cargo.toml"`
    - `requires-python = ">=3.11"`, dynamic version, MIT
    - runtime deps `typer>=0.12`, `polars>=1.44`
    - `[project.optional-dependencies] test = ["pytest>=8", "airr>=2.0"]`
    - script `abstar = "abstar.cli:main"`
    - pytest `testpaths = ["tests/python"]`
  - `python/abstar/__init__.py`, `python/abstar/_core.pyi`
  - `tests/python/test_smoke.py`
  - `.github/workflows/rewrite.yml`
- Modify: `AGENTS.md`. Prepend a "Rewrite branch" section:
  - new code lives in `crates/` and `python/`; `legacy/` is reference-only
  - the build/test commands below
  - point to the spec, the plan, and `python/abstar/decisions/`

  Keep the existing correctness rules below that section unchanged.
- Modify: `.gitignore`. Add `target/`.

**Interfaces:**
- Produces:
  - `abstar._core.version() -> str`: the crate version
  - `abstar.__version__: str`, from `importlib.metadata.version("abstar")`
  - `abstar/__init__.py` imports nothing heavy (no polars, no typer, no `_core`) at import time

- [ ] **Step 1: Do the moves and deletions; commit them alone**

```bash
git commit -m "chore: move legacy implementation under legacy/ for the rewrite"
```

- [ ] **Step 2: Write the failing smoke test** in `tests/python/test_smoke.py`

```python
import subprocess, sys
import abstar

def test_version_matches_core():
    from abstar import _core
    assert abstar.__version__ == "1.0.0.dev0"
    assert _core.version() == "1.0.0-dev.0"

def test_import_is_light():
    code = "import sys, abstar; print(any(m in sys.modules for m in ('polars', 'typer', 'abstar._core')))"
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True)
    assert out.stdout.strip() == "False"
```

- [ ] **Step 3: Run it and confirm it fails.** Run `python -m pytest tests/python/test_smoke.py -q`. Expected: an import or collection failure.
- [ ] **Step 4: Create the workspace, crates, `pyproject.toml` and package skeleton.** Then run `maturin develop --release`.
- [ ] **Step 5: Run the tests and confirm they pass.** Run `cargo test --workspace && python -m pytest tests/python -q`. Expected: all pass.
- [ ] **Step 6: Write `.github/workflows/rewrite.yml`**
  - Trigger: pushes and PRs on `rewrite`.
  - **test job.** Matrix: `ubuntu-latest`, `ubuntu-24.04-arm`, `macos-14`. Python 3.11 and 3.13. Steps:
    - Rust stable
    - `pip install maturin`
    - `maturin develop --release --extras test`
    - `cargo test --workspace`
    - `python -m pytest -q`
  - **wheels job.** `PyO3/maturin-action` for linux x86_64, linux aarch64, macOS arm64, and windows x86_64 (Windows is build only), plus an sdist. Then install the sdist in a fresh venv and run `python -c "import abstar._core"`. This proves the sdist builds.
- [ ] **Step 7: Commit**

```bash
git add -A && git commit -m "feat: Rust workspace, maturin wheel, and CI skeleton"
```

---

### Task 3: Decision records

**Files:**
- Create: `python/abstar/decisions/NNNN-<slug>.md`, one per decision below
- Create: `python/abstar/decisions/README.md`, the format description
- Test: `tests/python/test_decisions.py`

**Interfaces:**
- Produces:
  - Record files named `^\d{4}-[a-z0-9-]+\.md$`, numbered contiguously from `0001`.
  - Each file has a first line `# NNNN: <Title>` and the sections `## Status`, `## Context`, `## Decision`, `## Alternatives considered`, `## Consequences`.
  - `## Status` holds one of `Accepted (YYYY-MM-DD)` or `Superseded by NNNN (YYYY-MM-DD)`.
  - Later tasks add records here.

Each record's content comes from the spec section named. Write the rationale in your own words, and keep each record to one screen.

| ID | Title | Spec section |
|---|---|---|
| 0001 | Rust engine behind a Python shell; Typer CLI | Approaches considered |
| 0002 | No external binaries (MMseqs, Parasail dropped) | Intent, Goals |
| 0003 | `run()` returns a DataFrame or a RunSummary, chosen by the caller | Python API |
| 0004 | One output row per input record; status columns replace failure logs | Record status and diagnostics |
| 0005 | Germline database format: `manifest.toml` + `alleles.tsv`; a database is a set of loci | Germline databases |
| 0006 | Germline sources: OGRDB first, best curated source otherwise | Constraints |
| 0007 | Mixed-species IDs: internal prefix, plain names in output, per-segment species | Germline databases → Mixed species |
| 0008 | Single primary gene call plus `*_call_ties`; tie-break order | Gene calls and ties |
| 0009 | Donor keep-lists restrict only the locus/segment pairs they name | Build modes |
| 0010 | User databases may not shadow built-in names | Locations |
| 0011 | Agents as first-class users (manual, stdout/stderr, `--json`, codes) | CLI, `abstar manual` |
| 0012 | Output fields single-sourced in `fields.toml` | Output schema |
| 0013 | Multiple inputs: per-input output files, flat layout, `source_file` column | Multiple inputs |
| 0014 | The generator may not depend on annotator code | Repository and build layout |
| 0015 | Python ≥ 3.11 | Constraints |
| 0016 | Plan-level conventions | This plan's Global Constraints |

Record 0016 covers these plan-level conventions:
- `sequence_id` is the first header token;
- the valid nucleotide alphabet;
- `sequence` is not uppercased;
- the exit-code table;
- germline anchors are 0-based ungapped offsets;
- functionality values `F`/`ORF`/`P`/`NF`/`unknown`.

- [ ] **Step 1: Write the failing test** in `tests/python/test_decisions.py`

```python
def test_records_are_well_formed():
    # every file in python/abstar/decisions except README.md:
    #   - its name matches ^\d{4}-[a-z0-9-]+\.md$
    #   - numbers run 0001..N with no gaps
    #   - its first line is f"# {num}: <title>"
    #   - it contains all five required headings in order
    #   - its Status line matches Accepted (YYYY-MM-DD) or Superseded by NNNN (YYYY-MM-DD)
    ...

def test_seed_records_present():
    # ids 0001..0016 exist
    ...
```

Fill in the bodies; the assertions are fixed by the comments.
- [ ] **Step 2: Run the test and confirm it fails.** Run `python -m pytest tests/python/test_decisions.py -q`. Expected: FAIL (no files).
- [ ] **Step 3: Write the 16 records and the README.**
- [ ] **Step 4: Run the test and confirm it passes.** Run `python -m pytest tests/python/test_decisions.py -q`.
- [ ] **Step 5: Commit**

```bash
git add python/abstar/decisions tests/python/test_decisions.py
git commit -m "docs: seed decision records from the foundation design"
```

---

### Task 4: Germline data model, file format, and validation (`abstar_germline`)

**Files:**
- Create: `crates/germline/src/{lib.rs,model.rs,io.rs,validate.rs,error.rs}`
- Test: `crates/germline/tests/format.rs` and fixtures in `crates/germline/tests/data/`

**Interfaces:**
- Produces (all `pub`, re-exported from the crate root):

```rust
pub enum Locus { IGH, IGK, IGL, TRA, TRB, TRG, TRD }          // Display/FromStr as those strings
pub enum Segment { V, D, J, C }                                // Display/FromStr "V","D","J","C"
pub enum Functionality { F, ORF, P, NF, Unknown }              // strings "F","ORF","P","NF","unknown"
pub struct Allele {
    pub id: String, pub name: String, pub gene: String,
    pub locus: Locus, pub segment: Segment, pub species: String,
    pub functionality: Functionality, pub sequence_gapped: String,
    pub anchor: Option<usize>, pub source: String, pub usage: Option<f64>,
}
impl Allele { pub fn sequence(&self) -> String }              // gapped with '.' removed
pub struct Manifest {
    pub format_version: u32, pub name: String, pub description: String,
    pub abstar_version: String, pub built: String /* YYYY-MM-DD */,
    pub alleles_sha256: String, pub species: Vec<String>, pub loci: Vec<Locus>,
    pub prefixes: BTreeMap<String, String>, pub provenance: Vec<Provenance>,
    pub warnings: Vec<BuildWarning>,
}
pub struct Provenance { pub species: String, pub files: Vec<String>, pub sha256: Vec<String>, pub description: String }
pub struct BuildWarning { pub allele: String, pub message: String }
pub struct Database { pub manifest: Manifest, pub alleles: Vec<Allele> }
pub fn load(dir: &Path) -> Result<Database, GermlineError>
pub fn write(db: &Database, dir: &Path, overwrite: bool) -> Result<(), GermlineError>
pub fn validate(alleles: &[Allele]) -> Result<Vec<BuildWarning>, GermlineError>
pub enum GermlineError { /* variants map 1:1 onto codes */ }
impl GermlineError { pub fn code(&self) -> &'static str }     // "DATABASE_NOT_FOUND", "DATABASE_INVALID", "DATABASE_EXISTS", "RECIPE_INVALID", "IO_ERROR"
```

- **`alleles.tsv` format:**
  - Header: `id name gene locus segment species functionality sequence_gapped anchor source usage`, tab-separated.
  - Empty cells mean `None`.
  - Rows are sorted by locus (enum order), then segment (enum order), then `id` (byte order).
- **`manifest.toml`** serializes `Manifest` with serde. `prefixes` is omitted when empty; `provenance` and `warnings` are arrays of tables.
- **`write`:**
  - Computes `alleles_sha256` over the exact `alleles.tsv` bytes and sets it in the manifest.
  - Writes both files into a sibling staging directory `<dir>.staging-<pid>`, then renames it into place.
  - If `dir` exists and `overwrite` is false → `DATABASE_EXISTS`. If `overwrite` is true, it swaps via rename and removes the old directory.
  - On any error the staging directory is removed.
- **`load`:**
  - Missing directory → `DATABASE_NOT_FOUND`.
  - Hash mismatch, bad header, unparseable cell, or anchor out of range → `DATABASE_INVALID`, with the row number in the message.
  - It runs `validate` and returns `DATABASE_INVALID` on its errors.
- **`validate`:**
  - Errors:
    - duplicate `id`
    - empty sequence
    - a `sequence_gapped` character outside `ACGTN.` (case-insensitive)
    - `anchor >= sequence().len()`
  - Warnings (returned, not errors):
    - a V without an anchor
    - a J without an anchor

- [ ] **Step 1: Write the failing tests** in `crates/germline/tests/format.rs`
  - `round_trip`: write a 4-allele database (one each of V/D/J/C) into a temp dir, `load` it, and assert equality. Assert the file has the header line exactly as above and that its rows come out sorted.
  - `hash_mismatch_is_invalid`: edit one byte of `alleles.tsv` after writing, then assert `load(..).unwrap_err().code() == "DATABASE_INVALID"`.
  - `duplicate_id_rejected` and `bad_character_rejected`: assert `validate` errors.
  - `missing_v_anchor_warns`: a V allele with `anchor: None` yields exactly one warning naming it.
  - `existing_dir_without_overwrite`: the second `write` returns `DATABASE_EXISTS`; with `overwrite=true` it succeeds and leaves no `*.staging-*` sibling.
  - `failed_write_leaves_nothing`: make the destination's parent read-only, then assert an error and no staging or final directory.
- [ ] **Step 2: Run the tests and confirm they fail.** Run `cargo test -p abstar_germline`.
- [ ] **Step 3: Implement** `model.rs`, `io.rs` (serde + toml + sha2), `validate.rs` and `error.rs` (thiserror).
- [ ] **Step 4: Run the tests and confirm they pass.** Run `cargo test -p abstar_germline`.
- [ ] **Step 5: Commit** with message `feat(germline): database model, format, and validation`.

---

### Task 5: Germline build from recipes (`abstar_germline::build`)

This covers build modes 1, 3 and 4 for gapped inputs (spec: *Build modes*). Gapping ungapped V genes waits for the annotator's aligner. If a V FASTA has no `.` in a sequence and that V has no anchor, the build reports it in a warning.

**Files:**
- Create: `crates/germline/src/{recipe.rs,sources.rs,anchors.rs,build.rs}`
- Test:
  - `crates/germline/tests/build.rs`
  - fixtures in `crates/germline/tests/data/`: a trimmed OGRDB JSON holding IGHV4-28*05, IGHJ6*03 and IGHJ4*02, copied verbatim from `legacy/abstar/germline_dbs/bcr/human/raw/Homo_sapiens_IGH_VDJ_rev_9_ex.json`; small gapped FASTAs; and recipes

**Interfaces:**
- Consumes: Task 4 types, `write`, `validate`, `load`.
- Produces:

```rust
pub fn find_v_anchor(gapped: &str) -> Option<usize>
pub fn find_j_anchor(ungapped: &str) -> Option<usize>
pub struct BuildOptions { pub overwrite: bool, pub abstar_version: String }
pub struct BuildReport {
    pub name: String, pub path: PathBuf, pub n_alleles: usize,
    pub counts: BTreeMap<String, usize>,          // key "LOCUS/SEGMENT/species", e.g. "IGH/V/human"
    pub restricted: Vec<String>, pub passed_through: Vec<String>,   // "LOCUS/SEGMENT" keys, for keep-lists
    pub warnings: Vec<BuildWarning>,
}
pub fn build(recipe: &Path, output_dir: &Path, database_paths: &BTreeMap<String, PathBuf>, opts: &BuildOptions) -> Result<BuildReport, GermlineError>
```

**Recipe format (TOML; paths are relative to the recipe file):**

```toml
name = "human+c57bl6"
description = "Human and C57BL/6 mouse IG"
keep = "genotype.txt"              # optional donor keep-list

[prefixes]                          # required iff the final database has >1 species; one per species
human = "hs"
mouse = "mm"

[[source]]
species = "human"                   # required unless `database` is set
description = "free text"           # optional; OGRDB sources auto-append germline_set_name, release_version, release_date
ogrdb_json = ["Homo_sapiens_IGH_VDJ_rev_9_ex.json"]
v_fasta = []                        # IMGT-gapped V
d_fasta = []
j_fasta = []
c_fasta = []

[[source]]
database = "human"                  # existing database: looked up in `database_paths` (filled by Python, Task 8); else used as a directory path; else RECIPE_INVALID
```

**Build rules:**

- **OGRDB JSON:**
  - Read `GermlineSet[*].allele_descriptions[*]`.
  - `name` = `label`; segment = `sequence_type`; locus = `locus`.
  - Functionality = `F` if `functional` is true, otherwise `NF`.
  - V: `sequence_gapped` = `v_gene_delineations[0].aligned_sequence` (the IMGT scheme).
  - D/J: `sequence_gapped` = `coding_sequence`.
  - J anchor = `j_cdr3_end - 1` when present. Otherwise use the motif.
- **FASTA:**
  - Segment comes from the key (`v_fasta` → V, and so on). It is never parsed from the name, because `IGHD*01` (constant) collides with D-gene names.
  - If the header contains `|`, then `name` = field 2 and functionality = field 4 with `()[]` stripped (`F`/`ORF`/`P`, else `unknown`). Otherwise `name` = the first token and functionality = `unknown`.
  - Locus = `name[..3]`, which must parse as a `Locus`; otherwise `RECIPE_INVALID`, naming the allele.
- **`gene`:** the part of `name` before `*`, or the whole name if there is no `*`.
- **`id`:**
  - `name` if the final database has one species.
  - Otherwise `"{prefix}|{name}"`. A missing prefix for any species → `RECIPE_INVALID`.
  - An allele from a `database` source that was already prefixed is re-prefixed from its `name`.
- **`source`:** the source file's basename, or `database:<name>`.
- **Anchors:**
  - V: if `gapped.len() >= 312` and `gapped[309..312]` (case-insensitive) ∈ {`TGT`, `TGC`}, the anchor is the number of non-`.` characters in `gapped[..309]`. Otherwise `None` plus a warning: `no Cys104 codon at IMGT positions 310-312`.
  - J: the first index `i` in ascending order (any frame) where `ungapped[i..i+12]` matches `(TGG|TTT|TTC)GG[ACGT]{4}GG[ACGT]`. If OGRDB supplies `j_cdr3_end` and the motif index differs from `j_cdr3_end - 1`, keep the OGRDB value and warn `motif anchor X disagrees with OGRDB j_cdr3_end Y`.
  - D and C: `None`, and no warning.
- **`keep`:**
  - One allele per line; `#` starts a comment; blank lines are ignored.
  - Each entry matches an `id`, else a unique `name`. An unknown or ambiguous entry → `RECIPE_INVALID`, listing every offender.
  - Let S be the set of (locus, segment) pairs among the matched alleles. Drop every allele whose (locus, segment) is in S and is not kept. Alleles outside S pass through.
  - Report S as `restricted` and the remaining pairs present as `passed_through`.
- **Provenance:** one `Provenance` per source, with file basenames, their sha256, and the description.
- **Output:** validate, then `write(db, output_dir, overwrite)`. Validation warnings and anchor warnings both go into `manifest.warnings` and the report.

- [ ] **Step 1: Write the failing tests** in `crates/germline/tests/build.rs`. These use exact values verified against the legacy data:

```rust
#[test] fn v_anchor_ighv4_28_05() { assert_eq!(find_v_anchor(IGHV4_28_05_GAPPED), Some(285)); }
#[test] fn j_anchor_ighj6_03() { assert_eq!(find_j_anchor("ATTACTACTACTACTACTACATGGACGTCTGGGGCAAAGGGACCACGGTCACCGTCTCCTCAG"), Some(29)); }
#[test] fn j_anchor_ighj4_02() { assert_eq!(find_j_anchor("ACTACTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG"), Some(14)); }
#[test] fn j_anchor_none_without_motif() { assert_eq!(find_j_anchor("ACGTACGTACGTACGTACGT"), None); }
```

Then add these build tests:
  - `ogrdb_build`: the fixture JSON plus one C FASTA builds. Assert:
    - 4 alleles
    - `IGHJ6*03` has anchor 29 and functionality `F`
    - no warnings
    - gene `IGHV4-28` for `IGHV4-28*05`
  - `mixed_species_prefixes`: two sources (human, mouse), each with an allele named `IGHV1-2*01`. Assert ids `hs|IGHV1-2*01` and `mm|IGHV1-2*01`, and both names are plain `IGHV1-2*01`.
  - `mixed_without_prefix_fails`: same setup with no `[prefixes]` → `RECIPE_INVALID`.
  - `keep_list_restricts_only_named_pairs`: a database with IGHV×3, IGHJ×2, IGKV×1 and a keep-list naming one IGHV. Assert:
    - the result has 1 IGHV, 2 IGHJ, 1 IGKV
    - `restricted == ["IGH/V"]`
    - `passed_through == ["IGH/J", "IGK/V"]`
  - `keep_list_unknown_entry_fails`: the error message names the unknown entry.
  - `fasta_header_imgt_style`: `>X62106|IGHV1-2*01|Homo sapiens|(F)|V-REGION|...` gives name `IGHV1-2*01`, functionality `F`.
  - `unknown_locus_fails`: an allele named `XYZV1*01` → `RECIPE_INVALID` naming it.
  - `failed_build_leaves_no_directory`: a recipe with a duplicate id → error, and `output_dir` does not exist afterwards.
- [ ] **Step 2: Run the tests and confirm they fail.** Run `cargo test -p abstar_germline`.
- [ ] **Step 3: Implement** `anchors.rs`, `recipe.rs` (serde), `sources.rs` (OGRDB via serde_json; FASTA via needletail) and `build.rs`.
- [ ] **Step 4: Run the tests and confirm they pass.** Run `cargo test -p abstar_germline`.
- [ ] **Step 5: Add record `0017-germline-recipes.md`.** It covers the recipe format, segment-from-key, the OGRDB J anchor rule, and the anchor-disagreement policy. Then run `python -m pytest tests/python/test_decisions.py -q`.
- [ ] **Step 6: Commit** with message `feat(germline): recipe-driven database builds`.

---

### Task 6: Output schema from `fields.toml`, plus TSV and Parquet writers (`abstar_annotate`)

**Files:**
- Create:
  - `crates/annotate/fields.toml`. This is the spec's `schema/fields.toml`, placed inside the owning crate so sdists build.
  - `crates/annotate/build.rs`
  - `crates/annotate/src/{schema.rs,writers.rs}`
- Test: `crates/annotate/tests/writers.rs`

**Interfaces:**
- **`fields.toml` format:**

```toml
airr_version = "2.0"

[[field]]
name = "sequence_id"
type = "string"            # string | bool | int64 | float64 | list<string>
nullable = false
origin = "airr-required"   # airr-required | airr-optional | abstar
coordinates = ""           # coordinate space, when the field holds positions
description = "..."
decision = ""              # required (4-digit id) when origin = "abstar"
```

- **Initial fields, in this order:**
  1. The 14 AIRR-required fields: `sequence_id sequence rev_comp productive v_call d_call j_call sequence_alignment germline_alignment junction junction_aa v_cigar d_cigar j_cigar`. Types are per AIRR: `rev_comp` and `productive` are bool; the rest are string.
  2. `c_call` (airr-optional, string).
  3. For each segment s in v, d, j, c: `{s}_gene`, `{s}_species` (string, abstar, decision 0008/0007), `{s}_call_ties`, `{s}_species_ties` (list<string>, abstar, 0008).
  4. `annotation_status` (string, non-null, abstar, 0004), `failure_reason` (string, abstar, 0004), `source_file` (string, abstar, 0013).

  Non-null fields: `sequence_id`, `sequence`, `annotation_status`.
- **`build.rs`** generates, into `OUT_DIR`, code that `schema.rs` includes:
  - `pub struct Record`, with one field per entry: non-null string → `String`; nullable string → `Option<String>`; `bool`/`int64`/`float64` → `Option<bool|i64|f64>`; `list<string>` → `Option<Vec<String>>`. It derives `Default`, `Clone`, `Debug` and `PartialEq`.
  - `pub fn arrow_schema() -> arrow::datatypes::SchemaRef`. The schema metadata holds `airr_version` and `abstar_version`.
  - `pub fn to_record_batch(records: &[Record]) -> arrow::record_batch::RecordBatch`
  - `pub const FIELDS_TOML: &str` (via `include_str!`)

  `build.rs` fails the build if an `abstar` field lacks `decision`, a name repeats, or a type is unknown.
- **`writers.rs`:**

```rust
pub enum Format { Airr, Parquet }                      // FromStr "airr" | "parquet"
pub struct TsvWriter<W: Write>;  impl TsvWriter<W> { pub fn new(w: W) -> io::Result<Self>; pub fn write(&mut self, b: &RecordBatch) -> io::Result<()>; pub fn finish(self) -> io::Result<W> }
pub struct ParquetSink;          impl ParquetSink { pub fn create(path: &Path) -> Result<Self, ParquetError>; pub fn write(&mut self, b: &RecordBatch) -> Result<(), ParquetError>; pub fn finish(self) -> Result<(), ParquetError> }   // Task 7 maps errors to IO_ERROR
```

- **TSV rules:**
  - The header is written on `new`, even when no batches follow.
  - `T`/`F` for bools; empty cell for null; lists joined with `,`.
  - A string value containing a tab or newline is a programming error: panic with the field name.
- **Parquet:** zstd compression. Key/value metadata copies the schema metadata.

- [ ] **Step 1: Write the failing tests** in `crates/annotate/tests/writers.rs`
  - `schema_matches_toml`: the `arrow_schema()` field names equal the toml names in order. Assert `annotation_status` is non-nullable and `v_call_ties` has type `List(Utf8)`.
  - `tsv_encoding`: for a record with `rev_comp=Some(true)`, `productive=None`, `v_call_ties=Some(vec!["IGHV1-2*02","IGHV1-2*04"])`, the TSV row contains the cells `T`, empty, and `IGHV1-2*02,IGHV1-2*04` in the right columns.
  - `tsv_header_only_for_zero_rows`: the output equals the header line plus `\n`.
  - `parquet_round_trip`: reading back with the parquet crate yields the identical `RecordBatch` and metadata `airr_version == "2.0"`.
  - `id_strings_preserved`: `sequence_id` values `10E8` and `0012` come back byte-identical in TSV and Parquet.
- [ ] **Step 2: Run the tests and confirm they fail.** Run `cargo test -p abstar_annotate`.
- [ ] **Step 3: Implement** `fields.toml`, `build.rs` (parse with `toml` as a build-dependency), `schema.rs` and `writers.rs`.
- [ ] **Step 4: Run the tests and confirm they pass.** Run `cargo test -p abstar_annotate`.
- [ ] **Step 5: Commit** with message `feat(annotate): single-source output schema and writers`.

---

### Task 7: Streaming run pipeline with stub annotator (`abstar_annotate`)

**Files:**
- Create: `crates/annotate/src/{input.rs,pipeline.rs,stub.rs,error.rs,runmeta.rs}`
- Test:
  - `crates/annotate/tests/pipeline.rs`
  - fixtures in `crates/annotate/tests/data/`: multi-line, CRLF, lowercase, blank-line, gzip, empty-sequence and duplicate-ID FASTA; one FASTQ

**Interfaces:**
- Consumes: Task 6 `Record`, `to_record_batch`, writers; Task 4 `Database`.
- Produces:

```rust
pub enum Input { File(PathBuf), Stdin, Records(Vec<(String, String)>) }
pub struct RunOptions { pub threads: Option<usize>, pub strict: bool, pub loci: Option<Vec<Locus>>, pub batch_size: usize /* default 4096 */ }
pub struct StatusCounts { pub annotated: u64, pub unassigned: u64, pub invalid_input: u64, pub error: u64 }
pub struct InputSummary { pub source: String, pub outputs: BTreeMap<String, PathBuf>, pub counts: StatusCounts }
pub struct RunSummary { pub inputs: Vec<InputSummary>, pub counts: StatusCounts, pub run_json: Option<PathBuf> }
pub struct RunMeta { pub abstar_version: String, pub germline_db: String, pub db_path: PathBuf }
pub enum RunError { /* codes: INPUT_NOT_FOUND, INPUT_UNREADABLE, INVALID_ARGUMENT, OUTPUT_NAME_CONFLICT, IO_ERROR, INTERNAL, STRICT_RECORD_FAILURE */ }
impl RunError { pub fn code(&self) -> &'static str }

pub fn annotate_stub(id: &str, seq: &str, db: &Database) -> Record      // stub.rs; the annotator project replaces this
pub fn run_to_files(inputs: &[Input], output_dir: &Path, db: &Database, formats: &[Format], opts: &RunOptions, meta: &RunMeta) -> Result<RunSummary, RunError>
pub fn run_to_writer<W: Write>(inputs: &[Input], db: &Database, opts: &RunOptions, out: W) -> Result<StatusCounts, RunError>   // TSV
pub fn run_in_memory(inputs: &[Input], db: &Database, opts: &RunOptions) -> Result<Vec<RecordBatch>, RunError>  // ≥1 batch, possibly 0 rows
```

`pipeline.rs` has an internal variant that takes `annotator: &(dyn Fn(&str, &str, &Database) -> Record + Sync)`, used by tests.

**Rules:**

- **Stub annotator (`annotate_stub`):**
  - Empty sequence → `invalid_input` / `empty_sequence`.
  - A character outside the alphabet → `invalid_input` / `invalid_characters`.
  - Otherwise → `unassigned` / `engine_not_implemented`.
  - `sequence_id` and `sequence` are copied verbatim. Every other field is `None`.
- **Errors per record:** each record's annotator call runs under `catch_unwind`. A panic → status `error`, `failure_reason` = `internal_error: <panic message>`. If `strict` is set, the first `error` aborts with `STRICT_RECORD_FAILURE`.
- **Batching:**
  - Read `batch_size` records.
  - Map them with a rayon `par_iter` inside a pool of `threads` (None = all cores), collecting in order.
  - Fill `source_file` (the input path as given, `"-"` for stdin, `None` for `Records`).
  - Convert, write, repeat.
  - Mark this with `// ponytail: sequential batches; overlap read/annotate/write if benches show I/O stalls`.
- **Inputs:** needletail reads FASTA and FASTQ, gzip-transparently. `sequence_id` = the header bytes up to the first ASCII whitespace.
  - A missing file → `INPUT_NOT_FOUND`.
  - A parse failure → `INPUT_UNREADABLE`, naming the file and record number.
  - An empty file → zero records, which is not an error.
- **`run_to_files`:**
  - Output stem = the file name minus a trailing `.gz`, then minus `.fasta|.fa|.fna|.fastq|.fq`. `Stdin` → `stdin`; `Records` → `records`.
  - Any stem collision → `OUTPUT_NAME_CONFLICT`, listing the colliding sources, before any file is written.
  - Writes `<output_dir>/<stem>.tsv` and/or `<stem>.parquet` (creating `output_dir`) and overwrites existing ones.
  - Then writes `<output_dir>/run.json` with the keys:
    - `abstar_version`, `airr_version`
    - `started`, `finished` (RFC 3339)
    - `parameters` {`germline_db`, `db_path`, `formats`, `threads`, `strict`, `loci`}
    - `database` {`name`, `alleles_sha256`}
    - `inputs` [{`source`, `outputs`, `counts`}]
    - `counts`
- **`loci`:** every requested locus must be in `db.manifest.loci`, else `INVALID_ARGUMENT`. Otherwise it is recorded only.

- [ ] **Step 1: Write the failing tests** in `crates/annotate/tests/pipeline.rs`
  - `messy_fasta_one_row_per_record`. The fixture has 7 records covering multi-line, CRLF, lowercase, a blank line between records, an empty sequence, a duplicate ID, and IDs `10E8`/`0012`. Assert:
    - 7 rows
    - IDs in file order and exact
    - the lowercase `sequence` is unchanged
    - status `invalid_input`/`empty_sequence` on the empty one
    - `unassigned`/`engine_not_implemented` on the others
  - `gzip_and_fastq_inputs`: the gzip copy of the same fixture produces identical TSV bytes. The FASTQ fixture yields its records.
  - `invalid_characters`: `ACGTX` → `invalid_input`/`invalid_characters`. `ACGTRYN` → `unassigned`.
  - `deterministic_across_threads_and_batches`: 10,000 records (the fixture replicated). TSV from (threads=1, batch=7) is byte-identical to TSV from (threads=4, batch=4096).
  - `output_stem_conflict`: inputs `a/s1.fasta` and `b/s1.fasta` → `OUTPUT_NAME_CONFLICT`, and `output_dir` contains no files.
  - `run_json_counts`: `run.json` parses, and `counts.unassigned` plus `counts.invalid_input` equals the input record count.
  - `panic_becomes_error_row` / `strict_aborts`: inject an annotator that panics on ID `boom`. Non-strict: that row has status `error` and the other rows are present. Strict: `STRICT_RECORD_FAILURE`.
  - `zero_records`: an empty file gives a header-only TSV and `run_in_memory` gives one batch with 0 rows and the full schema.
  - `unknown_locus_filter`: `loci = [TRG]` against a database without TRG → `INVALID_ARGUMENT`.
- [ ] **Step 2: Run the tests and confirm they fail.** Run `cargo test -p abstar_annotate`.
- [ ] **Step 3: Implement** `input.rs`, `stub.rs`, `pipeline.rs`, `runmeta.rs` and `error.rs`.
- [ ] **Step 4: Run the tests and confirm they pass.** Run `cargo test -p abstar_annotate`.
- [ ] **Step 5: Add record `0018-run-outputs.md`.** It covers output stems, the conflict policy, overwrite-on-rerun, the `run.json` keys, and panics becoming `error` rows.
- [ ] **Step 6: Commit** with message `feat(annotate): streaming run pipeline with stub annotator`.

---

### Task 8: Python bindings, errors, and `abstar.germline`

**Files:**
- Modify: `crates/py/src/lib.rs`, `python/abstar/_core.pyi`
- Create: `python/abstar/errors.py`, `python/abstar/germline.py`
- Test: `tests/python/test_germline_api.py`, `tests/python/conftest.py`

**Interfaces:**
- **Consumes:** Tasks 4–7.
- **Produces (`_core`; every function releases the GIL during Rust work):**

```python
class CoreError(Exception): ...            # args == (code: str, message: str)
def version() -> str: ...
def fields_toml() -> str: ...
def germline_build(recipe: str, output_dir: str, overwrite: bool, database_paths: dict[str, str]) -> dict: ...   # BuildReport fields; path as str
def germline_summary(path: str) -> dict: ...  # manifest fields + "counts": {"IGH/V/human": n, ...} + "n_alleles"
def run_to_files(inputs: list[str], output_dir: str, db_path: str, formats: list[str], threads: int | None, strict: bool, loci: list[str] | None) -> dict: ...
def run_to_stdout(inputs: list[str], db_path: str, threads: int | None, strict: bool, loci: list[str] | None) -> dict: ...  # counts; BrokenPipe -> returns normally
def run_in_memory(inputs: list[str], records: list[tuple[str, str]], db_path: str, threads: int | None, strict: bool, loci: list[str] | None) -> object: ...  # has __arrow_c_stream__ (pyo3-arrow)
```

  Input strings `"-"` mean stdin.
- **`python/abstar/errors.py`:**
  - `class AbstarError(Exception)` with `.code: str`, `.message: str` and `.exit_code: int`.
  - `ERROR_CODES: dict[str, tuple[int, str]]`, mapping code → (exit code, one-line description), with exactly the codes and exit codes in Global Constraints.
  - `def from_core(e: CoreError) -> AbstarError`.
  - `EXIT_CODES: dict[int, str]`.
- **`python/abstar/germline.py`:**

```python
@dataclass(frozen=True)
class DatabaseInfo: name: str; kind: Literal["builtin", "user", "path"]; path: Path; species: list[str]; loci: list[str]; n_alleles: int; counts: dict[str, int]; description: str
@dataclass(frozen=True)
class BuildReport: name: str; path: Path; n_alleles: int; counts: dict[str, int]; restricted: list[str]; passed_through: list[str]; warnings: list[dict]
def abstar_home() -> Path                      # $ABSTAR_HOME or ~/.abstar
def builtin_names() -> list[str]               # dirs under abstar/germline_dbs containing manifest.toml
def resolve(name_or_path: str | Path) -> tuple[Path, Literal["builtin", "user", "path"]]
def build(recipe: str | Path, *, output_dir: str | Path | None = None, overwrite: bool = False) -> BuildReport
def list() -> builtins.list[DatabaseInfo]      # builtins first, then user; each sorted by name
def info(name_or_path: str | Path) -> DatabaseInfo
```

**Rules:**

- `resolve` checks, in order:
  1. an existing directory containing `manifest.toml`
  2. a built-in name
  3. a user name

  Otherwise it raises `DATABASE_NOT_FOUND`, with the available names in the message.
- `build` with `output_dir=None` targets `abstar_home()/germline/<recipe name>`. If that name is a built-in, it raises `BUILTIN_NAME_CONFLICT` before any work.
- `build` reads the recipe with `tomllib`, resolves every `[[source]] database` value with `resolve`, and passes the resulting `{value: path}` map as `database_paths`.
- `germline.py` must not import polars or typer.

- [ ] **Step 1: Write the failing tests.** In `tests/python/conftest.py`, add an autouse fixture that sets `ABSTAR_HOME` to `tmp_path / "home"`. In `tests/python/test_germline_api.py`, using the Task 5 fixture recipe copied into `tests/python/data/`:
  - `test_build_to_user_home`: the database lands at `tmp_path/home/germline/<name>`, and `info(name).kind == "user"`.
  - `test_builtin_name_conflict`: a recipe named `human` → `AbstarError` with code `BUILTIN_NAME_CONFLICT`. Assert that nothing is created.
  - `test_exists_then_overwrite`: the second build raises `DATABASE_EXISTS`; `overwrite=True` succeeds.
  - `test_resolve_order_and_not_found`: an unknown name raises `DATABASE_NOT_FOUND`, and the message lists the user database name.
  - `test_error_table`: `ERROR_CODES` equals the Global Constraints table exactly.
  - `test_database_source_by_name`: a recipe with `database = "<user db name>"` plus a keep-list builds, and `restricted` is as expected.
- [ ] **Step 2: Run the tests and confirm they fail.** Run `maturin develop --release && python -m pytest tests/python/test_germline_api.py -q`.
- [ ] **Step 3: Implement** the bindings (pyo3, pyo3-arrow), `errors.py` and `germline.py`. Map Rust error `.code()` → `CoreError(code, message)`.
- [ ] **Step 4: Run the tests and confirm they pass.** Run `cargo test --workspace && python -m pytest tests/python -q`.
- [ ] **Step 5: Commit** with message `feat: Python bindings, error codes, and germline API`.

---

### Task 9: Built-in germline databases

**Files:**
- Create:
  - `germline-sources/recipes/{human,c57bl6,balbc,macaque,human+c57bl6}.toml`
  - `germline-sources/<db>/…` (source files copied from legacy, as below)
  - `scripts/build_builtin_germlines.py`
  - `python/abstar/germline_dbs/<db>/{manifest.toml,alleles.tsv}` (generated, committed)
- Test: `tests/python/test_builtin_databases.py`

**Interfaces:**
- **Consumes:** `abstar.germline.build(recipe, output_dir=..., overwrite=True)`.
- **Produces:** the built-in names `human`, `c57bl6`, `balbc`, `macaque`, `human+c57bl6`.
- **Sources** (copied from `legacy/abstar/germline_dbs/`):

  | Database | Sources |
  |---|---|
  | `human` | OGRDB JSON ×3 from `bcr/human/raw/` (species `human`); C from `bcr/human/imgt_gapped/c.fasta`; TR from `tcr/human/imgt_gapped/{v,d,j,c}.fasta`, copied with the `__homo_sapiens` header suffix stripped (species `human`) |
  | `c57bl6` | OGRDB JSON ×5 + `gapped_constants.fasta` from `bcr/c57bl6/raw/` (species `mouse`) |
  | `balbc` | OGRDB JSON ×5 + `gapped_constants.fasta` from `bcr/balbc/raw/` (species `mouse`) |
  | `macaque` | `bcr/macaque/imgt_gapped/{v,d,j,c}.fasta` (species `macaque`) |
  | `human+c57bl6` | the human IG sources + the c57bl6 sources; prefixes `human = "hs"`, `mouse = "mm"` |

- Each recipe's `description` cites the true per-locus source. For `balbc` and `c57bl6`, cite the *Mus musculus* OGRDB sets (the spec fixes this provenance error).
- **`scripts/build_builtin_germlines.py [NAME ...]`** rebuilds the named databases (default all) into `python/abstar/germline_dbs/` with `overwrite=True`, and prints each report's warnings.

- [ ] **Step 1: Write the failing tests** in `tests/python/test_builtin_databases.py`
  - `test_builtin_names`: `germline.builtin_names() == ["balbc", "c57bl6", "human", "human+c57bl6", "macaque"]`.
  - `test_segment_counts`: `info(db).counts` summed per segment equals:

    | Database | V | D | J | C |
    |---|---|---|---|---|
    | `human` | 588 | 37 | 109 | 124 |
    | `c57bl6` | 196 | 8 | 22 | 35 |
    | `balbc` | 265 | 10 | 22 | 35 |
    | `macaque` | 1071 | 52 | 28 | 46 |
    | `human+c57bl6` | 538 | 39 | 45 | 138 |

  - `test_sequences_match_legacy`. For each database and segment, the set of (name, ungapped sequence) equals the legacy `ungapped/<seg>.fasta` set with any `__<species>` suffix stripped. Legacy `human` is the union of `bcr/human` and `tcr/human`. For `human+c57bl6`, compare (species, name, sequence), with legacy species taken from the suffix. **If this fails, report the difference; do not adjust the expectation.**
  - `test_mixed_ids_prefixed_names_plain`. In `human+c57bl6`:
    - every id starts with `hs|` or `mm|`
    - no `name` contains `|`
    - `hs|IGHV1-2*02` exists with name `IGHV1-2*02`
  - `test_human_tr_j_anchor_warnings`: the `human` manifest warnings include the alleles `TRAJ16*01`, `TRAJ16*02`, `TRAJ35*01`, `TRBJ2-7*02`, `TRGJP1*01` and `TRGJP2*01`, each with a no-motif message.
  - `test_human_v_anchor`: `IGHV3-23*01` has anchor 285.
  - `test_builtins_load`: every built-in database loads cleanly through `germline_summary`.
- [ ] **Step 2: Run the tests and confirm they fail.** Run `python -m pytest tests/python/test_builtin_databases.py -q`.
- [ ] **Step 3: Copy the sources, write the recipes and the script, then run `python scripts/build_builtin_germlines.py`.** Read every printed warning. Any warning on an allele that legacy marked functional, or any OGRDB anchor disagreement, goes into the commit message.
- [ ] **Step 4: Run the tests and confirm they pass.** Run `python -m pytest tests/python/test_builtin_databases.py -q`.
- [ ] **Step 5: Add record `0019-builtin-databases.md`.** It lists each database's sources and notes that TR and macaque functionality is `unknown`, because their source FASTAs carry no functionality field.
- [ ] **Step 6: Commit** with message `feat: built-in germline databases from recipes`.

---

### Task 10: Python `abstar.run()`

**Files:**
- Create: `python/abstar/api.py`
- Modify: `python/abstar/__init__.py`. Export `run`, `RunSummary`, `InputSummary`, `AbstarError`, and lazy `germline` (module `__getattr__`) so that `import abstar` stays light.
- Test: `tests/python/test_run_api.py`

**Interfaces:**
- **Consumes:** `_core.run_to_files`, `_core.run_in_memory`, `germline.resolve`, `errors.from_core`.
- **Produces:**

```python
@dataclass(frozen=True)
class InputSummary: source: str; outputs: dict[str, Path]; counts: dict[str, int]
@dataclass(frozen=True)
class RunSummary: inputs: list[InputSummary]; counts: dict[str, int]; run_json: Path
def run(sequences, output_dir: str | Path | None = None, *, germline_db: str = "human", locus: str | list[str] | None = None,
        formats: str | Iterable[str] = ("airr",), threads: int | None = None, strict: bool = False) -> "pl.DataFrame | RunSummary"
```

`counts` keys are exactly `annotated`, `unassigned`, `invalid_input` and `error`.

**Input normalization:**

- **Paths:** a `str` or `Path` that exists is a path; a directory expands recursively to files ending `.fasta .fa .fna .fastq .fq`, optionally with `.gz`, sorted by relative path. A string that does not exist but contains `/` or `\`, or ends with one of those extensions, raises `INPUT_NOT_FOUND`.
- **Sequences:** any other `str` is one sequence with ID `"1"`. An iterable of strings that are not paths becomes sequences with IDs `"1".."N"`. An iterable of `(id, seq)` tuples becomes records. Mixing paths and sequences raises `INVALID_ARGUMENT`.
- **Iterables:** an iterable is consumed exactly once.
- **Validation:**
  - `formats` ⊆ {`airr`, `parquet`}, and not empty.
  - `threads` is None or ≥ 1.
  - `output_dir=None` with `formats` other than the default is allowed but ignored. Document this in the docstring.
- **Return values:**
  - `output_dir=None` → `pl.DataFrame(_core.run_in_memory(...))`, importing polars inside the function.
  - Otherwise → `RunSummary`.

- [ ] **Step 1: Write the failing tests** in `tests/python/test_run_api.py`. Use the built-in `human` database.
  - `test_single_sequence_dataframe`: `run("CAGGTGCAG")` returns a `pl.DataFrame` with one row, `sequence_id == "1"`, `annotation_status == "unassigned"`, and columns equal to the `fields.toml` names in order.
  - `test_zero_inputs`: `run([])` returns 0 rows with the full schema.
  - `test_missing_path_not_sequence`: `run("sample1.fasta")` and `run("data/x")` both raise `AbstarError`, code `INPUT_NOT_FOUND`.
  - `test_tuples_and_generator`: `run((p for p in [("a","ACGT"),("b","")]))` gives IDs `["a","b"]` and statuses `["unassigned","invalid_input"]`.
  - `test_directory_input_and_source_file`. A dir with `x.fasta` and `sub/y.fa.gz`:
    - in-memory rows come in `sub/y.fa.gz`-then-`x.fasta` order (sorted relative paths)
    - `source_file` holds each path
  - `test_files_summary`. `run(dir, out, formats=["airr","parquet"])` returns a `RunSummary` with:
    - `x.tsv`, `x.parquet`, `y.tsv`, `y.parquet`
    - `run.json`
    - `counts` summed equal to the record total
  - `test_airr_validation`: `airr.validate_rearrangement` on `x.tsv` returns True. If it rejects nulls in required fields, **stop and report**; do not alter the data to pass.
  - `test_ids_preserved_parquet`: IDs `10E8` and `0012` read back as Utf8, exactly.
  - `test_mixed_inputs_rejected`: `run([existing_path, "ACGT"])` raises `INVALID_ARGUMENT`.
- [ ] **Step 2: Run the tests and confirm they fail.** Run `python -m pytest tests/python/test_run_api.py -q`.
- [ ] **Step 3: Implement** `api.py` and the `__init__.py` exports.
- [ ] **Step 4: Run the tests and confirm they pass.** Run `python -m pytest tests/python -q`.
- [ ] **Step 5: Add record `0020-run-input-normalization.md`** covering the path-vs-sequence rule and default IDs. Then commit with message `feat: abstar.run() Python API`.

---

### Task 11: Typer CLI

**Files:**
- Create: `python/abstar/cli.py`
- Test: `tests/python/test_cli.py` (subprocess calls to the installed `abstar` entry point)

**Interfaces:**
- **Consumes:** `api.run`, `germline.*`, `errors.*`, `_core.run_to_stdout`.
- **Produces:** `app: typer.Typer`, `germline_app: typer.Typer`, and `def main() -> None` (the entry point). `main` catches `AbstarError` and exits with its `exit_code`.
- **Commands** (spec: *CLI (Typer)*):

  ```
  abstar --version
  abstar run [SEQUENCES] [OUTPUT_DIR] [--seq TEXT] [--germline-db/-g TEXT=human] [--locus TEXT]... [--format [airr|parquet]]... [--threads INT] [--strict] [--json] [--quiet/-q] [--verbose/-v]
  abstar germline build RECIPE [--output DIR] [--overwrite] [--json]
  abstar germline list [--json]
  abstar germline info NAME [--json]
  abstar manual [TOPIC] [SUBTOPIC]        # implemented in Task 12; register the command here as a stub that Task 12 fills
  ```

**Rules for `run`:**

- **Inputs:** exactly one of SEQUENCES or `--seq` is required (`INVALID_ARGUMENT` otherwise). SEQUENCES `-` means stdin.
- **No OUTPUT_DIR:**
  - TSV goes to stdout via `_core.run_to_stdout`.
  - `--format parquet` → `INVALID_ARGUMENT`.
  - A one-line summary, or the JSON summary with `--json`, goes to stderr.
- **With OUTPUT_DIR:** call `api.run`, then print the summary to stdout (human-readable, or the `RunSummary` as JSON with `--json`).
- **Errors:**
  - Without `--json`: `error[CODE]: message`.
  - With `--json`: `{"error": {"code": ..., "message": ...}}` on stderr.
  - Exit with the code's exit code. Typer usage errors exit 2.
- **`--quiet`** suppresses the summary, though `--json` output is still printed. **`--verbose`** also prints the resolved database path and the input list to stderr before running.
- **Imports:** heavy imports (`polars`, `abstar._core`, `rich`) happen inside command bodies only.

- [ ] **Step 1: Write the failing tests** in `tests/python/test_cli.py`
  - `test_seq_to_stdout`: `abstar run --seq ACGT` exits 0. Its stdout is exactly the TSV header plus one row (2 lines, first cell of row 2 is `1`), and nothing but TSV appears on stdout.
  - `test_stdin`: piping a 3-record FASTA into `abstar run -` gives 4 stdout lines.
  - `test_broken_pipe`: run `abstar run <10k-record fasta>` with stdout read for 1 line and then closed. Exit code is 0, and stderr contains no `Traceback`.
  - `test_output_dir_json_summary`: `abstar run in.fasta out --json --format airr --format parquet` prints JSON whose `counts` sum equals the record count. The files exist.
  - `test_parquet_without_output_dir`: exit 2, and stderr contains `error[INVALID_ARGUMENT]`.
  - `test_json_error_shape`: `abstar run missing.fasta out --json` exits 2. Stderr parses as JSON with `error.code == "INPUT_NOT_FOUND"`.
  - `test_germline_commands`:
    - `abstar germline list --json` includes `human` with kind `builtin`
    - `abstar germline info human --json` has `counts["IGH/V/human"] > 0`
    - `abstar germline build <recipe> --json` (with `ABSTAR_HOME` set to tmp) prints a report
  - `test_no_seq_and_no_input`: `abstar run` exits 2.
- [ ] **Step 2: Run the tests and confirm they fail.** Run `python -m pytest tests/python/test_cli.py -q`.
- [ ] **Step 3: Implement** `cli.py`. In `run_to_stdout`, Rust treats `BrokenPipe` as end of output. In Python, after an early-closed pipe, redirect stdout to devnull before exit, so interpreter shutdown doesn't raise.
- [ ] **Step 4: Run the tests and confirm they pass.** Run `python -m pytest tests/python -q`.
- [ ] **Step 5: Commit** with message `feat: Typer CLI for run and germline`.

---

### Task 12: `abstar manual`, plus generated docs

**Files:**
- Create:
  - `python/abstar/manual/__init__.py`
  - `python/abstar/manual/topics/{conventions,ambiguity}.md` (hand-written prose, the only hand-written part)
  - `docs/conf.py`, `docs/index.md`, `docs/requirements.txt` (`sphinx`, `myst-parser`, `furo`)
- Modify: `python/abstar/cli.py` (the `manual` command), `.github/workflows/rewrite.yml` (add a docs job: `sphinx-build -W -b html docs docs/_build`)
- Test: `tests/python/test_manual.py`

**Interfaces:**
- **Consumes:** `cli.app` (via `typer.main.get_command`), `_core.fields_toml()`, the decision files, `errors.ERROR_CODES`, `errors.EXIT_CODES`.
- **Produces:**

```python
TOPICS: dict[str, str]                       # topic -> one-line summary: commands, run, germline, manual, fields, conventions, ambiguity, errors, decisions
def render(topic: str | None = None, subtopic: str | None = None) -> str   # Markdown
```

**Rules:**

- **`render(None)`** concatenates every topic, in `TOPICS` order, under `#` headings.
- **Command topics** walk the Click command tree. For each command they list the usage line, the help text, and every parameter: names, type, default, help. Every command and option in `cli.app` must appear.
- **`fields`** lists every field: name, type, nullable, origin, coordinates, description, and its decision ID and title.
- **`fields NAME`** shows one field plus the full text of its decision record.
- **`decisions`** lists ID, title and status. **`decisions ID`** prints the record.
- **`errors`** renders both tables from `errors.py`.
- **Conventions and ambiguity:**
  - `conventions` covers coordinate spaces, TSV encoding, record statuses and `sequence_id` rules.
  - `ambiguity` covers the primary call and ties, and cross-references decision 0008.
- **Unknown topic or subtopic** → `INVALID_ARGUMENT`. The message lists the valid topics or names.
- **`docs/conf.py`** writes `render()` to `docs/_generated/manual.md` at build time. `docs/index.md` includes it with a MyST `{include}`, and the decision records are included the same way. Add `docs/_generated/` and `docs/_build/` to `.gitignore`.

- [ ] **Step 1: Write the failing tests** in `tests/python/test_manual.py`
  - `test_every_command_and_option_documented`: for every command and every param `opts` string in the Click tree, it appears in `render("commands")`.
  - `test_every_field_documented`: every `fields.toml` name appears in `render("fields")`.
  - `test_field_with_decision`: `render("fields", "v_call_ties")` contains `0008` and the record's `## Decision` text.
  - `test_errors_topic`: every code in `ERROR_CODES` and every exit code appears.
  - `test_unknown_topic`: `abstar manual bogus` exits 2, and stderr lists `fields`.
  - `test_cli_manual_stdout`: `abstar manual` exits 0 and stdout starts with `# abstar`.
- [ ] **Step 2: Run the tests and confirm they fail.** Run `python -m pytest tests/python/test_manual.py -q`.
- [ ] **Step 3: Implement** `manual/__init__.py`, the two prose topics, the CLI command and the docs files.
- [ ] **Step 4: Run the tests and confirm they pass.** Run `python -m pytest tests/python -q && pip install -r docs/requirements.txt && sphinx-build -W -b html docs docs/_build`. Expected: the tests pass and Sphinx builds with no warnings.
- [ ] **Step 5: Commit** with message `feat: generated abstar manual and docs`.

---

### Task 13: Benchmark harness

**Files:**
- Create: `crates/annotate/benches/pipeline.rs` (criterion; `[[bench]] harness = false`), `scripts/bench.py`
- Test: none beyond running both (they report; they are not gates)

**Interfaces:**
- **Consumes:** `run_to_writer` (with `io::sink()`), the installed `abstar` CLI, and `legacy/test_data/bcr_corpus/sequences.fasta.gz`.
- **Produces:**
  - **Criterion bench `pipeline_10k`:** parse, stub-annotate and write TSV for 10k corpus records into a sink.
  - **`scripts/bench.py [--replicate N=20] [--out PATH]`** prints JSON with:
    - `help_ms`: median of 10 `abstar --help` runs
    - `single_seq_ms`: median of 10 `abstar run --seq <first corpus sequence>` runs
    - `throughput_rps`: records/sec for `abstar run <corpus × N> <tmp out>`
    - `peak_rss_mb`: via `resource.getrusage(RUSAGE_CHILDREN)`
    - machine info (platform, CPU count)

    It also prints a pass/fail line against the spec targets (`help_ms ≤ 150`, `single_seq_ms ≤ 300`), and exits 0 either way.

- [ ] **Step 1: Implement both, then run them.** Run `cargo bench -p abstar_annotate --bench pipeline` and `python scripts/bench.py`. Expected: both complete. Record the numbers in the commit message. If `help_ms` or `single_seq_ms` misses its target, profile with `python -X importtime` and fix any eager import before committing.
- [ ] **Step 2: Commit** with message `feat: benchmark harness with startup and throughput targets`.
