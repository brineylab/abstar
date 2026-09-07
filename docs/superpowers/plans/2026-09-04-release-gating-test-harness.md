# Release-Gating Test Harness Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a release-gating test harness that conserves every input record, checks exact biological annotations on adjudicated real sequences, emits AIRR 2.0-compliant TSV, and blocks publication of untested artifacts.

**Architecture:** Keep fast deterministic logic in the default pytest scope, with explicit `integration`, `e2e`, and `slow` boundaries. A deterministic discovery tool nominates records from the external Figshare corpus; only small attributed and adjudicated fixtures enter the repository. Annotation keeps Python-native 0-based half-open coordinates internally, while a single serializer converts to AIRR 2.0's 1-based closed coordinates and textual values.

**Tech Stack:** Python 3.10-3.13, pytest, pytest-cov, Hypothesis, AIRR Python 2.0.0, Polars, PyArrow, MMseqs2 through abutils, Click CliRunner, GitHub Actions, Sphinx.

**Spec:** `docs/superpowers/specs/2026-09-04-release-gating-test-harness-design.md`

## Global Constraints

- Support Python 3.10 and newer; the CI matrix covers 3.10, 3.11, 3.12, and 3.13.
- Target the AIRR Data Standards 2.0 Rearrangement schema exactly for AIRR TSV output.
- Keep internal query coordinates 0-based and half-open; convert once at the AIRR serialization boundary to 1-based closed intervals.
- Preserve every external sequence identifier as a string and use a separate immutable `row_id` for joins, ordering, and duplicate identifiers.
- Every input record ends as an `annotated` or `unassigned` public row, or as a structured failure attached to a raised run-level exception.
- Internal exceptions and external-tool failures must never become a successful empty or partial result.
- Preserve deterministic row ordering across input shapes, chunk sizes, and process counts.
- Use only pytest temporary directories for generated databases, inputs, outputs, and logs; redirect user germline paths away from the developer's real home.
- The packaged BCR databases are `human`, `macaque`, `c57bl6`, `balbc`, and `human+c57bl6`; the packaged TCR database is `human`.
- Published BCR fixtures come from Figshare record `19617633`, version 3, under the MIT license; the 2 GB source corpus stays outside the repository.
- Current `abstar` output and Cell Ranger output are evidence for adjudication, never self-generated golden truth.
- Do not add unconditional skips or xfails for supported behavior. A supported matrix cell finishes without failures, skips, xfails, xpasses, or project-originated warnings.
- Any production correction begins with a focused failing regression test and lands in a separate reviewable commit.

## File and Responsibility Map

### New files

- `requirements-test.txt`: test-only dependency ranges shared by local development and CI.
- `abstar/annotation/airr.py`: AIRR 2.0 field selection, coordinate conversion, value encoding, CIGAR construction, and TSV writing.
- `abstar/core/results.py`: structured per-record failures and the run-level exception type.
- `abstar/tests/helpers.py`: immutable fixture loading, normalized row comparison, and cardinality assertions.
- `abstar/tests/corpus.py`: deterministic hashing, candidate models, real-case loading, and derived-sequence transformations.
- `abstar/tests/test_corpus.py`: discovery selection, input-schema, and report reproducibility tests.
- `scripts/discover_bcr_cases.py`: optional command-line entry point over the external published corpus.
- `abstar/tests/data/germline_counts.json`: reviewed packaged-database gene counts.
- `abstar/tests/data/real_bcr/README.md`: source, license, attribution, selection, and adjudication notes.
- `abstar/tests/data/real_bcr/sequences.fasta`: bounded original published sequences.
- `abstar/tests/data/real_bcr/cases.json`: source evidence and accepted exact or explicitly ambiguous expectations.
- `abstar/tests/data/real_bcr/derived_cases.json`: reproducible mutations of accepted parent sequences.
- `abstar/tests/data/tcr/README.md`, `sequences.fasta`, and `cases.json`: human TCR golden provenance and expectations.
- `abstar/tests/test_database_integrity.py`: packaged FASTA/index/manifest invariants.
- `abstar/tests/test_properties.py`: Hypothesis invariants for coordinates, alignments, identifiers, ranking, and transformations.
- `abstar/tests/test_real_bcr.py`: original and derived BCR biological regressions.
- `abstar/tests/test_tcr_e2e.py`: TRA, TRB, TRD, and TRG public-pipeline goldens.
- `abstar/tests/test_airr.py`: AIRR 2.0 structural and semantic validation.
- `abstar/tests/test_output_parity.py`: AIRR/Parquet normalization and equality.
- `abstar/tests/test_failure_contracts.py`: invalid input, external failure, worker failure, and partial-run behavior.
- `abstar/tests/test_public_e2e.py`: API, installed CLI, input-shape, chunking, and multiprocessing matrix.
- `abstar/tests/test_coverage_contract.py`: critical-module coverage floor enforcement.
- `scripts/check_coverage.py`: validate an existing coverage JSON report against recorded per-module floors.
- `coverage-floors.json`: measured nondecreasing critical-module coverage floors.
- `.coveragerc`: measured package-wide coverage floor and branch-coverage settings.
- `.github/workflows/nightly-corpus.yml`: optional larger sentinel cohort using an explicitly provisioned corpus.

### Existing files changed

- `pyproject.toml`: pytest strictness, marker registration, and warning policy.
- `requirements.txt`: remove pytest from runtime dependencies.
- `abstar/assigners/mmseqs.py`: carry `row_id`, retain unassigned rows, preserve ordering, and expose checked external failures.
- `abstar/annotation/annotator.py`: return structured chunk outcomes, repair real-sequence regressions, assemble AIRR-correct alignments, and populate CIGAR fields.
- `abstar/annotation/antibody.py`: add public annotation status fields without exposing the internal row key.
- `abstar/annotation/schema.py`: separate internal and public schemas and add explicit status/AIRR fields.
- `abstar/annotation/mask.py`: derive masks from actual assembled segment spans.
- `abstar/annotation/productivity.py`: evaluate V-to-J frame in the correct coordinate origin.
- `abstar/core/abstar.py`: validate public arguments, aggregate outcomes, enforce conservation, sort deterministically, call format-specific writers, and raise on internal failures.
- `abstar/core/germline.py`: only changes required by focused custom-database failure tests.
- `abstar/preprocess/merging.py`: only changes required by checked-process and cleanup tests.
- `abstar/scripts/abstar.py`: map run failures to a nonzero Click exit with a concise diagnostic.
- Existing modules under `abstar/tests/`: remove stale markers, consolidate duplicate fixtures, and replace field-existence assertions with contract assertions.
- `.github/workflows/pytest.yml`: fast and supported-version jobs.
- `.github/workflows/python-publish.yml`: publish the exact artifacts produced by the successful release gate.
- `README.md`, `AGENTS.md`, and `docs/source/{installation,python_api,cli,output_formats}.rst`: supported commands, return/failure behavior, AIRR 2.0 semantics, and migration notes.

## Spec Traceability

| Design requirement | Implemented by |
| --- | --- |
| Fast/integration/e2e/slow scopes and strict markers | Tasks 1, 2, 14, 17 |
| Property and invariant testing | Tasks 3, 4, 8 |
| Published BCR discovery, provenance, adjudication, and derived cases | Tasks 5-8 |
| Record identity, conservation, ordering, and visible failures | Tasks 9-11, 14, 17 |
| Exact BCR biology and productivity | Tasks 11-12 |
| TCR chains and all packaged BCR databases | Tasks 3, 13 |
| AIRR 2.0 semantics and reference validation | Task 15 |
| Parquet/AIRR parity | Task 16 |
| Dependency, tool, worker, merging, UMI, and custom-database failures | Tasks 17-19 |
| Coverage floors | Task 20 |
| Supported-version CI, installed artifacts, docs, and publication gate | Tasks 21-22 |

---

## Phase 1: Harness Foundation

### Task 1: Strict pytest configuration and test-only dependencies

**Files:**
- Create: `requirements-test.txt`
- Modify: `requirements.txt`
- Modify: `pyproject.toml`
- Test: `abstar/tests/test_schema.py`

**Interfaces:**
- Consumes: Python `>=3.10` from project metadata.
- Produces: registered markers `integration`, `e2e`, and `slow`; reproducible test dependency installation with `python -m pip install -r requirements-test.txt`.

- [ ] **Step 1: Record the untouched coverage baseline**

Before editing configuration or tests, run:

```bash
python -m pip install "pytest-cov>=6,<8"
python -VV
python -m pip show abstar abutils polars pyarrow parasail pytest
python -m pytest --cov=abstar --cov-branch --cov-report=term-missing --cov-report=json:/tmp/abstar-coverage-before.json -q
```

Record exact totals, module values, test counts, warnings, interpreter, dependencies, revision, and wall time in the execution notes. Preserve `/tmp/abstar-coverage-before.json` through Task 20.

- [ ] **Step 2: Add a configuration test that asserts pytest owns the expected test root and markers**

```python
from pathlib import Path

import pytest


def test_pytest_configuration_registers_test_scopes(pytestconfig):
    assert Path(pytestconfig.getini("testpaths")[0]) == Path("abstar/tests")
    markers = "\n".join(pytestconfig.getini("markers"))
    for name in ("integration", "e2e", "slow"):
        assert f"{name}:" in markers
```

- [ ] **Step 3: Run the test and record the missing configuration failure**

Run: `python -m pytest abstar/tests/test_schema.py::test_pytest_configuration_registers_test_scopes -q`

Expected: FAIL because `testpaths` and the three project markers are not configured.

- [ ] **Step 4: Add the test dependency file and strict pytest configuration**

```text
# requirements-test.txt
-e .
pytest>=8.4,<10
pytest-cov>=6,<8
hypothesis>=6.138,<7
airr==2.0.0
build>=1.2,<2
```

```toml
[tool.pytest.ini_options]
testpaths = ["abstar/tests"]
strict_config = true
strict_markers = true
xfail_strict = true
markers = [
  "integration: exercises a real component or external-tool boundary",
  "e2e: exercises abstar.run or the installed abstar command",
  "slow: bounded checks excluded from the ordinary pull-request matrix",
]
filterwarnings = [
  "error::Warning:abstar.*",
]
```

Remove `pytest` from `requirements.txt`; keep every runtime import declared there.

- [ ] **Step 5: Install and verify the fast and complete command selectors**

Run:

```bash
python -m pip install -r requirements-test.txt
python -m pytest --collect-only -q
python -m pytest -m "not integration and not e2e and not slow" --collect-only -q
```

Expected: both collections succeed, and the second excludes every explicitly marked integration, end-to-end, or slow test.

- [ ] **Step 6: Commit the harness configuration**

```bash
git add pyproject.toml requirements.txt requirements-test.txt abstar/tests/test_schema.py
git commit -m "test: configure strict pytest scopes"
```

### Task 2: Remove stale marker debt and obsolete fixtures

**Files:**
- Modify: `abstar/tests/test_as_dataframe.py`
- Modify: `abstar/tests/test_germline.py`
- Modify: `abstar/tests/test_mmseqs.py`
- Modify: `abstar/tests/test_pipeline.py`
- Modify: `abstar/tests/conftest.py`
- Delete: `abstar/tests/reference/bnab_hc_antibodies.pkl`
- Delete: `abstar/tests/reference/bnab_lc_antibodies.pkl`
- Delete: unused files under `abstar/tests/reference/OLD/` if present at execution time
- Modify: `abstar/core/abstar.py:450`

**Interfaces:**
- Consumes: strict xfail and warning policy from Task 1.
- Produces: a baseline with zero stale `xfail`/`skip` markers and no `pl.count()` project warning.

- [ ] **Step 1: Turn all five stale Polars xfails into ordinary assertions**

Remove the decorators from the multi-sequence dataframe and pipeline tests. Strengthen each assertion to require the exact input count and ordered IDs rather than `len(result) >= 1`.

```python
assert result.height == 3
assert result.get_column("sequence_id").to_list() == ["10E8", "10J4", "10M6"]
```

- [ ] **Step 2: Replace exception xfails with direct exception tests**

```python
with pytest.raises(FileNotFoundError, match="NotAGermlineDatabase"):
    get_germline_database_path("NotAGermlineDatabase", receptor="bcr")

with pytest.raises(ValueError, match="receptor"):
    get_germline_database_path("human", receptor="abc")
```

Apply the same direct `pytest.raises` pattern to the two missing/nonmatching gene tests.

- [ ] **Step 3: Replace nonexistent `mouse` skips with supported database assertions**

Parameterize initialization over `c57bl6` and `balbc`. The pipeline smoke moves to Task 9, where it uses a database-derived compatible sequence and checks an exact locus.

```python
@pytest.mark.parametrize("database", ["c57bl6", "balbc"])
def test_mmseqs_initialization_mouse_database(database, temp_directories):
    output_dir, log_dir = temp_directories
    assigner = MMseqs(output_dir, log_dir, database, receptor="bcr")
    assert Path(assigner.germdb_path).name == database
```

- [ ] **Step 4: Make the existing baseline warning-free**

Change `output_df.select(pl.count())` to `output_df.select(pl.len())` in `abstar/core/abstar.py` and add an output-count assertion to the existing file-output test.

- [ ] **Step 5: Mark existing tests by the boundary they exercise**

Keep pure input/schema/ranking tests unmarked. Mark tests that launch MMseqs, fastp, or a database builder as `integration`; mark tests that call `abstar.run()` or Click as `e2e`; add `slow` only when the measured runtime exceeds the ordinary bounded cohort. Verify collection with:

```bash
python -m pytest --collect-only -m "not integration and not e2e and not slow" -q
python -m pytest --collect-only -m integration -q
python -m pytest --collect-only -m e2e -q
```

- [ ] **Step 6: Confirm the binary references are unused and remove them**

Run: `rg -n "bnab_hc_antibodies|bnab_lc_antibodies|tests/reference" abstar`

Expected: no active Python reference. Delete the two 72 MB aggregate pickle fixtures and unused `OLD` copies; keep small text FASTA inputs until Task 5 replaces them where useful.

- [ ] **Step 7: Run the affected tests**

Run:

```bash
python -m pytest abstar/tests/test_as_dataframe.py abstar/tests/test_germline.py abstar/tests/test_mmseqs.py abstar/tests/test_pipeline.py -q
```

Expected: PASS with zero skipped, xfailed, xpassed, or project warnings in these modules.

- [ ] **Step 8: Commit the baseline cleanup**

```bash
git add abstar/core/abstar.py abstar/tests
git commit -m "test: remove stale baseline debt"
```

### Task 3: Automate packaged germline database integrity

**Files:**
- Create: `abstar/tests/data/germline_counts.json`
- Create: `abstar/tests/test_database_integrity.py`
- Create: `abstar/tests/helpers.py`

**Interfaces:**
- Consumes: packaged database paths resolved by `abstar.core.germline.get_germline_database_path`.
- Produces: `read_fasta_records(path: Path) -> dict[str, str]` and a reviewed count snapshot used by CI.

- [ ] **Step 1: Add the reviewed count snapshot**

```json
{
  "bcr": {
    "human": {"v": 342, "d": 31, "j": 23, "c": 103},
    "macaque": {"v": 1071, "d": 52, "j": 28, "c": 46},
    "c57bl6": {"v": 196, "d": 8, "j": 22, "c": 35},
    "balbc": {"v": 265, "d": 10, "j": 22, "c": 35},
    "human+c57bl6": {"v": 538, "d": 39, "j": 45, "c": 138}
  },
  "tcr": {
    "human": {"v": 246, "d": 6, "j": 86, "c": 21}
  }
}
```

- [ ] **Step 2: Write failing database invariants**

Parameterize over every receptor/database/segment in the JSON and assert:

```python
assert len(ungapped) == expected_count
assert set(gapped) == set(ungapped)
assert len(gapped) == len(set(gapped))
assert {key: value.replace(".", "") for key, value in gapped.items()} == ungapped
assert (database_path / "manifest.txt").is_file()
for suffix in ("", ".dbtype", ".index", ".lookup", ".source", "_h", "_h.dbtype", "_h.index"):
    assert (database_path / "mmseqs" / f"{segment}{suffix}").is_file()
```

Also assert BCR IDs start with `IG`, TCR IDs start with `TR`, D loci are only `IGH`, `TRB`, or `TRD`, and every FASTA ID is unique within its segment.

- [ ] **Step 3: Run the invariants and investigate every mismatch**

Run: `python -m pytest abstar/tests/test_database_integrity.py -q`

Expected: PASS. A mismatch is a database defect to correct from source FASTA and regenerated indexes; never edit binary indexes directly or weaken the expected count.

- [ ] **Step 4: Commit database integrity coverage**

```bash
git add abstar/tests/helpers.py abstar/tests/data/germline_counts.json abstar/tests/test_database_integrity.py
git commit -m "test: verify packaged germline integrity"
```

### Task 4: Add property tests for pure invariants

**Files:**
- Create: `abstar/tests/test_properties.py`
- Modify: `abstar/annotation/positions.py` only if a minimized counterexample proves a defect
- Modify: `abstar/annotation/indels.py` only if a minimized counterexample proves a defect
- Modify: `abstar/assigners/mmseqs.py` only if ranking order changes output

**Interfaces:**
- Consumes: `get_gapped_position_from_raw`, `get_raw_position_from_gapped`, `select_best_hits`, and the indel annotation functions.
- Produces: reusable Hypothesis strategies `dna`, `gapped_dna`, `external_ids`, and `alignment_pairs` inside the test module.

- [ ] **Step 1: Define bounded deterministic strategies**

```python
dna = st.text(alphabet="ACGT", min_size=1, max_size=300)
external_ids = st.text(
    alphabet=st.characters(blacklist_categories=("Cs",), blacklist_characters="\t\n\r"),
    min_size=1,
    max_size=128,
)
```

Use `@settings(max_examples=200, deadline=None)` for pure functions and print the Hypothesis seed in CI on failure.

- [ ] **Step 2: Add reverse-complement and position round trips**

```python
@given(dna)
def test_reverse_complement_is_involutive(sequence):
    assert reverse_complement(reverse_complement(sequence)) == sequence

@given(dna, st.lists(st.integers(min_value=0, max_value=299), unique=True))
def test_imgt_gapped_round_trip(sequence, gap_offsets):
    gaps = {offset for offset in gap_offsets if offset < len(sequence)}
    gapped = "".join(("." if i in gaps else "") + base for i, base in enumerate(sequence))
    for raw_position in range(1, len(sequence) + 1):
        imgt = get_gapped_position_from_raw(raw_position, gapped)
        assert get_raw_position_from_gapped(imgt, gapped) == raw_position
```

- [ ] **Step 3: Add interval, indel, and ranking properties**

Assert segment intervals are ordered and within the query; applying insertion/deletion descriptions reconstructs the aligned pair; and shuffling hit rows cannot change `select_best_hits` calls or tie ordering.

```python
expected = select_best_hits(frame, "v").sort("v_query")
for permutation in permutations:
    actual = select_best_hits(frame[permutation], "v").sort("v_query")
    assert actual.equals(expected)
```

- [ ] **Step 4: Add identifier and iterable properties around `_process_inputs`**

Generate numeric, leading-zero, scientific-notation-like, Unicode, whitespace-bearing, long, and duplicate IDs. Assert `_process_inputs` writes the exact header bytes and sequence count for list, iterator, and generator inputs. The full MMseqs round trip for these IDs is covered in Task 10, where the immutable `row_id` prevents FASTA tools from interpreting external identifiers.

- [ ] **Step 5: Run and repair only minimized proven defects**

Run: `python -m pytest abstar/tests/test_properties.py -q`

Expected: PASS for 200 examples per property. Commit a focused production fix with its failing minimized example before continuing if Hypothesis finds a defect.

- [ ] **Step 6: Commit invariant coverage**

```bash
git add abstar/tests/test_properties.py abstar/annotation/positions.py abstar/annotation/indels.py abstar/assigners/mmseqs.py
git commit -m "test: add annotation invariant properties"
```

## Phase 2: Published-Data Discovery and Curation

### Task 5: Implement deterministic BCR corpus discovery primitives

**Files:**
- Create: `abstar/tests/corpus.py`
- Create: `abstar/tests/test_corpus.py`

**Interfaces:**
- Produces: `CorpusRecord`, `CandidateEvidence`, `stable_rank`, `select_sweep`, `normalize_gene`, and `load_cellranger_annotations`.

```python
@dataclass(frozen=True, slots=True)
class CorpusRecord:
    dataset: str
    sequence_id: str
    sequence: str
    donor: str
    flow_class: str
    chain: str | None

    @property
    def row_key(self) -> str:
        return f"{self.dataset}\0{self.sequence_id}"


def stable_rank(record: CorpusRecord, *, seed: str, algorithm_version: int = 1) -> str:
    payload = f"{algorithm_version}\0{seed}\0{record.row_key}".encode()
    return hashlib.sha256(payload).hexdigest()


def select_sweep(records: Iterable[CorpusRecord], *, per_dataset: int = 200, seed: str = "abstar-real-bcr-v1") -> list[CorpusRecord]:
    materialized = list(records)
    by_dataset: dict[str, list[CorpusRecord]] = defaultdict(list)
    for record in materialized:
        by_dataset[record.dataset].append(record)

    selected: list[CorpusRecord] = []
    for dataset in sorted(by_dataset):
        candidates = by_dataset[dataset]
        by_stratum: dict[tuple[str, str, str | None], list[CorpusRecord]] = defaultdict(list)
        for record in candidates:
            by_stratum[(record.donor, record.flow_class, record.chain)].append(record)

        chosen: dict[str, CorpusRecord] = {}
        for stratum in sorted(by_stratum, key=lambda value: tuple("" if item is None else item for item in value)):
            record = min(by_stratum[stratum], key=lambda value: stable_rank(value, seed=seed))
            chosen[record.row_key] = record

        for record in sorted(candidates, key=lambda value: stable_rank(value, seed=seed)):
            if len(chosen) >= min(per_dataset, len(candidates)):
                break
            chosen.setdefault(record.row_key, record)

        selected.extend(sorted(chosen.values(), key=lambda value: stable_rank(value, seed=seed)))
    return selected
```

- [ ] **Step 1: Test stable selection independent of input order**

Create two donors, two flow classes, all three BCR chains, and two datasets in memory. Assert identical ordered row keys after reversing and shuffling inputs, a maximum of 200 per dataset, and representation of every available stratum before filling remaining slots by stable hash.

- [ ] **Step 2: Test explicit CSV schemas and identifiers**

Write a temporary manifest and Cell Ranger CSV containing `10E8`, `00123`, and duplicate contig IDs in different datasets. Assert string preservation and uniqueness of `dataset + sequence_id`; reject a duplicate within one dataset with a cardinality error.

- [ ] **Step 3: Implement immutable models and deterministic selection**

Materialize each input iterable exactly once. Sort strata by `(dataset, donor, flow_class, chain)`, select the lowest stable hash from each stratum, then fill each dataset quota by stable hash. Include `algorithm_version=1` and the seed in returned report metadata.

- [ ] **Step 4: Run and commit the discovery primitives**

Run: `python -m pytest abstar/tests/test_corpus.py -q`

```bash
git add abstar/tests/corpus.py abstar/tests/test_corpus.py
git commit -m "test: add deterministic corpus selection"
```

### Task 6: Add the optional discovery command and candidate report

**Files:**
- Create: `scripts/discover_bcr_cases.py`
- Modify: `abstar/tests/test_corpus.py`

**Interfaces:**
- Consumes: Task 5 discovery primitives and `abstar.run()`.
- Produces: `main(argv: Sequence[str] | None = None) -> int` and JSON Lines candidate reports with one metadata header followed by one record per selected contig.

- [ ] **Step 1: Write a CLI test using a two-dataset temporary corpus**

Invoke:

```python
exit_code = main([
    "--fasta-dir", str(fasta_dir),
    "--manifest", str(manifest),
    "--cellranger-root", str(cellranger_root),
    "--output", str(report),
    "--per-dataset", "2",
    "--seed", "fixture-seed",
    "--n-processes", "1",
])
assert exit_code == 0
```

Assert that no path under `/home/` appears in the report and rerunning with reversed file discovery produces byte-identical JSON Lines.

- [ ] **Step 2: Implement checked path validation and report metadata**

The header records schema version, selection algorithm/seed, abstar git revision/version, Python/dependency versions, receptor/database, and the SHA-256 of each germline manifest. Each record carries source fields, raw/normalized calls, status, exception category, junction/CDR3 comparison, productivity comparison, ties, indels, no-D state, lengths, and selection reasons.

- [ ] **Step 3: Run the default external sweep**

Run:

```bash
python scripts/discover_bcr_cases.py \
  --fasta-dir /home/bryanbriney/Projects/lc_coherence/data/bcr_fastas \
  --manifest /home/bryanbriney/Projects/lc_coherence/data/sample_manifest.csv \
  --cellranger-root /home/bryanbriney/Projects/lc_coherence/data/raw_download \
  --output /tmp/abstar-bcr-candidates.jsonl \
  --per-dataset 200 \
  --seed abstar-real-bcr-v1
```

Expected: at most 18,800 selected records, no source-tree writes, and an explicit outcome for every selected record. Keep the report in `/tmp` or another external review location; do not commit it.

- [ ] **Step 4: Commit the discovery command**

```bash
git add scripts/discover_bcr_cases.py abstar/tests/test_corpus.py
git commit -m "test: add published BCR discovery command"
```

### Task 7: Curate and commit the original real-BCR fixture set

**Files:**
- Create: `abstar/tests/data/real_bcr/README.md`
- Create: `abstar/tests/data/real_bcr/sequences.fasta`
- Create: `abstar/tests/data/real_bcr/cases.json`
- Modify: `abstar/tests/corpus.py`
- Create: `abstar/tests/test_real_bcr.py`
- Modify: `abstar/tests/conftest.py`

**Interfaces:**
- Produces: `load_real_bcr_cases() -> tuple[RealBCRCase, ...]`, where every case verifies its sequence SHA-256 at load time.

- [ ] **Step 1: Define and test the fixture schema**

```python
@dataclass(frozen=True, slots=True)
class RealBCRCase:
    dataset: str
    sequence_id: str
    sequence: str
    sequence_sha256: str
    selection_reasons: tuple[str, ...]
    source: Mapping[str, object]
    expected: Mapping[str, object]
    evidence: tuple[str, ...]

    def as_sequence(self) -> Sequence:
        return Sequence(self.sequence, id=self.sequence_id)
```

Load the FASTA as an ordered record list and zip it to the ordered JSON case list; this permits the same external contig ID to occur in different datasets without rewriting either ID. The loader rejects an order/ID mismatch, missing FASTA record, duplicate `(dataset, sequence_id)`, hash mismatch, unknown expected field, absent evidence, or a gene call that is neither an exact string nor a sorted allowed set. Add `pilot_loss_cases` and `real_bcr_cases` fixtures to `conftest.py`; both return new immutable tuples on each request.

- [ ] **Step 2: Select a bounded 25-40-case cohort deterministically**

Start with all eight mandatory dataset `1279068` loss records:

```text
ACGATACCAGGTTTCA-1_contig_2
CAAGATCAGAGCTTCT-1_contig_2
CCATGTCCAGTCTTCC-1_contig_1
CTAAGACAGCAATCTC-1_contig_2
CTGTTTACAGGTGCCT-1_contig_1
GCTGCTTAGAAACGAG-1_contig_2
GGTGCGTAGGGAAACA-1_contig_1
GTAACTGAGTATTGGA-1_contig_2
```

Then take the lowest stable-hash adjudication candidate in each available bucket until the cohort contains no more than 40 unique records: two concordant controls per IGH/IGK/IGL; two productivity disagreements per locus; both junction disagreements; two tied calls per locus; two accepted insertions; two accepted deletions; two no-D IGH records; and shortest/longest junction candidates. When a candidate fails adjudication, record the rejection in the external report and select the next stable-hash candidate in that bucket.

- [ ] **Step 3: Adjudicate each selected record**

For every accepted case, manually review source sequence, germline alignments, coordinate spaces, Cell Ranger evidence, and the biological rule under test. Store exact locus, strand, accepted call or sorted ambiguity set, zero-based half-open internal segment coordinates, junction/CDR3, productivity, reason codes, status, and selection evidence. Exact allele calls require direct evidence; otherwise store the accepted gene-level ambiguity.

- [ ] **Step 4: Add source and licensing documentation**

The README names “Functional antibodies exhibit light chain coherence,” Figshare record `19617633`, version `3`, URL `https://figshare.com/articles/preprint/Functional_antibodies_exhibit_light_chain_coherence/19617633/3`, and MIT license. State that copied fixtures are redistributed under the repository MIT license and retain original dataset/contig IDs.

- [ ] **Step 5: Test fixture integrity without running MMseqs**

Run: `python -m pytest abstar/tests/test_real_bcr.py -m "not e2e" -q`

Expected: PASS for count bounds, unique source keys, sequence hashes, provenance, accepted expectation types, and required inclusion of all eight loss IDs.

- [ ] **Step 6: Commit the adjudicated original fixtures**

```bash
git add abstar/tests/corpus.py abstar/tests/conftest.py abstar/tests/test_real_bcr.py abstar/tests/data/real_bcr
git commit -m "test: add adjudicated published BCR fixtures"
```

### Task 8: Add reproducible real-background transformations

**Files:**
- Modify: `abstar/tests/corpus.py`
- Create: `abstar/tests/data/real_bcr/derived_cases.json`
- Modify: `abstar/tests/test_real_bcr.py`
- Modify: `abstar/tests/test_properties.py`

**Interfaces:**
- Produces: `derive_sequence(parent: str, operation: DerivedOperation) -> str` and frozen operation records for reverse complement, truncation, substitution/ambiguity, insertion, and deletion.

- [ ] **Step 1: Write exact transformation examples and hash checks**

```python
assert derive_sequence("AACCGG", DerivedOperation("substitute", 2, "C", "T")) == "AATCGG"
assert derive_sequence("AACCGG", DerivedOperation("insert", 3, "", "GGA")) == "AACGGACGG"
assert derive_sequence("AACCGG", DerivedOperation("delete", 1, "ACC", "")) == "AAG"
```

Reject an out-of-range offset or a parent slice that differs from `affected_bases`. Verify the generated SHA-256 against `derived_cases.json` before annotation.

- [ ] **Step 2: Add properties for reversible transformations**

Assert reverse complement twice restores the parent; deleting a just-inserted payload restores the parent; a substitution changes only the named offset; and derived operations never mutate their parent fixture.

- [ ] **Step 3: Create the derived-case matrix**

For clean adjudicated parents, add forward/reverse complements, 5′/3′ truncations, one ambiguous base, one substitution, 1/2/3-base insertions and deletions in V and junction contexts. Each entry includes parent key, operation, input-query offset, affected/replacement bases, resulting hash, and exact invariant or exact expected change.

- [ ] **Step 4: Run and commit transformations**

Run:

```bash
python -m pytest abstar/tests/test_properties.py abstar/tests/test_real_bcr.py -m "not e2e" -q
```

```bash
git add abstar/tests/corpus.py abstar/tests/test_properties.py abstar/tests/test_real_bcr.py abstar/tests/data/real_bcr/derived_cases.json
git commit -m "test: add real-background sequence variants"
```

## Phase 3: Biological and Public Contracts

### Task 9: Introduce structured record outcomes and run failures

**Files:**
- Create: `abstar/core/results.py`
- Create: `abstar/tests/test_failure_contracts.py`
- Modify: `abstar/annotation/antibody.py`
- Modify: `abstar/annotation/schema.py`
- Modify: `abstar/__init__.py`

**Interfaces:**
- Produces: `RecordFailure`, `AnnotationChunkResult`, and `AnnotationRunError`; public status fields `annotation_status` and `failure_reason`.

```python
FailureStage = Literal["preprocess", "assignment", "annotation", "output"]
FailureCategory = Literal["invalid_input", "unassigned", "external_tool", "internal_error"]


@dataclass(frozen=True, slots=True)
class RecordFailure:
    row_id: str
    sequence_id: str
    stage: FailureStage
    category: FailureCategory
    message: str
    traceback_text: str | None = None


@dataclass(frozen=True, slots=True)
class AnnotationChunkResult:
    output_path: str
    failures: tuple[RecordFailure, ...]
    failed_log_path: str | None
    succeeded_log_path: str | None


class AnnotationRunError(RuntimeError):
    def __init__(self, failures: Iterable[RecordFailure], partial_output_paths: Iterable[str] = ()):
        self.failures = tuple(failures)
        self.partial_output_paths = tuple(partial_output_paths)
        counts = Counter((failure.stage, failure.category) for failure in self.failures)
        summary = ", ".join(
            f"{stage}/{category}={count}"
            for (stage, category), count in sorted(counts.items())
        )
        super().__init__(f"abstar run failed: {summary}")
```

- [ ] **Step 1: Test immutability, pickling, and diagnostic formatting**

Round-trip each dataclass through `pickle`, assert mutation raises `FrozenInstanceError`, and assert `str(AnnotationRunError)` reports counts by stage/category without dumping sequence data or tracebacks.

- [ ] **Step 2: Add public status fields and separate the internal schema**

`OUTPUT_SCHEMA` gains `annotation_status: pl.String` and `failure_reason: pl.String`. Create `ANNOTATION_WORK_SCHEMA = {"row_id": pl.String, **OUTPUT_SCHEMA}`. Add `row_id` as an internal `Antibody` field and exclude it when `__post_init__` constructs `airr_fields`. `Antibody.to_dict()` must return a fresh field list and never mutate `self.airr_fields` when `include` is used.

- [ ] **Step 3: Export only the public exception**

Expose `AnnotationRunError` from `abstar.__init__`; keep chunk plumbing internal.

- [ ] **Step 4: Run and commit the outcome model**

Run: `python -m pytest abstar/tests/test_failure_contracts.py abstar/tests/test_schema.py abstar/tests/test_annotator.py -q`

```bash
git add abstar/core/results.py abstar/annotation/antibody.py abstar/annotation/schema.py abstar/tests/test_failure_contracts.py abstar/tests/test_schema.py abstar/tests/test_annotator.py abstar/__init__.py
git commit -m "feat: define structured annotation outcomes"
```

### Task 10: Enforce row identity, conservation, and deterministic order

**Files:**
- Modify: `abstar/assigners/mmseqs.py`
- Modify: `abstar/annotation/annotator.py`
- Modify: `abstar/core/abstar.py`
- Modify: `abstar/tests/test_failure_contracts.py`
- Modify: `abstar/tests/test_pipeline.py`

**Interfaces:**
- Consumes: Task 9 result types.
- Produces: `_assert_record_conservation(input_count: int, annotated_count: int, unassigned_count: int, failures: Sequence[RecordFailure]) -> None`; `annotate(input_file: str, output_directory: str, germline_database: str, log_directory: str | None = None, umi_pattern: str | None = None, umi_length: int | None = None, debug: bool = False) -> AnnotationChunkResult`.

- [ ] **Step 1: Add failing duplicate-ID and conservation tests**

Use input IDs `10E8`, `00123`, `duplicate`, and `duplicate`. Assert output preserves those four values in order while internal `row_id` values are unique. Add one unassignable `N` record and assert it appears as `annotation_status="unassigned"` with a nonempty `failure_reason` and null biological calls.

- [ ] **Step 2: Carry `row_id` through every assignment join**

Generate `row_id` from sample ordinal and record ordinal before MMseqs. Write it as the MMseqs query ID, retain the external `sequence_id` separately in an explicitly typed TSV, join every V/D/J/C result on `row_id`, validate one-to-one cardinality, and left-join against original input so records with no V hit remain present.

- [ ] **Step 3: Return unassigned rows and structured annotation failures**

For `v_call is None`, emit a public row with status `unassigned`; do not call `annotate_single_sequence`. For a caught internal exception, add `RecordFailure(category="internal_error", stage="annotation")` to `AnnotationChunkResult` and omit it from ordinary public rows.

Use `annotation_status="annotated"` and `failure_reason=None` for successful rows. Supported biological ambiguity such as `N` is not invalid input; reject only characters outside the documented IUPAC nucleotide alphabet and FASTA/FASTQ structural errors.

- [ ] **Step 4: Aggregate in submission order and fail on internal errors**

Keep the existing indexed `chunk_results` pattern, concatenate frames in sample/chunk order, sort by `row_id`, and drop `row_id` only at the public boundary. Raise `AnnotationRunError` after logs/failure artifacts are assembled if any internal or external failure exists. Preserve the input-shape return contract: one input returns one `Sequence`; multiple inputs return a list even when only one annotates; `as_dataframe=True` always returns a DataFrame.

- [ ] **Step 5: Prove the conservation equation**

```python
def _assert_record_conservation(input_count, annotated_count, unassigned_count, failures):
    accounted = annotated_count + unassigned_count + len(failures)
    if accounted != input_count:
        raise AnnotationRunError([
            RecordFailure("run", "<run>", "output", "internal_error", f"record conservation failed: input={input_count}, accounted={accounted}")
        ])
```

Call this once per sample before serialization and once for the complete API return.

- [ ] **Step 6: Run serial and multiprocessing tests**

Run:

```bash
python -m pytest abstar/tests/test_failure_contracts.py abstar/tests/test_pipeline.py abstar/tests/test_mmseqs.py -q
```

Expected: exact IDs, cardinality, statuses, and order match for `n_processes=1` and `n_processes=2`, with `chunksize=1` and `chunksize=3`.

- [ ] **Step 7: Commit record conservation**

```bash
git add abstar/assigners/mmseqs.py abstar/annotation/annotator.py abstar/core/abstar.py abstar/tests/test_failure_contracts.py abstar/tests/test_pipeline.py
git commit -m "fix: conserve records across annotation"
```

### Task 11: Repair the two observed real-sequence exception classes

**Files:**
- Modify: `abstar/annotation/annotator.py:478-530`
- Modify: `abstar/annotation/mask.py:179-245`
- Modify: `abstar/tests/test_real_bcr.py`

**Interfaces:**
- Consumes: the eight original cases and structured failures.
- Produces: `_clear_empty_d_alignment(ab: Antibody) -> Antibody` and masks whose length equals the corresponding ungapped assembled sequence.

- [ ] **Step 1: Add the eight-record end-to-end regression**

```python
@pytest.mark.e2e
def test_pilot_loss_cohort_has_no_internal_failures(pilot_loss_cases):
    result = abstar.run([case.as_sequence() for case in pilot_loss_cases], n_processes=1)
    assert [row.id for row in result] == [case.sequence_id for case in pilot_loss_cases]
    assert all(row["annotation_status"] in {"annotated", "unassigned"} for row in result)
```

Expected before fixes: `AnnotationRunError` containing four empty D global-alignment failures from `annotator.py:_segment_identities` and four gene-segment-mask length failures from `mask.py:generate_nongermline_mask`.

- [ ] **Step 2: Discard a D alignment that contributes zero residues**

After `process_dgene_alignment`, check both `ab.d_sequence` and `ab.d_germline` before identity alignment. If either is empty, clear every D call/score/identity/coordinate/sequence field and treat the full V-to-J interval as `np1`; never ask the alignment library to align an empty value.

```python
if not ab.d_sequence or not ab.d_germline:
    ab = _clear_empty_d_alignment(ab)
    d_loc = None
```

- [ ] **Step 3: Build gene-segment masks from assembled spans**

For nucleotides, concatenate `V * len(v_sequence)`, `N * len(np1)`, optional `D * len(d_sequence)`, optional `N * len(np2)`, and `J * len(j_sequence)`. For amino acids, derive the mask from codons in the aligned nucleotide mask using `ab.frame`, assigning a codon its segment letter only when all contributing nucleotides agree and `N` otherwise. Assert mask length equals ungapped assembled sequence length before returning.

- [ ] **Step 4: Assert exact adjudicated outcomes for the eight cases**

Use `cases.json` to assert locus, strand, accepted call sets, segment coordinates, junction/CDR3, productivity, and reasons for each record. Do not accept “did not crash” as sufficient.

- [ ] **Step 5: Run and commit each exception-class fix separately**

Run after the D fix: `python -m pytest abstar/tests/test_real_bcr.py -k "pilot_loss and empty_d" -q`

```bash
git add abstar/annotation/annotator.py abstar/tests/test_real_bcr.py
git commit -m "fix: ignore empty D realignments"
```

Run after the mask fix: `python -m pytest abstar/tests/test_real_bcr.py abstar/tests/test_mask.py -k "pilot_loss or mask" -q`

```bash
git add abstar/annotation/mask.py abstar/tests/test_real_bcr.py abstar/tests/test_mask.py
git commit -m "fix: derive masks from assembled segments"
```

### Task 12: Enforce real-BCR biological goldens and correct frame origin

**Files:**
- Modify: `abstar/tests/test_real_bcr.py`
- Modify: `abstar/annotation/productivity.py:76-116`
- Modify: `abstar/tests/test_productivity.py`

**Interfaces:**
- Consumes: all accepted original and derived cases.
- Produces: `junction_is_in_frame(junction_start: int, v_sequence_start: int, frame: int) -> bool`.

- [ ] **Step 1: Add parameterized exact golden assertions**

For every original case, assert exact status, original ID, row count, locus, strand, gene-level or allele ambiguity set, retained score/support evidence, internal V/D/J/C coordinates, junction, CDR3, productivity, and reason codes. When C is assigned, also assert C identity and equal-length gapped C query/germline output. For every derived case, assert the named invariant and exact changed fields, including mutation/indel representation and reconstruction.

- [ ] **Step 2: Reproduce the frame-origin defect**

Add a unit case where `v_sequence_start=137`, `junction_start=437`, and `frame=1`. It is in frame because `(437 - 137 - 0) % 3 == 0`; the existing absolute-query calculation rejects it.

```python
def junction_is_in_frame(junction_start, v_sequence_start, frame):
    return (junction_start - v_sequence_start - (frame - 1)) % 3 == 0
```

- [ ] **Step 3: Use the V-region origin for productivity**

Replace the absolute-query modulo calculation with `junction_is_in_frame`. Keep the existing checks for valid frame, minimum junction length, length modulo three, ambiguity, stop codons, V/J locus agreement, and receptor-specific terminal motif.

- [ ] **Step 4: Compare adjudicated productivity disagreements**

Run: `python -m pytest abstar/tests/test_productivity.py abstar/tests/test_real_bcr.py -k "productiv or frame" -q`

Expected: the accepted real cases match the reviewed productive/reason values; any remaining disagreement is resolved by evidence, not by copying Cell Ranger.

- [ ] **Step 5: Commit productivity and golden coverage**

```bash
git add abstar/annotation/productivity.py abstar/tests/test_productivity.py abstar/tests/test_real_bcr.py
git commit -m "fix: evaluate productivity in V frame"
```

### Task 13: Add human TCR goldens and packaged-database smoke coverage

**Files:**
- Create: `abstar/tests/data/tcr/README.md`
- Create: `abstar/tests/data/tcr/sequences.fasta`
- Create: `abstar/tests/data/tcr/cases.json`
- Create: `abstar/tests/test_tcr_e2e.py`
- Modify: `abstar/tests/test_pipeline.py`

**Interfaces:**
- Produces: one exact golden each for TRA, TRB, TRD, and TRG; a smoke input for every packaged BCR database.

```python
@dataclass(frozen=True, slots=True)
class TCRCase:
    sequence_id: str
    sequence: str
    locus: str
    allowed_v_calls: tuple[str, ...]
    allowed_d_calls: tuple[str | None, ...]
    allowed_j_calls: tuple[str, ...]
    junction: str
    productive: bool


def load_tcr_cases() -> tuple[TCRCase, ...]:
    records = read_fasta_records(TCR_DATA_DIR / "sequences.fasta")
    definitions = json.loads((TCR_DATA_DIR / "cases.json").read_text())
    return tuple(TCRCase(sequence=records[item["sequence_id"]], **item) for item in definitions)
```

- [ ] **Step 1: Construct transparent TCR goldens from packaged germlines**

For each locus, assemble one sequence from named packaged human V and J alleles, plus a D allele for TRB/TRD, with a recorded junction payload and no somatic mutations. Record every source allele, concatenation boundary, and SHA-256 in the README and JSON so the synthetic construction is independently reproducible.

- [ ] **Step 2: Assert receptor and locus propagation end to end**

```python
@pytest.mark.e2e
@pytest.mark.parametrize("case", load_tcr_cases(), ids=lambda case: case.locus)
def test_tcr_goldens(case):
    row = abstar.run(case.as_sequence(), receptor="tcr", germline_database="human", n_processes=1)
    assert row["locus"] == case.locus
    assert row["v_call"] in case.allowed_v_calls
    assert row["j_call"] in case.allowed_j_calls
    assert row["d_call"] in case.allowed_d_calls
    assert row["junction"] == case.junction
    assert row["productive"] is case.productive
```

TRA/TRG require a null D call; TRB/TRD require a compatible D call. Reverse-complement each case and assert identical biological values with `rev_comp=True`.

- [ ] **Step 3: Generate database-compatible BCR smoke inputs**

For each BCR database, concatenate the first lexicographically sorted compatible V and J germline sequences with a documented in-frame junction. The test asserts an annotated row, database identity, expected species/locus, and a V/J call belonging to that database; it does not claim cross-species biological goldens.

- [ ] **Step 4: Run TCR and database smoke tests**

Run:

```bash
python -m pytest abstar/tests/test_tcr_e2e.py abstar/tests/test_pipeline.py -m e2e -q
```

Expected: TRA/TRB/TRD/TRG and all five BCR databases annotate without skips or internal failures.

- [ ] **Step 5: Commit receptor/database coverage**

```bash
git add abstar/tests/data/tcr abstar/tests/test_tcr_e2e.py abstar/tests/test_pipeline.py
git commit -m "test: add TCR and packaged database goldens"
```

### Task 14: Exercise public input shapes, chunking, multiprocessing, and CLI

**Files:**
- Create: `abstar/tests/test_public_e2e.py`
- Modify: `abstar/tests/conftest.py`
- Modify: `abstar/core/abstar.py` only for proven contract defects
- Modify: `abstar/scripts/abstar.py` only for CLI/API argument mismatches

**Interfaces:**
- Consumes: adjudicated BCR controls, structured failures, and stable output status fields.
- Produces: `assert_same_annotations(left: Sequence[Mapping[str, object]], right: Sequence[Mapping[str, object]], fields: Sequence[str]) -> None` in `abstar/tests/helpers.py`.

- [ ] **Step 1: Build one shared three-record BCR fixture**

Use one accepted IGH, IGK, and IGL case. Write equivalent FASTA, FASTQ, flat-directory, and nested-directory inputs under `tmp_path`; expose list, iterator, and generator factories that create fresh iterables per call.

- [ ] **Step 2: Parameterize the API matrix**

Run each input through `abstar.run` with `(n_processes, chunksize)` values `(1,1)`, `(1,2)`, `(2,1)`, and `(2,3)`. Assert exact row count, IDs, order, status, calls, coordinates, junction/CDR3, and productivity equal the serial list baseline.

- [ ] **Step 3: Test API return shapes and empty/invalid boundaries**

Assert a single input returns `Sequence`, multiple inputs return `list[Sequence]`, `as_dataframe=True` returns a Polars DataFrame even for one record, project mode returns `None` after writing, and empty list/directory, nonpositive process/chunk values, and unsupported formats raise before creating a project directory. Import `abstar.gl`, `abstar.pp`, and `abstar.tl` and assert their documented public functions resolve to the active implementations. With `copy_inputs_to_project=True`, assert nested relative paths are preserved.

- [ ] **Step 4: Test Click through its public command**

Use `click.testing.CliRunner` against `abstar.scripts.abstar.cli` for `--help`, `run --help`, one AIRR run, and a forced `AnnotationRunError`. Success is exit 0 with exact output rows; failure is nonzero with the structured summary and a failure artifact path.

- [ ] **Step 5: Run and commit the public matrix**

Run: `python -m pytest abstar/tests/test_public_e2e.py -m e2e -q`

```bash
git add abstar/tests/test_public_e2e.py abstar/tests/conftest.py abstar/tests/helpers.py abstar/core/abstar.py abstar/scripts/abstar.py
git commit -m "test: cover public annotation entry points"
```

### Task 15: Emit exact AIRR 2.0 Rearrangement TSV

**Files:**
- Create: `abstar/annotation/airr.py`
- Create: `abstar/tests/test_airr.py`
- Modify: `abstar/annotation/annotator.py`
- Modify: `abstar/annotation/schema.py`
- Modify: `abstar/core/abstar.py`

**Interfaces:**
- Produces: `to_airr_interval`, `build_cigar`, `to_airr_row`, and `write_airr_tsv`.

Define `AIRR_REQUIRED_FIELDS` as `sequence_id`, `sequence`, `rev_comp`, `productive`, `v_call`, `d_call`, `j_call`, `sequence_alignment`, `germline_alignment`, `junction`, `junction_aa`, `v_cigar`, `d_cigar`, and `j_cigar`. Define `AIRR_FIELDS` as those fields in that order followed by every public `OUTPUT_SCHEMA` field not already present; omit only `row_id`, which belongs to `ANNOTATION_WORK_SCHEMA`.

```python
AIRR_SCHEMA_VERSION = "2.0"


def to_airr_interval(start: int | None, end: int | None) -> tuple[int | None, int | None]:
    if start is None or end is None:
        return None, None
    if start < 0 or end <= start:
        raise ValueError(f"invalid half-open interval: [{start}, {end})")
    return start + 1, end


def build_cigar(aligned_query: str, aligned_germline: str, *, query_start: int, germline_start: int) -> str:
    if len(aligned_query) != len(aligned_germline):
        raise ValueError("aligned query and germline must have equal lengths")
    operations: list[str] = ["S"] * query_start + ["N"] * germline_start
    for query_base, germline_base in zip(aligned_query, aligned_germline):
        if query_base == germline_base == "-":
            raise ValueError("an alignment column cannot contain two gaps")
        if query_base == "-":
            operations.append("D")
        elif germline_base == "-":
            operations.append("I")
        else:
            operations.append("M")
    return "".join(
        f"{sum(1 for _ in group)}{operation}"
        for operation, group in groupby(operations)
    )


def to_airr_row(row: Mapping[str, object]) -> dict[str, object]:
    output = {field: row.get(field) for field in AIRR_FIELDS}
    output["sequence"] = row.get("sequence_input")
    for segment in ("v", "d", "j", "c"):
        for axis in ("sequence", "germline"):
            start_field = f"{segment}_{axis}_start"
            end_field = f"{segment}_{axis}_end"
            output[start_field], output[end_field] = to_airr_interval(
                row.get(start_field), row.get(end_field)
            )
    for region in ("fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4"):
        start_field = f"{region}_start"
        end_field = f"{region}_end"
        output[start_field], output[end_field] = to_airr_interval(
            row.get(start_field), row.get(end_field)
        )
    return output


def write_airr_tsv(frame: pl.DataFrame, path: str | Path) -> None:
    rows = [to_airr_row(row) for row in frame.iter_rows(named=True)]
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=AIRR_FIELDS,
            delimiter="\t",
            lineterminator="\n",
            quoting=csv.QUOTE_NONE,
        )
        writer.writeheader()
        for row in rows:
            encoded = {}
            for field, value in row.items():
                if isinstance(value, bool):
                    value = "T" if value else "F"
                elif value is None:
                    value = ""
                elif "\t" in str(value) or "\n" in str(value) or "\r" in str(value):
                    raise ValueError(f"AIRR field {field!r} contains a forbidden delimiter")
                encoded[field] = value
            writer.writerow(encoded)
```

- [ ] **Step 1: Test the serialization boundary before changing annotation**

Assert `[0, 1) -> [1, 1]`, `[137, 439) -> [138, 439]`, null pairs remain null, and empty/reversed intervals raise. Test both query and germline coordinate fields. Booleans encode only `T`/`F`, nulls encode as empty strings, and tabs/newlines in values are rejected. `to_airr_row` must not mutate its source mapping.

- [ ] **Step 2: Test CIGAR construction and reconstruction**

Use AIRR-compatible `M`, `I`, `D`, leading `S`, and leading `N` operations. Assert equal runs collapse, `aligned_query`/`aligned_germline` lengths match, and replaying CIGAR operations reconstructs the segment alignment and the published coordinate spans.

```python
assert build_cigar("AC-GT", "ACCGT", query_start=2, germline_start=1) == "2S1N2M1D2M"
```

- [ ] **Step 3: Populate segment CIGAR and alignment fields from actual alignments**

Store V/D/J/C CIGAR immediately after each accepted segment realignment. Assemble `sequence_alignment` from observed aligned segment strings plus observed NP bases; assemble `germline_alignment` over the same columns with gaps for NP bases. Never copy NP nucleotides into germline. Keep existing assembled sequence fields as documented abstar extensions.

- [ ] **Step 4: Map official fields explicitly**

`sequence` is the unmodified query (`sequence_input`); when `rev_comp=True`, all calls, alignments, and coordinates refer to its reverse complement as AIRR requires. Add D segment coordinates plus `fwr1_start/end`, `cdr1_start/end`, `fwr2_start/end`, `cdr2_start/end`, `fwr3_start/end`, `cdr3_start/end`, and `fwr4_start/end` to `Antibody` and `OUTPUT_SCHEMA`. Convert region-local boundaries to oriented-query coordinates when annotation assigns them. Map official calls, identities, support, coordinates, junction/CDR3, regions, productivity, and CIGAR fields by name. Emit official required fields first in AIRR schema order, then public schema fields in their declared order.

- [ ] **Step 5: Validate with AIRR Python 2.0.0 and semantic assertions**

```python
reader = airr.read_rearrangement(output_path, validate=True)
rows = list(reader)
assert len(rows) == expected_count
```

Also inspect raw TSV text for `T`/`F` and empty nulls. Assert every non-null coordinate slices the correctly oriented query after conversion back to `[start - 1:end]`; junction includes conserved endpoints while CDR3 excludes them; alignment strings have equal length; and germline NP columns are gaps.

- [ ] **Step 6: Run AIRR tests against forward, reverse, indel, truncated, productive, nonproductive, and unassigned cases**

Run: `python -m pytest abstar/tests/test_airr.py -q`

Expected: AIRR Python 2.0.0 validation succeeds for every file and all semantic assertions pass.

- [ ] **Step 7: Commit AIRR compliance**

```bash
git add abstar/annotation/airr.py abstar/annotation/annotator.py abstar/annotation/schema.py abstar/core/abstar.py abstar/tests/test_airr.py
git commit -m "fix: emit AIRR 2.0 rearrangements"
```

### Task 16: Prove AIRR and Parquet output parity

**Files:**
- Create: `abstar/tests/test_output_parity.py`
- Modify: `abstar/tests/helpers.py`
- Modify: `abstar/core/abstar.py` only for proven divergence

**Interfaces:**
- Produces: `normalize_airr_row(row) -> dict[str, object]` and `normalize_parquet_row(row) -> dict[str, object]` for shared official fields.

- [ ] **Step 1: Run one annotation once and write both formats**

Use `output_format=["airr", "parquet"]` on the three-locus BCR fixture. At the file-output boundary, map Parquet's official `sequence` field to the unmodified `sequence_input`, while keeping native booleans, nulls, and 0-based half-open coordinates. Assert each file has exactly three ordered rows and original identifiers.

- [ ] **Step 2: Normalize only representation differences**

Convert AIRR `T`/`F` to booleans, empty strings to null, and 1-based closed coordinates back to 0-based half-open. Keep gene ambiguity strings, unmodified query sequence, aligned sequences, coordinates, identities, junction/CDR3, productivity, status, and reason fields otherwise unchanged. Do not normalize away a genuine logical mismatch.

- [ ] **Step 3: Compare every shared field**

```python
for airr_row, parquet_row in zip(airr_rows, parquet_rows):
    assert normalize_airr_row(airr_row) == normalize_parquet_row(parquet_row)
```

Repeat for reverse-complement, indel, unassigned, and duplicate-ID inputs.

- [ ] **Step 4: Run and commit output parity**

Run: `python -m pytest abstar/tests/test_output_parity.py -q`

```bash
git add abstar/tests/test_output_parity.py abstar/tests/helpers.py abstar/core/abstar.py
git commit -m "test: enforce AIRR and Parquet parity"
```

### Task 17: Inject dependency, external-tool, worker, and output failures

**Files:**
- Modify: `abstar/tests/test_failure_contracts.py`
- Modify: `abstar/core/abstar.py` only for failure propagation/cleanup defects
- Modify: `abstar/assigners/mmseqs.py` only for checked-process diagnostics

**Interfaces:**
- Consumes: `AnnotationRunError` and `RecordFailure`.
- Produces: checked MMseqs errors with captured command, exit status, stdout, and stderr; deterministic cleanup after worker/output failures.

- [ ] **Step 1: Inject a missing dependency capability**

Monkeypatch `abutils.tl.translate` away and assert startup or first use raises an actionable run error naming `abutils` and `translate`; no ordinary output file may be reported as successful.

- [ ] **Step 2: Inject MMseqs and worker failures**

Monkeypatch the MMseqs wrapper to raise a checked-process error containing sentinel stdout/stderr, then assert a nonzero CLI exit, `external_tool` category, retained diagnostic text, and cleaned owned temporary paths. Monkeypatch one annotation worker to raise and assert an `internal_error`; a mixed two-record worker test must expose the failed ID rather than returning one silent survivor.

- [ ] **Step 3: Cover malformed and unwritable inputs**

Assert empty FASTA, malformed FASTQ, a character outside the documented IUPAC nucleotide alphabet, incomplete germline database, and unwritable output target produce distinct `invalid_input`, `assignment`, or `output` diagnostics without creating misleading AIRR/Parquet success artifacts.

- [ ] **Step 4: Run and commit core failure coverage**

Run: `python -m pytest abstar/tests/test_failure_contracts.py -q`

```bash
git add abstar/tests/test_failure_contracts.py abstar/core/abstar.py abstar/assigners/mmseqs.py
git commit -m "test: expose pipeline failures"
```

### Task 18: Enforce merging and UMI process contracts

**Files:**
- Modify: `abstar/tests/test_merging.py`
- Modify: `abstar/tests/test_umi.py`
- Modify: `abstar/preprocess/merging.py` only for checked-process defects
- Modify: `abstar/core/abstar.py` only for preprocessing propagation defects

**Interfaces:**
- Consumes: `AnnotationRunError` and checked fastp wrapper behavior.
- Produces: paired/interleaved record-conservation tests and visible fastp failures.

- [ ] **Step 1: Add paired and interleaved FASTQ fixtures under `tmp_path`**

Use two read pairs with distinct IDs and known overlaps. Assert the fixture reader sees exactly two pairs before invoking merge code.

- [ ] **Step 2: Test checked fastp arguments and failure diagnostics**

Monkeypatch the subprocess boundary to capture the argument list and return sentinel stdout/stderr with a nonzero exit. Assert no `shell=True` string is used, the run raises an `external_tool` failure containing exit status and captured diagnostics, and no merged file is reported.

- [ ] **Step 3: Prove merge and UMI record conservation**

Run real fastp tests under `@pytest.mark.integration`; assert merged record count/IDs and cleanup for paired and interleaved files. Keep pure UMI parsing tests unmarked. Through `abstar.run`, assert both UMI-present and UMI-absent records survive in original order.

- [ ] **Step 4: Run and commit preprocessing contracts**

Run: `python -m pytest abstar/tests/test_merging.py abstar/tests/test_umi.py -q`

```bash
git add abstar/tests/test_merging.py abstar/tests/test_umi.py abstar/preprocess/merging.py abstar/core/abstar.py
git commit -m "test: conserve preprocessing records"
```

### Task 19: Make custom-database construction transactional

**Files:**
- Modify: `abstar/tests/test_custom_germline.py`
- Modify: `abstar/core/germline.py` only for transactional-build defects

**Interfaces:**
- Produces: staged custom database builds that become discoverable only after all FASTA, manifest, and MMseqs index validations pass.

- [ ] **Step 1: Redirect all user database paths to `tmp_path`**

Snapshot the developer database path before each test, monkeypatch the effective user database root, and assert the real `~/.abstar/germline_dbs` tree has identical metadata after the test.

- [ ] **Step 2: Inject failure at every external build boundary**

Fail gapping and each MMseqs segment-index step in turn. Assert the final database path is absent, the staging directory is removed, the exception includes the failed argument list/status/stdout/stderr, and an existing database remains byte-identical.

- [ ] **Step 3: Prove successful validation and atomic visibility**

Build a minimal valid database under the redirected root. Assert V/J required files, optional D/C handling, manifest, gapped/ungapped identity, and MMseqs indexes are all validated inside a sibling staging directory before one atomic rename exposes the final database path.

- [ ] **Step 4: Run and commit custom-database contracts**

Run: `python -m pytest abstar/tests/test_custom_germline.py -q`

```bash
git add abstar/tests/test_custom_germline.py abstar/core/germline.py
git commit -m "fix: stage custom germline builds"
```

## Phase 4: Measured Release Gates

### Task 20: Measure and ratchet coverage

**Files:**
- Modify: `pyproject.toml`
- Create: `.coveragerc`
- Create: `coverage-floors.json`
- Create: `scripts/check_coverage.py`
- Create: `abstar/tests/test_coverage_contract.py`
- Modify: `AGENTS.md`

**Interfaces:**
- Consumes: stable fast/integration/e2e suite from Phases 1-3.
- Produces: package-wide and critical-module integer coverage floors equal to the measured post-suite baseline and never below the pre-change baseline.

- [ ] **Step 1: Load the untouched baseline recorded in Task 1**

Read `/tmp/abstar-coverage-before.json` and the Task 1 execution notes. If either is unavailable, check out Task 1's parent in a temporary worktree, install the same recorded dependency versions, rerun the baseline command there, and preserve the result outside both worktrees. Never measure the modified tree and label it “before.”

- [ ] **Step 2: Measure the stabilized suite**

Run the same command after Tasks 1-19 and extract integer branch+statement coverage for the package and these files: `annotation/annotator.py`, `annotation/germline.py`, `annotation/positions.py`, `annotation/regions.py`, `core/abstar.py`, `core/germline.py`, `assigners/mmseqs.py`, `preprocess/merging.py`, and `annotation/umi.py`.

- [ ] **Step 3: Test and implement the critical-module coverage checker**

Write a temporary coverage JSON with two file entries and a temporary floor JSON. Assert `check_coverage(report_path, floors_path) -> list[str]` returns an empty list at or above each floor and returns one diagnostic per file below its floor or absent from the report. The CLI exits 0 for no diagnostics and 1 after printing sorted diagnostics otherwise.

- [ ] **Step 4: Set nondecreasing floors**

Set the package `fail_under` in `.coveragerc` to the stabilized total rounded down, provided it is not below the recorded pre-change total. Store each critical file and its stabilized integer percentage in `coverage-floors.json`, again never below its pre-change value. The dedicated coverage job writes `/tmp/abstar-coverage.json` and then runs `scripts/check_coverage.py /tmp/abstar-coverage.json coverage-floors.json` after pytest creates the report.

- [ ] **Step 5: Verify both gates fail when coverage is below the floor**

Run once with `--cov-fail-under` one point above the measured value and observe a nonzero exit. Lower one temporary report entry below its per-file floor and observe `check_coverage.py` exit 1. Restore the real report and committed floors and observe both commands succeed.

- [ ] **Step 6: Commit measured coverage gates**

```bash
git add .coveragerc coverage-floors.json scripts/check_coverage.py pyproject.toml AGENTS.md abstar/tests/test_coverage_contract.py
git commit -m "test: ratchet critical coverage"
```

### Task 21: Split CI, test installed artifacts, and gate publication

**Files:**
- Modify: `.github/workflows/pytest.yml`
- Create: `.github/workflows/nightly-corpus.yml`
- Modify: `.github/workflows/python-publish.yml`
- Modify: `requirements.txt`

**Interfaces:**
- Consumes: marker commands, test requirements, coverage floors, database checks, AIRR tests, and docs commands.
- Produces: named required jobs `fast`, `integration-e2e`, `dependency-bounds`, `database-airr`, `coverage`, `build-artifacts`, `non-linux-install`, and `docs`.

- [ ] **Step 1: Add the fast supported-Python matrix**

Run `python -m pytest -m "not integration and not e2e and not slow" -q` on Ubuntu with Python 3.10-3.13 after installing `requirements-test.txt`. Upload the pytest report on failure.

- [ ] **Step 2: Add bounded integration and end-to-end jobs**

Run integration/e2e on Python 3.10 and 3.13 with real MMseqs and packaged databases. Run database integrity and AIRR conformance as separate focused jobs so failures are attributable.

- [ ] **Step 3: Add dependency-bound jobs**

Before encoding the job, test candidate bounds for `abutils`, Polars, and PyArrow on Python 3.10 using the fast suite plus one BCR/TCR sentinel. Raise `abutils` above the known-broken 0.5 series and declare the lowest passing versions plus reviewed exclusive upper bounds in `requirements.txt`; document the evidence in the commit. Install the lowest declared compatible versions in one job and newest permitted versions in another, then run the fast suite and one BCR/TCR end-to-end sentinel. The job prints installed versions before testing and fails if pip resolves outside the declared ranges.

- [ ] **Step 4: Build once and test the exact artifacts**

Build wheel and sdist using `python -m build`, upload them as a named artifact, install each into a clean virtual environment, and run imports, `abstar --help`, `abstar run --help`, one CLI annotation, and one Python API annotation. The CLI test asserts biological output, not only exit status.

- [ ] **Step 5: Add non-Linux install and documentation jobs**

On macOS, install the wheel and verify imports/CLI discovery without claiming external-binary annotation support. Install `docs/doc_requirements.txt`, then build Sphinx with `python -m sphinx -W --keep-going -b html docs/source docs/_build/html`.

- [ ] **Step 6: Add the optional nightly corpus workflow**

Make it `workflow_dispatch` plus scheduled, require an explicit configured corpus path/artifact, run a bounded sentinel cohort, and upload the candidate/outcome report. If the corpus is not provisioned, the job fails with an explicit configuration message rather than silently skipping.

- [ ] **Step 7: Make publishing consume the tested artifact**

On release, call or depend on the complete test workflow, download the `build-artifacts` output associated with the release commit, verify its SHA-256 manifest, and pass those exact files to `pypa/gh-action-pypi-publish`. Remove the independent rebuild from `python-publish.yml`.

- [ ] **Step 8: Validate workflow syntax and commit CI gates**

Run any repository-available YAML parser plus `python -m pytest --collect-only -q`. Review every `needs:` edge manually: publishing must be unreachable until all required jobs succeed.

```bash
git add .github/workflows/pytest.yml .github/workflows/nightly-corpus.yml .github/workflows/python-publish.yml requirements.txt
git commit -m "ci: gate releases on tested artifacts"
```

### Task 22: Document contracts and run final verification

**Files:**
- Modify: `README.md`
- Modify: `AGENTS.md`
- Modify: `docs/source/installation.rst`
- Modify: `docs/source/python_api.rst`
- Modify: `docs/source/cli.rst`
- Modify: `docs/source/output_formats.rst`

**Interfaces:**
- Consumes: final public behavior and exact commands from all previous tasks.
- Produces: contributor and user documentation aligned with implementation.

- [ ] **Step 1: Document AIRR 2.0 semantics and the compatibility boundary**

State that Python internals and Parquet use documented 0-based half-open coordinates, AIRR TSV uses 1-based closed coordinates, booleans are `T`/`F`, nulls are empty, `sequence` is the unmodified query, `rev_comp` controls the orientation of annotations, and non-templated bases are gaps in `germline_alignment`.

- [ ] **Step 2: Document record outcomes and migration behavior**

Explain `annotation_status`, `failure_reason`, unassigned rows, `AnnotationRunError`, partial output artifacts, stable input-shape return types, duplicate identifier preservation, and deterministic ordering. Include before/after examples for a mixed annotated/unassigned input and an internal failure.

- [ ] **Step 3: Document the test scopes and commands**

Add exact commands for fast, integration, end-to-end, slow, complete, coverage, AIRR, database integrity, docs, and optional corpus discovery. State that corpus discovery requires explicit local paths and never runs in ordinary CI.

- [ ] **Step 4: Run focused verification**

```bash
python -m pytest abstar/tests/test_database_integrity.py abstar/tests/test_properties.py -q
python -m pytest abstar/tests/test_real_bcr.py abstar/tests/test_tcr_e2e.py -q
python -m pytest abstar/tests/test_airr.py abstar/tests/test_output_parity.py -q
python -m pytest abstar/tests/test_failure_contracts.py abstar/tests/test_public_e2e.py -q
```

Expected: all pass with zero skipped, xfailed, xpassed, or project warnings.

- [ ] **Step 5: Run the complete release verification**

```bash
python -VV
python -m pip show abstar abutils polars pyarrow parasail pytest pytest-cov hypothesis airr
python -m pytest --cov=abstar --cov-branch --cov-report=term-missing --cov-report=json:/tmp/abstar-coverage-final.json -q
python scripts/check_coverage.py /tmp/abstar-coverage-final.json coverage-floors.json
python -m build
python -m sphinx -W --keep-going -b html docs/source docs/_build/html
```

Install the built wheel and sdist in separate clean environments and run the installed CLI/API sentinel commands defined in Task 21. Record exact pass/fail/skip/xfail/xpass/warning counts, coverage, artifact hashes, interpreter, dependencies, and wall time.

- [ ] **Step 6: Check repository hygiene**

Run:

```bash
git diff --check
git status --short
find abstar/tests -type f -size +1M -print
```

Expected: no whitespace errors, only intended tracked changes, and no unexplained test fixture larger than 1 MB.

- [ ] **Step 7: Commit documentation and the verification record**

```bash
git add README.md AGENTS.md docs/source
git commit -m "docs: describe release-gating contracts"
```

The branch is ready for code review only after the complete verification output has been attached to the pull request and every required CI job is green.
