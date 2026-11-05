# Unified Public Output Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Every public abstar output (dataframe return, `Sequence` return, final
Parquet, AIRR TSV) carries identical field names, order, and meanings, built by
one conversion, and gains VDJ-only `sequence_vdj` / `sequence_vdj_aa` fields.

**Architecture:** One pure function, `to_public_row`, converts an `Antibody.to_dict()`
row (internal meanings: `sequence` is the assembled V(D)J, coordinates 0-based
half-open) into the public row (AIRR meanings, 1-based closed coordinates). It runs
once per record inside the annotation workers, so work Parquets already hold public
rows. `run()` only drops `row_id`; the Parquet writer writes the frame as-is; the
TSV writer only encodes text. `Antibody` objects keep internal conventions.

**Tech Stack:** Python 3.10+, Polars, abutils, pytest.

**Spec:** Approved in conversation on 2026-11-05 (no separate design doc). Its
decisions are restated verbatim in Global Constraints.

## Global Constraints

- New field names exactly `sequence_vdj` and `sequence_vdj_aa`, for all loci (no `vj` variant).
- `sequence_vdj` = internal assembled V(D)J (`Antibody.sequence`); `sequence_vdj_aa` = internal `Antibody.sequence_aa`. Both null unless `annotation_status == "annotated"`.
- Public `sequence` = `sequence_input`; public `sequence_aa`, `sequence_alignment_aa`, `germline_alignment_aa` keep the meanings final Parquet/TSV already use (`translate_airr_query` / `translate_airr_alignment`).
- Public coordinates are 1-based closed for every `X_start`/`X_end` pair in `OUTPUT_SCHEMA`: public start = internal start + 1, public end = internal end (numerically unchanged).
- `OUTPUT_SCHEMA` order is the single public field order: `AIRR_REQUIRED_FIELDS` first, then the remaining fields in existing schema order, `sequence_vdj` immediately after `sequence` and `sequence_vdj_aa` immediately after `sequence_aa` within the remainder.
- Booleans/nulls stay native in Parquet/dataframes and `T`/`F`/empty in TSV.
- `Antibody` attributes and `annotate_single_sequence` results are unchanged (internal, 0-based, VDJ `sequence`). Masks (`*_mask`, `*_mask_aa`) align to `sequence_vdj` / `sequence_vdj_aa`.
- Do not bump `abstar/version.py`. Docs label the change "0.9.0 (unreleased)".
- No new dependencies. Coverage floors in `.coveragerc` and `coverage-floors.json` must still pass; never lower them.
- The corpus baseline is only regenerated through the documented `--record-baseline` review procedure (`test_data/bcr_corpus/README.md`).
- Commits end with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.

## Review Focus

- Unassigned records (no V or J): public `sequence` is the original input, `sequence_vdj`/`sequence_vdj_aa` and all coordinates null — pinned in Task 1.
- Reverse-complemented input: public `sequence` is the original (un-reoriented) input while `sequence_vdj` is in the oriented sense — pinned in Task 2.
- Light chains / TCR without D: `d_*` coordinates stay null, `sequence_vdj` still populated — pinned in Task 1.
- A record whose public conversion raises (invalid interval) becomes an isolated `RecordFailure`; other records in the chunk are still written — pinned in Task 2.
- Running the audit script on 0.9 Parquet must not report off-by-one inconsistencies — pinned in Task 3.

---

### Task 1: Public schema order and `to_public_row`

**Files:**
- Modify: `abstar/annotation/schema.py` (add fields, define `AIRR_REQUIRED_FIELDS`, public order)
- Modify: `abstar/annotation/airr.py` (replace `to_airr_row` with `to_public_row`; `AIRR_FIELDS = tuple(OUTPUT_SCHEMA)`)
- Test: `abstar/tests/test_airr.py`, `abstar/tests/test_schema.py`

**Interfaces:**
- Produces: `abstar.annotation.schema.AIRR_REQUIRED_FIELDS: tuple[str, ...]` (moved from `airr.py`, same 14 names/order); `OUTPUT_SCHEMA` in public order including `sequence_vdj` and `sequence_vdj_aa` (both `pl.String`).
- Produces: `abstar.annotation.airr.to_public_row(row: Mapping[str, object]) -> dict[str, object]` — input is an `Antibody.to_dict()` row (internal meanings); output has exactly the keys of `OUTPUT_SCHEMA`, in that order. Pure; raises `ValueError` on invalid intervals (via `to_airr_interval`).
- Removes: `to_airr_row` (all callers move to `to_public_row` in Task 2).

- [ ] **Step 1: Write failing tests** in `abstar/tests/test_airr.py` (update `test_exact_public_field_order` and `test_row_conversion_is_pure_and_converts_every_coordinate` to target `to_public_row`), adding:

```python
def test_public_order_starts_with_required_and_places_vdj_fields():
    fields = list(OUTPUT_SCHEMA)
    assert tuple(fields[:len(AIRR_REQUIRED_FIELDS)]) == AIRR_REQUIRED_FIELDS
    assert fields.index("sequence_vdj_aa") > fields.index("sequence_aa")
    assert AIRR_FIELDS == tuple(OUTPUT_SCHEMA)

def test_public_row_maps_sequences_and_every_coordinate():
    internal = annotated_internal_row()  # helper: realistic annotated row, VDJ sequence, 0-based coords
    public = to_public_row(internal)
    assert list(public) == list(OUTPUT_SCHEMA)
    assert public["sequence"] == internal["sequence_input"]
    assert public["sequence_vdj"] == internal["sequence"]
    assert public["sequence_vdj_aa"] == internal["sequence_aa"]
    for field in OUTPUT_SCHEMA:
        if field.endswith("_start") and internal[field] is not None:
            assert public[field] == internal[field] + 1
            assert public[field[:-6] + "_end"] == internal[field[:-6] + "_end"]

def test_unassigned_public_row_has_no_vdj_or_coordinates():
    public = to_public_row(unassigned_internal_row())  # annotation_status="unassigned", sequence=oriented input
    assert public["sequence"] == public["sequence_input"]
    assert public["sequence_vdj"] is None and public["sequence_vdj_aa"] is None
    assert all(public[f] is None for f in OUTPUT_SCHEMA if f.endswith(("_start", "_end")))

def test_light_chain_public_row_keeps_null_d_coordinates():
    public = to_public_row(light_chain_internal_row())  # d_* fields None
    assert public["sequence_vdj"] and public["d_sequence_start"] is None
```

Build the three row helpers from an existing fixture in `abstar/test_data/` via `annotate_single_sequence(...).to_dict()` where possible; synthetic dicts only for the unassigned case.

- [ ] **Step 2: Run to verify failure**

Run: `python -m pytest abstar/tests/test_airr.py abstar/tests/test_schema.py -q -m "not integration and not e2e and not slow"`
Expected: FAIL (`to_public_row` / `AIRR_REQUIRED_FIELDS` import errors).

- [ ] **Step 3: Implement** schema changes and `to_public_row`. Coordinates: convert every `X_start`/`X_end` pair present in `OUTPUT_SCHEMA` (derive the list from the schema, not a hand list). Sequence fields per Global Constraints; `frame` selection as in the current `to_airr_row` (`v_frame`, else `frame`).

- [ ] **Step 4: Run to verify pass** (same command). Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add abstar/annotation/schema.py abstar/annotation/airr.py abstar/tests/test_airr.py abstar/tests/test_schema.py
git commit -m "Build public output rows in one conversion with VDJ-only fields"
```

---

### Task 2: Convert once in the workers; writers and returns only read

**Files:**
- Modify: `abstar/annotation/annotator.py:263-272` (`annotate()` serializes public rows)
- Modify: `abstar/annotation/airr.py` (`write_airr_tsv` encodes rows from the frame; no conversion)
- Modify: `abstar/core/abstar.py:180-210` (`_write_sample_outputs` writes `public_df` directly; delete the per-row loop and now-unused imports)
- Test: `abstar/tests/test_output_parity.py`, `abstar/tests/test_annotator.py`, then every suite the full run flags

**Interfaces:**
- Consumes: `to_public_row`, `OUTPUT_SCHEMA` order from Task 1.
- Produces: work Parquets (and `ANNOTATION_WORK_SCHEMA = {"row_id": ..., **OUTPUT_SCHEMA}`) hold public rows. `run()` dataframe/`Sequence` returns, final Parquet, and TSV are all `work_df.drop("row_id")`. `Sequence` returns keep `seq_key="sequence"` (now the input sequence).

- [ ] **Step 1: Write failing tests** in `abstar/tests/test_output_parity.py`:
  - Replace `test_no_project_api_keeps_internal_sequence_meanings_and_return_shapes` with `test_no_project_api_returns_public_rows`: same parametrization; assert `row["sequence"] == case.sequence`, `row["sequence_vdj"] == case.sequence[s:e]` (fixture 0-based `s`, `e`), `row["v_sequence_start"] == s + 1`, `row["j_sequence_end"] == e`, `len(row["gene_segment_mask"]) == len(row["sequence_vdj"])`, `len(row["gene_segment_mask_aa"]) == len(row["sequence_vdj_aa"])`.
  - Add `test_dataframe_return_equals_final_parquet`: run the same cases once with `as_dataframe=True` and once to a `project_path` with `output_format="parquet"`; `assert_frame_equal` on the two frames (same columns, order, dtypes, values).
  - Simplify the TSV↔Parquet normalizers so they convert representation only (`T`/`F`, empty→null, int parsing); the existing parity tests must then pass with no coordinate arithmetic.
  - Add `test_reverse_complement_public_sequence_is_original_input` using a reverse-complemented fixture case: `sequence == input as given`, `sequence_vdj` equals the oriented slice.
  - In `abstar/tests/test_annotator.py` add `test_public_conversion_failure_is_isolated`: monkeypatch `abstar.annotation.annotator.to_public_row` to raise `ValueError` for one `sequence_id`; `annotate()` returns one `RecordFailure` (`stage="annotation"`, `category="internal_error"`) and the output Parquet contains every other record.

- [ ] **Step 2: Run to verify failure**

Run: `python -m pytest abstar/tests/test_output_parity.py abstar/tests/test_annotator.py -q`
Expected: FAIL.

- [ ] **Step 3: Implement.** In `annotate()`, call `to_public_row(ab.to_dict())` inside the existing per-record `try` for assigned records (so a conversion error follows the existing failure path) and for unassigned records directly; build the work frame from those rows plus `row_id`. Remove the conversion from `write_airr_tsv` and `_write_sample_outputs`.

- [ ] **Step 4: Run focused tests to verify pass** (same command). Expected: PASS.

- [ ] **Step 5: Run the full suite and update public-surface assertions**

Run: `python -m pytest -q`
Tests that read `abstar.run()` returns, `annotate()` work Parquets, or final Parquet/TSV and compare coordinates against 0-based fixture values now fail (expected in `test_real_bcr`, `test_tcr_e2e`, `test_public_e2e`, `test_pipeline`, the boundary/endpoint suites, `test_umi`). Fix each by asserting the public form (`start + 1`, `sequence_vdj` where VDJ was meant). Never change fixture JSON coordinates or assertions on `annotate_single_sequence` / `Antibody` results — those stay internal. Expected after fixes: all pass.

- [ ] **Step 6: Commit**

```bash
git add abstar/annotation/annotator.py abstar/annotation/airr.py abstar/core/abstar.py abstar/tests
git commit -m "Emit identical public rows from every output surface"
```

---

### Task 3: Parquet-reading scripts use public coordinates

**Files:**
- Modify: `scripts/audit_annotation_consistency.py` (reference slice becomes `sequence_oriented[v_sequence_start - 1:j_sequence_end]`; region checks use `start - 1`; adjacency check `left_end + 1 == right_start` → mismatch when `left_end != right_start - 1`; update the report `reference` string)
- Modify: `scripts/build_ci_corpus.py:119` (comment only: it reads pre-0.9 0-based source Parquets; threshold unchanged — selection is frozen)
- Test: `abstar/tests/test_audit_annotation_consistency.py`

**Interfaces:**
- Consumes: final Parquet in the Task 2 public form.
- Produces: audit report schema unchanged except the `reference` string.

- [ ] **Step 1: Update the audit test fixtures to 1-based starts and add** `test_audit_clean_on_current_run_output`: annotate two fixture sequences with `abstar.run(..., output_format="parquet")`, audit the output, assert zero inconsistencies.
- [ ] **Step 2: Run** `python -m pytest abstar/tests/test_audit_annotation_consistency.py -q`. Expected: FAIL.
- [ ] **Step 3: Implement** the audit changes listed above.
- [ ] **Step 4: Run** the same command. Expected: PASS. Also run `grep -n "_start\|_end" scripts/compare_annotation_runs.py scripts/discover_bcr_cases.py scripts/run_corpus.py` and confirm none slice sequences by coordinates (if any do, apply the same `- 1`).
- [ ] **Step 5: Commit** `scripts/` and the test: `git commit -m "Audit public Parquet with 1-based coordinates"`.

---

### Task 4: Regenerate and review the corpus baseline

**Files:**
- Modify: `test_data/bcr_corpus/baseline.parquet`, `test_data/bcr_corpus/baseline-metadata.json`, `test_data/bcr_corpus/README.md` (168 → 170 public fields)

- [ ] **Step 1: Record a proposal** per `test_data/bcr_corpus/README.md` ("Review a proposed baseline or cohort change"): `python scripts/run_corpus.py --output <scratchpad>/corpus-proposal --record-baseline`. Expected: 49,999 annotated, the one documented expected failure, zero audit inconsistencies.
- [ ] **Step 2: Verify the only differences** with a throwaway script in the scratchpad (not committed), joining old and new baselines on `corpus_ordinal`:
  - new columns are exactly `{"sequence_vdj", "sequence_vdj_aa"}`; no columns removed;
  - every `*_start` equals old + 1 (nulls stay null); every `*_end` and every other shared column is identical;
  - for annotated rows, `sequence_vdj == sequence_oriented[v_sequence_start - 1:j_sequence_end]`.
  Expected: all three hold for all 50,000 rows. Any other difference is a stop-and-report, not a baseline update.
- [ ] **Step 3: Copy** the accepted baseline and metadata into `test_data/bcr_corpus/`, update the README field count, and run `python scripts/run_corpus.py --output <scratchpad>/corpus-check`. Expected: exact baseline match.
- [ ] **Step 4: Commit** `git commit -m "Refresh corpus baseline for unified public output"`.

---

### Task 5: Document the breaking change

**Files:**
- Modify: `docs/source/output_formats.rst` (coordinate section near line 20, field meanings near 61 and 104-110, field table near 450-456, native-Parquet sections near 514-571, region semantics near 661)
- Create: `docs/source/changelog.rst`; Modify: `docs/source/index.rst` (toctree entry)
- Modify: `AGENTS.md` (line ~161 reference slice; any 0-based public-output statements), `docs/EDGE-CASES.md` (any 0-based public-output statements), `README.md` (link to changelog if it documents outputs)

- [ ] **Step 1: Write the docs.** `output_formats.rst` gets a prominent `.. warning::` "Breaking change in 0.9.0" at the top stating: all public outputs (Python returns, Parquet, TSV) now share field names, order, and meanings; Parquet/dataframe coordinates are 1-based closed (starts +1, ends unchanged); dataframe and `Sequence` returns' `sequence`/`sequence_aa` are now the input sequence and its translation (use `sequence_vdj`/`sequence_vdj_aa` for the V(D)J); masks align to `sequence_vdj`; migration for old Parquet: `start - 1` restores 0-based half-open. Include one before/after example record (a few fields). Remove every "internal meanings in Python returns" and "Parquet stays 0-based" statement. `changelog.rst` opens with "0.9.0 (unreleased)" — Breaking changes (same bullets, linking the warning) and Added (`sequence_vdj`, `sequence_vdj_aa`).
- [ ] **Step 2: Verify** `grep -rn -i "0-based\|zero-based\|half-open" docs/source AGENTS.md docs/EDGE-CASES.md README.md` shows only internal (`Antibody`) or migration mentions; `python -m sphinx -W --keep-going -b html docs/source docs/_build/html` passes.
- [ ] **Step 3: Full gate**

```bash
python -m pytest --cov=abstar --cov-branch --cov-report=term-missing --cov-report=json:<scratchpad>/cov.json -q
python scripts/check_coverage.py <scratchpad>/cov.json coverage-floors.json
```
Expected: all pass; floors met (raise floors if measured coverage rose).
- [ ] **Step 4: Commit** `git commit -m "docs: document unified public output breaking change"`.
