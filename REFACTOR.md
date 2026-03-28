# abstar Refactoring Plan

Comprehensive analysis and refactoring plan for the `abstar` codebase. Findings are organized
by category, with issues ranked by severity within each section.

---

## Table of Contents

1. [Repository Structure & Packaging](#1-repository-structure--packaging)
2. [Dead Code & Legacy Artifacts](#2-dead-code--legacy-artifacts)
3. [Testing Framework](#3-testing-framework)
4. [Code Organization & Naming](#4-code-organization--naming)
5. [Syntax Modernization](#5-syntax-modernization-python-311)
6. [Implementation Plan](#6-implementation-plan)

---

## 1. Repository Structure & Packaging

### 1.1 Migrate to `src/` Layout (HIGH)

**Current state:** Flat layout with `abstar/` at repository root.

**Target state:**
```
abstar/
├── src/
│   └── abstar/
│       ├── __init__.py
│       ├── version.py
│       ├── gl.py
│       ├── pp.py
│       ├── tl.py
│       ├── annotation/
│       ├── assigners/
│       ├── core/
│       ├── preprocess/
│       ├── scripts/
│       ├── utils/
│       └── germline_dbs/
├── tests/                      # move out of package
├── pyproject.toml
├── .github/workflows/
└── docs/
```

**Plan:**
1. Create `src/abstar/` and move all package code from `abstar/` into it.
2. Move `abstar/tests/` to top-level `tests/`.
3. Move `abstar/test_data/` to `tests/data/` (or `tests/test_data/`).
4. Update `pyproject.toml` package discovery (see 1.2).
5. Update all relative imports if needed (src layout shouldn't change internal imports).
6. Update CI workflows and pytest configuration.

### 1.2 Switch Build Backend to Hatchling (HIGH)

**Current state:** `setuptools` with dynamic version/dependency loading from separate files.

**Current `pyproject.toml`:**
```toml
[build-system]
requires = ["setuptools>=68", "wheel"]
build-backend = "setuptools.build_meta"
```

**Target state:**
```toml
[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"

[project]
name = "abstar"
version = "0.8.0"  # or use hatch-vcs for git-tag-based versioning
requires-python = ">=3.11"
dependencies = [
    "abutils>=0.5.1",
    "click",
    "matplotlib",
    "numpy",
    "pandas",
    "parasail",
    "polars",
    "pyarrow",
]

[project.optional-dependencies]
dev = ["pytest", "ruff"]

[tool.hatch.build.targets.wheel]
packages = ["src/abstar"]
```

**Plan:**
1. Replace `setuptools` build backend with `hatchling`.
2. Inline dependencies directly in `pyproject.toml` (eliminate `requirements.txt` for runtime deps).
3. Move `pytest` to `[project.optional-dependencies] dev`.
4. Bump `requires-python` to `>=3.11`.
5. Remove `MANIFEST.in` — hatchling handles package data via config.
6. Remove `abstar/version.py` if using `hatch-vcs`, or keep with hatchling's `version` field.
7. Configure package data inclusion for `germline_dbs/` and related data files.

### 1.3 Update CI Workflows (HIGH)

#### pytest.yml
- Update Python version matrix: drop 3.10, keep 3.11-3.13.
- Install with `pip install -e ".[dev]"` instead of separate `requirements.txt` + `pytest` installs.
- Add `ruff check` and `ruff format --check` steps.
- Add pytest configuration in `pyproject.toml` (markers, coverage, etc.).

#### python-publish.yml
- Update to use hatchling: `pip install build` remains valid (`python -m build` works with any
  PEP 517 backend).
- Remove commented-out legacy `twine`/`setuptools` workflow (lines 39-61).

### 1.4 Remove Legacy CI and Packaging Artifacts (MEDIUM)

| File/Directory | Status | Action |
|---|---|---|
| `.travis.yml` | Obsolete (GitHub Actions in use) | Delete |
| `bin/` (5 scripts) | Legacy scripts, not wired into packaging | Delete |
| `requirements.txt` | Replaced by inline deps in `pyproject.toml` | Delete |
| `MANIFEST.in` | Not needed with hatchling | Delete |
| `python-publish.yml` lines 39-61 | Commented-out old workflow | Remove |

### 1.5 Add Tool Configuration to `pyproject.toml` (MEDIUM)

```toml
[tool.pytest.ini_options]
testpaths = ["tests"]
markers = [
    "slow: marks tests as slow (deselect with '-m \"not slow\"')",
    "integration: marks integration tests",
]

[tool.ruff]
line-length = 100
target-version = "py311"

[tool.ruff.lint]
select = ["E", "F", "I", "UP"]

[tool.ruff.lint.isort]
known-first-party = ["abstar"]
```

---

## 2. Dead Code & Legacy Artifacts

### 2.1 Files to Delete (HIGH)

#### `abstar/assigners/blastn.py` (679 lines)
- **Why:** Imports non-existent modules (`..core.vdj.VDJ`, `..core.germline.GermlineSegment`).
  Cannot be imported or executed. Fully replaced by `mmseqs.py`.
- **Action:** Delete entirely. Retained in git history.

#### `abstar/utils/build_germline_dbs_OLD.py` (424 lines)
- **Why:** Explicitly marked as old (`_OLD` suffix). Not imported anywhere. Replaced by
  `core/germline.py`.
- **Action:** Delete entirely.

### 2.2 Large Commented-Out Code Blocks (HIGH)

These are all substantial blocks of legacy code preserved as comments. They should be removed
(git history serves as the archive).

| File | Lines | Size | Description |
|---|---|---|---|
| `assigners/assigner.py` | 68-442 | ~375 lines | Old `BaseAssigner` class with extensive documentation |
| `annotation/regions.py` | ~208-308 | ~100 lines | Previous `get_region_sequence()` implementation |
| `scripts/abstar.py` | ~304-377 | ~73 lines | Old `build_germdb_from_igdiscover()` CLI command |
| `annotation/umi.py` | ~601-644 | ~44 lines | Old `UMI` class implementation |
| `utils/mixins.py` | ~121-217 | ~97 lines | Old `LoggingMixin` implementation |
| `tests/test_integration.py` | 1-47 | 47 lines | Entirely commented-out integration tests |

**Total: ~736 lines of commented-out code to remove.**

### 2.3 Commented-Out Imports (MEDIUM)

| File | Line(s) | Import |
|---|---|---|
| `annotation/annotator.py` | 11 | `import pandas as pd` |
| `annotation/antibody.py` | 9 | `from ..utils.mixins import LoggingMixin` |
| `core/germline.py` | 35, 40-42 | `weakref`, `click`, `natsort`, `networkx` |
| `scripts/abstar.py` | 12-14 | `build_germdb_from_igdiscover` |
| `tests/__init__.py` | 7-11 | Python 2/3 pickle compatibility code |

**Action:** Remove all commented imports.

### 2.4 Scattered Commented-Out Code in Active Files (MEDIUM)

| File | Description |
|---|---|
| `annotation/positions.py` | Line 137: commented variable assignment |
| `annotation/regions.py` | Lines 117-142: old position calculation code |
| `annotation/umi.py` | Lines 134-138, 373-438, 492, 545: `extra_length_for_alignment` parameters |
| `annotation/annotator.py` | Lines 847, 934, 946: commented function parameters |

**Action:** Remove all. If any represent TODOs, convert to proper `# TODO:` comments with context.

### 2.5 Misleading `__all__` in `annotation/germline.py` (MEDIUM)

Functions used by `annotator.py` are commented out of `__all__`:
- `realign_germline`
- `reassign_dgene`
- `process_vgene_alignment`
- `process_jgene_alignment`
- `process_dgene_alignment`

**Action:** Uncomment these entries in `__all__` or remove `__all__` and rely on explicit imports.

### 2.6 Unused Import in Active Code (LOW)

- `annotation/regions.py` line 11: `get_position_from_aligned_reference` imported but only used
  in commented-out code.

**Action:** Remove from import statement.

---

## 3. Testing Framework

### 3.1 Test Suite Summary

| Metric | Value |
|---|---|
| Test files | 14 (1 entirely commented out) |
| Active test cases | ~168 |
| Total fixtures | ~69 |
| Tests marked `xfail` | 9 (5.4%) |
| Tests marked `skip` | 2 (1.2%) |
| Modules with zero test coverage | 11 |
| Pytest configuration | None |

### 3.2 Known Bugs Hiding Behind `xfail` (HIGH)

9 tests are marked `xfail`, indicating known bugs rather than expected failures:

**Polars schema errors (5 tests):**
- `test_as_dataframe.py::test_dataframe_row_count_multiple_sequences` — Polars schema mismatch
  when D gene assignment is missing
- `test_as_dataframe.py::test_multiple_sequences_as_dataframe` — Same root cause
- `test_as_dataframe.py::test_multiple_sequences_default_returns_list` — Same root cause
- `test_mmseqs.py::test_mmseqs_multiple_sequences` — Same root cause
- `test_pipeline.py::test_run_returns_sequence_list` — Same root cause

**Root cause:** When a sequence lacks a D gene assignment (e.g., light chains), the schema
produced by MMseqs assignment has a missing/null column that causes Polars `concat` or
`read_parquet` to fail with a schema mismatch.

**Germline retrieval errors (4 tests):**
- `test_germline.py::test_get_germline_database_name_invalid` — expects `FileNotFoundError`
- `test_germline.py::test_get_germline_database_receptor_invalid` — expects `ValueError`
- `test_germline.py::test_get_single_germline_nonexistent` — expects `ValueError`
- `test_germline.py::test_get_single_germline_nonunique` — expects `ValueError`

**Action:**
1. Fix the Polars schema bug (ensure consistent schema regardless of D gene presence).
2. Investigate the germline `xfail` tests — they may be testing for error behavior that isn't
   implemented yet. Convert to proper tests once behavior is defined.

### 3.3 Modules with Zero Test Coverage (HIGH)

Listed by criticality:

| Module | Why It Matters | Priority |
|---|---|---|
| `annotation/antibody.py` | Core data model (170+ fields); used everywhere | **Critical** |
| `core/germline.py` | Germline database building; user-facing API via `abstar.gl` | **High** |
| `scripts/abstar.py` | CLI entry point; user-facing | **High** |
| `preprocess/merging.py` | Paired-end read merging; user-facing via `abstar.pp` | **High** |
| `assigners/assigner.py` | Base class for assigners | Medium |
| `assigners/blastn.py` | Legacy — delete instead of testing | N/A |
| `utils/callbacks.py` | CLI argument parsing utilities | Medium |
| `utils/mixins.py` | LoggingMixin (may be redundant with `abutils` version) | Low |
| `gl.py`, `pp.py`, `tl.py` | Namespace re-export modules | Low |

**Plan:**
1. Write tests for `annotation/antibody.py`: construction, field defaults, serialization.
2. Write tests for `core/germline.py`: database building with real IMGT data.
3. Write CLI tests for `scripts/abstar.py` using Click's `CliRunner`.
4. Write tests for `preprocess/merging.py` if fastp is available in CI.

### 3.4 Integration Tests Disabled (HIGH)

`test_integration.py` is entirely commented out. It previously tested:
- 1000 BCR sequences (`test_1k.fasta`)
- 221 HIV bnAb heavy chains
- 207 HIV bnAb light chains
- Chunking with different core counts

**Plan:**
1. Rewrite integration tests as proper pytest functions.
2. Mark with `@pytest.mark.slow` and `@pytest.mark.integration`.
3. Use the existing test data files (`test_data/test_1k.fasta`, etc.).
4. Run in CI with a separate job or `pytest -m integration`.

### 3.5 Weak Assertions (MEDIUM)

Many tests only check field existence or non-null rather than validating actual values.

**Examples of weak patterns found:**
```python
# Too weak — only checks existence
assert result["v_call"] is not None
assert "IGHV" in result["v_call"]

# Better — validates specific expected result
assert result["v_call"] == "IGHV4-59*01"
```

**Affected files:** `test_annotator.py` (field-existence checks), `test_as_dataframe.py`
(type/null checks), `test_productivity.py` (issue presence only), `test_mask.py` (string
equality without edge cases), `test_schema.py` (minimal coverage).

**Plan:** For each test file, add value-level assertions using known reference antibody
sequences (e.g., 10E8, VRC01) where expected annotation values are well-characterized.

### 3.6 Redundant/Overlapping Tests (MEDIUM)

1. **`test_as_dataframe.py` overlaps with `test_pipeline.py`**: Both test `run()` with
   DataFrame output, single sequence input, and return types. Consolidate into `test_pipeline.py`.

2. **MMseqs initialization tests**: 5 tests with slight parameter variations could be
   consolidated with `@pytest.mark.parametrize`.

3. **Annotation field existence tests**: `test_annotator.py` and `test_pipeline.py` both verify
   that annotation fields (v_call, j_call, cdr3, etc.) are populated. Deduplicate.

**Plan:** Merge `test_as_dataframe.py` into `test_pipeline.py`. Parametrize initialization tests.

### 3.7 Missing Edge Case Coverage (MEDIUM)

| Area | Missing Edge Cases |
|---|---|
| Mutations | Empty sequences, mismatched lengths, ambiguous bases |
| Productivity | Combined failures, None inputs, ambiguous codons |
| Masks | Gaps, missing regions, partial CDR3, sequences with D genes |
| Schema | Schema violations, missing required fields, type mismatches |
| UMI | Mismatch tolerance, invalid patterns, empty input |

### 3.8 Add Pytest Configuration (MEDIUM)

No pytest configuration exists. Add to `pyproject.toml`:
- `testpaths`
- Custom markers (`slow`, `integration`)
- Coverage configuration
- Default options (`-v`, `--tb=short`)

---

## 4. Code Organization & Naming

### 4.1 Confusing Duplicate Module Names (HIGH)

Two modules named `germline.py` with different purposes:

| Module | Purpose |
|---|---|
| `annotation/germline.py` | Process germline alignments (realign, extract V/D/J/C) |
| `core/germline.py` | Build germline databases from FASTA/JSON |

**Plan:** Rename for clarity:
- `annotation/germline.py` → `annotation/germline_alignment.py`
- `core/germline.py` → `core/germline_builder.py`

Update all imports and the `gl.py` re-export module accordingly.

### 4.2 Cryptic Namespace Aliases (MEDIUM)

The `gl`, `pp`, `tl` aliases follow the scanpy/AnnData convention but are opaque outside that
ecosystem.

| Current | Full Name | Contents |
|---|---|---|
| `abstar.gl` | germline | `get_germline()`, `get_germline_database_path()`, `build_germline_database()` |
| `abstar.pp` | preprocess | `merge_fastqs()`, `group_paired_fastqs()` |
| `abstar.tl` | tools | `parse_umis()`, `build_germline_database()` |

**Issues:**
- `tl` is not intuitive (usually "tools" but could mean many things).
- `build_germline_database` appears in both `gl` and `tl` (via `import *`).
- These modules use `from X import *` without `__all__` in the re-export files.

**Plan:** This is a stylistic choice inherited from the scanpy ecosystem. If the team wants to
keep this convention (likely, given the domain), then at minimum:
1. Add `__all__` to `gl.py`, `pp.py`, and `tl.py` to make exports explicit.
2. Remove `build_germline_database` from `tl.py` (it belongs in `gl`).
3. Add docstrings to each file explaining the namespace.

### 4.3 Missing `__all__` Declarations (MEDIUM)

Most modules lack `__all__`, making the public API ambiguous. Modules missing `__all__`:

- `annotation/__init__.py` (empty)
- `core/__init__.py` (empty)
- `preprocess/__init__.py` (empty)
- `gl.py`, `pp.py`, `tl.py`
- `annotation/mask.py`
- `annotation/mutations.py`
- `annotation/regions.py`
- `annotation/schema.py`
- `annotation/antibody.py`
- `annotation/productivity.py`

**Plan:** Add `__all__` to every module that exports public symbols.

### 4.4 Abbreviated Variable Names (LOW)

Some modules use heavily abbreviated names that reduce readability:

| Abbreviation | Full Name | Where Used |
|---|---|---|
| `ab` | `antibody` | `germline.py`, `mask.py`, `indels.py`, `mutations.py`, `productivity.py`, `regions.py` |
| `aln` | `alignment` | `germline.py`, `annotator.py` |
| `sg` | `semiglobal` | `germline.py` |
| `loc` | `local` | `germline.py` |
| `germ` | `germline` | `germline.py`, `annotator.py` |

**Note:** `ab` for `antibody` is used pervasively (dozens of occurrences across 8+ files). This
is a domain convention and arguably fine for a bioinformatics package. However, function
*signatures* should use the full name for clarity:

```python
# Current
def assess_productivity(ab: Antibody) -> Antibody:

# Better (parameter name matches type)
def assess_productivity(antibody: Antibody) -> Antibody:
```

**Plan:** Rename parameters in function signatures to use full names. Internal variable names
can keep abbreviations where the scope is small and the meaning is obvious.

### 4.5 Empty `__init__.py` Files (LOW)

`annotation/__init__.py`, `core/__init__.py`, `preprocess/__init__.py` are all empty. These
could either:
- Export key symbols for convenient imports, or
- Be left empty (which is fine for namespace packages).

**Plan:** Add exports if there are commonly-used symbols users should access directly. Otherwise
leave as-is.

### 4.6 Fragile Dynamic Loading in `assigners/__init__.py` (LOW)

The assigners `__init__.py` uses glob-based dynamic module loading to populate `__all__`. This
is fragile and makes the import graph hard to trace.

**Plan:** Replace with explicit imports:
```python
from .mmseqs import MMseqs
from .assigner import AssignerBase

__all__ = ["MMseqs", "AssignerBase"]
```

---

## 5. Syntax Modernization (Python 3.11+)

### 5.1 `os.path` → `pathlib.Path` (HIGH — large scope)

Extensive `os.path` usage across the codebase. This is the largest single modernization task.

**Affected files (by approximate count of `os.path` calls):**

| File | `os.path` calls |
|---|---|
| `assigners/mmseqs.py` | ~50+ |
| `preprocess/merging.py` | ~30+ |
| `core/germline.py` | ~25+ |
| `core/abstar.py` | ~20+ |
| `assigners/assigner.py` | 2 |

**Common patterns to replace:**
```python
# Before
os.path.join(path1, path2)      →  Path(path1) / path2
os.path.isfile(f)               →  Path(f).is_file()
os.path.isdir(d)                →  Path(d).is_dir()
os.path.exists(p)               →  Path(p).exists()
os.path.abspath(p)              →  Path(p).resolve()
os.path.basename(p)             →  Path(p).name
os.path.dirname(p)              →  Path(p).parent
os.path.expanduser(p)           →  Path(p).expanduser()
os.makedirs(p, exist_ok=True)   →  Path(p).mkdir(parents=True, exist_ok=True)
```

**Plan:** Convert file-by-file, starting with `core/abstar.py` (the main entry point), then
`assigners/mmseqs.py`, then the rest. Update type hints from `str` to `str | Path` where
paths are accepted.

### 5.2 Add `from __future__ import annotations` (HIGH)

No source files include this import. Adding it enables:
- PEP 604 union syntax (`X | Y`) in all contexts without runtime overhead.
- Forward references without quotes.
- Consistent behavior across Python 3.11+.

**Action:** Add to every `.py` file in `src/abstar/`.

### 5.3 `Tuple[]` → `tuple[]` and Remove `typing` Import (MEDIUM)

`annotation/germline.py` line 33 imports `Tuple` from `typing` and uses `Tuple[...]` in a
return type annotation (line 214).

**Action:** Replace `Tuple[...]` with `tuple[...]` and remove the `from typing import Tuple`.

### 5.4 `Optional[X]` / `Union[X, Y]` in Docstrings (MEDIUM)

Several files use old-style type annotations in docstrings:

| File | Lines | Pattern |
|---|---|---|
| `core/abstar.py` | 124, 135, 151, 155, 158, 165, 179, 460 | `Optional[str]`, `Union[str, Sequence, ...]` |
| `core/germline.py` | 130-136 | `Optional[str]`, `Union[str, Iterable[str]]` |
| `annotation/annotator.py` | 77, 80, 83, 94, 98 | `Optional[str]`, `Optional[int]` |
| `annotation/germline.py` | 274 | `Tuple[Optional[...], Optional[...]]` |
| `utils/callbacks.py` | 23 | `Optional[str]` |

**Action:** Update all docstrings to use modern `X | Y` and `X | None` syntax.

### 5.5 `.format()` → f-strings (LOW)

Remaining `.format()` calls are mostly in legacy files that will be deleted (`blastn.py`,
`build_germline_dbs_OLD.py`). Only one instance found in active code:

- `core/germline.py` line 842 (commented-out code — will be removed with dead code cleanup).

**Status:** Active codebase already uses f-strings extensively (~284+ instances). No action
needed beyond dead code removal.

### 5.6 `== None` → `is None` (LOW)

One instance found in `core/germline.py` line 420 (in commented-out code — will be removed).

**Status:** Active code already uses `is None` correctly. No action needed.

---

## 6. Implementation Plan

### Phase 1: Cleanup (Low Risk, High Value)

Remove dead code and legacy artifacts. These changes are purely subtractive and can't break
anything.

| Task | Files | Est. Lines Removed |
|---|---|---|
| Delete `assigners/blastn.py` | 1 file | 679 |
| Delete `utils/build_germline_dbs_OLD.py` | 1 file | 424 |
| Delete `.travis.yml` | 1 file | 36 |
| Delete `bin/` directory | 5 files | ~40 |
| Remove commented code in `assigners/assigner.py` | 1 file | 375 |
| Remove commented code in `annotation/regions.py` | 1 file | 100 |
| Remove commented code in `scripts/abstar.py` | 1 file | 73 |
| Remove commented code in `annotation/umi.py` | 1 file | 44 |
| Remove commented code in `utils/mixins.py` | 1 file | 97 |
| Remove commented code in `tests/test_integration.py` | 1 file | 47 |
| Remove commented imports across codebase | 5 files | ~20 |
| Remove scattered commented-out code in active files | 4 files | ~30 |
| Remove commented-out workflow in `python-publish.yml` | 1 file | 23 |
| Fix `__all__` in `annotation/germline.py` | 1 file | — |
| **Total** | | **~1,988 lines** |

### Phase 2: Repository Restructure (Medium Risk)

Reorganize to `src/` layout and modernize packaging.

1. Create `src/abstar/` directory and move package code.
2. Move `tests/` to top level.
3. Move `test_data/` into `tests/`.
4. Switch to hatchling build backend.
5. Inline dependencies in `pyproject.toml`.
6. Delete `requirements.txt`, `MANIFEST.in`.
7. Add tool configuration (pytest, ruff) to `pyproject.toml`.
8. Update CI workflows.
9. Bump `requires-python` to `>=3.11`.
10. Update Python version classifiers.

### Phase 3: Code Organization (Medium Risk)

Naming and structural improvements.

1. Rename `annotation/germline.py` → `annotation/germline_alignment.py`.
2. Rename `core/germline.py` → `core/germline_builder.py`.
3. Update all imports referencing renamed modules.
4. Add `__all__` to all public modules.
5. Replace dynamic loading in `assigners/__init__.py` with explicit imports.
6. Add docstrings to `gl.py`, `pp.py`, `tl.py`.
7. Remove `build_germline_database` from `tl.py` (keep only in `gl.py`).

### Phase 4: Syntax Modernization (Low Risk)

Systematic Python 3.11+ updates.

1. Add `from __future__ import annotations` to all source files.
2. Replace `Tuple[...]` → `tuple[...]` and remove `typing.Tuple` import.
3. Update docstring type annotations (`Optional` → `| None`, `Union` → `|`).
4. Convert `os.path` to `pathlib.Path` (file by file, largest task in this phase):
   - `core/abstar.py`
   - `assigners/mmseqs.py`
   - `assigners/assigner.py`
   - `core/germline.py` (now `germline_builder.py`)
   - `preprocess/merging.py`

### Phase 5: Testing (Medium Risk)

Fix and expand test coverage.

1. **Fix the Polars schema bug** causing 5 `xfail` tests to fail.
2. **Investigate and resolve** the 4 germline `xfail` tests.
3. **Write new tests** for untested modules:
   - `annotation/antibody.py` (critical — core data model)
   - `core/germline.py` (database building)
   - `scripts/abstar.py` (CLI via `CliRunner`)
   - `preprocess/merging.py` (if fastp available)
4. **Rewrite integration tests** (from commented-out `test_integration.py`):
   - Mark with `@pytest.mark.slow` / `@pytest.mark.integration`.
   - Use existing test data files.
5. **Strengthen assertions** in existing tests:
   - Replace field-existence checks with value-level validation.
   - Add edge case tests for mutations, productivity, masks, and UMI modules.
6. **Consolidate redundant tests**:
   - Merge `test_as_dataframe.py` into `test_pipeline.py`.
   - Parametrize MMseqs initialization tests.
7. **Add pytest configuration** to `pyproject.toml`.

---

## Appendix: Quick Reference

### Files to Delete
```
.travis.yml
bin/abstar
bin/basespace_downloader
bin/batch_mongoimport
bin/build_abstar_germline_db
bin/make_basespace_credfile
abstar/assigners/blastn.py
abstar/utils/build_germline_dbs_OLD.py
```

### Files to Rename
```
abstar/annotation/germline.py  →  abstar/annotation/germline_alignment.py
abstar/core/germline.py        →  abstar/core/germline_builder.py
```

### Files with Most Commented-Out Code (lines to remove)
```
assigners/assigner.py      ~375 lines
annotation/regions.py      ~100 lines
utils/mixins.py            ~97 lines
scripts/abstar.py          ~73 lines
annotation/umi.py          ~44 lines
tests/test_integration.py  ~47 lines
```

### xfail Tests to Fix
```
test_as_dataframe::test_dataframe_row_count_multiple_sequences
test_as_dataframe::test_multiple_sequences_as_dataframe
test_as_dataframe::test_multiple_sequences_default_returns_list
test_mmseqs::test_mmseqs_multiple_sequences
test_pipeline::test_run_returns_sequence_list
test_germline::test_get_germline_database_name_invalid
test_germline::test_get_germline_database_receptor_invalid
test_germline::test_get_single_germline_nonexistent
test_germline::test_get_single_germline_nonunique
```
