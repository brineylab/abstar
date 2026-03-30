# Technical Analysis of `abstar`

## Scope and Method

This analysis is based on:

- Source review of the core pipeline, CLI, annotation modules, preprocessing, and germline DB builder.
- Test review across the full `tests/` tree.
- CI workflow review in `.github/workflows/pytest.yml`.
- Local validation runs in the current environment:
  - `pytest -q` completed successfully with `234 passed, 2 skipped` in `144.71s`.
  - Bundled `abutils` binaries were present locally for `mmseqs`, `fastp`, and `vsearch`.
  - A real `fastp` merge succeeded through `merge_fastqs_fastp()`.
  - Several API edge cases were reproduced manually and are called out below.

## Executive Summary

`abstar` is in a better state than many scientific Python repositories in two important ways:

- The annotation core is genuinely exercised end-to-end, not just unit-tested in isolation.
- The current suite is still CI-viable on CPU-only Linux runners: the full suite completed locally in about 2 minutes 24 seconds.

The main issues are not a lack of tests in general, but uneven test coverage around the highest-risk operational paths:

- The V/D/J assignment and annotation path is covered reasonably well.
- Preprocessing and custom germline DB creation are not.
- The public API is convenient for simple notebook use, but it is brittle for production pipeline use because `run()` is highly polymorphic, side-effectful, and contains several correctness bugs in edge cases.

The most important conclusion is this:

- The current codebase already supports a strong CPU-only CI test harness, but it is not yet testing all of the external-tool paths that matter in real-world use.
- Several important extension and preprocessing paths are currently broken or under-specified enough that future flexibility will be impaired unless they are fixed and covered now.

## 1. Test Harness Assessment

### Current strengths

- The repository has substantial breadth: `236` collected tests across annotation, masking, positions, mutations, productivity, MMseqs assignment, pipeline behavior, CLI behavior, UMI parsing, and integration runs.
- The suite is not purely synthetic. It uses real antibody sequences repeatedly, especially the 10E8 heavy chain and two HIV bnAb datasets in `tests/test_integration.py`.
- The core MMseqs2 assignment path is exercised for real, not mocked, in `tests/test_mmseqs.py:286-357`.
- The full pipeline is exercised through the public `run()` API in `tests/test_pipeline.py:82-357`.
- The CLI is exercised through `CliRunner` in `tests/test_cli.py:80-175`.
- Integration tests are marked as `slow` and `integration`, which is the right idea for suite stratification (`pyproject.toml:49-54`).

### Where the current suite is strong enough for CI

- The existing suite already fits a CPU-only CI model.
- The current GitHub Actions workflow runs `pytest` across `2 x 3` Linux jobs and is still reasonable in wall-clock terms (`.github/workflows/pytest.yml:28-47`).
- The integration tests use small enough datasets to remain practical.

### Major test gaps

#### 1. FASTQ merging is not actually tested

- `tests/test_merging.py` explicitly says it only tests filename parsing and grouping and does **not** test real merging (`tests/test_merging.py:5-10`).
- This is a major gap because preprocessing is a real user-facing workflow, especially for the CLI.
- The bundled `fastp` binary is available locally and a real merge succeeded during validation, so these tests are practical on CPU-only CI.

What should be added:

- A true binary-backed `fastp` test through `merge_fastqs_fastp()`.
- A true public-API test through `merge_fastqs()`.
- An interleaved FASTQ test.
- A `run(..., merge=True)` CLI/API e2e test with tiny paired FASTQs.

#### 2. Custom germline database creation is almost completely untested

- The only germline-builder CLI test is `--help` (`tests/test_cli.py:65-72`).
- There are no tests for:
  - FASTA ingestion
  - JSON ingestion
  - IMGT gap addition
  - manifest transfer
  - MMseqs DB creation
  - overwrite behavior
  - addon DB resolution

This is not just a coverage hole; it is already masking real bugs described later.

#### 3. TCR support lacks meaningful end-to-end coverage

- TCR is covered only at the database-access level (`tests/test_germline.py:60-80`) and a CLI smoke test that uses a BCR sequence with `--receptor tcr` and only checks exit status (`tests/test_cli.py:160-175`).
- There is no true TCR sequence through MMseqs assignment + annotation + output validation.

Given the stated product goal of supporting both antibody and TCR annotation, this is a material gap.

#### 4. Negative-path and operational robustness tests are thin

The current suite is mostly happy-path. It does not cover several failure modes that matter in real notebooks and pipelines:

- Non-assignable sequences.
- Empty J-query / C-query FASTA intermediates.
- Generator inputs.
- Directory inputs with duplicate basenames.
- `run(..., merge=True)` via the Python API.
- `copy_inputs_to_project`.
- Output/log cleanup behavior.
- Real unassigned-sequence log contents.
- Mixed valid/invalid input batches.

This matters because some of these paths are already broken in practice.

#### 5. “Slow/integration” markers exist, but CI does not use them

- The project defines `slow` and `integration` markers (`pyproject.toml:49-54`).
- The CI workflow still runs plain `pytest` with no marker selection (`.github/workflows/pytest.yml:46-47`).

That is acceptable today because the suite is still small enough, but it will become a problem as soon as preprocessing/custom-DB/TCR integration coverage is added.

Recommended CI layout:

- `unit`: pure Python logic, run on the full OS/Python matrix.
- `binary`: real `mmseqs` / `fastp` / DB-build smoke tests, run on the full Linux matrix or one Linux matrix.
- `integration/e2e`: real datasets and CLI end-to-end tests, run on one Linux/Python target.

#### 6. Regression sensitivity is weaker than it looks

The integration tests are useful, but many assertions are generic:

- “has V gene”
- “is heavy chain”
- “row count > 0”

That is good for smoke coverage, but weaker for regression detection than a small golden dataset with exact expected outputs.

Recommended addition:

- A tiny golden-set test of 3-5 sequences asserting exact AIRR/parquet outputs for core fields.

## 2. API Structure and Usability

### What is good about the current API

- A single top-level `run()` is easy to discover.
- Returning `Sequence` objects is convenient for notebooks.
- The `gl` / `pp` / `tl` namespaces are the right general idea for interactive use.
- The API already supports several realistic entry modes:
  - single sequence
  - iterable of sequences
  - FASTA/Q file
  - directory of FASTA/Q files

### Main red flags

#### High: `run()` is doing too many jobs at once

`src/abstar/core/abstar.py` mixes:

- input normalization
- temp/project directory management
- logging setup
- FASTQ merging
- MMseqs assignment
- multiprocessing
- output assembly
- return-type switching

This makes the function convenient in the short term, but difficult to extend safely. The current bug profile is consistent with that architecture.

The biggest usability issue is the return contract:

- `Sequence`
- `list[Sequence]`
- `polars.DataFrame`
- `None`

all come back from the same function depending on the input shape and flags (`src/abstar/core/abstar.py:217-224`, `382-445`).

That is easy for ad hoc notebook use, but awkward for typed pipeline code and very easy to misuse.

#### High: interactive temp directories leak

- In API mode, `run()` creates a `TemporaryDirectory` but immediately discards the object and keeps only `.name` (`src/abstar/core/abstar.py:220-224`).
- That means the automatic cleanup handle is lost.
- During local validation, a successful `run(seq)` call left new `/tmp/abstar*` directories behind.

This is a real resource leak for repeated notebook use.

#### High: non-assignable sequences fail hard instead of degrading gracefully

- J-gene assignment always runs after the V pass (`src/abstar/assigners/mmseqs.py:188-211`).
- Unlike D-gene assignment, there is no guard for an empty J-query FASTA.
- During validation, a short non-antibody sequence raised a `RuntimeError` from MMseqs because the generated `.jquery` FASTA had no entries.

For an interactive API intended for notebooks and pipelines, this is the wrong failure mode. A sequence with no recognizable rearrangement should be reported as unassigned, not crash the whole run.

#### High: generator inputs are broken

- `_process_inputs()` checks iterables using `all(isinstance(s, Sequence) for s in sequences)` and then iterates over `sequences` again to write FASTA (`src/abstar/core/abstar.py:483-491`).
- That exhausts generators.
- This was reproduced locally: a generator of `Sequence` objects resulted in an empty temp FASTA.

This is particularly relevant for larger pipelines, where generator-style streaming inputs are common.

#### High: `as_dataframe=True` loses earlier files for multi-file inputs

- In return-by-dataframe mode, `sequence_df` is overwritten inside the per-file loop (`src/abstar/core/abstar.py:383-386`).
- The function returns the last `sequence_df` only (`src/abstar/core/abstar.py:439-441`).

So for a directory input with multiple files, the current code returns only the final sample’s DataFrame instead of the aggregate result. This is a correctness bug in the public API.

#### High: custom germline DB support is currently broken

There are two separate bugs here.

1. `build_germline_database()` calls `get_database_directory()` with swapped arguments (`src/abstar/core/germline_builder.py:163`).

- Default usage reproduced locally as:
  - `AttributeError: 'NoneType' object has no attribute 'lower'`

2. addon DB lookup in `get_germline_database_path()` returns the wrong path:

- It constructs `addon_dir / receptor / germdb_name` even though `addon_dir` already includes the receptor (`src/abstar/annotation/germline_alignment.py:77-80`).
- That makes the documented custom DB override path incorrect.

For a tool that explicitly advertises custom germline databases, this is a serious extensibility problem.

#### High: the Python API and CLI defaults have drifted

Example:

- The CLI normalizes `merge_kwargs` via a callback that returns `{}` when missing (`src/abstar/utils/callbacks.py:12-47`).
- The Python API leaves `merge_kwargs=None`, but `run()` unpacks it with `**merge_kwargs` (`src/abstar/core/abstar.py:271-290`).

That means the public Python API can fail on `merge=True` even though the CLI path is safe.

This is a classic sign that CLI and library concerns are too tightly coupled in one implementation.

#### Medium: preprocessing command execution is brittle and partially broken

The preprocessing/builder layer uses `shell=True` in multiple places:

- `merge_fastqs_vsearch()` (`src/abstar/preprocess/merging.py:542-560`)
- `merge_fastqs_fastp()` (`src/abstar/preprocess/merging.py:681-741`)
- `_make_mmseqs_db()` (`src/abstar/core/germline_builder.py:639-652`)

Problems:

- It is more fragile with paths/quoting.
- It increases command-injection risk if user-supplied values reach CLI args.
- It is harder to test precisely.

There is also a concrete bug in `merge_fastqs_vsearch()`:

- `--fastaout {merged_file}` and `--fastqout {merged_file}` are missing `f`-string interpolation (`src/abstar/preprocess/merging.py:543-546`).
- That function is effectively broken and currently untested.

`MergeGroup.merge()` only accepts `fastp` (`src/abstar/preprocess/merging.py:184-188`), so `vsearch` is also dead/stranded code from the public workflow perspective.

#### Medium: the data model is too large and duplicated

- `Antibody` is a very large mutable dataclass with well over 100 fields (`src/abstar/annotation/antibody.py:16-178`).
- The output schema is separately hand-maintained in `schema.py` (`src/abstar/annotation/schema.py:22-177`).

This creates multiple maintainability problems:

- adding/removing fields requires synchronized changes in multiple places
- the real source of truth is unclear
- accidental output drift is easy

There is no automated consistency check between the dataclass and `OUTPUT_SCHEMA`.

#### Medium: field typing is inconsistent over an object’s lifecycle

- `Antibody.productivity_issues` starts as a `list` (`src/abstar/annotation/antibody.py:76`)
- `assess_productivity()` converts it into a pipe-delimited string (`src/abstar/annotation/productivity.py:60-63`)

That kind of type mutation is a usability smell. It makes downstream code harder to reason about and complicates serialization logic.

#### Medium: performance defaults are not well matched to notebook use

- `run()` defaults to `mp.cpu_count()` processes (`src/abstar/core/abstar.py:297-299`).
- It always uses a `ProcessPoolExecutor` for annotation (`src/abstar/core/abstar.py:367-379`), even for tiny jobs.

That is okay for batch CLI use, but it is not ideal for the interactive notebook use case described in the prompt.

Related issue:

- `MMseqs.__call__()` mutates `self.threads` when the first sample is small (`src/abstar/assigners/mmseqs.py:97-98`).
- Because one assigner instance is reused across files, a small first sample can pin later larger samples to a single MMseqs thread.

#### Medium: directory recursion can create output-name collisions

- `_process_inputs()` recursively lists files from nested directories (`src/abstar/core/abstar.py:470-480`).
- Output naming uses only the basename-derived sample name (`src/abstar/core/abstar.py:326-329`, `405-409`).

If two different subdirectories contain `sample1.fasta`, their outputs/logs will collide. This is a realistic operational problem for project-organized sequencing runs.

#### Medium: package import has global side effects

- Importing `abstar` globally suppresses `BiopythonWarning` and all `FutureWarning`s (`src/abstar/__init__.py:3-8`).

Suppressing all `FutureWarning`s at package import is not maintainable. It can hide real deprecations in upstream libraries that would otherwise surface during development.

## 3. Documentation and Usability Drift

There is noticeable drift between the real API and the user-facing docs/examples.

### README drift

- `pyproject.toml` requires Python `>=3.11` (`pyproject.toml:10`)
- `README.md` still says `Python 3.10+` (`README.md:59-60`)

### Namespace module examples are stale or incorrect

- `abstar.pp` docstring shows `merge_fastqs(forward, reverse, merged_output)` even though the exported public function signature is `merge_fastqs(files, output_directory, ...)` (`src/abstar/pp.py:17-21`).
- `abstar.tl` docstring shows `parse_umis(..., umi_pattern=...)`, but the real parameter is `pattern` (`src/abstar/tl.py:14-17`; real function in `annotation/umi.py`).
- `abstar.gl` docstring shows `build_germline_database(..., fasta_files=[...])`, but the real argument is `fastas` (`src/abstar/gl.py:22-23`; real function in `core/germline_builder.py`).

For notebook users who rely on docstrings and tab-completion, this matters more than it would in a purely CLI-first tool.

## 4. Other Relevant Observations

### The current suite is already a good base to build on

This is not a “rewrite the tests from scratch” situation.

The right next step is:

- keep the current annotation/MMseqs coverage
- add missing binary-backed preprocessing/custom-DB tests
- add a small number of targeted regression tests for the bugs identified here

### `abutils` coupling is deep

This codebase depends heavily on `abutils` for:

- sequence types
- FASTA/Q IO
- logging
- alignments
- MMseqs search wrappers
- packaged binary resolution

That is fine architecturally if intentional, but it means:

- compatibility with `abutils` changes should be treated as a first-class concern
- binary-backed integration tests are essential, not optional

## 5. Prioritized Recommendations

### Immediate fixes

1. Fix `build_germline_database()` argument ordering and addon DB path resolution.
2. Fix `run()` API bugs:
   - temp directory leak
   - generator exhaustion
   - multi-file DataFrame overwrite
   - `merge_kwargs=None` unpacking
3. Guard MMseqs J/C searches against empty query FASTAs and treat those cases as unassigned rather than exceptions.
4. Replace `pl.count()` with `pl.len()` (`src/abstar/core/abstar.py:419`).
5. Remove or fix dead/broken `vsearch` code.

### Test-harness additions

1. Add real binary-backed tests for:
   - `merge_fastqs_fastp()`
   - `merge_fastqs()` via the public preprocessing API
   - `_make_mmseqs_db()` / custom DB build path on tiny germline inputs
2. Add true TCR end-to-end tests with real TCR sequences.
3. Add regression tests for:
   - generator inputs
   - unassignable sequences
   - multi-file DataFrame returns
   - temp cleanup
   - duplicate-basename directory inputs
4. Add small golden-output tests for AIRR/parquet content, not just existence.

### Structural refactor

1. Introduce a typed `RunConfig` object instead of continuing to grow the `run()` parameter list.
2. Split `run()` into explicit entry points or explicit phases:
   - normalize inputs
   - preprocess
   - assign
   - annotate
   - emit outputs
3. Return a stable result type for the Python API instead of switching between `Sequence`, `list`, `DataFrame`, and `None`.
4. Make the output model a single source of truth instead of maintaining `Antibody` fields and `OUTPUT_SCHEMA` separately.

## Bottom Line

The repository already has a respectable test foundation and a real end-to-end core. The most important technical problem is not lack of ambition in the test suite; it is that coverage falls off sharply exactly where the code becomes most operationally complex: preprocessing, custom database creation, negative-path robustness, and Python API edge cases.

If the goal is a durable CLI for larger jobs and a reliable notebook/pipeline API for interactive use, the next round of work should focus on:

- fixing the broken extensibility paths
- adding real binary-backed preprocessing/custom-DB tests
- simplifying and stabilizing the public API surface

That will materially improve both immediate reliability and long-term maintainability.
