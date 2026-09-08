# Committed BCR Corpus Implementation Plan

> Execute the independent selection and runner tasks with parallel agents;
> integrate and review their public boundary together before completion.

**Goal:** Add a fixed, difficult-case-enriched BCR corpus and a self-contained
push/PR regression job targeting five minutes, with explicit failure diagnosis.

**Approved design:** The conversation approved a committed enriched sample,
approximately 50,000 records subject to measured runtime; one Linux/Python job;
record conservation, independent consistency checks, and strict biological
baseline comparison; reproducibility and baseline-update guidance for agents.

**Architecture:** Root-level `test_data/bcr_corpus/` stores compressed original
FASTA, per-record source metadata, a manifest, and a compressed native Parquet
baseline. `scripts/build_ci_corpus.py` selects from explicit local inputs only.
`scripts/run_corpus.py` annotates into a fresh external workspace, audits output,
and compares per-record results. Existing small biological goldens remain in
ordinary pytest scopes. No corpus execution is added to those scopes.

**Global constraints:** No production annotation changes; no full-corpus rerun;
original source files remain read-only; immutable `(source_file, record_ordinal)`
identities preserve duplicate external IDs; no implicit baseline refresh or
tolerances; expected failures are individual reviewed records; output workspaces
are new and outside source/corpus trees. Keep earlier nightly removal and guide.

## Shared on-disk contract

- `sequences.fasta.gz`: original IDs and sequences in fixed corpus order.
- `records.parquet`: `corpus_ordinal` (integer), `source_file` (string),
  `record_ordinal` (integer), `sequence_id` (string), `sequence_sha256` (string),
  `donor`, `flow_class`, `locus`, `panel` (all strings), and
  `selection_reasons` (list of strings). Source ordinals are zero-based.
- `manifest.json`: `schema_version: 1`, `corpus_version: "bcr-ci-v1"`,
  `record_count`, `files` (filename-to-SHA256), `selection` and `sources`
  provenance, and `parameters` containing receptor `bcr`, database `human`,
  `n_processes: 2`, `mmseqs_threads: 2`, `chunksize: 500`.
- Manifest `expected_failures`: list of objects with `corpus_ordinal`, `stage`,
  `category`, and `exception_type`; initially only independently reviewed
  known failures. New failures block baseline creation and regression checks.
- `baseline.parquet`: every public output field, plus `corpus_ordinal` for
  alignment to source identity. Only successful/unassigned rows live here.
- `baseline-metadata.json`: environment, MMseqs binary/version, reference hashes,
  implementation provenance, comparison policy, and baseline-run parameters.

## Tasks

- [x] Selection: test deterministic stratification, difficulty classification,
  duplicate IDs and exact sequence deduplication, mandatory original cases,
  safe output/provenance; implement explicit-path builder; select ~40k background
  plus ~10k difficult records using existing annotations only as sampling evidence.
- [x] Runner: test missing/duplicate/reordered records, changed biological fields,
  unexpected/expected failures, tampered corpus and invalid workspace, and a tiny
  real public-API run. Implement `--corpus` and required `--output`; default corpus
  is the committed directory. `--record-baseline` writes proposed baseline files
  only into the new output directory after conservation/auditing/expected-failure
  validation. Normal checking never changes corpus inputs. Use existing native
  Parquet auditor and strict field/schema comparison with actionable JSON output.
- [x] Workflow/docs: add fixed dependency profile, push/PR workflow on Python 3.12,
  fixed process/thread settings, report upload on failure/success, and generous
  hard timeout with five-minute target. Test shell behavior against a tiny corpus
  rather than running 50k in pytest. Exclude large fixture data from distributions.
- [x] Baseline: run the actual selected corpus, inspect all failures/consistency
  findings, record baseline, rerun independently with identical cohort/settings,
  and require exact agreement. Record timings and limits of local measurements.
- [x] Review and verify: focused modules, full pytest/coverage floors, warning-free
  docs, built distribution contents, full corpus runner, source-link/whitespace
  checks. Do not claim GitHub timing until a hosted run establishes it.

## Execution record

The current feature branch contains only earlier authorized workflow removal
and debugging-guide changes. Work continues there to preserve that task context.
Local timing probes (two-CPU affinity) measured 9,462 records in 29.93 seconds;
that was a sizing probe, not the final enriched selection or GitHub validation.

The installed selection contains 50,000 unique original sequences, with 40,000
background and 10,000 difficult records across 94 datasets, four donors, six
flow classes, and all 298 observed strata. It includes all 100 matching fixture
candidates (97 distinct sequences) and preserves 758 repeated external IDs.
The corpus, provenance, and compressed all-field baseline total about 43.5 MB.

Review caught and fixed an incomplete fixture inventory, reuse of mismatched
annotation evidence, loss of mandatory reasons during duplicate replacement,
and imported-package/checkout provenance mismatch. Each has a regression test.
No production assignment or annotation module changed.

The first two full runs agreed exactly across all 168 public fields. Each
conserved 49,999 annotated rows and the independently reviewed ambiguous-boundary
failure at corpus ordinal 11457, with zero internal-consistency findings. Local
runner totals were 151.18 and 151.26 seconds; hosted timing remains unverified.
The first run also exercised rejection of that failure before its reviewed
manifest expectation was loaded; no baseline was emitted by the rejected run.

Verification:

- Full statement/branch suite: 1,678 passed in 406.32 seconds, coverage 91.25%.
- Final focused selection/runner/workflow modules: 44 passed in 6.11 seconds,
  including the additional imported-package provenance regression added after
  the full suite collected tests.
- Package floor raised to 91, MMseqs to 97, and orchestration to 95, using the
  floored measured percentages; both coverage gates pass.
- Sphinx warning-as-error build passed; local Markdown links resolve.
- Wheel and sdist built with installed setuptools 84.0.0 / wheel 0.48.0 using
  `python -m build --no-isolation`; archive inspection confirms corpus exclusion.
  The isolated build could not fetch build dependencies in the network-restricted
  environment, so the installed build tools were used instead.
- Final independent review found no remaining actionable integration issues.

The final normal checking run, restricted to CPU affinity 0 and 1, passed in
147.74 seconds (annotation 145.55, audit 0.31, comparison 0.25). All 50,000 records
were conserved, all 168 public fields matched exactly, and no inconsistencies or
unexpected failures were found. This is the third agreeing full run; it does not
include dependency installation or establish GitHub-hosted timing. The final
coverage-contract module also passed all 23 tests after raising the floors.
