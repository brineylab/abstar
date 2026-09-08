# Repository instructions

This file is the source of truth for coding-agent guidance in this repository.
Keep it current when the package layout, supported workflows, or development
commands change.

## Project overview

`abstar` assigns immunoglobulin and T-cell receptor germline genes and produces
detailed V(D)J annotations. Its public surfaces are:

- the `abstar` Click CLI;
- `abstar.run()` for annotation from Python;
- `abstar.gl`, `abstar.pp`, and `abstar.tl` convenience namespaces;
- AIRR TSV and Parquet output;
- bundled BCR and TCR germline databases.

The package requires Python 3.10 or newer. Core dependencies include
`abutils`, Biopython, Polars, PyArrow, Parasail, and Click. MMseqs searches use
checked argument-list subprocesses with binaries resolved by the public
`abutils.bin.get_path` accessor; read merging uses binaries exposed by `abutils`.
Biopython's strict FASTQ iterator validates complete, potentially multiline
records before the normal input parser and assignment run.

## Repository map

- `abstar/core/abstar.py`: top-level `run()` orchestration, input discovery,
  chunking, multiprocessing, output assembly, logging, and cleanup.
- `abstar/assigners/mmseqs.py`: active V, J, D, and C assignment pipeline.
- `abstar/annotation/annotator.py`: per-sequence annotation orchestration.
- `abstar/annotation/junction.py`: conservative raw FWR3 anchor recovery and
  joint V-region boundary recovery constrained by retained alignment coordinates.
- `scripts/audit_annotation_consistency.py`: read-only native-Parquet audit of
  region assemblies and mask lengths against the direct oriented-query V(D)J
  slice, preserving file/row identity independently of external sequence IDs.
- `scripts/run_corpus.py`: fixed-corpus annotation, conservation, consistency,
  and exact per-record baseline gate; `scripts/build_ci_corpus.py` selects its
  inputs from explicit read-only source paths.
- `test_data/bcr_corpus/`: committed enriched BCR corpus and baseline, excluded
  from Python distributions and ordinary pytest matrix/coverage runs.
- `abstar/annotation/antibody.py`: annotation data model and logging state.
- `abstar/annotation/germline.py`: germline lookup and segment realignment.
- `abstar/annotation/{regions,positions,indels,mutations,mask,productivity}.py`:
  biological feature calculations.
- `abstar/annotation/schema.py`: Polars output schema.
- `abstar/preprocess/merging.py`: paired/interleaved FASTQ merging.
- `abstar/annotation/umi.py`: UMI parsing.
- `abstar/core/germline.py`: custom germline database construction.
- `abstar/scripts/abstar.py`: Click commands and CLI-to-API argument mapping.
- `abstar/germline_dbs/`: packaged databases and MMseqs indexes.
- `abstar/tests/`: pytest suite; shared fixtures are in `conftest.py` and input
  fixtures are in `abstar/test_data/`.
- `docs/source/`: editable Sphinx documentation sources.
- `.github/workflows/`: test and PyPI publishing workflows.

The main data flow is:

```text
input -> optional UMI/merge preprocessing -> MMseqs assignment
      -> temporary Parquet chunks -> parallel annotation
      -> AIRR TSV and/or Parquet output
```

## Environment and common commands

Create or use an isolated Python 3.10+ environment. Invoke Python tools through
the selected interpreter so that `pytest`, imports, and package metadata come
from the same environment.

On Apple Silicon, Parasail 1.3.4 builds from source. Install the Xcode Command
Line Tools and `brew install autoconf automake libtool m4`, then set
`M4="$(brew --prefix m4)/bin/m4"` in the installation environment. The macOS
wheel CI job provisions these tools and exercises a native Parasail alignment.

```bash
python -m pip install -e .
python -m pytest -q
python -m pytest abstar/tests/test_regions.py -q
python -m pytest abstar/tests/test_regions.py::test_get_region_sequence_fwr1 -q
```

Run the complete statement-and-branch coverage gate with:

```bash
python -m pytest --cov=abstar --cov-branch --cov-report=term-missing --cov-report=json:/tmp/abstar-coverage.json -q
python scripts/check_coverage.py /tmp/abstar-coverage.json coverage-floors.json
```

The first command enforces the package floor in `.coveragerc`; coverage.py
compares that total at the configured two-decimal precision. The second
enforces the critical-module floors in `coverage-floors.json` against raw JSON
percentages. Both stored floors are selected by mathematically rounding down
coverage.py's combined statement-and-branch `percent_covered` value. Raise
floors when measured coverage increases, and never lower them without treating
the change as a release-gate regression.

The supported pytest scopes are:

```bash
# fast unit/property scope used on every supported Python
python -m pytest -m "not integration and not e2e and not slow" -q

# real component boundaries, complete public entry points, and slow checks
python -m pytest -m "integration" -q
python -m pytest -m "e2e" -q
python -m pytest -m "slow" -q

# combined real integration/end-to-end CI gate and complete local suite
python -m pytest -m "integration or e2e" -q
python -m pytest -q
```

Focused release gates are:

```bash
python -m pytest abstar/tests/test_airr.py -q
python -m pytest abstar/tests/test_database_integrity.py -q
python -m sphinx -W --keep-going -b html docs/source docs/_build/html
```

The dedicated `.github/workflows/corpus.yml` job runs the committed enriched BCR
subset on pushes and pull requests, separately from the pytest version matrix.
Use its fixed Python 3.12 reproduction profile and a new external output path:

```bash
python -m pip install -r requirements-corpus.txt
POLARS_MAX_THREADS=2 OMP_NUM_THREADS=2 python scripts/run_corpus.py \
  --output /tmp/abstar-corpus-check
```

The runner checks immutable source identities, expected failures, independent
region/mask consistency, and exact public fields against the committed baseline.
It retains `report.json`, actual outputs, and diagnostics. Never refresh the
baseline to silence an unexplained difference. Follow
`abstar/tests/README.md` and `test_data/bcr_corpus/README.md` for reproduction
and reviewed baseline updates. Five minutes is the target for the dedicated
job; local timings do not establish hosted-runner performance.

Broader discovery against the full published BCR corpus is optional and never part
of ordinary push or pull-request CI. It requires all input paths explicitly and
must write outside the source and input trees:

```bash
python scripts/discover_bcr_cases.py \
  --fasta-dir /path/to/bcr_fastas \
  --manifest /path/to/sample_manifest.csv \
  --cellranger-root /path/to/cellranger \
  --output /tmp/abstar-bcr-candidates.jsonl \
  --per-dataset 25 --n-processes 2
```

Audit existing native Parquet annotations without rerunning assignment:

```bash
python scripts/audit_annotation_consistency.py /path/to/parquet \
  --report /tmp/abstar-consistency.json
```

The report must be new and outside the checkout and input trees. Exit 1 means
inconsistencies were found; exit 0 means the implemented checks passed. The
reference is `sequence_oriented[v_sequence_start:j_sequence_end]`, translated
using the one-based `frame` for AA checks, not the final-file `sequence` fields.
This checks internal consistency, not the biological correctness of boundaries.

For an exact regression comparison after an explicitly requested corpus rerun,
use the original FASTAs, both run roots, and reviewed expectations for every
baseline failure. Write a new report outside the checkout and run/input trees:

```bash
python scripts/compare_annotation_runs.py \
  --baseline /path/to/original_run --candidate /path/to/rerun \
  --fasta-dir /path/to/bcr_fastas \
  --expected abstar/test_data/lc_anchor_failures.json \
  --report /tmp/abstar-comparison.json
```

This compares all public Parquet fields and schemas for original successful
records, verifies recovered failures against literal expectations, and checks
record conservation using original FASTA ordinals even when IDs repeat.

Full external corpus discovery is a manual workflow; there is no scheduled corpus
job. The scheduled documentation linkcheck is independent of corpus discovery
and is not an ordinary push or pull-request gate.

Useful CLI checks after an editable install:

```bash
abstar --help
abstar run --help
abstar build_germline_database --help
abstar run path/to/sequences.fasta path/to/project_directory
```

Build the documentation with:

```bash
python -m pip install -r docs/doc_requirements.txt
make -C docs html
```

Build release artifacts, when relevant, with:

```bash
python -m pip install build
python -m build
```

There is currently no committed Ruff, Black, mypy, or other formatter/type
checker configuration. Do not imply that such a check is an official project
gate. Follow the surrounding style and report the command and version if an
ad-hoc checker is used.

## Current baseline caveat

Do not assume that the unmodified test suite is green under every dependency
version permitted by `requirements.txt`. At `v0.8.0`, one permitted environment
with `abutils==0.5.4` produced 48 failures, 107 passes, 2 skips, and 9 xfails;
the major failure class was missing `abutils.tl.translate` and
`abutils.tl.reverse_complement` APIs.

Before diagnosing a regression, record the interpreter and important package
versions:

```bash
python -VV
python -m pip show abstar abutils polars pyarrow parasail
```

Do not hide baseline or new failures with broad catches, relaxed assertions,
or new unconditional `xfail` markers. State which failures predated the change
and which are new.

## Implementation guardrails

### Protect scientific correctness

Plausible but incorrect annotations are more harmful than explicit failures.
For changes in assignment or annotation logic:

- Preserve the original sequence identifier exactly. Treat all external IDs as
  strings and use a separate immutable internal row key for joins.
- Propagate `receptor`, locus, strand, and germline-database identity explicitly
  through assignment, lookup, realignment, and serialization. Do not rely on a
  BCR default in shared BCR/TCR code.
- Define the coordinate space used by every start/end value: input query,
  oriented query, alignment, ungapped germline, IMGT-gapped germline, or AIRR.
  Convert at named boundaries rather than with scattered arithmetic.
- Preserve score evidence for gene calls. Ranking must use documented metrics
  and deterministic tie-breakers; retain biologically meaningful allele ties.
- Constrain J and D candidates by receptor and compatible chain locus.
- Verify productivity using frame, junction, stop-codon, ambiguity, and
  receptor-appropriate motif rules.
- Sequence annotation exceptions are recoverable by default: omit their incomplete
  rows, persist per-record diagnostics and a failure index, report counts (and an
  API warning), and continue other records/samples. `strict=True` / `--strict`
  restores abort-on-record-error behavior. An all-failed sample may have empty
  output only with explicit persistent failure accounting. Biological nonassignment,
  invalid input, and programming errors must remain distinguishable.
- Record diagnostics use `logs/<input-stem>/<encoded-id>__<row-key>.failed`;
  `logs/failures.tsv` indexes exact IDs and source paths. `logs/run.json` records
  parameters, versions, and source hashes. Preserve diagnostics across reruns and
  no-project API cleanup; worker/tool/output/storage failures remain fatal.
- Check record conservation: every input record must result in an annotation or
  an explicit, inspectable failure status.

Any change to gene calls, coordinates, junction/CDR3 construction, mutation or
indel representation, identity, productivity, or AIRR fields needs a focused
biological regression test with exact expected values.

### Dataframe and parallel-processing boundaries

- Supply explicit Polars schemas for identifier-bearing TSV/CSV reads; do not
  allow type inference to reinterpret IDs such as `10E8` or values with leading
  zeroes.
- Validate join cardinality and uniqueness before selecting the first match.
- Consume arbitrary iterables only once, or materialize them deliberately.
- Keep output ordering deterministic across process counts and chunk sizes.
- Keep the Python API return shape stable; do not make it depend on how many
  records happened to annotate successfully.
- Validate `chunksize`, process counts, output formats, and empty inputs near
  the public boundary.
- Treat worker and external-tool failures as run failures unless the public API
  explicitly returns their status.

### Germline databases

The currently packaged database names are:

- BCR: `human`, `macaque`, `c57bl6`, `balbc`, and `human+c57bl6`;
- TCR: `human`.

Each database may contain:

- `ungapped/` FASTA files used for alignment;
- `imgt_gapped/` FASTA files used for IMGT-aware annotation;
- `mmseqs/` generated search indexes;
- `manifest.txt` provenance metadata.

When changing a packaged database:

- keep gapped and ungapped IDs identical and unique;
- confirm that removing IMGT dots from each gapped sequence reproduces its
  ungapped counterpart;
- validate locus and segment naming, including BCR/TCR distinctions;
- update provenance and expected gene-count tests;
- regenerate MMseqs indexes from the source FASTA rather than editing index
  files directly;
- stage custom-database builds and expose them only after every required step
  succeeds.

User databases normally live under `~/.abstar/germline_dbs/<receptor>/`. Tests
must redirect this to a temporary location and must not modify a developer's
real user database.

### External commands and temporary files

- Prefer argument lists with `subprocess.run(..., check=True)` or an equivalent
  checked wrapper over `shell=True` command strings.
- Quote and validate paths and never interpolate untrusted additional arguments
  into a shell command.
- Capture enough stdout/stderr to diagnose MMseqs or read-merging failures.
- Own temporary directories with a context manager whose lifetime encloses all
  workers that use them.
- Tests must use pytest temporary directories and verify cleanup. Do not write
  fixtures or results into the source tree.

## Testing expectations

For assignment-dependent regression failures and corpus baseline changes, follow
[the test harness debugging guide](abstar/tests/README.md). Reproduce the exact
cohort and environment before reducing inputs or changing expected results.

Place tests in `abstar/tests/` and follow the existing `test_<module>.py`
layout. For a change:

1. Add the smallest focused regression test that fails for the original bug.
2. Run the affected test module.
3. Run the full suite when dependencies and runtime permit.
4. Report exact pass/fail/skip/xfail counts and any environmental blocker.

Tests should assert biological values and contracts, not only file existence or
column presence. In particular, assert exact or explicitly ambiguous gene
calls, coordinates, strand, junction/CDR3, productivity and reason codes,
sequence IDs, row counts, and error status.

End-to-end changes should cover both the CLI and Python API where relevant and
exercise real public entry points. Important matrix dimensions include:

- single and multi-record FASTA/FASTQ;
- file, directory, list, iterator, and generator inputs;
- BCR and supported TCR chains;
- forward and reverse-complement reads;
- truncated, ambiguous, indel-containing, nonproductive, and unassigned reads;
- serial and multiprocessing execution with multiple chunk sizes;
- AIRR TSV and Parquet output;
- paired/interleaved merging, UMI parsing, and custom databases when touched.

For expected exceptions, use `pytest.raises`. Reserve `xfail` for a documented,
tracked defect whose behavior is not the subject of the current change. Empty
dataframes and header-only files are not sufficient success criteria for a
nonempty input.

## Output and compatibility requirements

- Treat `abstar.run()` and the Click commands as public APIs.
- Keep CLI defaults and the Python API intentionally aligned; test both after
  changing arguments.
- AIRR coordinates, field meanings, nulls, booleans, and sequence semantics
  must follow the targeted AIRR schema, not only match a Polars dtype.
- Internal Python/dataframe coordinates and raw Parquet coordinates are
  zero-based half-open. AIRR TSV converts once to one-based closed intervals,
  writes booleans as `T`/`F` and nulls as empty cells, and retains the original
  query in `sequence`. `rev_comp` means annotation coordinates and alignments
  address its reverse complement. NP bases align to gaps in
  `germline_alignment`.
- Final Parquet uses the official `sequence`, `sequence_aa`,
  `sequence_alignment_aa`, and `germline_alignment_aa` meanings while retaining
  native booleans/nulls and internal coordinate offsets. No-project Python
  objects and dataframe returns retain the legacy assembled V(D)J meanings for
  those sequence and amino-acid fields.
- Changes to schemas must be reflected in serializers, dataframe return paths,
  documentation, and compatibility tests.
- FWR/CDR amino-acid regions project established nucleotide intervals into the
  continuous V(D)J query frame. A codon crossing a region boundary belongs to
  the downstream region; partial terminal codons do not translate. Keep
  `cdr_mask_aa` and `cdr3_length` consistent with this partition. `junction_aa`
  and CDR3 V/N/D/J subdivisions retain local junction-frame semantics and must
  not be used to reconstruct continuous-query regions on out-of-frame reads.
- FWR3 and CDR3 share the established junction-start-codon endpoint. Synchronize
  the FWR3 region to that boundary without altering retained V alignment or
  junction evidence; reject an invalid shared interval as a record error.
- Preserve backward compatibility deliberately. If a behavioral break is
  necessary, document it and add a migration note.

## Documentation and generated files

- Edit sources under `docs/source/`, not rendered HTML or doctrees.
- Do not hand-edit `docs/source/_build/`; it contains legacy generated output.
- Keep README examples, Click help, Python docstrings, and Sphinx pages in sync
  with actual signatures and packaged database names.
- Do not claim support, AIRR compliance, performance, or a green version matrix
  unless an automated test establishes it.
- Update documentation in the same change when a public argument, default,
  return type, output field, database, or workflow changes.

## Change discipline

- Inspect nearby code and tests before editing; several modules contain legacy
  commented implementations that are not active behavior.
- Make focused changes and preserve unrelated user work.
- Prefer deleting obsolete code through version control over retaining large
  commented copies.
- Avoid adding another compatibility wrapper when one clear shared
  implementation will suffice.
- Do not manually edit bundled binary/index artifacts unless the task
  explicitly requires regenerating them.
- End each task with the relevant verification commands, observed results, and
  any remaining risk or untested path.
