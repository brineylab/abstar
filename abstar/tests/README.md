# Debugging scientific regression tests

Read [AGENTS.md](../../AGENTS.md) for environment setup and supported test scopes.
This guide explains how to investigate small assignment differences without
weakening the scientific checks. It applies to the existing reviewed fixtures
and the [fixed, committed BCR corpus](../../test_data/bcr_corpus/README.md).
The larger corpus runs in its own push/PR workflow, outside the pytest matrix.

## What a baseline establishes

A frozen abstar output establishes previous behavior, not biological truth.
Unexpected changes to gene calls, coordinates, junction/CDR3, productivity, or
record outcomes require investigation even if the new annotation is internally
consistent. A one-base boundary change can alter region translation and mutation
numbering; a different allele can change many dependent fields.

Keep three kinds of evidence distinct:

- **Reviewed biological fixtures:** exact source sequences, reference evidence,
  and independently reviewed expectations. See the [real BCR fixture
  guide](data/real_bcr/README.md) and [TCR fixture guide](data/tcr/README.md).
- **Regression baselines:** outputs from a recorded implementation, environment,
  database, and input cohort. A mismatch identifies a behavior change that
  needs explanation.
- **Internal consistency checks:** record conservation, coordinate validity,
  region reconstruction, and mask/alignment agreement. These detect inconsistent
  output but cannot establish that a gene call or boundary is biologically right.

Cell Ranger agreement is comparison evidence, not an automatic expected answer.
Likewise, a productive-looking junction is not grounds for preferring one
alignment over another.

## First preserve the failing run

Before updating dependencies, changing parameters, or reducing the input:

1. Preserve the test report, actual output, failure index, per-record diagnostics,
   and `logs/run.json`. Keep temporary investigations outside the checkout and
   source-data trees.
2. Record the tested commit and any local changes. For CI, use the actual checked
   out SHA from the run, which may differ from a branch's current head.
3. Identify records by source file/dataset and original record ordinal. Preserve
   external IDs exactly; duplicate IDs must not collapse into one comparison row.
4. Check the original sequence, cohort membership, order, file grouping, and
   database against the recorded hashes. A changed input is a different test.

Useful environment evidence, using the same interpreter as pytest:

```bash
git rev-parse HEAD
git status --short
python -VV
python -m pip show abstar abutils biopython polars pyarrow parasail pytest
```

Resolve MMseqs through `abutils.bin.get_path("mmseqs")`; record its version and
binary checksum, not just the name of a program on `PATH`. Confirm the resolved
germline database too. A user database under `~/.abstar/germline_dbs/` can shadow
a packaged database. Automated tests must isolate user-database lookup in their
temporary environment, as required by `AGENTS.md`.

## Reproduce the complete assignment context first

Use the same input files and order, receptor, germline database, annotation
process count, MMseqs thread count, assignment batch boundaries, and annotation
chunk size. Distinguish assignment batches from annotation chunks: they are
different stages. Keep the dependency versions and binary/database contents
fixed while comparing the baseline and candidate implementations.

Do not begin by annotating only the failing sequence and assume this is the
same experiment. Changing the cohort, splitting one input into many files, or
altering threading can change the search context. A one-record reproduction
that passes does not invalidate the original failure. After reproducing the
full context, a reduced input is useful for isolating the cause.

There is a concrete precedent in the [junction-anchor validation
report](../../docs/superpowers/reports/2026-09-07-junction-anchor-validation.md#separate-assignment-reproducibility-finding):
five original successful records changed during a corpus rerun, including four
support-value changes and one V-allele change. Controls with anchor recovery
disabled reproduced those differences. That separated the assignment variation
from the recovery change; it did **not** prove a universal explanation for
MMseqs variability or justify ignoring future assignment differences.

Repeat the unchanged baseline and the candidate under the same complete
context. Interpret the results explicitly:

| Observation | Next investigation |
| --- | --- |
| Baseline stable; candidate consistently differs | Trace the candidate change and its biological consequences |
| Repeated runs of the same implementation differ | Isolate assignment reproducibility or another uncontrolled input |
| Difference disappears only after changing grouping or threads | Preserve the original reproduction and investigate that parameter separately |
| Both now agree, but neither matches the stored baseline | Check the recorded environment, cohort, binary, and reference provenance |

A retry that happens to pass does not repair nondeterminism. Retain both outcomes
and explain what differed; do not add an automatic retry-until-green policy.

## Find the first stage that differs

Inspect assignment calls and search evidence before comparing downstream
annotations. In [select_best_hits](../assigners/mmseqs.py), ranking uses bit
score, E-value, identity, target/query coverage, and alignment length. Exact
evidence ties retain multiple calls, with deterministic representative details.

- A changed support value alone is not proof of a changed biological annotation,
  but it still needs a reproducibility explanation. Compare the hit ranking and
  retained query/reference traces before treating it as numerical noise.
- A changed ambiguity set is meaningful. Removing allele suffixes, comparing only
  the first call, or accepting any member of a broadened set can hide a regression.
- A different primary reference can legitimately change its associated alignment,
  mutation, or boundary evidence. Validate those fields together against that
  exact reference; do not reuse coordinates authenticated for a different allele.

Use `debug=True` in a temporary project to retain assignment work files and
alignment logs when needed. Replaying the same saved assignment through
annotation helps distinguish an MMseqs change from an annotation change. Existing
examples are in [test_junction_anchor.py](test_junction_anchor.py),
[test_region_boundary_recovery.py](test_region_boundary_recovery.py), and
[test_real_bcr.py](test_real_bcr.py). Preserve receptor/database identity when
reconstructing an annotation object.

Relevant focused commands are:

```bash
python -m pytest abstar/tests/test_mmseqs.py -q
python -m pytest abstar/tests/test_germline.py -q
python -m pytest abstar/tests/test_real_bcr.py -q
```

For native Parquet output, [audit_annotation_consistency.py](../../scripts/audit_annotation_consistency.py)
can check assemblies against the oriented-query V(D)J slice without rerunning
assignment. Pass a fresh report path outside the checkout and input trees.
Consult its `--help` and `AGENTS.md`; AIRR coordinates and legacy in-memory
sequence fields are not interchangeable with native file fields.

## When an expectation may change

An intentional biological change needs independently supported new values,
a focused regression test, and a reviewable explanation of the old/new result.
An intentional environment or reference update needs its own provenance and
reproducibility assessment. Neither is a reason to regenerate every expected
answer without examining the differences.

For a baseline update, retain:

- Exact affected source identities and changed fields, including linked changes
  to calls, coordinates, junctions, mutations, productivity, and record status.
- The cause, supporting query/reference evidence, and a reproducible comparison
  between the old and new implementations or environments.
- An explicit explanation of any newly accepted allele ambiguity, expected
  failure, or narrowly justified numeric comparison policy.
- Confirmation that unaffected records retain their expected results and that
  record conservation and independent consistency checks still pass.

Do not repair a failure with blanket float tolerances, gene-name normalization,
aggregate-only thresholds, silent record filtering, broad exception catches, or
unconditional `xfail`. An allowed alternate reference must have its own reviewed
evidence and any reference-dependent expectations. The narrow J4/J5 treatment
in `expected_for_selected_reference` in [test_real_bcr.py](test_real_bcr.py) is
an example of an explicit evidence-based distinction, not a general permission
to accept whichever gene the latest run selects.

## The committed corpus gate

The fixed sample is enriched in difficult cases, committed with the testing
harness, and run in one dedicated push/PR job targeting about five minutes.
Local timings are sizing evidence, not proof of GitHub runtime.

The corpus README documents the reproduction and proposed-baseline commands;
the manifest records input grouping, selection, and process/thread settings.
`requirements-corpus.txt` pins the critical reproduction dependencies. Start with
`report.json` in the runner output directory, which includes comparison evidence
and a pointer to this guide. The baseline's metadata records its environment and
reference provenance.

Generate the baseline from the exact committed cohort and verify repeatability
before freezing it. Record per-sequence selection reasons and source provenance;
separate representative-sample summaries from the deliberately enriched cases.
Never regenerate expectations in the CI checking path.

The corpus gate must fail on unexpected biological-output changes and new record
failures. Expected unassigned or ambiguous-failure records must be listed
individually and conserved, rather than covered by an allowed failure percentage.
Retain their diagnostic categories; do not require incidental filesystem paths
or timestamps to match a biological baseline.

Keep the selected cohort fixed across commits. Do not silently drop slow,
difficult, or failing records to meet the runtime target. Changing the cohort
or comparison policy is a reviewed harness change. Keep the large corpus out
of repeated version-matrix and coverage runs, and retain the smaller reviewed
fixtures as independent checks of biological correctness.
