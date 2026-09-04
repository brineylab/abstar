# Release-Gating Test Harness Design

**Date:** 2026-09-04

**Branch:** `test/release-gating-suite`

## Purpose

`abstar` needs a test harness that detects biological misannotation, silent
record loss, serialization errors, dependency incompatibilities, and failures
in real public entry points before a release is published. The existing suite
has useful unit coverage and several audit-driven regression tests, but it is
not yet a trustworthy release gate: many tests assert only that a field or file
exists, biological expectations are concentrated around one human heavy-chain
sequence, integration coverage is commented out, and CI does not validate
installed artifacts or AIRR semantics.

This design establishes a layered test system whose highest-value cases use
published human BCR sequences. It also defines how test-driven production fixes
may be included when a new regression case exposes a broken public contract.

## Current Baseline

On the current project interpreter (Python 3.12.14, `abutils==0.6.0`,
`polars==1.44.1`, `pyarrow==25.0.0`, `parasail==1.3.4`, and `pytest==9.1.1`),
the suite collects 241 tests and reports:

```text
230 passed, 2 skipped, 4 xfailed, 5 xpassed, 6 warnings in 162.93s
```

The five unexpected passes are stale Polars workarounds. The four expected
failures misuse `xfail` around assertions that already use `pytest.raises`.
The two skips refer to a packaged database named `mouse`, although the
supported mouse database names are `c57bl6` and `balbc`. The six warnings are
an `abstar` use of deprecated `pl.count()`.

The suite has no committed pytest configuration, no property-testing or
coverage dependency, no active integration module, no installed-CLI test, no
AIRR reference validation, and no automated release dependency on a successful
test run.

## Scope

The work will:

- make pytest configuration and dependency boundaries explicit;
- remove stale skips and expected-failure markers from supported behavior;
- consolidate duplicated fixtures and remove unused opaque reference pickles;
- add invariant, property, database-integrity, integration, end-to-end, and
  failure-injection tests;
- create a small, attributed fixture corpus from published human BCR data;
- create deterministic variants of adjudicated real sequences for controlled
  coordinate, indel, frame, ambiguity, truncation, and strand tests;
- validate AIRR TSV with the AIRR reference library plus semantic assertions;
- compare Parquet and AIRR values from the same run;
- test the installed `abstar` command and `abstar.run()` public API;
- add measured coverage floors for critical modules and the package overall;
- make tests, artifact installation, database integrity, AIRR validation, and
  documentation build prerequisites for publishing.

Focused production fixes are in scope when an in-scope regression test proves
that the current implementation violates an intended public or biological
contract. Each such fix must be introduced test-first and kept in a separate,
reviewable commit. Refactoring unrelated production code, restoring the legacy
BLAST assigner, or redesigning the complete annotation data model is outside
this effort.

`python -m abstar` is not currently an advertised public interface and will not
be added solely to satisfy the old audit matrix. The installed Click command,
Python API, convenience namespaces, AIRR TSV, and Parquet are the public
surfaces tested here.

## Test Architecture

### Fast tests

Unmarked tests are the fast default. They exercise pure or narrowly isolated
logic without launching MMseqs, fastp, process pools, or full annotation runs.
They include existing unit regressions, schema behavior, and new property tests.
They must be deterministic and suitable for every supported Python version.

Hypothesis will cover invariants with broad input domains and minimized failure
examples:

- raw, aligned, germline, IMGT-gapped, and AIRR coordinate conversions round
  trip at their named boundaries;
- reverse complementation is involutive;
- gapping and ungapping preserve biological positions;
- segment intervals stay ordered and within sequence bounds;
- insertion and deletion representations reconstruct aligned sequences;
- hit ranking is invariant to input row order and retains genuine ties;
- list, iterator, and generator inputs conserve the same records;
- identifier-bearing serialization preserves strings such as `10E8`, leading
  zeroes, Unicode, whitespace, long values, and duplicates;
- serializers preserve values and do not mutate input annotation objects.

### Integration tests

Tests marked `integration` exercise one real component boundary, such as an
MMseqs search against a packaged database, fastp merging, a custom database
build, an AIRR validator call, or a worker process. They use temporary
directories and real public wrappers but avoid unnecessarily repeating a full
pipeline run.

### End-to-end tests

Tests marked `e2e` execute `abstar.run()` or the installed `abstar` command.
They cover:

- FASTA, FASTQ, `Sequence`, list, iterator, and generator inputs;
- single files, flat directories, and nested directories;
- serial and multiprocessing annotation with multiple chunk sizes;
- human IGH, IGK, IGL, TRA, TRB, TRD, and TRG;
- smoke cases for every packaged BCR database;
- productive, nonproductive, ambiguous, truncated, reverse-complemented,
  indel-containing, and no-D records;
- AIRR TSV, Parquet, and Python dataframe/object return paths.

These tests assert exact cardinality, original identifiers, deterministic
ordering, strand, calls or documented allele ambiguity, coordinates,
junction/CDR3, productivity and reason codes, and absence of unexpected
internal errors or temporary artifacts. File-output tests must inspect rows and
biological values; a file's existence or a header-only file is never sufficient
for nonempty valid input.

### Slow tests

The `slow` marker is reserved for larger corpus checks or combinations whose
runtime makes them unsuitable for every pull-request matrix cell. A nightly or
explicit workflow may run a broader published-data sentinel cohort. The
ordinary release gate remains self-contained and bounded.

Pytest will enable strict configuration, strict registered markers, and strict
xfails. Project-originated warnings are errors. Any narrow third-party warning
filter must name the exact warning, source, and reason.

## Published BCR Corpus

### Source and license

The discovery corpus is the version 3 Figshare record:

```text
Functional antibodies exhibit light chain coherence
https://figshare.com/articles/preprint/Functional_antibodies_exhibit_light_chain_coherence/19617633/3
Version: 3
License: MIT
```

The local material comprises 3,441,852 filtered BCR contigs in 94 FASTA files
(approximately 2.0 GB). All 94 datasets occur in `sample_manifest.csv`, which
describes four donors and six flow-class labels. Matching Cell Ranger
`filtered_contig_annotations.csv` files provide chain, V/D/J/C gene calls,
productivity, CDR3 sequence, read count, and UMI count.

The bulk corpus remains external to the repository. CI never relies on a home
directory path, network download, or the presence of the source project.
Selected sequences are copied into the repository under the source MIT license
with the Figshare URL, version, title, and attribution recorded beside them.

### Discovery utility

An optional developer utility accepts explicit paths for the FASTA directory,
sample manifest, and Cell Ranger output root. It does not hardcode a developer
home directory. It produces a candidate report rather than test expectations.

Selection is deterministic. Each record is keyed by an immutable internal
combination of dataset and original contig identifier; the external identifier
is never rewritten. Stable hashing chooses a configurable number of records per
dataset while stratification ensures representation across donor, flow class,
and IGH/IGK/IGL. The default routine sweep selects 200 records from each
dataset, yielding at most 18,800 records before deduplication. A fixed selection
algorithm version and seed are written into the report, so the same inputs and
configuration reproduce the same cohort.

The report records:

- source dataset, original identifier, donor, flow class, chain, and sequence
  SHA-256;
- Cell Ranger calls, productivity, CDR3, read count, and UMI count;
- `abstar` revision and version, Python and dependency versions, receptor,
  germline database, and germline manifest checksum;
- annotation status and exception category;
- normalized gene-level agreement and raw allele-level output;
- junction/CDR3 agreement, productivity agreement, call ties, putative indels,
  no-D status, and sequence/CDR3 length extremes;
- a machine-readable selection reason for every nominated candidate.

Cell Ranger output is independent comparison evidence, not ground truth.
Database and allele nomenclature differences are normalized only for candidate
ranking; raw values are retained.

### Initial observed regression cohort

A pilot run on dataset `1279068` processed 1,935 real contigs. `abstar` exited
successfully while returning 1,927 rows. Eight full-length productive Cell
Ranger contigs disappeared after internal exceptions: four empty global
alignment errors and four gene-segment-mask length errors.

Among the 1,927 surviving records:

- all chain/locus classifications agreed;
- gene-level agreement was 86.1% for V, 96.7% for J, and 96.1% for C;
- 1,925 Cell Ranger CDR3 nucleotide sequences matched the `abstar` junction;
- Cell Ranger marked every filtered contig productive, while `abstar` marked
  551 productive and rejected 1,362 as out of frame;
- `abstar` emitted 176 tied V/J calls, 42 putative V insertions, 27 putative V
  deletions, six IGH records without a D call, and CDR3 lengths from 5 to 30
  amino acids.

The eight lost records are the first mandatory record-conservation regression
cohort. Representative frame disagreements, the two junction disagreements,
ties, putative indels, no-D records, and CDR3 extremes are priority candidates
for adjudication.

### Adjudication and oracle policy

Current `abstar` output must never generate its own golden expectations.
Candidate cases are promoted only after review of the input sequence, relevant
germline alignments, coordinate spaces, Cell Ranger evidence, and targeted
biological rules. A third annotator may contribute evidence, but agreement
between tools alone does not establish correctness.

Each committed case records:

- Figshare record URL and version;
- dataset and exact original contig identifier;
- original sequence and SHA-256;
- donor and flow class where useful to selection provenance;
- source Cell Ranger fields;
- selection reason;
- accepted exact values or an explicit allowed ambiguity set;
- the evidence and reviewer decision supporting accepted values;
- whether the case is original or derived.

Gene-level expectations are used when source databases cannot independently
support an exact allele. Exact allele expectations require direct evidence.
Coordinates name their space and convention explicitly.

### Committed fixture set

The initial repository fixture set will contain approximately 25 to 40
adjudicated cases:

- concordant IGH, IGK, and IGL controls;
- all eight silent-loss examples from the pilot;
- representative productivity/frame disagreements across loci;
- the two observed junction/CDR3 boundary disagreements;
- biologically meaningful allele ties;
- putative insertion and deletion cases that survive adjudication;
- no-D, short-CDR3, and long-CDR3 cases.

The intended layout is:

```text
abstar/tests/data/real_bcr/
├── README.md
├── sequences.fasta
├── cases.json
└── derived_cases.json
```

`README.md` carries the Figshare attribution, version, MIT license statement,
selection method, and field definitions. `cases.json` is human-reviewable and
contains original cases and their accepted expectations. `derived_cases.json`
describes transformations; it does not duplicate generated sequence data when
the transformation can be reproduced exactly.

### Derived real-background cases

Derived cases start from clean, adjudicated parents. Each transformation names
the parent, operation, input-query offset, affected bases, resulting sequence
hash, and expected invariant or exact change. Supported operations are:

- reverse complement;
- five-prime and three-prime truncation;
- single-nucleotide ambiguity or substitution;
- one- and two-nucleotide frameshift insertions/deletions;
- three-nucleotide in-frame insertions/deletions;
- V-region and junction-region edits.

Tests assert the consequences dictated by the transformation: strand,
coordinate shifts, indel representation, reconstruction, identity change,
frame, productivity, and stable unaffected calls. Mutation offsets are defined
in the original input-query coordinate space and converted only at named
boundaries.

The real BCR corpus does not replace TCR or non-human fixture sources. Human
TCR goldens and packaged non-human database smoke cases remain separate
requirements.

## Failure Injection

The harness will inject:

- missing or incompatible dependency capabilities;
- MMseqs and fastp nonzero exits with captured diagnostics;
- malformed, empty, ambiguous, and duplicate-ID inputs;
- missing or incomplete germline databases;
- one broken record among valid records;
- unwritable output targets;
- worker-process exceptions;
- custom-database build failure and rollback.

Expected biological non-assignment, invalid input, external-tool failure, and
internal programming errors must remain distinguishable. Every input record
must produce an annotation or an explicit inspectable failure result. A run
with nonempty input and only internal failures must fail. Tests must prove that
an internal exception cannot become an ordinary successful empty output.

## AIRR and Parquet Contracts

`abstar` targets the AIRR Data Standards 2.0 Rearrangement schema. AIRR TSV
output will be read and validated by version 2.0.0 of the official AIRR Python
reference library. Internal annotation code may keep Python-native 0-based,
half-open query intervals; a single named serialization boundary converts them
to AIRR's required 1-based, closed intervals. Additional semantic tests will
assert behavior that structural validation alone cannot prove:

- one-based AIRR coordinate conventions and correct endpoints;
- `T`/`F` boolean encoding in TSV;
- original/oriented query semantics for `sequence`;
- separation of templated germline and non-templated bases;
- null and empty-value behavior;
- coherent calls, alignments, CIGAR values, identity, junction/CDR3, and
  productivity.

For a shared annotation run, normalized AIRR TSV values and Parquet values must
agree field by field. Serialization cannot alter identifiers, ordering,
booleans, coordinates, sequence values, or ambiguity sets.

## Fixture and Runtime Discipline

Shared immutable source data lives in `abstar/tests/data/`; reusable fixture
loaders and assertion helpers live under `abstar/tests/`. Expensive pipeline
runs are performed once per logically shared fixture and their results are read
afresh by tests, avoiding the present pattern of re-annotating 10E8 for many
independent field-existence assertions. Mutable annotation objects are not
shared between tests.

All generated inputs, databases, logs, and outputs use pytest temporary
directories. Tests that redirect the user database also redirect the effective
home/database root and prove that the developer's real `~/.abstar` tree is not
modified. Cleanup tests snapshot owned temporary paths before and after a run.

The unused binary pickle references and their `OLD/` copies are removed after
confirming no active test consumes them. Test data should be compact,
diff-reviewable, reproducible, and traceable to source.

## Dependencies and Configuration

Test tools move out of runtime requirements into a dedicated test dependency
file or project optional dependency group. Direct runtime imports remain
declared as runtime dependencies. The test environment includes pytest,
pytest-cov, Hypothesis, and the AIRR Python reference library with versions
selected across supported Python 3.10 through 3.13.

Pytest configuration in `pyproject.toml` defines test paths, marker
descriptions, strict configuration, strict markers, strict xfails, and warning
handling. Official commands in `AGENTS.md` and contributor documentation are
updated when the fast, integration, end-to-end, and complete commands change.

Coverage is a regression gate rather than a vanity target. Once the new suite
is in place, coverage is measured on the supported reference environment.
Package-wide and critical-module floors are set to the achieved integer
percentages and cannot be lower than the current baseline. Critical modules
include annotation orchestration, germline realignment, region/coordinate
logic, top-level run orchestration, custom germline construction, MMseqs
assignment, and merging/UMI preprocessing. Future changes may raise but not
silently lower these floors.

## Implementation Decomposition

The work stays on `test/release-gating-suite`, but it is implemented as four
sequential, independently reviewable phases. Every phase ends with its affected
tests green and a recorded baseline; later phases build on those results.

1. **Harness foundation:** add pytest configuration and test dependencies,
   repair stale marker/skip/warning debt, centralize reusable fixtures, remove
   unused binary references, and establish database-integrity checks.
2. **Published-data discovery and curation:** add the optional deterministic
   discovery utility, generate the first candidate report from the external
   corpus, adjudicate the bounded committed BCR cohort, record provenance and
   licensing, and add deterministic derived-case machinery.
3. **Biological and public-contract enforcement:** add property, integration,
   end-to-end, failure-injection, AIRR, and Parquet tests; introduce focused
   test-first production fixes for demonstrated contract violations; add the
   separately sourced TCR and non-human coverage.
4. **Measured release gates:** measure and set coverage floors after the suite
   stabilizes, test installed artifacts, divide the supported CI matrix by test
   scope, build documentation, and make publishing consume the exact artifacts
   that passed every required job.

The implementation plan will name the files, test command, expected initial
failure, minimal change, and verification evidence for each task within these
phases. If a phase uncovers a production redesign outside this design's scope,
that redesign receives its own follow-up plan and does not block completion of
unrelated phase work.

## CI and Release Gates

The workflow is divided into clear jobs:

1. Fast tests run on Ubuntu for Python 3.10, 3.11, 3.12, and 3.13.
2. Integration and end-to-end tests run on the minimum and maximum supported
   Python versions with real MMseqs and packaged databases.
3. Minimum- and maximum-supported dependency jobs verify explicit compatibility
   bounds, including `abutils`, Polars, and PyArrow.
4. A wheel and source distribution are built, installed into a clean
   environment, imported, and exercised through `abstar --help`, `abstar run
   --help`, a small installed-CLI annotation, and the Python API.
5. A non-Linux job validates installation, imports, and CLI discovery without
   claiming unsupported external-binary behavior.
6. Database-integrity and AIRR-conformance jobs run their focused gates.
7. Sphinx builds from `docs/source/` with warnings as errors. External link
   checking runs on a scheduled workflow so transient network failures do not
   make ordinary pull requests nondeterministic.

The publishing workflow depends on successful tests and validates the exact
artifacts it uploads. It may not publish a separately rebuilt, untested wheel.
The final supported matrix must have no unexpected failures, skips, xfails,
xpasses, or project warnings.

## Success Criteria

The effort is complete when:

- the current stale xfail/skip/warning debt is removed;
- fast, integration, end-to-end, and slow scopes are documented and selectable;
- all packaged database invariants are automated;
- the committed real BCR fixtures are traceable, licensed, adjudicated, and
  exercise the observed failure and edge-case classes;
- human BCR and TCR goldens assert biological values across all supported loci;
- public API and installed CLI matrices conserve records and order;
- AIRR output passes reference validation and semantic assertions;
- Parquet and AIRR serialization agree;
- failure injection proves that internal and external failures remain visible;
- coverage floors protect critical modules;
- documentation and artifact-install checks run in CI;
- publishing is mechanically blocked unless the tested artifacts pass all
  release gates;
- the full verification report states exact pass, failure, skip, xfail, xpass,
  warning, and coverage results with interpreter and dependency versions.

## Risks and Mitigations

Biological adjudication is the primary schedule risk. The fixture set starts
small, prioritizes already observed failures, and permits gene-level ambiguity
where allele evidence is insufficient. It never promotes bulk current output
to expected output.

Runtime growth is controlled through markers, shared immutable run artifacts,
and a bounded committed cohort. Broader corpus discovery is optional and
external to CI.

Differences between Cell Ranger and `abstar` may reflect database versions or
nomenclature rather than a defect. Candidate reports retain raw evidence and
normalize only for ranking; manual review decides the contract.

Large production changes uncovered by a test are split into independently
reviewable test-first commits. If satisfying a case requires redesign outside
the stated scope, the candidate report preserves the reproduction and evidence
for a dedicated follow-up. Such a case is not added as an expected failure or
represented as success in the release-gating suite.
