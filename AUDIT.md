# Deep critical review of `abstar` v0.8.0

Audit date: 2026-09-03
Audited revision: `c0b7c29` (`v0.8.0`)

## Scope and bottom line

This review covered all production Python modules, tests, workflows, bundled
germline databases, packaging, CLI/API paths, and repository and hosted
documentation. It included the complete test suite, production-code coverage,
static analysis, targeted runtime reproductions, and integrity checks against
the bundled FASTA pairs.

Four release-blocking correctness problems were identified, together with
several additional high-risk annotation defects and major gaps in the tests and
documentation. The most serious recurring pattern is that internal annotation
exceptions are converted into apparently successful, empty output files. This
makes genuine software defects indistinguishable from "no rearrangements
found."

## Executive summary

| Priority | Finding | Consequence |
| --- | --- | --- |
| P0 | Allowed `abutils` releases break annotation, while exceptions are swallowed | Complete silent loss of annotations |
| P0 | TCR annotation looks up genes in the BCR database | Effectively all TCR records fail and disappear |
| P0 | Polars infers sequence-ID types | IDs such as `10E8` are silently changed |
| P0 | MMseqs hit ranking uses an uninformative `nident` column | Wrong gene/allele can be selected without warning |
| P1 | AIRR output violates important semantic and serialization requirements | Files appear valid but are not reliably interoperable |
| P1 | Productivity logic is incomplete and has incorrect TCR motif rules | False productive/unproductive calls |
| P1 | TCR D-segment logic uses the wrong loci | TRA is searched for D, TRB is not |
| P1 | Custom database construction and discovery are broken | Documented feature is largely unusable |
| P1 | Generators, multiple files, nested samples, partial reads, merging, and UMI paths have data-loss defects | Real-world workflows fail or silently omit records |
| P2 | Dead implementations, large commented blocks, broad mutable state, and weak dependency boundaries | High maintenance and regression risk |
| P2 | The current test suite is red and does not test biological correctness end-to-end | Releases cannot be trusted from CI results |
| P2 | README, RTD, database lists, and API examples disagree with implementation | Users are led into broken or misleading workflows |

## P0 — release blockers

### 1. A permitted dependency version causes complete, silent annotation loss

`abstar` permits any `abutils>=0.5.1` in
[`requirements.txt`](requirements.txt#L1). The installed, permitted version
`0.5.4` does not expose either `abutils.tl.translate` or
`abutils.tl.reverse_complement`, both of which `abstar` calls.

A canonical forward-strand 10E8 assignment completed successfully, then
annotation raised:

```text
AttributeError: module 'abutils.tools' has no attribute 'translate'
```

The exception originates in
[`germline.py`](abstar/annotation/germline.py#L422), but
[`annotator.py`](abstar/annotation/annotator.py#L110) catches every
`Exception`, records it as a failed annotation, and continues. The top-level
run then:

- writes a syntactically valid AIRR TSV containing only its header;
- can write a valid zero-row Parquet file;
- prints a zero-annotation summary;
- returns successfully;
- exposes the failure only in a secondary log, sometimes under a leaked
  temporary directory.

The full test suite reproduces the dependency incompatibility. PyPI currently
lists newer `abutils` releases as well, so the unbounded upper range is an
ongoing compatibility risk, not only a local-environment peculiarity. See the
[`abutils` release history](https://pypi.org/project/abutils/).

Recommended action:

- Replace the broad compatibility range with a tested range or an exact
  compatible release.
- Import the required operations from their stable documented locations.
- Add a startup capability check.
- Treat "input count > 0, output count = 0 because every record raised" as a
  run failure.
- Return explicit per-record failure status from the Python API.

### 2. TCR annotations are resolved against the BCR database

MMseqs assignment selects the correct TCR database, but downstream realignment
does not preserve the receptor type. `realign_germline()` and its callers
ultimately use the default BCR receptor in
[`germline.py`](abstar/annotation/germline.py#L90), including V and J
processing around [`germline.py`](abstar/annotation/germline.py#L202).

A direct lookup demonstrates the defect:

- `TRAV1-1*01__homo_sapiens` with the default receptor fails.
- The same lookup with `receptor="tcr"` succeeds.

Because the exception is swallowed by the annotation loop, valid TCR records
become empty "successful" output rather than a failed run. This contradicts
the hosted documentation's claim of broad TCR support in the
[TCR annotation guide](https://abstar.readthedocs.io/en/latest/tcr.html).

Recommended action:

- Make `receptor` required in all annotation and germline lookup calls.
- Store it in an immutable run context rather than relying on defaults.
- Add end-to-end golden tests for TRA, TRB, TRD, and TRG before continuing to
  advertise full TCR support.

### 3. Sequence identifiers are silently corrupted by Polars inference

MMseqs query data is written as TSV in
[`mmseqs.py`](abstar/assigners/mmseqs.py#L518), then read without an explicit
string schema around [`mmseqs.py`](abstar/assigners/mmseqs.py#L134) and
[`mmseqs.py`](abstar/assigners/mmseqs.py#L399).

The real-world ID `10E8` is inferred as scientific notation and converted to:

```text
1000000000.0
```

The original V output retained `10E8`, while subsequent D/J queries and the
final annotation used `1000000000.0`. Numeric IDs, leading-zero IDs, and
scientific-notation-like sample names are all at risk.

This can cause:

- loss of provenance;
- incorrect joins;
- collisions between originally distinct IDs;
- incorrect D/J association;
- user-visible annotations attached to the wrong identifier.

Duplicate user IDs are also not rejected, so dataframe joins can form multiple
matches.

Recommended action:

- Introduce an internal immutable row key independent of `sequence_id`.
- Declare all identifier columns as strings at every read boundary.
- Preserve the original identifier byte-for-byte.
- Reject or explicitly disambiguate duplicate IDs.
- Add round-trip property tests using numeric, Unicode, whitespace, very long,
  duplicated, and scientific-notation-like IDs.

### 4. MMseqs hit selection can choose the wrong germline hit

The MMseqs output requests `nident` and uses `--alignment-mode 3` in
[`mmseqs.py`](abstar/assigners/mmseqs.py#L132). In inspected V, D, and J result
files, every hit had `nident=0`, despite strongly differentiated E-values.

Nevertheless, V, D, J, and C selection sort only on `nident` before retaining
the first hit, for example
[`mmseqs.py`](abstar/assigners/mmseqs.py#L175) and
[`mmseqs.py`](abstar/assigners/mmseqs.py#L291). E-value, bit score, alignment
coverage, and normalized identity are not tie-breakers.

Consequences:

- the selected gene can be determined by database/output order rather than
  alignment quality;
- raw identical-base count would favor longer alignments even if it were
  populated;
- genuinely tied alleles are collapsed to one arbitrary allele instead of
  producing an ambiguity set;
- J searches are not constrained by the V-chain locus, allowing biologically
  incompatible cross-locus J calls;
- the D threshold permits E-values as high as `1e6` with alignments as short as
  five bases, making weak D overcalling likely.

In the 10E8 reproduction, IGHV3-15 alleles `*07` and `*01` tied exactly but
only one was emitted without ambiguity.

Recommended action:

- Verify the meaning and population of every requested MMseqs field.
- Rank by a documented tuple such as normalized identity, coverage, E-value,
  and bit score.
- Preserve tied allele sets.
- Constrain J/D candidates by receptor and chain locus.
- Calibrate D thresholds using negative controls and curated rearrangements.

## P1 — high-priority correctness defects

### 5. AIRR output is structurally plausible but semantically non-compliant

Several important discrepancies exist:

- The code acknowledges AIRR's 1-based convention but emits 0-based,
  half-open coordinates in
  [`germline.py`](abstar/annotation/germline.py#L488).
- Polars serializes booleans as lowercase `true`/`false`; AIRR TSV uses `T`/`F`.
- `sequence` is constructed as an assembled V(D)J slice rather than preserving
  the unmodified/oriented query in
  [`annotator.py`](abstar/annotation/annotator.py#L644).
- Observed NP-region nucleotides are copied directly into the assembled
  `germline`, making non-templated bases appear germline-derived.
- V/D/J CIGAR fields are declared but never populated.
- D and J identity are not calculated.
- Constant-region gapped output appends ungapped values in
  [`annotator.py`](abstar/annotation/annotator.py#L663).

These conflict with the coordinate, boolean, and sequence semantics in the
official
[AIRR Rearrangement schema](https://docs.airr-community.org/en/v2.0.0/datarep/rearrangements.html).

Calling the output "fully AIRR compatible" should stop until it passes a
schema validator and semantic golden tests.

### 6. Productivity classification can silently be wrong

[`productivity.py`](abstar/annotation/productivity.py#L48) primarily checks
stop codons, uppercase `N`, locus agreement, and expected junction motifs. It
does not reliably establish:

- an open reading frame;
- junction length modulo three;
- correct V-to-J frame;
- validity of an empty or truncated junction;
- other ambiguous IUPAC bases or lowercase ambiguity characters.

The terminal motif rule also groups `IGH`, `TRA`, and `TRD` as expecting
tryptophan. IMGT distinguishes the immunoglobulin-heavy J-TRP motif from the
T-cell-receptor J-PHE motif; see the
[IMGT feature definitions](https://www.imgt.org/download/LIGM-DB/ftable_doc.html).

This creates both false productive and false unproductive calls.

### 7. TCR D-segment logic uses the wrong loci

D-query construction and annotation fallback use:

```python
["IGH", "TRA", "TRD"]
```

in [`mmseqs.py`](abstar/assigners/mmseqs.py#L661) and
[`annotator.py`](abstar/annotation/annotator.py#L431).

That:

- searches TRA, which has no D segment;
- omits TRB, which does have D segments;
- conflicts with the bundled human TCR D database, which contains `TRBD` and
  `TRDD`;
- is actively encoded as expected behavior in
  [`test_mmseqs.py`](abstar/tests/test_mmseqs.py#L245).

The D reassignment path also hardcodes `IGHD` in
[`germline.py`](abstar/annotation/germline.py#L334), making it unsuitable for
TRBD/TRDD.

There are additional boundary defects: the minimum-length test examines the J
alignment span instead of the extracted pre-J query, and forward slicing
includes the first J base. Similar one-base contamination is present in reverse
J/C query construction.

### 8. Custom germline database support is broken and unsafe

The documented database builder calls:

```python
get_database_directory(location, receptor)
```

in [`core/germline.py`](abstar/core/germline.py#L170), but the function
signature is `(receptor, db_location)` at
[`core/germline.py`](abstar/core/germline.py#L486). The default call therefore
fails while trying to lowercase `None`.

Additional defects include:

- custom lookup repeats the receptor directory component in
  [`annotation/germline.py`](abstar/annotation/germline.py#L75);
- if any V sequence already contains an IMGT gap, all sequences are returned
  unchanged, leaving mixed inputs partially ungapped;
- all loci, including IGK/IGL and TCR, are aligned against an IGHV reference in
  [`core/germline.py`](abstar/core/germline.py#L307);
- overwrite does not clear the previous database, allowing stale files to
  survive;
- V, D, and J builds are attempted unconditionally;
- external commands use `shell=True` and their return codes are not enforced;
- IDs, sequence alphabets, duplicates, loci, and minimum segment sets are not
  validated;
- database creation is not staged transactionally, so failures can leave a
  discoverable partial database.

The builder should be considered unsupported until it has a clean-room
end-to-end test.

### 9. Several input shapes silently lose data

Distinct defects affect common Python and directory workflows:

- `_process_inputs()` consumes an iterable using `all(...)`, then attempts to
  iterate over it again in [`core/abstar.py`](abstar/core/abstar.py#L484). A
  two-sequence generator produced a zero-byte FASTA.
- For multiple input files, `sequence_df` is overwritten on every iteration in
  [`core/abstar.py`](abstar/core/abstar.py#L383), so only the last file's
  dataframe is returned.
- An empty directory can leave `sequence_df` unbound.
- Recursive input discovery uses only basenames as sample keys in
  [`core/abstar.py`](abstar/core/abstar.py#L318). Nested files with the same
  name, or files with different extensions but the same stem, can overwrite
  one another.
- Input copying flattens the directory hierarchy.
- Completed chunks are accumulated in `as_completed` order, making output
  order nondeterministic.
- The return type depends on the number of successful records: a multi-record
  run with one survivor can unexpectedly return a scalar `Sequence`, concealing
  the other failures.
- Unsupported `output_format` values are not rejected consistently and can
  produce output directories without output.

These APIs need deterministic, stable return types and explicit input/output
cardinality contracts.

### 10. Truncated reads can raise instead of being annotated as partial

`get_region_sequence()` returns bare `""` for a missing region in
[`regions.py`](abstar/annotation/regions.py#L91), while callers unpack three
return values in [`annotator.py`](abstar/annotation/annotator.py#L930). A
legitimate 5-prime-truncated read can therefore raise and disappear through the
broad exception handler.

Missing biological regions should be represented as a typed, nullable result,
not by changing the function's return shape.

### 11. Several coordinate formulas are demonstrably unsafe

Notable examples:

- J germline end adds the germline start twice in
  [`germline.py`](abstar/annotation/germline.py#L499).
- C sequence end double-counts one query offset in
  [`germline.py`](abstar/annotation/germline.py#L717).
- D germline end uses the query alignment span, which is wrong when D contains
  indels.
- CDR3 and segment masking recover coordinates by calling `.find()` on sequence
  motifs in [`regions.py`](abstar/annotation/regions.py#L340) and
  [`mask.py`](abstar/annotation/mask.py#L181). Repeated motifs or `find() == -1`
  can produce plausible but incorrect slices.
- IMGT deletion endpoints are derived by adding raw length, ignoring
  intervening IMGT gap characters in
  [`indels.py`](abstar/annotation/indels.py#L127).
- Nongermline masking indexes before validating bounds and can terminate early
  with a shortened mask.

The region-boundary adjustment code also advances by codon-sized steps while
testing nucleotide gaps. That is a high-risk area requiring biological golden
cases before its intended behavior can be trusted.

The underlying design problem is the absence of one canonical coordinate
model. Raw query, oriented query, local alignment, germline, IMGT-gapped, and
AIRR coordinates are transformed ad hoc.

### 12. Reported identity is incomplete

V and C identity are calculated approximately as:

```python
1 - mutation_count / len(germline)
```

in [`annotator.py`](abstar/annotation/annotator.py#L905). Insertions and
deletions are excluded, and the denominator does not necessarily correspond to
the aligned span expected by AIRR. D and J identity are left null.

This can produce identity values that look precise but are not comparable with
standard aligner output.

### 13. Paired-read merging has multiple broken paths

- `run(merge=True)` expands `**merge_kwargs` even when its default is `None` in
  [`core/abstar.py`](abstar/core/abstar.py#L279).
- Interleaved merging opens a temporary file in binary mode and writes strings
  in [`merging.py`](abstar/preprocessing/merging.py#L366).
- Pair discovery merely destructures two sorted files without confirming that
  one is R1 and one R2 in
  [`merging.py`](abstar/preprocessing/merging.py#L191).
- Filename construction uses chained `rstrip()` character sets rather than
  suffix removal, allowing stem corruption and collisions.
- `fastp` execution is assembled through `shell=True`, with incompletely quoted
  paths and user-supplied additional arguments.

### 14. UMI handling has contradictory and invalid edge behavior

In [`umi.py`](abstar/preprocessing/umi.py#L133):

- `length=None` is documented as inferable but is compared numerically and
  raises.
- Negative-length patterns reverse-complement the pattern but then add the
  negative length to the match endpoint.
- Built-in patterns document two mismatches, but the function's default of one
  prevents that pattern-specific default from taking effect.
- Pattern matching searches the full sequence even though documentation
  describes end-limited matching, permitting coincidental internal matches.
- `pattern=None, length=None` eventually calls `len(None)`.
- Standalone file/iterable modes drop records without a detected UMI, while the
  main annotation pipeline may retain them; the behavior differs by entry
  point.
- File mode can replace the input when no separate output is supplied.

## P2 — design, simplification, and maintenance

### 15. Exception handling conflates biological outcomes with software defects

The central annotation loop treats all exceptions alike. "No D gene assigned,"
"truncated read," dependency `AttributeError`, invalid coordinate math, and
programming bugs all become failed records.

Introduce a small result hierarchy:

- `Annotated`
- `BiologicallyUnassigned`
- `InvalidInput`
- `InternalError`

Internal errors should fail the run by default or trip a configurable
error-rate threshold. They should never silently produce a normal empty result.

### 16. Temporary-directory ownership is incorrect

[`core/abstar.py`](abstar/core/abstar.py#L220) creates
`TemporaryDirectory(...).name` without retaining the manager. The object can
immediately delete its directory, after which subdirectories are recreated
manually and no longer cleaned up.

After the test runs, 119 `/tmp/abstar*` directories remained. Project-local
`tmp` directories also do not consistently satisfy the cleanup behavior
promised by documentation.

Use one run-scoped context manager whose lifetime encloses all worker jobs and
output consolidation.

### 17. There is substantial dead and duplicated implementation

Examples:

- [`blastn.py`](abstar/assigners/blastn.py#L40) defines `Blastn` twice, imports
  removed modules, and cannot be imported successfully.
- `mmseqs.py` is 1,630 lines and contains a very large commented legacy
  implementation after the active code.
- Similar commented legacy blocks remain in `assigner.py`, `umi.py`, regions,
  and indel code.
- `utils/build_germline_dbs_OLD.py` contains Python 2 `raw_input`.
- Several `bin/` scripts reference nonexistent modules or obsolete APIs.
- The Dockerfile is based on Ubuntu 14.04, Anaconda 2, and MongoDB 3.2 even
  though the package now requires Python 3.10+.
- Built Sphinx output is committed under `docs/_build`, adding 82 generated
  files.

Delete obsolete implementations through version control rather than retaining
them as comments. If BLAST remains a supported alternative, rebuild it behind
the same tested assigner interface; otherwise remove it.

### 18. Dependency and package boundaries are weak

The production requirements include testing and plotting packages that are not
directly required by core runtime paths, while directly imported packages such
as Biopython, `natsort`, and `tqdm` are inherited transitively through
`abutils`.

There is also circular package metadata: current `abstar` requires `abutils`,
while the installed `abutils` declares an `abstar` dependency. That makes
independent compatibility and installation behavior difficult to reason about.

Recommended separation:

- direct runtime dependencies;
- optional visualization dependencies;
- test dependencies;
- documentation dependencies;
- explicit, tested `abutils` compatibility bounds.

### 19. Process-wide warning suppression hides important signals

[`abstar/__init__.py`](abstar/__init__.py#L5) suppresses `BiopythonWarning` and
all `FutureWarning` globally. Importing `abstar` therefore changes warning
behavior for the host application and hides upstream compatibility warnings.

Suppress only narrowly understood warnings inside the smallest relevant
context.

### 20. The annotation model is overly broad and weakly typed

The primary annotation object exposes roughly 147 fields, with `airr_fields`
derived dynamically from `self.__dict__.keys()`. Many fields are declared but
never populated, and some declared types disagree with schema expectations.

A better division would be:

- immutable input/read identity;
- assignment candidates and score evidence;
- canonical coordinate/alignment model;
- biological interpretation;
- format-specific serializers.

AIRR, Parquet, and internal representations should not all depend on one
mutable object's incidental attributes.

## Test assessment

### Current state

Running the release tag with the project interpreter produced:

```text
48 failed, 107 passed, 2 skipped, 9 xfailed
```

Production-code coverage was approximately **40%**:

| Module | Coverage |
| --- | ---: |
| `annotation/annotator.py` | 18% |
| `annotation/germline.py` | 33% |
| `annotation/regions.py` | 28% |
| `core/germline.py` | 13% |
| `preprocessing/merging.py` | 18% |
| `core/abstar.py` | 67% |
| `assigners/mmseqs.py` | 92% |

The high MMseqs line coverage is misleading because the tests generally verify
that files or columns exist, not that the selected gene, coordinates, scores,
or record associations are correct.

Other test problems:

- [`test_integration.py`](abstar/tests/test_integration.py#L1) is entirely
  commented out.
- Five xfails acknowledge existing Polars/multi-sequence annotation defects.
- Four other xfails are tests for expected exceptions and should be ordinary
  `pytest.raises` tests.
- Mouse tests skip because the expected `mouse` database does not exist under
  that name.
- Output tests accept header-only TSVs and zero-row dataframes.
- No test asserts conservation of input records.
- No test validates AIRR semantics.
- No test builds and then uses a custom database.
- No test exercises the public CLI as an installed package.
- No test detects temporary-directory leakage.

### Required end-to-end suite

A release-gating E2E suite should include the following matrix.

#### Entry points

- installed `abstar` CLI;
- `python -m abstar`;
- Python API;
- FASTA, FASTQ, `Sequence`, list, iterator, and generator inputs;
- single file, flat directory, and nested directory;
- serial and multiprocessing modes.

#### Biological coverage

- human IGH, IGK, and IGL;
- each advertised non-human BCR database;
- human TRA, TRB, TRD, and TRG;
- productive, out-of-frame, stop-containing, ambiguous, truncated,
  reverse-complemented, indel-containing, and no-D cases;
- reads with known tied allele calls;
- weak D-like random sequence negative controls.

#### Workflow coverage

- paired and interleaved merging;
- UMI extraction from both ends and missing-UMI behavior;
- custom database build, discovery, overwrite, failure rollback, and
  annotation;
- TSV and Parquet serialization;
- supported database discovery.

#### Assertions

- exact input/output cardinality and identifiers;
- exact V/D/J/C calls or an explicit allowed ambiguity set;
- exact raw, germline, IMGT, and AIRR coordinates;
- junction/CDR3 sequence and translation;
- productivity with reason codes;
- strand/orientation;
- deterministic record ordering;
- no temporary artifacts;
- no unexpected per-record internal errors;
- AIRR schema and semantic validation.

#### Failure injection

- unavailable or incompatible dependency API;
- MMseqs nonzero exit;
- malformed and empty input;
- duplicate and pathological IDs;
- missing database files;
- one broken record among otherwise valid records;
- unwritable output;
- worker-process failure.

A run with nonempty input and only internal failures must itself fail.

### Unit and property tests

Unit coverage should concentrate on invariants rather than getters:

- coordinate transforms round-trip between every coordinate space;
- gapping/ungapping preserves biological positions;
- reverse-complement transformations are involutive;
- segment intervals remain ordered and within bounds;
- insertion/deletion representations reconstruct the aligned sequences;
- hit ranking is invariant to input row order;
- tied candidates remain tied;
- iterator and list inputs produce identical output;
- serializers never mutate values or identifiers;
- all bundled database gapped/ungapped pairs agree after dot removal.

The last invariant was checked manually across the active bundled databases:
their paired FASTAs had matching IDs, no duplicate IDs, and matching sequences
after gap removal. That is a positive result, but it should be automated and
extended with manifests, versions, checksums, locus validation, and expected
gene counts.

## CI and release process

The workflow has linting commented out in
[`pytest.yml`](.github/workflows/pytest.yml#L35), no coverage threshold, no docs
build, no AIRR validation, and no package-install smoke test.

The release workflow builds and uploads on a release event using a generic
Python version, without an explicit successful test dependency in
[`python-publish.yml`](.github/workflows/python-publish.yml#L24).

Static analysis reported 165 findings. Many are cosmetic, but definite defects
include:

- duplicate `Blastn` definitions;
- undefined `raw_input`;
- unused/dead imports;
- bare exception handlers.

Minimum release gates should be:

- clean unit and E2E runs;
- coverage thresholds on critical modules, not only repository-wide coverage;
- Ruff and type checking;
- package build plus installation in a fresh environment;
- minimum and maximum supported dependency jobs;
- Linux plus at least one other supported OS;
- Sphinx build with warnings as errors and link checking;
- database integrity checks;
- AIRR conformance checks.

## Documentation defects

The documentation currently overstates support and includes invalid
instructions.

Most important discrepancies:

- The README CLI example omits the required `run` subcommand in
  [`README.md`](README.md#L23). The same stale usage appears on the
  [current PyPI page](https://pypi.org/project/abstar/).
- Installation documentation says Python 3.9, while
  [`pyproject.toml`](pyproject.toml#L9) requires Python 3.10+.
- [`installation.rst`](docs/source/installation.rst#L23) documents a
  nonexistent `mmseqs_binary` argument.
- The documented built-in databases do not match the directories actually
  shipped. BCR currently contains `human`, `macaque`, `c57bl6`, `balbc`, and
  `human+c57bl6`; TCR contains only `human`.
- TCR documentation claims mouse support, all-chain functionality, and correct
  beta/delta D handling despite the receptor lookup and locus defects.
- TCR documentation says constant genes are not assigned, but a TCR constant
  MMseqs database is bundled and C search is executed.
- Python examples iterate the return value as though it always were a list,
  while the implementation collapses one successful result to a scalar.
- "Full AIRR compatibility" is not justified.
- Log documentation mentions only `abstar.log`, omitting failed-annotation and
  MMseqs logs.
- The public API page exposes only a small fraction of the actual supported
  surface.
- No CI job verifies that RTD can build against the current source.

The
[hosted installation guide](https://abstar.readthedocs.io/en/latest/installation.html)
and [hosted TCR guide](https://abstar.readthedocs.io/en/latest/tcr.html) should
be corrected alongside code, then versioned so "stable" does not mix promises
from incompatible releases.

## Recommended remediation order

1. Stop silent success:
   - remove blanket exception swallowing;
   - add explicit result/error types;
   - fail runs dominated by internal errors;
   - fix and constrain `abutils`.
2. Repair correctness foundations:
   - propagate receptor/locus everywhere;
   - introduce canonical coordinate types;
   - repair MMseqs scoring and ambiguity handling;
   - preserve immutable row IDs;
   - rewrite productivity against AIRR/IMGT rules.
3. Establish biological golden E2E tests:
   - BCR and all supported TCR chains;
   - exact calls, coordinates, junctions, productivity, and AIRR output;
   - CLI and Python API.
4. Repair real-world I/O workflows:
   - generators, multi-file results, nested samples, merging, UMI, and
     truncated reads;
   - stable return type and deterministic ordering.
5. Rebuild custom database support transactionally.
6. Split the mutable annotation object and serializers into smaller
   components.
7. Remove dead BLAST/legacy code, generated docs, obsolete scripts, and ancient
   container definitions.
8. Correct documentation only after the tested support matrix is known.
