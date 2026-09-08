# Fixed BCR regression corpus

This versioned sample contains 50,000 distinct original sequences from the
3,441,852-record BCR corpus: 40,000 stratified background records and 10,000
deliberately enriched difficult records. It runs once per push and pull request
in [the corpus workflow](../../.github/workflows/corpus.yml), outside the ordinary
pytest version matrix and coverage gate. These files are excluded from source
and wheel distributions; running the gate requires a repository checkout.

The source is **“Functional antibodies exhibit light chain coherence,” Figshare
record 19617633, version 3**, with the same MIT attribution as the
[reviewed original BCR fixtures](../../abstar/tests/data/real_bcr/README.md).
Original dataset and contig identifiers and sequence contents are retained.
The baseline records abstar behavior; it is not independently adjudicated
biological truth. Small reviewed biological fixtures remain separate tests.

## Run the gate

Use Python 3.12 in an isolated environment, from the checkout root:

```bash
python -m pip install -r requirements-corpus.txt
python -m pip check
POLARS_MAX_THREADS=2 OMP_NUM_THREADS=2 python scripts/run_corpus.py \
  --output /tmp/abstar-corpus-check
```

The output directory must be new and outside the checkout and corpus trees.
The runner uses the editable checkout and packaged human BCR reference; it
rejects a user database that shadows that reference. One input FASTA preserves
the committed order. Assignment uses two MMseqs threads and a 1,000,000-record
assignment batch size; annotation uses two processes and 500-record chunks.
The dependency profile pins critical packages, with actual versions and binary
hashes recorded in each report. It is not a complete transitive lockfile.

The runner authenticates input bytes and metadata, runs the public Python API,
checks record conservation and individual failure outcomes, independently audits
region and mask consistency, then compares the ordered schema and every public
Parquet field exactly. It applies no numeric tolerance or gene normalization.
Environment differences are reported; reference, input, and parameter drift are
rejected. A conserved row can still fail the biological regression comparison.

Start debugging with `report.json`, then native output in `parquet/`,
`logs/run.json`, `logs/failures.tsv`, and per-record diagnostics. The workflow
uploads the run directory even when checking fails. Read the
[debugging guide](../../abstar/tests/README.md) before reducing an input, changing
threading, or replacing expectations: small assignment differences can depend
on the complete search context.

## Selection and identity

All 94 source datasets, four donors, six flow classes, and 298 observed
donor/flow/locus/dataset strata are represented. Background selection balances
these strata rather than reproducing their original abundance. Difficult
selection additionally balances mechanism labels. Thus neither panel estimates
the prevalence of biological events in the original corpus.

Candidates rank by SHA-256 of `bcr-ci-v1`, source filename, and zero-based record
ordinal, separated by NUL characters. Bounded candidate pools and sorted
round-robin strata make selection deterministic. Mandatory originals are chosen
first, then difficult records, then background records. Deduplication uses the
exact case-sensitive sequence hash across both panels. Background records may
also have difficulty labels; panel membership is exclusive, labels can overlap.

The 100 candidate cases in the original BCR, anchor-recovery, region-boundary,
and FWR4 endpoint fixtures matched 97 distinct original sequences. All are
included, with fixture reasons unioned when several fixtures share a sequence.
Synthetic/edited probes are never substituted for unmatched original reads.

Examples of overlapping selection labels in this version are:

| Label | Selected records | Sampling criterion |
| --- | ---: | --- |
| V deletion / insertion | 1,289 / 1,037 | Saved annotation reports a V indel |
| C deletion / insertion | 3 / 3 | Saved annotation reports a C indel |
| High SHM | 2,261 | V identity below 0.90 |
| Truncated V | 120 | V germline start at least 15 nt |
| Short / long CDR3 | 1,379 / 1,717 | At most 5 aa / at least 25 aa for IGH or 13 aa for light chains |
| Frameshift / stop codon | 122 / 2 | Saved productivity issue label |
| Nonproductive | 1,716 | Saved productive value is false |
| V / J allele ties | 1,974 / 1,794 | Multiple retained calls |
| Reverse complement | 26 | Saved reverse-complement flag |
| Missing annotation | 19 | No saved annotation for the original record |

Labels come from saved annotations and only guide sampling; new outputs were
generated from the fixed cohort for the baseline. The source scan found no
ambiguous-base originals, so that mechanism remains covered by focused tests.
The manifest records full label counts, thresholds, selection rules, and hashes
of source FASTAs, saved annotations, the source manifest, fixtures, and builder.

External IDs are not unique: 758 records repeat an earlier selected ID.
`records.parquet` preserves `(source_file, record_ordinal)` as the source identity
and `corpus_ordinal` as the fixed run order. Both ordinals are zero-based.
`sequence_id` remains a string and `sequence_sha256` authenticates its exact
sequence. No comparison joins on the external ID alone.

## Files and expected failure

- `sequences.fasta.gz`: deterministic gzip of the original two-line FASTA.
- `records.parquet`: identity, source stratum, panel, and selection reasons.
- `manifest.json`: version, authenticated artifact hashes, selection provenance,
  fixed parameters, and explicit expected failures.
- `baseline.parquet`: every native public output field plus `corpus_ordinal`.
- `baseline-metadata.json`: parameters, comparison policy, code/environment,
  reference and MMseqs provenance, and outcome counts.

Corpus ordinal **11457**, source `1279072.fasta` ordinal **2587**, external ID
`CTAATGGTCCAGAGGA-1_contig_3`, has an independently tested ambiguous region
boundary. Its sequence SHA-256 is
`d382fe17805110581bd08a48bce3baef535bd51a4f2dad57fb7909ce7598e912`.
[The original fixture](../../abstar/test_data/region_boundary_recovery.json)
and [focused tests](../../abstar/tests/test_region_boundary_recovery.py) establish
that competing optimal boundaries must fail explicitly. The manifest requires
stage `annotation`, category `internal_error`, and exception type
`abstar.annotation.junction.JunctionAnchorError` for this exact record.
An unexpected recovery also requires review. Message wording and diagnostic
filesystem paths are not compared as biological fields.

## Review a proposed baseline or cohort change

First run the existing comparison and inspect every changed record using the
debugging guide. For a justified update, generate proposed files externally:

```bash
POLARS_MAX_THREADS=2 OMP_NUM_THREADS=2 python scripts/run_corpus.py \
  --output /tmp/abstar-corpus-proposal --record-baseline
```

This still enforces conservation, consistency, and the manifest's individual
expected failures. It cannot overwrite committed expectations. Review the
old/new field differences and their biological evidence, retain the report,
copy only the accepted baseline and metadata into this directory, and update
their SHA-256 entries in `manifest.json`. Repeat a normal full-cohort check
before committing. CI never records or refreshes a baseline.

Changing the cohort is a separate reviewed harness change. The builder requires
explicit local sources, never downloads data, and never reruns annotation:

```bash
python scripts/build_ci_corpus.py \
  --fasta-dir /path/to/bcr_fastas \
  --manifest /path/to/sample_manifest.csv \
  --annotations /path/to/saved/native/parquet \
  --output /tmp/abstar-corpus-selection \
  --size 50000 --difficult-fraction 0.2
```

Authenticate any expected failures against independent fixture evidence before
baseline generation. Never remove difficult records or relax comparisons simply
to obtain a passing job or meet the runtime target. The job targets five minutes
with a ten-minute hard timeout for cold installation and diagnostic headroom;
hosted timing must be established by an actual GitHub run.

Initial local validation produced three agreeing full runs. The final normal
check, restricted to two CPU cores, took **147.74 seconds** and matched all
**168 public fields**: 49,999 annotated records, the one expected failure, and
zero consistency findings. This excludes dependency installation and is not a
GitHub-hosted benchmark. The complete corpus and baseline occupy about 43.5 MB.
