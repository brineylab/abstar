# Junction-anchor recovery validation — 2026-09-07

The conservative fallback recovers all 20 previously failed LC-coherence records
(19 unique sequences), covering all three observed alignment variations. It runs
only when the original junction-start endpoint has no ungapped query coordinate.
It fits raw FWR3 through IMGT codon 104 from a retained, independently mapped
upstream boundary. All optimal tracebacks must agree on a complete, contiguous
anchor, and retained anchor bases must agree with the recovered coordinates.
No cysteine search or productive-frame preference is used. Unresolved evidence
remains an explicit per-record failure.

## Automated validation

Environment: Python 3.12.14, abstar 0.8.0 (editable), abutils 0.6.0,
Biopython 1.88, Polars 1.44.1, PyArrow 25.0.0, Parasail 1.3.4.

- The initial 20 saved-assignment regression cases failed before implementation.
- Final full suite: **1,498 passed, 0 failed, 0 skipped, 0 xfailed**.
- Combined statement/branch coverage: **90.49%**. Critical-module gates pass;
  positions floor raised to 76%, new junction module floor set to 93%
  (measured 93.75%).
- Sphinx warning-as-error build and `git diff --check` pass.

Commands:

```bash
MPLCONFIGDIR=/tmp/abstar-mpl python -m pytest --cov=abstar --cov-branch \
  --cov-report=term-missing --cov-report=json:/tmp/abstar-coverage.json -q --tb=short
python scripts/check_coverage.py /tmp/abstar-coverage.json coverage-floors.json
python -m sphinx -W --keep-going -b html docs/source /tmp/abstar-anchor-docs
git diff --check
```

The 73 focused anchor tests cover literal reviewed coordinates and biological
values for all 20 cases, reverse complements and query offsets, one-/two-/three-
base upstream indels, competing versus equivalent optimal alignments, retained
alignment conflicts, incomplete references, truncated codons, stop/non-C/
ambiguous codons, nonproductive frameshifts, downstream cysteine decoys, and
TRA/TRB/TRD/TRG controls. Real CLI and Python API tests exercise AIRR and Parquet,
duplicate IDs, input order, workers, and chunk sizes. Failure-recovery contract
tests now inject an exception independently of this biological defect.

## Full corpus

Baseline: `/home/bryanbriney/Projects/lc_coherence/data/abstar_2026-09-07`.
Source: `/home/bryanbriney/Projects/lc_coherence/data/bcr_fastas`.
Candidate: `/tmp/abstar-lc-anchor-validation-w1j0j8tp`.

All 94 FASTAs were reprocessed under strict mode, in eight isolated partitions,
each with four annotation workers, two MMseqs threads, and chunks of 500.
Original inputs and baseline outputs were preserved. Per-partition logs and
metadata remain inside the candidate root; final sample Parquets are linked
under its `parquet/` directory.

| Measure | Count |
| --- | ---: |
| Input records | 3,441,852 |
| Baseline successful records | 3,441,832 |
| Candidate successful records | 3,441,852 |
| Recovered baseline failures | 20 |
| Candidate failures | 0 |
| Original successful rows identical in every public field | 3,441,827 |
| Original successful rows with differences | 5 |

The comparator checks schemas, every public field without float tolerance,
nulls, ordering, source sequences/IDs, record conservation by original FASTA
ordinal, and recovered fields against reviewed literal fixtures. All recovered
records match their expectations. All original successful junctions, CDR3s,
productivity values and coordinates match.

```bash
python scripts/compare_annotation_runs.py \
  --baseline /home/bryanbriney/Projects/lc_coherence/data/abstar_2026-09-07 \
  --candidate /tmp/abstar-lc-anchor-validation-w1j0j8tp \
  --fasta-dir /home/bryanbriney/Projects/lc_coherence/data/bcr_fastas \
  --expected abstar/test_data/lc_anchor_failures.json \
  --report /tmp/abstar-lc-anchor-comparison.json
```

This exact comparison deliberately exits **1** for the five differences below;
it is not an unconditional corpus parity pass. Reports require a fresh path.

## Separate assignment reproducibility finding

| Sample | Sequence ID | Difference |
| --- | --- | --- |
| 1279059 | CCTTACGAGACTGGGT-1_contig_2 | d_support: 4.994e-07 → 1.784e-05 |
| 1279075 | AAAGATGCATGCATGT-1_contig_2 | v_support: 1.014e-128 → 4.459e-130 |
| 1287165 | AGAATAGCATCGGTTA-1_contig_2 | IGHV3-23*03 → IGHV3-23*01; associated V/germline/mutation fields |
| 1287173 | GATCGTAGTTGACGTT-1_contig_1 | v_support: 9.394e-118 → 2.122e-116 |
| 1287188 | CTGTGCTAGGAGTTGC-1_contig_1 | v_support: 1.078e-140 → 2.463e-139 |

All five complete samples were independently rerun with recovery disabled:
the fallback was replaced with the original TypeError when reached. These
controls reproduced **every** difference above. Every original successful row
in each control sample matched the recovery-enabled rerun in every public
field. The assigner source hash also matches the baseline run metadata.
Thus the observed changes are independent of junction recovery. The precise
upstream source of assignment variation remains a separate investigation;
the baseline used automatic MMseqs threading and a different sample grouping.
No assignment logic was changed to suppress these discrepancies.

Detailed local artifacts:

- `/tmp/abstar-lc-anchor-comparison.json`: complete strict corpus comparison.
- `/tmp/abstar-anchor-discrepancy-details.json`: exact changed values.
- `/tmp/abstar-anchor-search-controls-final.json`: five-sample control comparison.
- `/tmp/abstar-anchor-search-control-7e0vyi5t`: first two control samples.
- `/tmp/abstar-anchor-search-control-p0kc2w7g`: remaining three control samples.
- `/tmp/abstar-anchor-full-tests.log`: complete test and coverage output.

This validates the observed failure modes and tested perturbations; it does not
establish correctness for every possible repertoire edge case or detect all
pre-existing errors in the unchanged successful annotation path.
