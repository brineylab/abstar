# Spike: V search and alignment throughput

Date: 2026-11-20. Task 1 of the foundation plan (spec: *Milestones and order*,
item 1). The code is in [`spikes/vsearch/`](../../../spikes/vsearch/). It is
throwaway and never joins the workspace.

## Result: PASS

| Gate | Target | Measured |
|---|---|---|
| Single-threaded V step | ≥ 3,000 records/sec/core | **4,282 records/sec** |
| Scaling at all physical cores (16) | ≥ 70% of linear | **82%** (56,338 records/sec) |

Both gates pass on this machine. One of five 16-thread runs was an outlier at
63% scaling (see [Variance](#variance)). The numbers come from an x86_64 AMD
machine, not the spec's M-series Mac, so re-check them on an M-series Mac
before relying on them there.

## Machine

- AMD Ryzen AI Max+ 395 (Radeon 8060S), x86_64, 16 physical cores and 32
  threads (SMT 2), 5.19 GHz maximum clock, 16 MiB L2, 64 MiB L3, 123 GiB RAM.
- SIMD used: AVX2 (block-aligner's `simd_avx2` backend; the CPU also has
  AVX-512BW, but block-aligner 0.5.1 has no AVX-512 backend).
- Linux 7.2.5, `performance` CPU governor, rustc 1.100.0.
- The machine was shared during the runs: one unrelated process held about one
  core for the whole session, and another user process started near the end.
- Clocks: about 5.0 GHz with one busy core, and 3.8–4.1 GHz across all cores
  under the 16-thread load (sampled from `/proc/cpuinfo`). Some of the gap
  from linear scaling is lower all-core clock, not software.

## What the prototype does

For each record, in [`src/main.rs`](../../../spikes/vsearch/src/main.rs):

1. **K-mer prefilter.** A CSR index maps every k-mer to the V alleles that
   contain it (each allele counted once per k-mer). The read and its reverse
   complement both count shared k-mers against all 342 alleles. The top N
   (allele, strand) pairs across both strands become candidates.
2. **Score-only SIMD alignment of each candidate.** This is semi-global affine
   alignment with block-aligner 0.5.1 (AVX2). The germline aligns end to end;
   read sequence before and after it is free (`FREE_QUERY_START_GAPS` and
   `FREE_QUERY_END_GAPS`, germline as query). Scoring: match +2, mismatch −3,
   gap open −12, gap extend −2, N scores 0. The best score wins, over both
   strands. Alleles tied at the best score on the winning strand are kept as a
   tie list.
3. **Traceback for the winner.** One more alignment with traceback produces the
   CIGAR, because a real engine needs V coordinates. This cost is in the timing.

Records are spread across threads with rayon `map_init`. Rayon calls the
scratch constructor once per work split, not once per thread, and the split
happens inside the timer. Each scratch allocates two 512 × 2048 alignment
blocks (score and traceback), two padded read buffers and a CIGAR, so these
allocations are repeated more often than a real engine would repeat them. The
bias is conservative: the spike understates throughput, possibly more at high
thread counts, where rayon splits work more finely. It may account for part
of the gap from linear scaling. The benchmarks were not re-run after this was
found.

The V database is `abstar/germline_dbs/bcr/human/ungapped/v.fasta`. It holds
342 human V alleles: 198 IGHV, 64 IGKV and 80 IGLV (the plan calls them
"IGHV", but the file holds all three loci). The corpus is also mixed: 16,796
IGH, 16,682 IGK and 16,521 IGL records.

## Choosing k and the candidate count

Sweep: 200,000 records (corpus ×4), 16 threads, agreement on the 49,999
baseline records.

| k | N=1 | N=2 | N=3 | N=5 | N=8 |
|---|---|---|---|---|---|
| 7 | 95,984 / 99.51% | 63,417 / 99.79% | **49,971 / 99.87%** | 36,060 / 99.91% | 23,264 / 99.93% |
| 9 | 89,737 / 99.32% | 62,385 / 99.74% | 51,947 / 99.83% | 35,755 / 99.88% | 25,973 / 99.91% |
| 11 | 86,964 / 99.15% | 66,715 / 99.65% | 50,742 / 99.79% | 39,230 / 99.87% | 27,065 / 99.90% |
| 12 | 98,972 / 99.03% | 71,402 / 99.56% | 56,334 / 99.74% | 39,322 / 99.85% | 26,975 / 99.89% |

Cells are records/sec and gene-level agreement. (These sweep runs were shorter
and noisier than the final runs below.)

- Alignment dominates the cost: throughput scales with N, and k hardly
  matters. Each extra candidate costs about 75 µs of one thread.
- Smaller k recovers the legacy call slightly more often, because the
  prefilter is more sensitive on highly mutated reads.
- **Chosen: k = 7, N = 3.** Going to N = 5 gains 0.04 points of agreement and
  costs about 28% of throughput; N = 1 is 0.36 points worse.

## Throughput (k = 7, N = 3)

1,000,000 records (the 50,000-record corpus cycled 20 times in memory; the
timer covers the V step only, not FASTA parsing). Runs were interleaved with
20-second pauses.

| Threads | Runs (records/sec) | Median | Per thread |
|---|---|---|---|
| 1 | 4,282; 4,245; 4,284 | **4,282** | 4,282 |
| 16 | 56,779; 56,338; 55,474 | **56,338** | 3,521 |
| 32 (SMT) | 64,087; 60,506; 64,238 | **64,087** | 2,003 |

- **Scaling efficiency**, against 16 physical cores:
  56,338 ÷ (16 × 4,282) = **82.2%** (81.0% at the slowest run).
- 32 threads (for information): 64,087 records/sec, 14% over 16 threads from
  SMT, or 93.5% of 16 × single-thread.
- 1M records take about 17.8 s at 16 threads and 3.9 minutes single-threaded.
- The checksum (the sum of scores plus CIGAR lengths) matched across all
  thread counts, so the parallel runs did the same work.

### Variance

The first 1M-record pass ran 1 → 16 → 32 threads back to back. Its 16-thread
run measured **42,935 records/sec (62.7% scaling)**, just after the
4-minute single-thread run. Repeated 16-thread runs on their own gave 56,106,
56,779, 56,338 and 55,474, and the 200k-record diagnosis runs gave 56.7k
unpinned and 58.6–59.8k pinned to one thread per core. The outlier is most
likely heat, power or other load on this shared machine. It is not a property
of the code, but it shows that sustained all-core throughput depends on the
machine. The spec's benchmark harness should repeat runs and report medians.

## Accuracy sanity check (informational)

Gene-level agreement with the legacy `v_call` in `baseline.parquet`, matched on
`corpus_ordinal` = FASTA record order. All 49,999 baseline sequences were
checked byte-identical to the FASTA at that ordinal; ordinal 11457 is the
corpus's expected annotation failure and has no baseline row. Allele suffixes
(`*NN`) are stripped. A record agrees if the spike's tie list and the legacy
comma-separated call share any gene, so ties on either side count as
agreement. A record with no spike call (no shared k-mer) counts as
disagreement; there were none in this corpus. The spike's tie list covers only
alleles among the top N = 3 candidates on the winning strand. It therefore
under-counts ties, which makes this agreement figure conservative.

| | Records | Agreement |
|---|---|---|
| All | 49,999 | **99.87%** |
| IGH | 16,796 | 99.75% |
| IGK | 16,682 | 99.95% |
| IGL | 16,521 | 99.90% |

The disagreements (67 records) are mostly close paralogs (IGHV3-23 vs
IGHV3-15, IGKV1-39 vs IGKV1-8, IGLV1-44 vs IGLV10-54 at low score). There are
also a few odd contigs, such as a legacy IGHV call where the spike's best hit
is an IGKV allele. The legacy call is not ground truth.

Orientation, also informational: legacy abstar marks 26 records
`rev_comp = true`. The spike calls 12 of them on the minus strand. For the
other 14, it finds the same gene on the plus strand with a high score (for
example, ordinal 42425 scores 584 against IGLV2-8 on the plus strand, where
the legacy V identity is 1.0). The spike calls no other record minus. These
records are worth a look when the real engine handles orientation.

## Headroom and caveats

- **Full DP band.** block-aligner's `FREE_QUERY_END_GAPS` requires a minimum
  block size larger than the germline (284–317 nt), so the block is 512 and
  each alignment fills essentially the whole germline × read matrix (about
  300 × 576 cells). A real engine could align only the read window suggested
  by the k-mer diagonal (about 340 nt instead of about 600), or use a banded
  or hand-written kernel. That is an untested speedup of roughly 40%.
- **Orientation.** Both strands are tried in the prefilter, and candidates
  from either strand compete on alignment score. In this corpus only 26
  records are reverse-complemented, so a "top N per strand" policy would cost
  almost nothing extra here, but up to 2× alignment work on mixed-orientation
  input.
- **Budget against the spec target** (≥ 10k records/sec end to end on a
  10-core M-series laptop): at this machine's per-core rate and scaling, 10
  cores would give about 4,282 × 10 × 0.82 ≈ 35k records/sec for V. That leaves
  about 70% of the per-record time for J, D and annotation. This is only an
  extrapolation until it is re-measured on an M-series Mac.
- Corpus replication cycles the same 50,000 records in memory rather than
  holding 20 physical copies, so cache behavior is slightly kinder than a
  true 1M-record input.

## Reproduce

From the repository root:

```bash
cd spikes/vsearch && RUSTFLAGS="-C target-cpu=native" cargo build --release && cd ../..
BIN=spikes/vsearch/target/release/vsearch
V=abstar/germline_dbs/bcr/human/ungapped/v.fasta
SEQ=test_data/bcr_corpus/sequences.fasta.gz

# Throughput: k=7, N=3, corpus x20 (1M records), 1/16/32 threads
$BIN $V $SEQ 7 3 20 16
$BIN $V $SEQ 7 3 20 32
$BIN $V $SEQ 7 3 20 1

# Predictions for one corpus pass, then agreement with the legacy calls
$BIN $V $SEQ 7 3 1 16 /tmp/preds.tsv
python spikes/vsearch/agreement.py /tmp/preds.tsv test_data/bcr_corpus/baseline.parquet
```

Usage: `vsearch <v.fasta> <seqs.fasta[.gz]> <k> <n_candidates> <reps>
<threads,csv> [preds.tsv]`. `agreement.py` needs Polars. On Apple Silicon, the
`Cargo.toml` selects block-aligner's `simd_neon` backend automatically.
