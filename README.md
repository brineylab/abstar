![](https://img.shields.io/pypi/v/abstar.svg?colorB=blue)
![](https://img.shields.io/badge/license-MIT-blue.svg)

# abstar

`abstar` annotates antibody and T-cell receptor sequences. For each read it calls
V, D, J and constant-region alleles (with ties), finds the junction and the IMGT
regions, and reports mutations, as AIRR Rearrangement TSV or Parquet. The
annotator is written in Rust; the command line and the Python API are a thin
Python layer over it.

Each topic below is covered in full by the built-in manual, `abstar manual`.

## Install

```bash
pip install abstar
```

Python 3.11 or later. Wheels are built for Linux (x86_64 and aarch64),
macOS (Apple silicon) and Windows (x64); they include the germline databases,
so there is nothing else to install and no network access at run time. On
other platforms pip builds from source, which needs a
[Rust toolchain](https://rustup.rs).

### On a cluster (HPC)

No root access, compiler or external tools are needed:

```bash
module load python/3.11        # or any Python >= 3.11 (conda, micromamba, uv)
python -m venv ~/venvs/abstar
~/venvs/abstar/bin/pip install abstar
```

- The Linux wheels need glibc 2.17 or later, which covers CentOS/RHEL 7 and
  newer.
- `abstar run` uses every core it is allowed. Under a scheduler, pass the
  allocation explicitly, e.g. `--threads $SLURM_CPUS_PER_TASK`.
- On CPUs without AVX2 (roughly pre-2013), abstar switches to a slower
  aligner by itself. If polars warns about missing CPU features, install its
  compatibility build: `pip install "polars[rtcompat]"`.
- User-built germline databases go in `~/.abstar`. Set `ABSTAR_HOME` to put
  them somewhere else, such as project storage when home quotas are small.

## Annotate

```bash
abstar run reads.fasta out/                    # out/reads.tsv and out/run.json
abstar run fastqs/ out/ --format parquet       # a directory, searched recursively
abstar run reads.fastq.gz > reads.tsv          # no OUTPUT_DIR: TSV on stdout
abstar run reads.fasta out/ -g c57bl6 --locus IGH
abstar run --seq CAGGTGCAGCTGGTGCAGTCTGGG       # one sequence
```

FASTA and FASTQ are accepted, gzipped or not, in either orientation (reverse
reads are detected per read). `abstar run --help` lists every option.

From Python:

```python
import abstar

df = abstar.run("reads.fasta")                 # a polars DataFrame
df = abstar.run([("id1", seq1), ("id2", seq2)])  # or a list of strings
summary = abstar.run("fastqs/", "out/", formats=["airr", "parquet"])
```

Every record gets a row. One that can't be annotated is `unassigned`, with a
`failure_reason`; `abstar manual fields` defines every column.

## Germline databases

Built in: `human` (IG and TR), `c57bl6`, `balbc` and `human+c57bl6` (for
humanized mice). Every V allele must be IMGT-gapped.

```bash
abstar germline list
abstar germline info human
abstar germline build donor.toml               # e.g. human restricted to a donor's genotype
abstar run reads.fasta out/ -g donor
```

`abstar manual germline` describes the recipe format.

## Synthetic reads

`abstar generate` simulates reads with known truth for benchmarking, and
`abstar.grade()` scores calls against it. See `abstar manual generate`.

## The manual

```bash
abstar manual              # everything, as Markdown
abstar manual examples     # common invocations
abstar manual migration    # moving from abstar 0.x
abstar manual errors       # error codes and exit codes
```

The documentation site built from `docs/` has the same content.

## Development

Inside a virtualenv with a Rust toolchain:

```bash
pip install maturin
maturin develop --release --extras test
cargo test --workspace
python -m pytest -q
```

`AGENTS.md` covers the repository layout and conventions, and design decisions
are recorded in `python/abstar/decisions/`.
