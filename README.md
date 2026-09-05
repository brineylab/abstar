![](https://img.shields.io/pypi/v/abstar.svg?colorB=blue)
[![tests](https://github.com/brineylab/abstar/actions/workflows/pytest.yml/badge.svg)](https://github.com/brineylab/abstar/actions/workflows/pytest.yml)
[![Documentation Status](https://readthedocs.org/projects/abstar/badge/?version=latest)](https://abstar.readthedocs.io/en/latest/?badge=latest)
![](https://img.shields.io/pypi/pyversions/abstar.svg)
![](https://img.shields.io/badge/license-MIT-blue.svg)

# abstar

`abstar` assigns immunoglobulin and T-cell receptor germline genes and produces
detailed V(D)J annotations. It supports the command line, a Python API, and
AIRR TSV and Parquet output.

- Source: [github.com/brineylab/abstar](https://github.com/brineylab/abstar)
- Documentation: [abstar.readthedocs.org](https://abstar.readthedocs.org)
- Package: [pypi.org/project/abstar](https://pypi.org/project/abstar/)
- Container: [hub.docker.com/r/brineylab/datascience](https://hub.docker.com/r/brineylab/datascience/)

## Install

`abstar` supports Python 3.10 through 3.13. The tested dependency ranges include
`abutils>=0.6,<0.7`, `polars>=1.5,<2`, and `pyarrow>=16.1,<26`.

```bash
pip install abstar
```

MMseqs2 and fastp executables are provided through `abutils`; no separate
system installation is required for the normal packaged workflow.

## Use

Annotate one FASTA or FASTQ file into a project directory:

```bash
abstar run path/to/sequences.fasta path/to/project_directory
```

Directories are discovered recursively:

```bash
abstar run path/to/input_directory path/to/project_directory
```

Use `-o parquet`, or repeat `-o` to write both formats. Use `--receptor tcr`
for TCR annotation. Run `abstar run --help` for merging, UMI, germline database,
and performance options.

The API preserves the input category:

```python
import abstar
from abutils import Sequence

one = abstar.run(Sequence("N", id="one"))
many = abstar.run([Sequence("N", id="first"), Sequence("N", id="second")])
frame = abstar.run("sequences.fasta", as_dataframe=True)
```

One input record returns one `abutils.Sequence`; multiple records return a list,
including when only one is assigned; `as_dataframe=True` always returns a Polars
DataFrame. Every input record has `annotation_status`. Biological non-assignment
produces an `unassigned` row with `failure_reason`; internal or external-tool
failures raise `abstar.AnnotationRunError` with structured failures and retained
artifact paths. Duplicate visible IDs and input ordering are preserved.

AIRR TSV targets the AIRR 2.0 Rearrangement schema. It writes 1-based closed
coordinates, `T`/`F` booleans, and empty null cells. Python returns and Parquet
coordinates are 0-based half-open. Both final file formats write the original,
unmodified query as `sequence`; `rev_comp` indicates that calls, alignments,
and coordinates address its reverse complement. See the
[output-format documentation](https://abstar.readthedocs.io/en/latest/output_formats.html)
for the official sequence and amino-acid compatibility boundary.

## Test

Install the development test environment:

```bash
python -m pip install -r requirements-test.txt
```

The ordinary fast gate and the complete suite are:

```bash
python -m pytest -m "not integration and not e2e and not slow" -q
python -m pytest -q
```

The integration and end-to-end gate uses real bundled executables and databases:

```bash
python -m pytest -m "integration or e2e" -q
```

Contributor commands, coverage floors, focused AIRR/database gates, and the
optional external-corpus discovery command are documented in `AGENTS.md`. The
bulk published corpus is never required by ordinary CI.
