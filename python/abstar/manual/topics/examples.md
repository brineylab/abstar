Common invocations. `abstar manual run`, `abstar manual germline` and `abstar manual commands` list every option.

## Annotate

```sh
# One FASTA/FASTQ file (gzip is detected automatically) to AIRR TSV files
abstar run reads.fasta out/
# -> out/reads.tsv and out/run.json; a summary on stdout

# A directory of FASTA/FASTQ files, searched recursively; TSV and Parquet
abstar run fastqs/ out/ --format airr --format parquet

# Stream TSV to stdout (no OUTPUT_DIR); the summary goes to stderr
abstar run reads.fasta > reads.tsv
abstar run reads.fasta | head

# Read stdin
zcat reads.fasta.gz | abstar run - > reads.tsv

# One sequence: TSV on stdout, or files in out/ (out/records.tsv)
abstar run --seq CAGGTGCAGCTGGTGCAGTCTGGG
abstar run --seq CAGGTGCAGCTGGTGCAGTCTGGG out/

# Another database, only some loci, a fixed thread count
abstar run reads.fasta out/ -g macaque --locus IGH --locus IGK --threads 8

# Machine-readable summary; fail with exit 3 on any record error
abstar run reads.fasta out/ --json --strict
```

SEQUENCES is always a file, a directory or `-`; use `--seq` for a sequence.

## Germline databases

```sh
abstar germline list                 # built-in, then user databases
abstar germline list --json
abstar germline info human           # manifest and allele counts
abstar germline build donor.toml     # into ~/.abstar/germline/<name>
abstar germline build donor.toml --output dbs/donor --overwrite
abstar run reads.fasta out/ -g donor # use it by name (or by path)
```

A donor recipe restricts a base database to the donor's genotype and adds novel alleles:

```toml
name = "donor"
keep = "genotype.txt"        # one allele name per line; restricts only the locus/segment pairs it names

[[source]]
database = "human"

[[source]]
species = "human"
v_fasta = ["novel_v.fasta"]  # also d_fasta, j_fasta, c_fasta, ogrdb_json
```

## The manual

```sh
abstar manual                    # everything, as Markdown
abstar manual fields             # every output field
abstar manual fields v_call_ties # one field and its decision record
abstar manual errors             # error codes and exit codes
abstar manual decisions 0008     # one decision record
```

## Python

```python
import abstar

df = abstar.run("reads.fasta")               # a Polars DataFrame
df = abstar.run("CAGGTGCAGCTGGTGCAG")        # one sequence, ID "1"
df = abstar.run([("a", "CAGGTG"), ("b", "GGCC")])  # (id, sequence) pairs keep their IDs
row = df.row(0, named=True)

summary = abstar.run("fastqs/", "out/", formats=["airr", "parquet"], germline_db="human", locus="IGH", threads=8)
summary.counts            # {"annotated": ..., "unassigned": ..., "invalid_input": ..., "error": ...}
summary.inputs[0].outputs # {"airr": Path("out/x.tsv"), "parquet": Path("out/x.parquet")}
summary.run_json          # Path("out/run.json")

try:
    abstar.run("missing.fasta")
except abstar.AbstarError as e:
    print(e.code, e.exit_code)  # INPUT_NOT_FOUND 2

abstar.germline.list()
abstar.germline.info("human")
abstar.germline.build("donor.toml")
```
