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
abstar run reads.fasta out/ -g c57bl6 --locus IGH --locus IGK --threads 8

# Machine-readable summary; fail with exit 3 on any record error
abstar run reads.fasta out/ --json --strict

# Protein input (one residue per letter); sequence_aa holds the input, sequence is null
abstar run antibodies.fasta out/ --sequence-type aa
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

Gap ungapped V alleles (e.g. novel alleles) into IMGT numbering by transfer from the best template:

```sh
abstar germline gap novel_v.fasta --output novel_v.gapped.fasta --report gap.tsv
abstar germline gap novel_v.fasta --templates human --templates c57bl6 > gapped.fasta
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

Add `gap_v = true` to a `species` source to IMGT-gap its ungapped V alleles during the build
(templates default to human, c57bl6, balbc; override with `gap_templates = ["human"]`). V alleles
that cannot be gapped are dropped with a build warning:

```toml
[[source]]
species = "human"
gap_v = true
v_fasta = ["novel_v.fasta"]
```

## Synthetic reads with truth

```sh
abstar generate --n 1000 --locus IGH --seed 1 > reads.fasta    # FASTA on stdout, summary on stderr
abstar generate bench/ --n 10000 --locus IGH --preset bnab --format airr --format parquet
# -> bench/generated.fasta, bench/generated.truth.tsv/.parquet, bench/generate.json
abstar generate bench/ --n 1000 --locus IGK --mutation uniform:0-0.1 --truncate-5p uniform:0-30 --productive-only
abstar model list                    # built-in, then user recombination models
abstar model info human-igh          # description and cited sources
```

`abstar manual generate` covers presets and grading; `abstar manual truth-fields` defines every truth column.

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
gapped, report = abstar.germline.gap("ungapped_v.fasta", templates=["human"])  # IMGT-gap V alleles; report is a DataFrame

truth = abstar.generate(1000, "IGH", seed=1)        # truth table as a Polars DataFrame
s = abstar.generate(1000, "IGH", seed=1, preset="bnab", output_dir="bench/", formats=["airr", "parquet"])
s.fasta, s.truth, s.generate_json                   # paths; s.seed, s.n, s.productive, s.attempts
abstar.model.list()                                 # recombination models, like abstar.germline.list()
abstar.model.info("human-igh")
```
