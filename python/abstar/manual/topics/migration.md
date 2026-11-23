How abstar 0.x options map to this version. Every option of the old `abstar run` and `abstar build_germline_database` commands and of the old `abstar.run()` is listed, with what replaces it or "removed". `abstar manual examples` shows the new commands in use.

## `abstar run`

| Old | Now |
|---|---|
| `INPUT_PATH` | `SEQUENCES` (still positional): a FASTA/FASTQ file, a directory, or `-` for stdin. Use `--seq` for a single sequence. |
| `PROJECT_PATH` | `OUTPUT_DIR` (still positional, now optional). Outputs go directly in it (`<stem>.tsv`, `<stem>.parquet`, `run.json`); there are no `tmp/` or log subdirectories. Without it, TSV streams to stdout. |
| `--germline_database` | `--germline-db` / `-g` (database name or path; default `human`). |
| `--receptor` | Removed. A database holds both IG and TR loci (decision 0005); `--locus` is an optional filter. |
| `-o`, `--output_format` | `--format` (repeatable; `airr` or `parquet`). `-o` is gone. |
| `--umi_pattern` | Removed (candidate for a later `abstar umi` command). |
| `--umi_length` | Removed (candidate for a later `abstar umi` command). |
| `-m`, `--merge` | Removed (candidate for a later `abstar merge` command). Merge reads before running abstar. |
| `--merge_kwargs` | Removed, with `--merge`. |
| `--interleaved_fastq` | Removed, with `--merge`. |
| `-c`, `--chunksize` | Removed: batching is internal. |
| `--mmseqs_chunksize` | Removed: there is no MMseqs2 (decision 0002). |
| `--mmseqs_threads` | Removed: there is no MMseqs2. Use `--threads`. |
| `--n_processes` | `--threads` (worker threads, at least 1; default all cores). |
| `--copy-inputs` / `--no-copy-inputs` | Removed: inputs are never copied. |
| `--verbose` / `--quiet` | `--verbose` / `-v` prints the database and inputs to stderr; `--quiet` / `-q` suppresses the summary. The default prints only the summary. |
| `--strict` | `--strict` (unchanged meaning; a record error now exits 3 with `STRICT_RECORD_FAILURE`). |
| `--debug` | Removed: there are no temporary files to keep. Errors carry stable codes (`abstar manual errors`). |
| `--started_from_cli` (hidden) | Removed. |

New: `--seq`, `--json`, `--locus`, `--include-nonfunctional` (non-functional alleles are no longer searched by default; decision 0038).

## `abstar build_germline_database`

Replaced by `abstar germline build RECIPE`. A TOML recipe file describes the sources (decision 0017); `abstar manual examples` has one.

| Old | Now |
|---|---|
| `NAME` | The recipe's `name`. |
| `-f`, `--fasta` | Recipe source `v_fasta`, `d_fasta`, `j_fasta` (the segment is given by the key). |
| `-j`, `--json` | Recipe source `ogrdb_json`. (`--json` on `abstar germline build` now means machine-readable output.) |
| `-c`, `--constant` | Recipe source `c_fasta`. |
| `-m`, `--manifest` | Recipe `description` and per-source `description`. The built `manifest.toml` records provenance, file hashes, build date and warnings. |
| `-r`, `--receptor` | Removed: loci are read from the alleles. |
| `-l`, `--location` | `--output DIR` (default: `~/.abstar/germline/<name>`, or `$ABSTAR_HOME/germline/<name>`). |
| `--reference` | Recipe source key `gap_v = true` (templates via `gap_templates`; default `human`, `c57bl6`, `balbc`), or `abstar germline gap` on a FASTA. V sequences are IMGT-gapped by protein alignment to the closest template allele (decision 0037). |
| `--include_species_in_name` / `--exclude_species_from_name` | Recipe `prefixes` (species to ID prefix), used when sources mix species (decision 0007). |
| `--verbose` / `--quiet` | Removed: the build prints a one-line summary, keep-list results and warnings; `--json` for machine-readable output. |
| `--debug` | Removed: no external tools are run. |

New: `--overwrite` (replaces an existing database only), `abstar germline list`, `abstar germline info`.

## `abstar.run()`

| Old | Now |
|---|---|
| `sequences` | `sequences`: a path, a directory, a sequence string, or an iterable of paths, sequence strings or `(id, sequence)` pairs. `abutils.Sequence` objects are not accepted; pass `(id, sequence)` pairs. |
| `project_path` | `output_dir`. |
| `germline_database` | `germline_db`. |
| `receptor` | Removed; `locus` is an optional filter. |
| `output_format` | `formats`. |
| `umi_pattern`, `umi_length` | Removed. |
| `merge`, `merge_kwargs`, `interleaved_fastq` | Removed. |
| `chunksize`, `mmseqs_chunksize`, `mmseqs_threads` | Removed. |
| `n_processes` | `threads`. |
| `copy_inputs_to_project` | Removed. |
| `verbose`, `concise_logging`, `debug` | Removed: `run()` prints nothing. |
| `as_dataframe` | Removed: without `output_dir`, `run()` always returns a Polars DataFrame; with it, a `RunSummary` (decision 0003). |
| `started_from_cli` | Removed. |
| `strict` | `strict` (unchanged meaning). |

## Output

- Failed records are no longer dropped into `logs/`: every input record is one output row, with `annotation_status` and `failure_reason` (`abstar manual conventions`).
- `logs/run.json` is now `<OUTPUT_DIR>/run.json`.
- Field definitions: `abstar manual fields`.
