`abstar generate` simulates V(D)J reads from a recombination model and a germline database. With each read it writes a truth record, so an annotator (abstar or any other) can be graded by joining on `sequence_id`. `abstar manual truth-fields` defines every truth column.

## Quick start

```sh
abstar generate --n 1000 --locus IGH --seed 1 > reads.fasta  # FASTA only; summary on stderr
abstar generate out/ --n 1000 --locus IGH --seed 1           # out/generated.fasta, .truth.tsv, generate.json
abstar generate out/ --n 1000 --locus IGH --preset bnab --format parquet
```

```python
df = abstar.generate(1000, "IGH", seed=1)                  # truth table as a Polars DataFrame
s = abstar.generate(1000, "IGH", seed=1, output_dir="out") # GenerateSummary
```

- Without `--seed`, a seed is drawn and printed in the summary (and recorded in `generate.json`). When FASTA streams to stdout, it is also printed to stderr first as `seed: N`, even with `-q`. The same seed gives identical output at any thread count.
- Rerunning into the same OUTPUT_DIR overwrites its files. `generate.json` is deleted first and written last, so a failed rerun cannot look finished.

## Models

A recombination model gives gene usage, P(D|J), trims, N lengths and composition, and whether reads mutate (`shm_model`). `--model` defaults to the built-in `<species>-<locus>` model for the database's species, such as `human-igh`. A database with more than one species (`human+c57bl6`) needs `--model`.

- Built-in: `human-igh`, `human-igk`, `human-igl`, `human-tra`, `human-trb` (decision 0027). Their values come from published figures, cited in each model's sources (`abstar model info human-igh`).
- User models: `~/.abstar/models/<name>/model.toml` (`$ABSTAR_HOME/models/`). A user model may not use a built-in name (`BUILTIN_NAME_CONFLICT`, decision 0010).
- `--model` also takes a path to a `model.toml` or a directory holding one. The format is in decision 0023.
- `abstar model list` and `abstar model info NAME` show them (Python: `abstar.model.list()`, `abstar.model.info(name)`). A user model that fails to load is still listed, with its error.

A model that does not fit the database (less than 90 % of its V, D or J mass maps to database alleles) is `MODEL_INVALID`. Genes it drops are listed as warnings.

## Presets

`--preset` sets starting parameters; any option given explicitly overrides the preset's value.

| Preset | Mutation | Long CDR3 weight | Indels | Productive only |
|---|---|---|---|---|
| (none), `naive`, `tcr` | 0 | 0 | 0 | no |
| `memory` | `uniform:0.02-0.12` | 0 | 0 | no |
| `bnab` | `uniform:0.15-0.35` | 0.08 | 0.110 | yes |

`bnab`'s long-CDR3 weight draws each junction's N length from the IGH model's long tail (20-60 nt) with probability 0.08. On `human-igh` it gives about 28% of reads with CDRH3 ≥ 24 aa, 2.4% with ≥ 35 aa and 10% with an indel, at about 6 attempts per read. These numbers are for productive output, which `bnab` turns on; `--no-productive-only` turns it off.

## Options

- `--mutation`, `--truncate-5p`, `--truncate-3p` take a number or a distribution drawn per read: `fixed:X`, `uniform:LO-HI` (`uniform:1e-3-0.1` works), or for mutation only `beta:A,B`. Mutation is the realized V divergence (0-1); truncations are nucleotides cut from each end of the read, so a beta draw (always under 1) is `INVALID_ARGUMENT` for them.
- `--long-cdr3` is the per-junction probability of drawing the N length from the model's long tail. `--indels` is the per-read probability of one indel (`--frameshift-indels` allows lengths that shift the frame).
- `--constant INT` appends that many nucleotides of a constant region, from the C allele's second base: its first is the J's last base, as in real transcripts, so truth's `c_germline_start` is 2 (Python: `isotype_weights` weights the C genes).
- `--productive-only` keeps only productive reads (`--no-productive-only` turns off a preset's); `attempts` in the summary counts the reads simulated. It also turns on selection during mutation, as in real antibodies (decision 0024). A rearrangement must be productive before it is mutated (mutation never repairs one), and a substitution is rejected and redrawn if it would create an in-frame stop codon, or make the Cys104 or Trp/Phe118 codon non-conserved. Without it, mutation is neutral.
- Mutation uses abstar's built-in motif model (decision 0024). `--shm-table DIR` uses S5F tables instead (`Mutability.csv`, `Substitution.csv`), which you download yourself; abstar does not bundle them.
- Mutating with a model whose `shm_model` is `none` (the TCR models) is `INVALID_ARGUMENT`. To mutate TCRs for testing, pass `shm_model="motif"` from Python, use a model file with `shm_model = "motif"`, or pass `--shm-table DIR`.

## Grading against truth

The truth table holds the true generative values under the same AIRR names as `abstar run` (`v_call`, `junction`, `v_sequence_end`, …). It also holds, for each ambiguous field, a `_valid` list of every equally correct answer (decision 0025):

- Allele and gene calls are correct when they are in `v_call_valid`, `d_call_valid`, `j_call_valid`, … These are the alleles tied for the fewest mismatches against the read. A null D call is also correct when `d_call_null_valid` is T.
- Boundaries (`v_sequence_end`, `d_sequence_start`, `d_sequence_end`, `j_sequence_start`) are correct when they are in their `_valid` list: every cut tied for the best matches-minus-mismatches score.
- `unrecoverable_v`, `unrecoverable_d` and `unrecoverable_j` mark reads where mutation made another allele strictly closer than the true one. No annotator can call these right, so they set the error floor.
- The difficulty flags (`high_mutation`, `long_cdr3`, `has_indel`, `truncated_5p`, `truncated_3p`, `ambiguous_v`, `ambiguous_d`, `ambiguous_j`) let you break accuracy down by case.
- Indels record their range of equivalent positions; a placement anywhere in it is correct.

`abstar.grade(truth, calls)` applies these rules and returns one row per field and subset (`field, subset, n, correct, accuracy`); the subsets are `all` and each difficulty flag. `calls` is any table with `sequence_id` and AIRR columns such as `abstar run` output; a field it lacks is reported with `n = 0`.

```python
import polars as pl
abstar.generate(10_000, "IGH", seed=1, preset="bnab", output_dir="bench", formats=["parquet"])
truth = pl.read_parquet("bench/generated.truth.parquet")
abstar.grade(truth, abstar.run("bench/generated.fasta"))
```

Grading rules:
- A segment's fields count only for reads where that segment is recoverable (`not unrecoverable_<seg>`). An empty `_valid` list means the read holds no base of that segment, so only a null call (or null boundary) is correct; a null `d_call` is also correct when `d_call_null_valid` is T.
- `junction`, `junction_aa`, `productive`, `vj_in_frame` and `stop_codon` are exact matches, scored on reads where V and J are recoverable.
- Regions, masks, `<seg>_mutations*` and indels are scored only where the primary call equals the true allele (regions: V, plus J for CDR3 and FWR4; masks: every segment present).
- Mutation and indel lists are compared after restricting both to the germline positions inside the intersection of the two `<seg>_germline_start`..`end` spans. Counts and frequencies are recomputed from the restricted lists. V positions are mapped through the IMGT gapping of `germline_db` (default `human`; pass `germline_db=` for another database).
- Masks cover only the V(D)J region (nt masks index `sequence_vdj`, aa masks `sequence_vdj_aa`; no constant-region bases). They are compared position by position where both `gene_segment_mask`s agree; an aa mask, at codons whose three bases have the same `gene_segment_mask` label in both. `gene_segment_mask` itself may differ only at read positions the `_valid` boundary lists leave open.

In the TSV truth, list columns are comma-joined strings; the Parquet truth keeps them as lists.

## Known limits

- V alleles are compared at their IMGT-gapped positions. A V allele with shifted gapping, such as human IGLV3-2\*03, is compared at the wrong positions and may be missing from or under-counted in truth sets (#138).
- Macaque is not shipped until its V genes are IMGT-gapped (#137). There are no built-in models for TRG, TRD or non-human species.
- The motif mutation model is a heuristic, not fit to data (decision 0024).
- The light-chain trims and every TCR and light-chain N composition are fallbacks borrowed from IGH (decision 0027).
- Models are generation (pre-selection) distributions; `--productive-only` filters reads but does not reweight genes toward a productive repertoire.
- Importing IGoR model files is not supported yet.
