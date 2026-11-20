# abstar rewrite — foundation design

Status: draft for approval · 2026-11-20

This is the first of four design documents for a ground-up rewrite of abstar.
It covers the **shared foundation**: repository and build layout, the
Rust/Python boundary, the germline database format and commands, the output
schema mechanism, the CLI and Python API shells, error/diagnostic behavior,
agent support, and performance targets. The three tools built on it
(`generate`, `run`, `number`) get their own design documents; what is already
decided about them is recorded in [Roadmap](#roadmap) so it is not lost.

## Intent

abstar keeps hitting accuracy failures and bugs that trace back to its
structure (Python orchestration around MMseqs and Parasail subprocesses,
temp-file chunking, multiple hand-synchronized schema/serializer/doc
surfaces). Patching will not get it to where it needs to be. The goal is a
first-principles rebuild that is **rock solid, completely reliable, and
hyper-efficient**.

Priorities, in order:

1. **Accuracy** — within the limit set by genuine biological ambiguity
   ("bf16 rounding error"), never chasing impossible perfection.
2. **Speed and memory efficiency** — a very close second.
3. **Maintainability** — one source of truth per concept.

abstar is also the first step of a wider abutils rework. abstar is designed on
its own terms; **abutils adapts to abstar, not the reverse.** Breaking current
callers (including abutils) is acceptable and encouraged when it yields a
cleaner abstar.

**Agents are first-class users.** Everything an agent needs to use abstar
correctly must be discoverable from the tool itself.

## Goals

- A Python package `abstar` with a CLI (`abstar run …`) and an importable API
  (`abstar.run(...)`), approximately compatible with today's surface.
- A Rust engine that owns every hot path; no external binaries, no temp files
  on the hot path.
- Output in AIRR-compliant TSV and Parquet with identical content. AIRR
  required fields plus abstar extension fields where AIRR has gaps (e.g.
  `sequence_vdj`).
- A germline database format and `abstar germline` commands that support
  built-in, mixed-species, new-species, and donor-specific databases.
- A single output-field definition that drives the Rust writers, the manual,
  and the docs.
- `abstar manual` and agent-friendly CLI behavior.
- A structure in which the synthetic generator **cannot** share annotation
  logic with the annotator.

## Non-goals

- Read merging and UMI parsing (dropped for now; candidates for later
  `abstar merge` / `abstar umi` commands).
- Error correction of long/noisy reads. They annotate as heavily mutated,
  indel-rich sequences.
- 100% accuracy on genuinely ambiguous cases.
- Compatibility with abutils or `abutils.Sequence` return types.
- The `abstar.gl` / `abstar.pp` / `abstar.tl` namespaces.
- A compiled/binary germline database format (add only if load time is
  measured to matter).
- A Python-free standalone binary (possible later on the same Rust crates).
- An MCP server or bundled agent skill (thin layers over `abstar manual`,
  addable later).
- Tested Windows support (wheels are built; it is not a CI test target).

## Anti-goals (failure even if it "works")

- Plausible but silently wrong annotations. An explicit failure is always
  better.
- A generator that reuses annotator logic, so that shared bugs score as
  correct.
- Accuracy benchmarks that only use our own generator (overfitting to our own
  idiosyncrasies).
- A manual, `--help`, or docs that drift from the code.
- Return shapes that depend on the input data or how many records succeeded.
- Fixed startup cost that makes small jobs feel slow.
- Output that changes with thread count or batch size.

## Constraints

- Python ≥ 3.11 (3.10 reached end of life in October 2026). *My call.*
- Rust core exposed via PyO3, built with maturin into a single `abstar` wheel.
- CLI built with **Typer** (user decision).
- Dependency versions: lower bounds only, unless a known incompatibility
  forces a cap.
- Germline sources: OGRDB for V/D/J wherever OGRDB covers the locus;
  otherwise the best curated alternative (e.g. IMGT GENE-DB, KIMDB). IMGT
  for constant regions. The source is recorded per locus in the manifest,
  and a locus switches to OGRDB when OGRDB adds it. *(User decision.)*
- Built-in databases: `human` (IG from OGRDB + TR from IMGT), `macaque`
  (IGH from KIMDB, IGK/IGL from IMGT), `c57bl6`, `balbc`, `human+c57bl6`.
  Today's `balbc` and `c57bl6` manifests cite the *Homo sapiens* OGRDB URL;
  fix the provenance when rebuilding.
- AIRR: target the latest released AIRR Rearrangement schema at
  implementation time; pin and record the version in the schema file.
- Platforms: prebuilt wheels and CI tests for Linux x86_64, Linux aarch64,
  and macOS arm64; Windows wheels built but untested. The Rust core uses the
  fast vector instructions each platform offers (AVX2 on x86_64, NEON on
  ARM), selected at runtime, with a portable fallback.
- License: MIT (unchanged).
- The correctness rules in the current `AGENTS.md` carry forward:
  - exact ID preservation
  - explicit coordinate spaces converted at named boundaries
  - codon-boundary region rules
  - receptor/locus propagation
  - deterministic tie-breaking with retained ties
  - record conservation

  The code does not carry forward.

## Approaches considered

| Approach | Summary | Outcome |
|---|---|---|
| **A. Python shell, Rust engine** | Typer CLI and Python API. Rust owns parsing, threading, annotation, and writing; Python passes paths/options and never touches per-sequence data on the bulk path. | **Chosen.** One language defines the CLI and API (avoiding the sync problems we have today). Python contributors can add commands, and Python-ecosystem commands fit naturally. The ~100 ms startup cost is negligible next to the MMseqs startup cost it replaces. Speed is identical because all heavy work is in Rust. |
| B. Rust all the way down | clap CLI in Rust shipped in the wheel and standalone; Python is only an API binding. | Rejected. Millisecond startup and a Python-free binary, but the CLI and API are defined twice in two languages. A higher bar for lab contributors; awkward for Python-ecosystem commands. Could be added later over the same crates. |
| C. Rust binary called as subprocess | Python shells out to a Rust executable. | Dropped. It is the MMseqs pattern being abandoned: IPC/serialization cost, no in-memory DataFrames, startup and temp-file overhead. |

## Design

### Repository and build layout

The rewrite happens **in this repository on a long-lived branch**. The current
implementation stays in place as a readable reference and is deleted when the
new version reaches parity. `AGENTS.md` is rewritten for the new structure at
that point, keeping the correctness rules.

```
Cargo.toml               workspace
crates/
  germline/              DB format, loading, validation (data only)
  generate/              recombination simulator      → depends on germline ONLY
  annotate/              search, alignment, annotation, writers → depends on germline
  number/                (later)
  py/                    PyO3 bindings: exposes coarse-grained entry points,
                         returns Arrow data to Polars without copying
python/abstar/           Typer CLI, Python API, manual rendering
schema/fields.toml       single source of truth for output fields
pyproject.toml           maturin mixed Rust/Python project → one wheel
```

**Hard rule:** `generate` must not depend on `annotate` or on any alignment
or annotation code. This is enforced by the crate dependency graph, not by
convention. Code that needs alignment, such as IMGT-gapping new-species
germlines during database builds, lives outside `germline` so that
`generate` stays alignment-free. *Exactly where it lives is the builder's
call*, as long as the rule holds.

### Rust/Python boundary

- Python calls a small number of coarse-grained Rust entry points, e.g.
  "annotate these inputs into these outputs with these options" and
  "annotate these in-memory sequences and return Arrow record batches".
- Rust releases the GIL (Python's global interpreter lock) and does
  everything: FASTA/FASTQ parsing, batching, multithreaded work, ordered
  reassembly, TSV/Parquet writing.
- Inputs are streamed in batches, so memory stays roughly flat from 100 to
  10⁹ sequences.
- Output order equals input order regardless of thread count or batch size.
- In-memory results cross into Python as Arrow data and become a Polars
  DataFrame without copying.
- The CLI path never imports Polars. Typer, the extension module, and Rich
  are imported lazily so that `abstar --help` and small runs start fast.

### Germline databases

**Format** — a database is a directory with:

- `manifest.toml` — database name, format version, receptor(s), species,
  provenance (source, release, date), and a content hash.
- `alleles.tsv` — one row per allele. Required columns:
  - allele ID (an opaque string, never parsed)
  - locus and segment (V/D/J/C)
  - species
  - functionality
  - IMGT-gapped sequence
  - V and J anchor positions (V Cys104; J Trp118 or Phe118)
  - source

  The ungapped sequence is derived by stripping gaps. Exact columns are the
  builder's call.

Plain text is deliberate: it is diffable, readable by agents, and the search
indexes build from it at load time in milliseconds. A cached or compiled
format is added only if load time is measured to matter.

**A database is a set of loci, not a receptor.** A database such as `human`
contains whichever loci it has, both IG and TR. `--receptor`/`--locus`
becomes an optional filter rather than a required selector. *My call*,
following the open TODO in today's `gl.py` about combined databases. It
removes a required argument and lets mixed BCR/TCR inputs work. IG and TR V
genes are far too divergent for cross-receptor confusion to be a practical
risk.

**Mixed species** *(user decision)*:

- **Inside a database, allele IDs are made unique with a short species
  prefix** (e.g. `hs|IGHV1-2*01`, `mm|IGHV1-2*01`) in mixed-species
  databases. Single-species databases use plain names.
- **Output never shows the prefix.** Call fields contain standard allele
  names only, so they stay compatible with AIRR tools.
- **Species is reported per segment** (`v_species`, `v_species_ties`, etc.;
  see [Gene calls and ties](#gene-calls-and-ties)). A tie across species with
  identical names reads `v_call_ties = [IGHV1-2*01, IGHV1-2*01]` and
  `v_species_ties = [human, mouse]`, with nothing deduplicated.
- Per-segment species is needed, rather than a single `species` column,
  because humanized-mouse sequences mix species within one rearrangement
  (e.g. human V/D/J with mouse C). *My call.*

**Gene grouping** — `alleles.tsv` has a required `gene` column, used for
gene-level calls and ties. At build time it is filled from IMGT-style names
when they parse, otherwise from a user-supplied mapping, and otherwise each
allele is its own gene. Allele IDs are parsed only at build time, never
during annotation. *My call.*

**Locations:**

- Built-in databases ship inside the package.
- User databases live under `~/.abstar/germline/<name>/`.
- Building a user database whose name collides with a built-in one is an
  error. Today's behavior, where a user database silently overrides a
  built-in one, is a footgun. *My call.*

**Build modes** (`abstar germline build`):

1. From OGRDB germline-set JSON (V/D/J) plus IMGT constant-region FASTA.
   This is how the built-in databases are made.
2. From FASTA for new species (e.g. bovine). Ungapped V sequences are
   IMGT-gapped by alignment to a reference species, as today; J anchors are
   found by motif.
3. Donor-specific: a base database, plus a keep-list of alleles (the donor's
   genotype), plus extra FASTA of novel alleles (e.g. from IgDiscover or
   TIgGER). The keep-list restricts only the locus/segment combinations it
   names, e.g. a genotype listing IGHV alleles restricts IGHV and leaves
   IGHD, IGHJ, light chains, and C untouched. The build prints what was
   restricted and what passed through. *My call.*
4. Mixed: any combination of the above sources in a single build.

Mode 1 and already-gapped inputs belong to the foundation. Modes that need
alignment (gapping ungapped V genes) land with the annotator's aligner.

**Build validation:**

- IDs are unique.
- Gapped and ungapped sequences are consistent.
- Locus and segment are present.
- Anchors are detectable. Every allele whose anchor cannot be found is
  reported by name; nothing fails silently later.
- Builds are staged and appear atomically, only when every step succeeds.

A donor database should build in seconds.

**Commands:**

- `abstar germline build`
- `abstar germline list`: built-in and user databases, readable with dozens
  of donor databases, `--json` available
- `abstar germline info <name>`: manifest and allele counts by
  locus/segment/species

**Separate from the database:** the generator's recombination model (allele
usage, trimming, N-addition, somatic mutation) lives in its own file. See the
Roadmap.

### Output schema

`schema/fields.toml` defines every output field: name, Arrow type,
nullability, description, coordinate space (where relevant), origin (`airr`
required, `airr` optional, or `abstar` extension), and the AIRR schema
version targeted. This one file drives:

- the column set and types in the Rust TSV and Parquet writers (via code
  generation or load-time validation — *builder's call*);
- `abstar manual fields`;
- the Sphinx field reference.

The foundation supplies the mechanism and the AIRR-required fields. The full
field list (including which current extension fields survive) is decided in
the annotator design.

### Gene calls and ties

*User decision.* Applies identically to V, D, J, and C (shown for V):

| Field | Contents |
|---|---|
| `v_call` | Exactly one allele, always. Searchable, countable, sortable. |
| `v_gene` | The gene of `v_call`. |
| `v_species` | The species of `v_call`. |
| `v_call_ties` | Every allele with the top score, including `v_call`. A native list in Parquet; comma-joined in TSV. |
| `v_species_ties` | Aligned entry by entry with `v_call_ties`. |

- Tied alleles share the top score by definition, so `v_score` covers them
  all. There is no score list.
- There is no ambiguity boolean: it is `len(v_call_ties) > 1`. Add one only
  if TSV users ask.
- The names deliberately avoid `v_calls`, which is one letter from `v_call`
  and easy to confuse.

**Primary-call rule.** `v_call` is the most probable allele among the ties,
chosen by these deterministic tie-breakers in order:

1. Best alignment score.
2. Functional before ORF before pseudogene.
3. Population usage frequency, if the database provides it (optional
   `usage` column in `alleles.tsv`).
4. A fixed natural sort order on allele name.

The same input and database always give the same call. The known bias is
small and runs in the least harmful direction: allele counts within tie sets
lean toward common, functional alleles. Rigorous analyses use `v_call_ties`.

**Evaluation consequences:**

- A call is correct if `v_call` is in the generator's set of equally valid
  answers. Coverage of that set by `v_call_ties` is scored separately.
- Because tie-break 3 shares usage data with the generator, evaluation must
  also include generated sets with **uniform** allele usage, so that
  agreement between the two cannot inflate accuracy.

Encoding rules:

- TSV follows AIRR conventions: `T`/`F` booleans, empty cells for nulls,
  1-based closed coordinates.
- Parquet carries identical content with native types (bool, int, null).
- Both writers read the same in-memory record, so there is no second
  conversion path.

### CLI (Typer)

```
abstar run SEQUENCES [OUTPUT_DIR]       # SEQUENCES: file, directory, or "-" (stdin)
abstar run --seq CAGGTGCAG...           # single inline sequence → TSV on stdout
abstar germline build|list|info
abstar manual [TOPIC]                   # e.g. `fields`, `fields junction_aa`, `run`
# later: abstar generate, abstar number, abstar merge, abstar umi
```

`run` options (initial set): `--germline-db/-g` (default `human`),
`--locus` (optional filter), `--format airr|parquet` (repeatable), `--threads`,
`--strict`, `--json` (summary output), `--quiet`/`--verbose`.

Changes from today's options:

| Today | New behavior |
|---|---|
| `--merge*`, `--umi*`, `--interleaved_fastq` | Removed |
| `--mmseqs_chunksize`, `--mmseqs_threads`, `--chunksize` | Removed: no MMseqs, and batching is internal |
| `--n_processes` | Replaced by `--threads` |
| `--copy_inputs_to_project`, `--concise_logging`, `--debug` | Removed |
| `PROJECT_PATH` | Becomes `OUTPUT_DIR` (still positional) |

A migration note lists every removed option.

**Agent-oriented behavior:**

- **Output streams:** results go to stdout only when writing to stdout;
  logs and progress always go to stderr. Progress bars appear only when
  stderr is a terminal.
- **No interactive prompts**, ever.
- **`--json`** on summaries, `germline list`/`info`, and errors.
- **Exit codes:** stable and documented (success, invalid usage/input, fatal
  run error, strict-mode record failure — exact values are the builder's
  call).
- **Error codes:** stable, machine-readable codes for each failure class.

### `abstar manual`

- `abstar manual` prints the complete reference as one Markdown document:
  - every command and option
  - output fields
  - coordinate conventions
  - ambiguity reporting
  - exit and error codes
  - examples
- `abstar manual <topic>` narrows it.
- It is generated from the Typer command definitions and
  `schema/fields.toml`, never hand-written, so it cannot drift.
- Prose that cannot be generated (e.g. coordinate conventions) lives in one
  Markdown source that both the manual and the Sphinx docs include.

### Python API

```python
abstar.run(
    sequences,              # path | dir | list of paths | str sequence | iterable of str or (id, seq)
    output_dir=None,
    *,
    germline_db="human",
    locus=None,
    formats=("airr",),
    threads=None,
    strict=False,
) -> pl.DataFrame | RunSummary
```

Exactly two return shapes, chosen by the caller, never by the data:

- `output_dir=None`: always a Polars DataFrame, including for one sequence
  or zero. A single record is `df.row(0, named=True)`.
- `output_dir` given: streams to disk and returns a small `RunSummary`
  (counts annotated/failed/unassigned, output paths, run metadata path).

A string argument is treated as a file path if it exists. Otherwise it is
validated as a nucleotide sequence, and anything else is an error. Germline
functions are exposed as `abstar.germline.build/list/info`.

### Record status and diagnostics

**Every input record produces exactly one output row.** Records are never
silently dropped.

- `annotation_status` and `failure_reason` columns distinguish:
  - annotated
  - biologically unassigned (e.g. no V found)
  - invalid input (e.g. non-nucleotide characters, empty)
  - per-record internal error (with message)
- Record conservation becomes a row count check.
- This **replaces** today's omit-failed-rows-plus-diagnostic-files approach.
  *My call*: it is simpler, it is directly visible to agents, and it removes
  the `logs/` per-record file tree.
- `--strict` turns any per-record internal error into a fatal run error.
- Fatal run errors (unreadable input, bad database, I/O failure) stop the
  run with a documented exit code and error code. Workers never fail
  silently.
- When writing to disk, a `run.json` beside the outputs records:
  - parameters
  - abstar version
  - database name and content hash
  - AIRR schema version
  - counts by status

### Performance targets

These are measured, not assumed. A throwaway search and alignment spike
confirms them at the **start** of the work (see
[Milestones](#milestones-and-order)):

| Workload | Target |
|---|---|
| `abstar --help` | ≤ 150 ms |
| `abstar run --seq …` (single sequence, end to end) | ≤ 300 ms |
| Tens to hundreds of sequences on a laptop | Effectively instant |
| ≥ 10k sequences/sec on an M-series MacBook Pro | ~1M sequences in about 2 minutes |
| Billions of sequences on a large server | Scales near-linearly with cores; flat memory |

A benchmark harness (Rust micro-benchmarks plus an end-to-end script over
fixed inputs) lands with the foundation, so regressions are visible from the
first commit.

### Multiple inputs

- Each input file produces its own output file per format, named after the
  input file's stem (as today).
- There is one `run.json` per run, with counts per input.
- In memory, multiple inputs come back as one DataFrame in input order. It
  always has a `source_file` extension column, null for sequences passed
  directly, so the columns never vary.
- *My call.*

### CI and release

- Maturin builds wheels for the platforms above; Rust tests and pytest run
  on the tested platforms.
- Existing workflows (tests, corpus, coverage floors, PyPI publish) keep
  running on `master`, unchanged. The rewrite branch gets its own CI.
- No pre-releases until the annotator produces real annotations.
- **Parity**, which triggers deleting the old code and releasing 1.0.0,
  means all three of:
  1. The new annotator meets the accuracy targets set in its design doc, on
     generator data and external benchmarks.
  2. Every difference from the old abstar on the committed BCR corpus has
     been reviewed, and the new call is judged at least as good.
  3. Every retained command works.
- *My calls.*

### Milestones and order

1. **Performance spike (throwaway):** a k-mer prefilter plus vectorized
   alignment against the human V set, to confirm the throughput targets
   before anything is built on them. If they are unreachable, this design
   is revisited.
2. **Foundation:**
   - the workspace and wheel build;
   - `germline build` for the gapped-input modes, plus `list` and `info`;
   - `fields.toml` and the writers;
   - the CLI, API, and manual;
   - `run.json`, status and exit codes.

   At this stage `abstar run` drives the full pipeline (parse → status →
   writers → `run.json` → exit codes), with every row marked unassigned and
   `failure_reason` set to `engine_not_implemented`. This is branch-only and
   never released. *My call*, so the pipeline is tested end to end from the
   start.
3. `generate` → 4. `run` (annotator, including alignment-dependent
   database builds and the corpus comparison runner) → 5. `number`.

## Roadmap

### Project 2 — `abstar generate` (next design doc)

- Built before the annotator, as the ground-truth source for test-driven
  development of annotation.
- It follows the rules of recombination:
  - allele choice
  - exonuclease trimming
  - P- and N-nucleotide addition
  - D-gene inversion where relevant
  - somatic hypermutation with hotspot bias
  - indels
  - truncation
- For each sequence it records the **set of equally valid answers**, not a
  single truth (e.g. every allele equidistant after mutation, ambiguous trim
  and N boundaries). An annotation counts as correct if it lands in that set.
- The recombination model is a separate file. Human ships as the default
  (drawn from literature / IGoR-style models). Other species may borrow it
  with a documented caveat, or supply their own.
- It is also a user-facing tool in its own right.

### Project 3 — `abstar run`, the annotator

- Purpose-built search and alignment in Rust, replacing MMseqs and Parasail.
  Germline sets are a few hundred alleles, so exhaustive SIMD (vectorized)
  alignment after a k-mer prefilter is feasible.
- Must handle large indels without falling over (keeps long reads viable).
- Ties are reported as in [Gene calls and ties](#gene-calls-and-ties): a
  single `v_call` plus a `v_call_ties` list.
- Accuracy is evaluated against our generator **and** external tools and
  datasets (e.g. IgBLAST, partis's simulator, IGoR/OLGA, real corpora) to
  avoid overfitting.
- Decides the full output field list.

### Project 4 — `abstar number`

- Antibody numbering. Its own design doc after the annotator.

### Later

- `abstar merge`, `abstar umi`.
- An optional MCP server or agent skill over `abstar manual`.
- A standalone Rust binary.

## Left to the builder

- Crate names and exact module boundaries, within the `generate` dependency
  rule.
- The library for passing Arrow data from Rust to Polars, the FASTA/FASTQ
  parser, and the threading library.
- Whether `fields.toml` drives the Rust code by code generation or by
  load-time validation.
- Exact `alleles.tsv` columns, exit code values, and error code names.
- Batch sizes and buffering.
