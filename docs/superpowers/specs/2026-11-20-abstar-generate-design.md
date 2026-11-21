# abstar generate — design

Status: approved · 2026-11-20

Second of the four design documents for the abstar rewrite (Project 2 in the
[foundation design](2026-11-20-abstar-foundation-design.md#project-2--abstar-generate-next-design-doc)).
It builds on the foundation: germline databases, the output-schema mechanism,
the writers, the CLI/API shells, error codes, decision records and `abstar
manual`.

## Intent

`abstar generate` produces synthetic antibody and TCR sequences by following
the rules of V(D)J recombination and somatic hypermutation, and records for
every sequence exactly what happened and which annotations are correct.

Its first and overriding job is **evaluating annotators**: above all the
rewritten `abstar run`, which will be built test-first against it. That
needs ground truth whose correct answers, and whose degree of genuine
ambiguity, are known exactly. Later it also serves **benchmarking other
annotators** and **generating ML training data**. Those uses eventually need
repertoire-scale features (clonal families, native-like heavy/light pairing,
fitted distributions), which come in later versions of the same tool.

## Goals

- Realistic single sequences for any locus that has both a germline database
  and a recombination model: IGH, IGK, IGL, TRA, TRB (TRG and TRD when a model
  exists).
- A truth record for every sequence that holds:
  - the **true generative values** (the allele actually used, the true trims,
    N/P bases, mutations and indels);
  - the **set of equally valid answers** for each ambiguous field, by the
    fewest-mismatches rule ([Truth](#truth-what-counts-as-correct)).
- Reach the hardest real cases: HIV bnAb-like sequences with about 35% V-gene
  mutation, CDRH3s of 35+ amino acids, and in-frame indels up to about
  12 codons.
- Difficulty labels on every sequence, so benchmark reports break accuracy
  down by difficulty instead of giving one blended number.
- Deterministic for a given seed; fast enough that millions of sequences take
  seconds on a laptop.
- A user-facing CLI command and Python API in the style of `abstar run`, with
  `abstar manual` coverage.

## Non-goals (this version)

- Repertoire-scale structure: clonal lineages, phylogenies, native-like
  heavy/light pairing, fitted repertoire-wide distributions. All of these are
  planned; none are in this version.
- Antigen selection on mutations. Mutation placement is neutral apart from
  context targeting.
- Sequencing-error models (Illumina, long-read). Truncation is supported; read
  errors are not.
- **Homologous-recombination grafts** (a segment of another gene inside a
  CDR). These are an annotator test case built by hand
  ([board #147](https://github.com/brineylab/abstar/issues/147)); the
  generator does not produce them.
- Fitting models from data. That is model source B, deferred until repertoire
  work.

## Anti-goals

- **Ground truth that is reconstructed after the fact**, or that is silently
  wrong. Every truth value comes from the generating events themselves.
- **Annotator logic leaking into the generator.** The foundation's hard rule
  stands: the `generate` crate depends on `germline` only, and contains no
  alignment code (see [Computing the sets](#computing-the-sets-without-alignment)).
- **Bundling GPL model files** (IGoR/OLGA) inside the MIT-licensed package.

## Constraints

Inherited from the foundation: Rust core + Python shell, Typer CLI, Python ≥
3.11, lower-bound-only dependencies, no external binaries or subprocesses, no
temp files on the hot path, AIRR 2.0 naming, the stable error/exit codes,
agent-friendly CLI behavior, and decision records for every consequential
choice.

## Approaches considered

| | Approach | Trade-off | Outcome |
|---|---|---|---|
| A | Rust `generate` crate behind `abstar generate` / `abstar.generate()`; our own model files plus an IGoR importer | Most code up front; exact truth, fast, shares the germline loader with the annotator, scales to repertoires later | **Chosen** |
| B | Pure Python/NumPy generator | Quickest to write, but ~100× slower, and a second germline-loading path that could drift from the annotator's | Rejected |
| C | Wrap IGoR / partis / immuneSIM, compute truth afterwards | Least code, but external binaries, GPL, and truth reconstructed instead of recorded | Rejected (violates anti-goals) |

## Architecture

- New crate `crates/generate` (`abstar_generate`), depending on
  `abstar_germline` and nothing from `abstar_annotate`. The crate graph
  enforces the rule.
- Output rows are built against a **truth schema** defined the same way as the
  annotation schema: a `truth.toml` in the crate, from which Rust types, the
  Arrow schema and the manual's field reference are generated. The TSV and
  Parquet writers are reused. The builder decides whether the writers move to
  a small shared crate that both `annotate` and `generate` depend on, or
  whether `generate` depends on a writers-only crate. Either way `generate`
  must not see annotation code.
- PyO3 bindings in `crates/py` expose `generate_to_files`,
  `generate_in_memory` (Arrow) and the model import; Python adds
  `abstar.generate()` and the Typer command.

## The recombination model

A model is a file (TOML or JSON; the builder chooses) describing one locus for
one species. Models live alongside germline databases: built-in models ship
in the package, user models in `$ABSTAR_HOME/models/<name>/`, with the same
no-shadowing rule as databases (record 0010).

**Contents:**

- **Gene usage:** P(V), P(J), and for loci with D, P(D | J). These are given
  per **gene**, with an optional per-allele split. When the split is absent,
  alleles of a gene are equally likely.
- **Trimming:** per-gene distributions of 3' V deletion, 5' and 3' D deletion,
  and 5' J deletion. A per-segment default is used when a gene has none.
  Negative values mean palindromic P-nucleotides, following the IGoR
  convention, so P-addition needs no separate mechanism.
- **N additions:** length distributions for each junction (VD and DJ, or VJ),
  with a long tail available. Composition is a first-order Markov chain over
  A/C/G/T, so N regions have realistic base content.
- **D inversion:** probability that the D is used in reverse complement.
  Defaults to 0 for human IGH unless the source gives a value.
- **Provenance:** every distribution names its source (paper, table, figure)
  inside the file.

**Name mapping.** Model gene names must resolve against the germline database
in use. A model gene with no allele in the database is reported in a warning
and dropped, and its probability mass is renormalized across the rest. A
database gene with no model entry is never generated. A run in which the
model and database share less than a set fraction of the model's mass fails
with `INVALID_ARGUMENT`; the threshold is the builder's call, recorded in a
decision.

**Default human models (source C).** These are built from published summary
distributions for IGH, IGK, IGL, TRA and TRB, citing each table used. They
are not fit to any dataset by us. Correlations between steps beyond P(D | J)
are not modelled. This is accepted, because grading accuracy needs coverage
more than joint realism.

**IGoR import.** `abstar model import-igor <igor_model_dir> --name
<name> --germline-db <db>` converts an IGoR model's marginals into our format.
Joint distributions are kept where our format can represent them and
marginalized where it can't, with a report of what was dropped. Users run it
on IGoR models they obtain themselves. We do not ship IGoR or OLGA files.

**Species.** Human ships first. A model for another species is a separate
file. Mouse models wait on name mapping, because OGRDB mouse names (`IGHV0-…`)
differ from IMGT/IGoR names. Macaque waits on V regapping
([board #137](https://github.com/brineylab/abstar/issues/137)), because its V
alleles lack anchors (see below).

## Somatic hypermutation

- **Targeting and substitution** *(revised 2026-11-20: no openly licensed S5F
  table exists. The Yale tables are CC BY-NC-SA 3.0, the paper's CC BY
  supplement has figures only, and SHazaM is AGPL-3. Decided with the user:
  option A)*:
  - **Default: abstar's own 5-mer motif model ("motif").** Mutability comes
    from published motif rules, with weights chosen by us and documented as a
    heuristic, not fit to data:
    - WRCY/RGYW and WRC/GYW hotspots;
    - WA/TW hotspots;
    - SYC/GRS coldspots.

    Substitution follows a simple transition-biased table. It ships in the
    package with a clean license. Grading annotators needs hotspot bias, not
    exact S5F rates.
  - **User-supplied S5F:** `--shm-table DIR` points at a directory holding the
    S5F `Mutability.csv` and `Substitution.csv` in the format the Kleinstein
    lab distributes. Users download them under their own license. Nothing
    S5F-derived ships with abstar.
  - **Later:** fit our own S5F-style model on data we may redistribute, with
    the repertoire work (board card).
- **Rate:** each sequence draws a target V mutation frequency from a
  user-chosen distribution (fixed, uniform range, or beta).
  - **Definition.** The target is the **realized** divergence of the final
    read's V-derived bases from the true V allele, so `high_mutation` and the
    `bnab` ranges measure the same quantity.
  - **Placement.** Mutations are placed across all of V(D)J by the context model, one at a
    time, updating 5-mer contexts as they go. Repeat hits and back-mutations
    are allowed. Placement continues until the realized V divergence reaches
    the target.
  - *(Builder-check answer; my call.)*
- **Where:** across the whole V(D)J region, including N/P bases, which are
  mutated like any other base. The constant region (if included) is not
  mutated.
- **TCR:** no mutation by default; it can be set explicitly for testing.
- **Neutral by default:** no selection, so stop codons can appear. Productivity
  is labelled; see below.
- **Selection when `productive_only` is set:**
  - These substitutions are rejected, and another position and base are drawn:
    - one that would create a stop codon in the V reading frame (V start
      through J end, in the read's V phase);
    - one that would make the Cys104 codon anything but TGT/TGC;
    - one that would make the Trp/Phe118 codon anything but TGG/TTT/TTC.
  - Synonymous changes within those sets are allowed.
  - A rejected draw does not count toward the target; the realized-divergence
    stopping rule and the safety cap still apply.
  - Mutation never repairs a read. The rearrangement (after the indel and
    truncation, before mutation) must already be productive: in-frame
    junction, no stop codon in the V frame from V start to J end, conserved
    Cys104 and Trp/Phe118 codons, and both anchors in the read. Otherwise the
    whole read is resampled before it is mutated. Letting SHM remove a stop
    would favour the more mutated reads and skew the mutation distribution.
  - Why: at `bnab` divergence nearly every neutral read has a stop, and most
    lose an anchor codon. Filtering neutral reads needed ~220 attempts per
    read and aborted large runs; with selection it needs ~6. Real mutated
    antibodies are productive because selection removes such reads
    (record 0024).

## Indels

- Optional, off by default, with a per-sequence probability.
- Placed in V or J, outside the conserved Cys104 and Trp/Phe118 codons, so
  the junction stays defined.
- Weighted toward CDRs and their flanks, as seen in VRC01-class and similar
  bnAbs.
- Length distribution up to about 12 codons. In-frame by default;
  frameshifting indels are allowed when explicitly requested.
- Every indel is recorded exactly: type, position in the read and in germline
  coordinates, length, and inserted bases. Inserted bases are random by
  default; duplication of neighbouring sequence is an option.
- **Grading indels** *(builder-check answer; my call)*:
  - **Equivalent placements.** An indel in a homopolymer or repeat can sit at
    any of several equivalent positions. Each indel records the full range of
    equivalent placements, from leftmost to rightmost (as with VCF
    normalization). An annotator's placement counts as correct anywhere in
    that range.
  - **Allele sets.** These are computed from the true indel placement. An
    annotator that picks a different but equally good placement could, rarely,
    rank alleles differently. This is accepted and documented.

## Truncation and constant region

- **Truncation:** optional 5' and 3' trimming of the finished read, drawn from
  user-set distributions, to imitate primer-based and partial reads. Recorded
  exactly.
- **Constant region:** optionally append the first N bases of a constant
  allele, chosen by isotype weights (for example IgM/IgG1/…/IgK). This gives
  `c_call` a truth to test against. Off by default.

## Productivity

- Every sequence gets the AIRR definitions of `productive`, `vj_in_frame` and
  `stop_codon`, computed from the final read.
- `--productive-only` resamples until the requested number of productive
  sequences exists, and reports the rejection rate.

## Presets

Named parameter bundles for common benchmark settings, for example `naive`,
`memory`, `bnab` and `tcr`. `bnab` produces productive output (it turns on
`productive_only`; `--no-productive-only` turns it off) and means:

- mutation drawn from 15–35%;
- the long-tail N lengths enabled, tuned so that for IGH at least 25% of
  sequences have CDRH3 ≥ 24 aa and at least 2% have CDRH3 ≥ 35 aa (my call;
  validation item 5 checks these fractions). The fractions are measured on
  productive sequences (`--productive-only`);
- indels on, at about 10% of productive sequences.

Presets are listed in `abstar manual generate`; explicit options override
them.

## Truth: what counts as correct

The rule was decided in conversation and is recorded as a decision record:
**the correct answer for each ambiguous field is the set of answers tied for
the fewest mismatches against the observed read.** The true generative value
is always recorded alongside it.

- **Allele calls (V, D, J, C):** the set of database alleles tied for fewest
  mismatches over the region of the read derived from that segment. Mutated N
  bases are not counted against a segment. If mutation made a different
  allele strictly closer than the true one, the true allele is **not** in the
  set. The sequence is then labelled `unrecoverable_<segment>`: no annotator
  could get it right, and it counts toward the reported error floor (the "bf16
  rounding error").
- **Gene calls:** the genes of the alleles in the allele set.
- **Boundaries:** `v_sequence_end`, `d_sequence_start` and `d_sequence_end`,
  and `j_sequence_start` (with the matching `np1`/`np2` lengths).
  - **Why mismatches alone can't decide this.** Counting mismatches alone is
    degenerate for cut points: N bases are never counted, so assigning every
    base to N would score zero. Boundaries are therefore scored by
    **matches minus mismatches**: each read base attributed to a germline
    segment scores +1 if it matches and −1 if it doesn't; N/P bases score 0.
  - **The valid set.** It is every joint configuration (V end, D start,
    D end, J start) that ties for the highest total score, with segments not
    overlapping. Two cases follow from this:
    - N bases that happen to match the adjacent germline are claimed by that segment (+1 beats 0), so the boundary moves past the true cut. The set widens only when two segments can claim the same bases equally, e.g. a base matching both V's continuation and J's start. The true cut is always recorded separately.
    - A mutated base at a segment edge can be given to N with no penalty.
  - **Which alleles.** Ties are computed for every allele in that segment's
    valid set, and their union is taken.
  - **How it is stored.** Each boundary stores the **list** of values that
    appear in any tied configuration. It is a list, not a `[min, max]` range,
    because tied cut points need not be contiguous. Each boundary is graded
    on its own against its list.
  - *(Builder-check answer; my call.)*
- **No D:** if the read's D-derived region is too short to tell apart from N
  sequence, "no D call" (null) is also a valid `d_call`. This applies when the
  best D match is shorter than a set length; the builder sets the length
  (proposal: 5 nt) and records it in a decision.
- **Junction and CDR3:** set by the conserved anchors, so there is exactly one
  answer. Indels never touch the anchors.
- **Derived fields** (mutation counts, `sequence_alignment`, etc.): given for
  the true allele, so grading can compare like with like on the recoverable
  subset. They are not given per allele in the set.

### Computing the sets without alignment

The generator knows the origin of every read base, so it never aligns:

- **V:** the true allele's bases map to IMGT gapped positions, and every
  candidate V allele is compared position by position in the same IMGT
  coordinates. Positions affected by indels are skipped.
- **J:** an ungapped comparison of the J-derived region against every J
  allele at every offset. This handles J alleles without an anchor.
- **D:** an ungapped comparison of the D-derived region against every D (and,
  where inversion is enabled, its reverse complement) at every offset.
- **C:** an ungapped prefix comparison.

**Coverage rule.** A read base derived from a segment that a candidate allele
doesn't cover counts as a mismatch for that candidate: a partial allele can't
explain bases it lacks. Candidate positions beyond the read's span are
ignored, which handles truncation. *(My call.)*

**What anchors and gapping are needed** *(builder-check answer; my call)*:

- **The model may only generate** V alleles that have an anchor and J alleles
  that have an anchor, because the junction truth depends on them. An
  anchorless allele in the model's gene usage is dropped, with a warning, as
  in name mapping.
- **Candidate sets** include every allele of the locus in the database,
  anchored or not:
  - V candidates are compared in IMGT coordinates, so they must be gapped. In
    the built-ins they all are.
  - J and D candidates use the sliding comparison, so they need no anchor.
- **Known limit:** a V allele whose IMGT gapping is shifted (for example
  human IGLV3-2\*03, [#138](https://github.com/brineylab/abstar/issues/138))
  is compared at the wrong positions. It is a candidate in the set but may be
  under-counted. The limit is documented in the decision record and fixed when
  #138 is fixed.
- **Macaque** waits on #137, because almost none of its V alleles are IMGT-gapped.

## Difficulty labels

Each truth row carries boolean flags. Benchmark reports group by them:

- `high_mutation`: V mutation ≥ 20%;
- `long_cdr3`: CDR3 ≥ 24 aa for heavy chains, with a locus-specific threshold
  otherwise;
- `has_indel`;
- `truncated_5p`, `truncated_3p`;
- `unrecoverable_v`, `unrecoverable_d`, `unrecoverable_j`;
- `ambiguous_v`, `ambiguous_d`, `ambiguous_j`: the set has more than one
  member.

Exact thresholds are the builder's call, written in the truth schema's field
descriptions.

## Output

`abstar generate OUTPUT_DIR` writes the following. With no OUTPUT_DIR, FASTA
goes to stdout and truth goes nowhere unless requested; the exact
stdout/stderr convention follows `abstar run`.

- `<name>.fasta`: the reads, with `sequence_id` as the header.
- `<name>.truth.tsv` and/or `<name>.truth.parquet`: one row per read.
  - True generative values use the **same AIRR column names as `abstar run`**
    (`v_call`, `junction`, `v_sequence_end`, …), so grading is a join on
    `sequence_id`.
  - Valid-answer sets are list columns with a `_valid` suffix
    (`v_call_valid`, `v_sequence_end_valid` as `[min, max]`, …).
  - Event columns describe what happened: trims, N/P bases, inversion,
    mutations as positions and substitutions, indels, truncation.
  - The difficulty flags.
  - `sequence`: the read.
- `generate.json`: provenance and parameters, mirroring `run.json`. It records:
  - the abstar version;
  - the model name and sha256;
  - the database name and its `alleles_sha256`;
  - the seed and every parameter;
  - counts, including the productive rejection rate.

Every truth field is defined in `truth.toml` with a description and, for
abstar-specific fields, a decision ID. That gives it the same `abstar manual`
coverage as `fields.toml`.

## CLI and Python API

```
abstar generate [OUTPUT_DIR] --n INT --locus TEXT [--germline-db/-g TEXT=human]
                [--model TEXT] [--preset TEXT] [--seed INT]
                [--mutation SPEC] [--indels FLOAT] [--truncate SPEC]
                [--constant] [--productive-only] [--format airr|parquet]...
                [--threads INT] [--json] [--quiet/-q]
abstar model list|info|import-igor        # top-level, like `abstar germline` (an optional OUTPUT_DIR on a group would swallow a subcommand name)
```

```python
abstar.generate(n, locus, *, germline_db="human", model=None, preset=None,
                seed=None, output_dir=None, **options)
    -> pl.DataFrame | GenerateSummary
```

- `--model` defaults to the built-in model for the database's species and the
  locus.
- `--seed` is optional. When absent, a seed is drawn and recorded in
  `generate.json`, so every run can be reproduced.
- Output is identical for a given seed regardless of thread count.

## Validation of the generator itself

Required before the generator grades anything:

1. **Distribution recovery:** in large samples, the empirical gene usage,
   trims, N lengths and composition, inversion rate and mutation-frequency
   distribution match the model within statistical tolerance.
2. **Brute-force truth:** on small random cases, every `_valid` set equals a
   slow exhaustive comparison against every allele and every cut point.
3. **Event replay:** every read is rebuilt exactly from its recorded events
   and the germline database.
4. **Determinism:** the same seed gives byte-identical output at any thread
   count.
5. **Hard cases:** the `bnab` preset produces the specified mutation,
   CDR3-length and indel ranges.
6. **Speed:** ≥ 1M sequences in ≤ 10 s on the M-series laptop target, with the
   truth table included.

## Decision records to add

- Fewest-mismatches truth sets plus true values; the unrecoverable label; boundary scoring (matches minus mismatches, joint ties, stored as lists); the coverage rule; anchor requirements for generated vs candidate alleles; indel placement ranges.
- Model file format and name mapping, including the renormalization threshold.
- Default model source C; IGoR import only; no bundled GPL files.
- The mutation model (own motif default; user-supplied S5F; why S5F isn't bundled).
- The no-D threshold.
- The truth schema's `_valid` column convention.
- The shared-writers crate arrangement.

## Left to the builder

- Model file serialization (TOML vs JSON) and exact layout.
- Whether writers move to a shared crate or another arrangement keeps
  `generate` free of annotation code.
- Difficulty-flag thresholds, the no-D threshold, and the name-mapping
  coverage threshold (each recorded in a decision).
- Exact CLI option spellings within the shape above, and the preset contents
  beyond `bnab`'s stated ranges.
- Which published tables the default human models draw from, cited in the
  model files.
