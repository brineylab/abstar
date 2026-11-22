# abstar run (the annotator) — design

Status: approved · 2026-11-20

Third of the four design documents for the abstar rewrite (Project 3 in the
[foundation design](2026-11-20-abstar-foundation-design.md#project-3--abstar-run-the-annotator)).

Built on:
- **the foundation:** germline databases, `fields.toml` and the writers in
  `abstar_table`, the streaming `run` pipeline, error codes, `abstar manual`,
  decision records;
- **`abstar generate`:** ground truth (spec:
  [2026-11-20-abstar-generate-design.md](2026-11-20-abstar-generate-design.md)).

## Intent

Replace the stub annotator in `abstar run` with real V(D)J assignment and
annotation. The goal is a version good enough to replace the current abstar
for the lab's highest-priority uses, and accurate within the limit set by
genuine biological ambiguity.

It is built test-first against `abstar generate`'s answer keys. It is checked
against external tools and real data so that it doesn't overfit to our own
generator.

## Goals

**Output for version 1**, as prioritized by the user:
1. **The AIRR minimum:** the 14 required fields, plus abstar's single
   primary call with `*_call_ties`, `*_species_ties`, `*_gene` and
   `*_species` (records 0007/0008).
2. **Regions:** FWR1–4 and CDR1–3 start/end positions, lengths, and nt and aa
   sequences, by IMGT definitions. CDR3 and junction too.
3. **Amino-acid sequences:** `sequence_aa`, `sequence_alignment_aa`,
   `germline_alignment_aa`, `junction_aa`.
4. **Assembled sequences:** `sequence_vdj`, `sequence_vdj_aa`,
   `sequence_vdjc`, `sequence_vdjc_aa`.
5. **Masks**, using legacy encodings and names for compatibility (nt and aa).
   Masks cover **only the V(D)J region** (user decision, 2026-11-20). Each nt
   mask indexes `sequence_vdj`, one character per base, and each aa mask
   indexes `sequence_vdj_aa`. They never include constant-region bases. This
   matches the legacy semantics, where the masks spanned the ungapped
   assembled V(D)J. They do **not** index the IMGT-gapped `sequence_alignment`.
   - `cdr_mask` (`0` FWR, `1`/`2`/`3` CDR1/2/3);
   - `gene_segment_mask` (`V`/`N`/`D`/`N`/`J`);
   - `nongermline_mask` (`1` for N bases and mismatches to germline, else
     `0`).
6. **Mutations for every segment, V, D, J and C (nt and aa):**
   - mutation calls, counts and frequencies;
   - insertions and deletions in their own fields.

   D and J mutation statistics are new; the old abstar did this only for V
   and C.

**Mutation positions use germline numbering**, so they compare across
sequences and lineages:
- **V nt:** IMGT-gapped position, e.g. `IMGT153:G>A`.
- **V aa:** IMGT amino-acid position.
- **D, J and C:** 1-based position along the called allele, e.g. `12:T>C`.
- **J aa:** in the frame set by the Trp/Phe118 codon.
- **D aa:** in the frame D actually occupies within the junction.
- **Frequency:** mutated positions divided by the aligned germline positions
  of that segment. Truncated reads are not penalized for bases they lack.

Mutation calls are made against each segment's **primary** call, so they are
a single value like `v_call`.

**Quality bars:**
- **Accuracy (agreed with the user, as a target).** The user expects 99.5%
  may be out of reach for some or all primary calls, possibly even V and J.
  It is pursued as a target, not assumed. The benchmark report states the
  achieved accuracy per field and subset, and explains each shortfall
  (genuine ambiguity the truth sets don't capture vs annotator error). The
  bar:
  - On `abstar generate` answer keys, every allele call, boundary, junction
    and mutation is correct on ≥ 99.5% of reads where the answer is
    recoverable. This is reported per difficulty subset: naive, heavily
    mutated, bnAb-like, indels, truncated.
  - It matches or beats the old abstar and IgBLAST in every subset.
  - Gene calls agree with the old abstar on the 50k-read real corpus, apart
    from explained differences.
- **Speed and memory (foundation targets):**
  - ≥ 10k sequences/s on a 10-core M-series laptop;
  - a single sequence in ≤ 300 ms end to end;
  - memory flat in input size.

## Scope of the quality bars (decided with the user)

- **The bars apply to:** human IGH, IGK, IGL, TRA and TRB, the loci with
  generator models.
- **Mouse** databases (c57bl6, balbc, human+c57bl6) are annotated. They are
  checked only for agreement with the old abstar on real data, because no
  generator model exists.
- **Macaque is deferred.** Until its V genes are regapped (#137), macaque reads
  report `unassigned` / `junction_undefined` (or partial annotation), stated
  in the manual. Regapping is not part of this project.
  (2026-11-21, record 0036: macaque is no longer shipped. A database whose V
  genes are not IMGT-gapped is refused at build and load; #137 regaps it.)
- **TRA/TRD cross-locus rearrangements** (TRAV/DV with TRDJ) are deferred. J
  candidates are constrained to the V's locus.

## Non-goals (this version)

- **The rest of the old abstar's ~170 columns.** Examples: leader-inclusive
  (`lvdjc`) sequences, gapped-sequence variants, CDR3 V/N/D/J subdivisions,
  `complete_vdj`, aa identity fields. All are planned; each follows when
  needed.
- **Antibody numbering schemes beyond IMGT** (Kabat, Chothia, etc.). Those
  belong to `abstar number`.
- **Error correction of long or noisy reads.** As in the foundation, they
  annotate as heavily mutated, indel-rich sequences.
- **Clonal clustering, lineage inference, pairing.**

## Anti-goals

- **Plausible but wrong annotations.** When a step can't be done reliably,
  the record says so. Examples: no J found, junction undefined.
- **Sharing code with the generator.** The crate rule stands: `generate`
  never depends on `annotate`. The annotator implements its own comparisons,
  even where the truth rules are the same.
- **Tuning to our generator alone.** Every accuracy claim is also checked
  against external tools and real data.

## Constraints

Inherited from the foundation and generate designs:
- Rust core, no external binaries or subprocesses in `crates/` or
  `python/abstar/`;
- one output row per input record, with `annotation_status` and
  `failure_reason`;
- AIRR 2.0 conventions (1-based closed coordinates in TSV/Parquet);
- identical TSV and Parquet content;
- deterministic output for any thread count;
- stable error codes;
- decision records for every consequential choice.

The scientific-correctness rules in `AGENTS.md` carry over in full:
- exact sequence IDs;
- every coordinate space named, with conversions at named boundaries;
- score evidence preserved for gene calls, with deterministic tie-breaks;
- J and D candidates constrained by locus;
- productivity by frame, junction, stop codons, ambiguity and motif rules;
- FWR3 and CDR3 share the junction-start endpoint;
- a codon that crosses a region boundary belongs to the downstream region;
- partial terminal codons don't translate.

## Approaches considered

| | Approach | Trade-off | Outcome |
|---|---|---|---|
| A | Staged pipeline: orient → V → J → D → C, then derive regions, mutations, masks | Fast and testable step by step; early boundary errors propagate into the junction, the hardest region | Rejected as-is |
| B | Joint probabilistic model of the whole V–N–D–N–J path (partis-style) | Best at junction ambiguity; 2–3 orders of magnitude slower; defines "equally good" by model probability, not fewest mismatches | Rejected |
| C | A's staged search for V and J, plus a **joint junction step** that settles V end, D and J start together, using the same matches-minus-mismatches scoring the answer key uses | Slightly more work than A; fixes A's main weakness at nearly A's speed; ties built the same way the answer key builds valid sets, without shared code | **Chosen** |

## Pipeline (per read)

Each step records why it stopped when it can't continue. A record is then
`unassigned` with a specific `failure_reason`, never an invented answer.

1. **Validate input** (foundation): empty or invalid characters give
   `invalid_input`.
2. **Orient.** Count k-mer hits against the database's V and J alleles on both
   strands, then pick the strand with the stronger V+J evidence. Ties between
   strands go to forward. Set `rev_comp`.
   - All coordinates from here on address the oriented sequence, as AIRR
     defines.
   - No V or J evidence on either strand gives `unassigned` with
     `failure_reason = no_v_hit`.
3. **V search.**
   - **Prefilter:** a k-mer prefilter ranks V **genes**. Every allele of every
     gene within a margin of the best gene's k-mer score is aligned, so tie
     sets are complete. A top-N allele cut, as in the speed spike, is not
     allowed for ties.
   - **Alignment:** SIMD affine-gap alignment, semi-global (free query ends;
     free germline 5' end for 5'-truncated reads; free germline 3' end).
     Gap penalties must allow in-frame indels of at least 12 codons
     (VRC01-class) without breaking the alignment.
   - **Tie set:** the best-scoring allele's alignment places the read on IMGT
     gapped positions. Every aligned candidate is then compared position by
     position in IMGT coordinates over the read's V-derived bases, skipping
     indel positions. A read base a candidate doesn't cover counts as a
     mismatch. The candidates with the fewest mismatches form the tie set.
     The primary call comes from the record 0008 tie-break: score, then
     F > ORF > P, then usage, then natural sort.
4. **J search.** Align all J alleles of the locus to the read downstream of the
   V alignment end, with a small overlap allowance. The tie set is the fewest
   mismatches over the J-derived bases, using the same coverage rule. No J
   found gives `failure_reason = no_j_hit`. V-only annotation is still
   reported for the fields that don't need J.
5. **Joint junction step** (loci with D: IGH, TRB, TRD; VJ loci skip D):
   - **Inputs:** the V tie set and J tie set, plus all D alleles of the locus
     **in both orientations**, at every offset in the gap between them. An
     inverted D match competes with forward matches on equal terms.
   - **Inverted D output:**
     - `d_inverted` (bool) is set when the primary D is used in reverse
       complement.
     - D coordinates and D mutation positions then refer to the
       reverse-complemented allele. This is the same convention the
       generator's truth uses (record 0025).
     - AIRR has no inversion field, so `d_inverted` is an abstar extension.
   - **Scoring:** cut points (V end, D start, D end, J start) are scored by
     matches minus mismatches over germline-attributed bases (N and P score
     0).
   - **D call and boundaries:** the D tie set holds the alleles reaching the
     best score. The boundaries come from the primary allele combination's
     best configuration, with deterministic tie-breaks. D is left uncalled
     unless the best D configuration scores at least 6 more than the best
     no-D configuration (in practice, an exact D run of at least 6 nt). Its
     abstar ties are then empty and `d_call` is null.
     - *Amended 2026-11-21 (user ruling):* originally "D is left uncalled
       when the best D match is shorter than 5 nt (the generator's no-D
       rule)". The margin of 6 (ruling B3, record 0030) removes chance D
       matches in N regions and raised true D accuracy on every generated
       set (TRB 85.7 → 91.9%). It is kept even though it lowers agreement
       with legacy (95.0 → 94.2%).
   - **Final V and J boundaries:** this step also settles the final V end and
     J start. The V/J body alignments from steps 3–4 fix everything upstream
     and downstream.
6. **C search** (when the read extends past J): a prefix alignment against the
   locus's C alleles gives `c_call`, its ties and the C boundaries.
7. **Derive annotations:**
   - **Regions:** IMGT nt positions on the primary V's gapped coordinates are
     FWR1 1–78, CDR1 79–114, FWR2 115–165, CDR2 166–195, FWR3 196–312.
     CDR3 runs from after the Cys104 codon to before the Trp/Phe118 codon,
     with the junction adding both anchors. FWR4 runs from the Trp/Phe118
     codon to J end. Truncated reads give partial or null regions, never
     invented ones.
   - **Alignment strings:**
     - `sequence_alignment` and `germline_alignment` use IMGT gaps for V, as
       AIRR describes.
     - N bases align to gaps in `germline_alignment`, as in the legacy
       convention.
     - The `_aa` versions translate in the V frame.
   - **Junction:** junction, `junction_aa`, CDR3, the AIRR coordinate fields
     (`v_sequence_start/end`, `v_germline_start/end`, … for D, J and C), np1
     and np2.
   - **Mutations:** the per-segment calls, counts, frequencies and indels
     described in Goals.
   - **Masks:** the three masks in Goals, nt and aa.
   - **Productivity:** AIRR `productive`, `vj_in_frame` and `stop_codon`,
     using the same definition as the generator.
   - **Assembled sequences:** `sequence_vdj(_aa)` and `sequence_vdjc(_aa)`.
8. **Record status.** Possible outcomes:
   - `annotated`;
   - `unassigned` with `no_v_hit`, `no_j_hit`, `junction_undefined`
     (anchors missing) or similar specific reasons;
   - `invalid_input`;
   - `error` (panics caught per record, as in the foundation).

## Output schema

- `fields.toml` grows to the v1 field list. The builder finalizes the exact
  names.
  - **Legacy names are reused** where a legacy field exists: `fwr1_start`,
    `cdr1_aa`, `v_mutations`, `v_mutation_count`, `nongermline_mask`, and so
    on. This eases migration for abutils.
  - **AIRR names are used** where AIRR defines the field.
  - **New fields follow the same pattern:** `d_mutations`,
    `d_mutation_count`, `d_mutation_frequency`, `j_…` (each with `_aa` too),
    and `v_insertions` / `v_deletions` for each segment.
- Every abstar-specific field names its decision record. `abstar manual
  fields` and the docs regenerate from the file.
- **Mutation call format:** `list<string>`.
  - V nt: `IMGT{pos}:{ref}>{alt}`.
  - D/J/C nt: `{pos}:{ref}>{alt}`.
  - aa: `{pos}:{ref}>{alt}` with IMGT aa positions for V.
  - Indels use `{pos}:ins:{bases}` / `{pos}:del:{bases}` in their own
    fields.
- **Counts and frequencies:** `int64` / `float64`.

## Generator additions (ground truth for the new fields)

`abstar generate`'s truth table must grow so every v1 field can be graded:
- FWR/CDR region start/end (true allele);
- the three masks (nt and aa);
- `sequence_vdj(_aa)` / `sequence_vdjc(_aa)`;
- the amino-acid fields;
- per-segment mutation calls in germline numbering (from the recorded events
  and true alleles);
- per-segment indels.

These go in the `generate` crate, implemented independently of the
annotator. They are valued for the true allele. Grading on fields that depend
on the allele call (regions, masks, mutation calls) is restricted to reads
whose primary call matches the true allele. Accuracy on the rest is already
graded through the call itself.

**Fields that depend on boundaries** *(builder-check answer; my call)*. The
affected fields are mutation calls/counts/frequencies, indels,
`nongermline_mask`, `gene_segment_mask`, np1/np2 and germline start/end.

- **The problem:** a valid boundary can differ from the true one. An N base
  that matches germline may join V; a mutated edge base may go to N. These
  fields then legitimately differ from truth even when the call and
  boundaries are correct.
- **The rule:** the grader compares them only over positions that **both**
  the annotator and the truth assign to the same segment. Positions in the
  disagreement zone are excluded, because they lie inside the `_valid`
  boundary range and are genuinely ambiguous. Counts and frequencies are
  compared over that shared span.
- **Why this rule:** no extra truth sets are needed, and the grader doesn't
  need a third implementation of the mutation logic.

## Evaluation harness

- **What the 99.5% bar means** *(builder-check answer; my call)*:
  - **Unit:** each field separately, per difficulty subset.
  - **Fields covered:**
    - the primary V, D, J and C calls (in their `_valid` sets);
    - the four boundaries (in their lists);
    - junction and `junction_aa` (exact);
    - productivity flags (exact);
    - mutation calls, counts and frequencies, and the masks, by the shared-span
      rule above.
  - **Not counted:** whether `*_call_ties` covers the whole valid set is
    reported but not part of the bar. Alignment strings and CIGARs are checked
    by unit tests, not graded.
  - **Mixed reads:** a read unrecoverable for one segment is excluded only for
    that segment's fields.
  - **CI gate:** it uses the agreed thresholds at a CI-friendly size with fixed
    seeds.
  - **If a subset misses:** v1 is not done until the user accepts the dated
    benchmark report. Any subset under the bar must be explained there, with
    the D call on bnAb-like reads the most likely candidate.
- **Grader** (`scripts/grade.py` or a small Python module outside the hot
  path):
  - joins `abstar run` output with a truth table on `sequence_id`;
  - grades each field by its rule: membership in `_valid` sets for calls and
    boundaries, exact match for the rest;
  - reports accuracy per difficulty flag and overall on the recoverable set.

  It is used in tests and for benchmark reports.
- **External tools:** a benchmark script outside the package runs IgBLAST and
  the old abstar on the same generated sets and grades them.
  - **Environment (approved by the user):** a separate micromamba/bioconda
    environment with `igblast`, plus `mmseqs2` and the legacy abstar
    dependencies. It is set up by a documented script, not a dependency of
    the package or CI.
  - **Same alleles for every tool:** IgBLAST runs on a database built with
    `makeblastdb` from abstar's built-in germlines, so all tools see the same
    alleles and names. Subprocesses are
  allowed here because it isn't in `crates/` or `python/abstar/`. Results go
  in a dated report under `docs/superpowers/` plus the commit message, not
  into CI. The tools aren't CI dependencies.
- **Real data:** compare gene calls with the legacy baseline
  (`legacy/test_data/bcr_corpus/baseline.parquet`, 50k reads). Every
  disagreement class is explained in the report.
- **CI:**
  - unit tests per pipeline step, each against hand-built cases with exact
    values (AGENTS.md: "focused biological regression test with exact
    expected values");
  - a generator-graded accuracy test with fixed seeds and thresholds:
    naive ≥ 99.5% and bnab ≥ 99.5% on recoverable reads, at a CI-friendly
    sample size;
  - the existing foundation tests.

## Decision records to add

- Approach C and the joint junction scoring.
- Tie-set construction for V/J/D/C, including the prefilter completeness rule.
- Mutation numbering, format and frequency denominator.
- Mask encodings (legacy-compatible).
- Region definitions and truncation behaviour.
- Status and failure reasons.
- D searched in both orientations; the `d_inverted` field and its coordinate convention.
- Field naming policy (legacy names reused).

## Left to the builder

- k-mer size and gene-margin rule for the prefilter (it must provably keep
  every allele that could tie the best).
- SIMD alignment implementation (block-aligner, used in the spike, is fine)
  and gap penalties, within the indel requirement.
- Exact field names for new fields, within the naming policy.
- How V-only records (no J) fill fields.
- Grader packaging and the external-benchmark script layout.
- Task order, provided the generator additions land before the fields they
  grade.
