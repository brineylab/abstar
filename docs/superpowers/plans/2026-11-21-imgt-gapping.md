# IMGT Gapping of V Germlines Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add IMGT gaps to ungapped or wrongly gapped V germline sequences (`abstar germline gap`, recipe `gap_v`, `abstar.germline.gap`), then use it to ship macaque again.

**Architecture:** A new `gap` module in `crates/germline` translates each V and aligns it (semi-global protein, BLOSUM62) to IMGT-gapped template V alleles of the same locus. It transfers framework IMGT positions from the best template and places CDR1/CDR2 gaps by IMGT's length rule. It checks the result and falls back to the next-ranked template, or reports a failure. Python exposes it through `crates/py`, `abstar.germline.gap` and a Typer command, and `abstar germline build` calls it for recipe sources with `gap_v = true`.

**Tech Stack:** Rust (abstar_germline, rayon), PyO3/maturin, Python 3.11+, Typer, polars, pytest.

**Spec:** `docs/superpowers/specs/2026-11-21-imgt-gapping-design.md` (read it first; this plan argues from it).

## Global Constraints

- Branch `rewrite`; never merge or open a PR. Commit after each task; end every commit message with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- `abstar_germline` must not depend on `abstar_annotate` or `abstar_generate` (`tests/python/test_crate_graph.py`).
- Only V alleles are gapped. D/J/C pass through unchanged. Any `.` in the input is stripped first; the input's own gaps are never trusted.
- IMGT position p occupies gapped nt `3(p−1)..3p`. FR1 1–26, CDR1 27–38, FR2 39–55, CDR2 56–65, FR3 66–104, Cys104 at nt 310–312 (0-based `309..312`). Positions after 104 are contiguous, with no gaps.
- CDR gaps open from the centre: CDR1 order `33, 32, 34, 31, 35, 30, 36, 29, 37, 28, 38, 27`; CDR2 order `61, 60, 62, 59, 63, 58, 64, 57, 65, 56`. A CDR1 longer than 12 or a CDR2 longer than 10 cannot be represented.
- Default templates: exactly `human`, `c57bl6`, `balbc`, never a database this tool gapped. `--templates`/`gap_templates`/`templates=` replace the default.
- Starting constants (their final values are set in Task 4 and recorded in decision 0037): BLOSUM62, gap open 11, extend 1; top-N retry 10; low-identity warning below 0.60; near-tie margin 5 score units.
- Warning codes (exact strings): `no_cys104`, `stop_codon`, `no_cys23`, `no_trp41`, `low_identity`, `template_disagreement`. Failure reasons (exact strings): `empty_sequence`, `no_alignment`, `framework_insertion`, `cdr_too_long`, `unplaceable_overhang`.
- New error code `NO_ALLELES_GAPPED` (exit 2, "no input V sequence could be gapped"). Some failed alleles are not an error.
- Every constant that changes calls lives in decision record 0037. Changing one needs a superseding record (the convention of record 0029).
- Verification before each commit:
  - `export PATH=$HOME/.cargo/bin:$PATH`
  - `cargo clippy --workspace --all-targets -- -D warnings`
  - `cargo test --workspace --release`
  - When Python changes: `unset CONDA_PREFIX; VIRTUAL_ENV=$PWD/.venv .venv/bin/maturin develop --release`, then `.venv/bin/python -m pytest tests/python -q --ignore tests/python/test_accuracy.py`.

## Review Focus

- Lowercase bases, `N` or other IUPAC codes inside a V are translated as `X`, scored as a mismatch, and the allele is still gapped (Task 3).
- An allele name without a known locus prefix (e.g. `seq1`), or a locus with no templates (TRBV against `--templates c57bl6`), falls back to all anchored templates instead of erroring (Task 3).
- A sequence whose length is not a multiple of 3 keeps its trailing bases contiguous after the last placed codon, and the dot-stripped output equals the input (Task 3).
- Duplicate names in the input FASTA are each gapped and reported, in input order, not de-duplicated (Task 5).
- Running `gap` on an already IMGT-gapped database reproduces it: the self round trip (Task 4) and an explicit human-IGHV idempotence test (Task 5).

---

### Task 1: IMGT position helpers and the CDR length rule

**Files:**
- Create: `crates/germline/src/gap/mod.rs` (module root; `pub mod imgt;` only for now), `crates/germline/src/gap/imgt.rs`
- Modify: `crates/germline/src/lib.rs` (add `pub mod gap;`)
- Test: `crates/germline/tests/gap_imgt.rs`

**Interfaces:**
- Produces, in `abstar_germline::gap::imgt`:
  - `pub enum Cdr { One, Two }`
  - `pub fn cdr_positions(cdr: Cdr, len: usize) -> Option<Vec<u16>>`: the IMGT positions a CDR of `len` residues occupies, ascending. `None` when `len` exceeds 12 (CDR1) or 10 (CDR2).
  - `pub fn translate(codon: &[u8]) -> u8`: standard code, uppercase; any codon with a non-ACGT base gives `b'X'`, stops give `b'*'`.
  - `pub struct Residue { pub pos: u16, pub aa: u8, pub nt: std::ops::Range<usize> }` (nt range in the gapped string).
  - `pub fn template_residues(gapped: &str) -> Vec<Residue>`: one residue per IMGT position whose codon is three bases (a `...` codon is absent; a codon holding partial dots, as in 5′-partial alleles, is skipped). Positions after 104 continue codon by codon to the end; a trailing partial codon is dropped.
  - `pub fn region(pos: u16) -> Region` with `pub enum Region { Fr1, Cdr1, Fr2, Cdr2, Fr3, Tail }` (Tail = 105+).

- [ ] **Step 1: Write the failing tests** in `crates/germline/tests/gap_imgt.rs`:
  - `cdr_rule_table`: `cdr_positions(Cdr::One, 8) == Some(vec![27,28,29,30,35,36,37,38])`, `(Cdr::One, 12)` is 27..=38, `(Cdr::One, 11)` excludes only 33, `(Cdr::Two, 3) == Some(vec![56,57,65])`, `(Cdr::Two, 8)` excludes 60 and 61, `(Cdr::One, 13) == None`, `(Cdr::Two, 11) == None`, `(Cdr::One, 0) == Some(vec![])`.
  - `cdr_rule_matches_shipped_databases`: load `python/abstar/germline_dbs/{human,c57bl6,balbc}` with `abstar_germline::load`. For every V allele with an anchor whose gapped string has no `.` before nt 78 (CDR1 test) or 165 (CDR2 test), i.e. a complete CDR, assert that the positions holding a codon equal `cdr_positions(cdr, n_codons)`. Assert at least 1,000 alleles were checked. (The only two shipped exceptions, TRBV5-4*04 and TRBV7-9*07, are 5′-partial and excluded by the completeness filter.)
  - `template_residues_human_ighv1_2`: for human `IGHV1-2*02`, residue 23 is `C`, 41 is `W`, 104 is `C` with `nt == 309..312`, and no residue has `pos` 10 or 73 (IGH framework gaps).
  - `template_residues_skip_partial_codons`: a gapped string of 16 dots followed by bases gives a first residue at position 7 (nt 18..21); the partial codon at nt 15..18 is skipped.

- [ ] **Step 2: Run to verify failure:** `cargo test -p abstar_germline --release --test gap_imgt`. Expected: compile errors (module missing).

- [ ] **Step 3: Implement** `crates/germline/src/gap/imgt.rs` per the Interfaces block. `cdr_positions` takes the CDR's full position range and removes positions in the Global Constraints order until `len` remain.

- [ ] **Step 4: Run to verify pass:** same command. Expected: all 4 pass.

- [ ] **Step 5: Commit:** `feat(germline): IMGT position helpers and CDR gap rule for V gapping`.

---

### Task 2: Semi-global protein aligner

**Files:**
- Create: `crates/germline/src/gap/align.rs`
- Modify: `crates/germline/src/gap/mod.rs` (`pub mod align;`)
- Test: `crates/germline/tests/gap_align.rs`

**Interfaces:**
- Produces, in `abstar_germline::gap::align`:
  - `pub struct Scoring { pub gap_open: i32, pub gap_extend: i32 }` with `pub const DEFAULT: Scoring = Scoring { gap_open: 11, gap_extend: 1 }` (a gap of length k costs `open + (k−1)·extend`).
  - `pub fn blosum62(a: u8, b: u8) -> i32`: the standard 24-symbol BLOSUM62 (`ARNDCQEGHILKMFPSTWYVBZX*`). Any symbol outside it scores as `X`.
  - `pub enum Col { Pair(usize, usize), QueryOnly(usize), TemplateOnly(usize) }`
  - `pub struct ProteinAlignment { pub score: i32, pub cols: Vec<Col> }`: columns from the first to the last aligned pair. Leading and trailing unaligned residues of either sequence are free and are not in `cols`.
  - `pub fn align(query: &[u8], template: &[u8], s: &Scoring) -> ProteinAlignment`: Gotoh affine-gap DP with free end gaps on both sequences. Ties prefer, in order: diagonal, then a gap in the query (`TemplateOnly`), then `QueryOnly`, and the leftmost placement.

- [ ] **Step 1: Write the failing tests** in `crates/germline/tests/gap_align.rs`:
  - `identical_sequences_align_fully`: `align(b"QVQLVQSGAEVKKPGASVKVSC", same)` has 22 `Pair(i,i)` columns and score equal to the sum of the diagonal BLOSUM62 scores.
  - `end_gaps_are_free`: query `b"VKKPGASVKV"` against template `b"QVQLVQSGAEVKKPGASVKVSCKAS"` aligns as 10 consecutive pairs starting at template index 10, with no gap columns.
  - `internal_deletion`: query = template with 2 residues removed from the middle of a 30-residue sequence. The alignment has exactly one run of 2 `TemplateOnly` columns, and the score is the pair sum minus 12.
  - `blosum62_spot_values`: `blosum62(b'W', b'W') == 11`, `(b'C', b'C') == 9`, `(b'A', b'R') == -1`, `(b'*', b'A') == -4`, `(b'X', b'A') == -1`.

- [ ] **Step 2: Run to verify failure:** `cargo test -p abstar_germline --release --test gap_align`. Expected: compile errors.

- [ ] **Step 3: Implement** `align.rs`. Scalar is fine.

- [ ] **Step 4: Run to verify pass:** same command. Expected: 4 pass.

- [ ] **Step 5: Commit:** `feat(germline): semi-global protein aligner (BLOSUM62) for V gapping`.

---

### Task 3: `gap_v`, the gapping algorithm

**Files:**
- Modify: `crates/germline/src/gap/mod.rs`, `crates/germline/Cargo.toml` (add `rayon = "1"`)
- Test: `crates/germline/tests/gap.rs`

**Interfaces:**
- Consumes: Task 1 (`imgt::*`) and Task 2 (`align::*`).
- Produces, in `abstar_germline::gap`:
  - `pub const DEFAULT_TEMPLATES: [&str; 3] = ["human", "c57bl6", "balbc"];`
  - `pub const TOP_N: usize = 10; pub const LOW_IDENTITY: f64 = 0.60; pub const TIE_MARGIN: i32 = 5;`
  - `pub struct GapResult { pub name: String, pub gapped: Option<String>, pub template: Option<String>, pub template_species: Option<String>, pub identity: Option<f64>, pub cdr1_length: Option<usize>, pub cdr2_length: Option<usize>, pub warnings: Vec<String>, pub failure: Option<String> }`. Here `template` is the template allele's `id`; `gapped` is `Some` exactly when `failure` is `None`.
  - `pub fn gap_v(queries: &[(String, String)], templates: &[Allele]) -> Vec<GapResult>`: one result per query, in input order (rayon over queries). Templates are the V alleles of `templates` with `anchor.is_some()`.

The algorithm is the spec's steps 1–8, with these decisions the spec leaves open:

- **Locus:** `name.get(..3)` parsed as `Locus` (as `sources::fasta` does). Unknown, or no templates of that locus, means all templates.
- **Query residues:** for each frame f in 0..3, the codons starting at f; a codon with a non-ACGT base is `X`. `cleaned = input with '.' and whitespace removed, uppercased`.
- **Ranking:** score every (frame, template) pair with `align(query_aa, template_aa, &Scoring::DEFAULT)`, sort by score (ties: lower template index, then lower frame), and try at most `TOP_N` pairs.
- **Placement**, for one pair:
  - A query residue in a `Pair` column takes the template residue's position, whatever its region (CDR residues get rearranged in the next bullet).
  - **Extra query residues** (`QueryOnly` in FR1/FR2/FR3, between template positions a < b): fill that many empty positions in (a, b), choosing the positions with the highest occupancy among the candidate templates of the locus (fraction of templates with a residue there). Ties go leftmost. Inserted residues are assigned in order to the chosen positions, sorted ascending. More inserted residues than empty positions gives `framework_insertion`.
  - **Missing query residues** (`TemplateOnly` runs in a framework): among placements of the run that give the same alignment score (slide the run left or right while the score is unchanged), choose the one whose emptied positions have the lowest total occupancy. Ties go leftmost.
  - **CDRs:** query residues between the last residue placed at ≤ 26 and the first at ≥ 39 are CDR1 (likewise ≤ 55 and ≥ 66 for CDR2), placed by `cdr_positions`. `None` gives `cdr_too_long`. A query that starts inside a CDR, with no FR residue before it, places its CDR residues on that CDR's last positions (right-aligned) *(plan ruling: 5′-partial alleles; the spec is silent)*.
  - **Overhangs:** query residues before the first `Pair` column take the positions just before it only if those positions are free and ≥ 1; otherwise `unplaceable_overhang`. Query residues after the last `Pair` continue contiguously only past position 104; before 104 it is `unplaceable_overhang`.
- **Nucleotides:** each placed residue writes its codon at `3(pos−1)`, and empty positions up to the last placed one are `...`. The f leading bases of frame f go directly before the first codon (leading dots fill the rest). If they would start before nt 0, it's `unplaceable_overhang`. Trailing bases after the last full codon are appended.
- **Checks (spec step 7):**
  - Required: `gapped.replace('.', "") == cleaned`, and every internal dot run is a multiple of 3 and codon-aligned. Otherwise try the next pair.
  - Cys104 preferred: walk the ranked pairs. The first pair passing the required checks and placing TGT/TGC at `309..312` wins. If none does, the first pair passing the required checks wins, with `no_cys104`.
  - Warnings: `stop_codon` (a `*` among placed residues), `no_cys23`, `no_trp41`, `low_identity` (identity = identical `Pair` residues / query residues in pairs, below `LOW_IDENTITY`), and `template_disagreement` (another pair passing the required checks, scoring within `TIE_MARGIN` of the winner, gives a different gapped string).
  - Failure: `empty_sequence` for an empty cleaned sequence. `no_alignment` when the best pair scores ≤ 0 or aligns fewer than 50 query residues *(plan ruling: below 50 aligned residues is not a V)*. Otherwise the reason from the best-ranked pair.

- [ ] **Step 1: Write the failing tests** in `crates/germline/tests/gap.rs`, using the loaded human, c57bl6 and balbc databases as needed:
  - `regaps_a_human_allele_from_human_templates`: strip `IGHV1-2*02` and gap it with human templates minus that allele. The output equals the shipped gapped string, `failure` is `None`, and `cdr1_length == Some(8)`.
  - `macaque_style_extra_gaps_are_ignored`: take the shipped IGHV1-2*02 gapped string, insert `...` after nt 27 and before nt 120 (macaque KIMDB style), and gap it with human templates minus the allele. The output equals the shipped string.
  - `five_prime_partial`: drop the first 16 nt of human IGHV3-23*01 (its shipped string has no `.` before nt 16), then gap. The output equals `"." * 16 + shipped[16:]`.
  - `trailing_partial_codon_kept`: append `"GA"` to a stripped human IGKV allele. The output ends with the shipped string + `"GA"`.
  - `iupac_bases_are_gapped`: replace one FR3 base with `N` and lowercase the whole input. Gapped, with the `N` at the same gapped index as the original base.
  - `unknown_locus_name_uses_all_templates`: rename IGHV1-2*02's query to `seq1`. Still gapped identically.
  - `missing_locus_falls_back`: gap human TRBV7-2*01 (stripped) with only c57bl6 templates. `failure` is `None`, and the result carries `low_identity`.
  - `stop_codon_is_a_warning`: put `TAG` at codon 50 of a human IGHV allele. Gapped, with `stop_codon` in warnings.
  - `prefers_cys104`: mutate IGHV1-2*02's Cys104 codon to `TCT` and gap. Gapped, `no_cys104`, and the codon at `309..312` is `TCT`.
  - `framework_insertion_fails`: insert 9 nt (3 codons) into FR2 of human IGHV1-2*02, between the codons of positions 45 and 46. IGH templates have a residue at every FR2 position, so there are no empty positions to fill. `failure == Some("framework_insertion")`.
  - `cdr_too_long_fails`: insert 15 nt inside CDR1 of IGHV1-2*02, making a 13-residue CDR1. `failure == Some("cdr_too_long")`.
  - `non_v_fails`: human IGHJ4*02 and a 300-nt random sequence (fixed seed) both fail with `no_alignment`. An empty string fails with `empty_sequence`.
  - `order_preserved`: 5 queries in, 5 results out, names in input order.

- [ ] **Step 2: Run to verify failure:** `cargo test -p abstar_germline --release --test gap`. Expected: compile errors.

- [ ] **Step 3: Implement** `gap_v` and the constants in `gap/mod.rs`, as specified above.

- [ ] **Step 4: Run to verify pass:** same command. Expected: 13 pass. If a test exposes a spec ambiguity that isn't resolved above, stop and report it rather than choosing.

- [ ] **Step 5: Commit:** `feat(germline): gap_v transfers IMGT positions from the best protein template`.

---

### Task 4: Round-trip accuracy, tuning, and decision record 0037

**Files:**
- Create: `crates/germline/tests/gap_roundtrip.rs`, `python/abstar/decisions/0037-imgt-gapping.md`
- Modify: `crates/germline/src/gap/mod.rs` and `align.rs` constants (only if tuning changes them)

**Interfaces:**
- Consumes: `gap_v`, `DEFAULT_TEMPLATES`, `TOP_N`, `LOW_IDENTITY`, `TIE_MARGIN`, `Scoring::DEFAULT`.

- [ ] **Step 1: Write the measurement tests** in `gap_roundtrip.rs`. Each test prints a per-locus table (n, exact gapped-string matches, Cys104 placed where IMGT places it) with `eprintln!` and asserts floors:
  - `self_round_trip_leave_one_out`: for each of human, c57bl6 and balbc, gap every anchored V allele with that database minus the allele. Assert exact matches ≥ 99% per database.
  - `mouse_from_human`: c57bl6 and balbc anchored V alleles gapped with human templates only.
  - `human_ig_from_mouse`: human IGH/IGK/IGL anchored V alleles gapped with c57bl6 + balbc templates only.
  - Cross-species floors: Cys104 placement ≥ 99% wherever IMGT's string has Cys104. Set the exact-match floor at the measured value rounded down to the next whole percent.
  - Running time: the whole file runs in under 120 s with `--release`. If it doesn't, add a k-mer prefilter (leave the floors unchanged).

- [ ] **Step 2: Run and tune:** `cargo test -p abstar_germline --release --test gap_roundtrip -- --nocapture`.
  - If a floor fails, sweep gap open (9–13), gap extend (1–2) and `TIE_MARGIN` (0–10). Keep the setting with the highest summed cross-species exact-match rate that still holds the self round trip at ≥ 99%.
  - Pin the chosen values in code and the measured floors in the test.

- [ ] **Step 3: Write decision record 0037** `python/abstar/decisions/0037-imgt-gapping.md` in the format of `python/abstar/decisions/README.md`. Cover:
  - the method;
  - the CDR rule and its orders;
  - the framework occupancy rules;
  - default templates and why gapped-by-abstar databases are excluded;
  - warning and failure codes;
  - every constant (gap open and extend, `TOP_N`, `LOW_IDENTITY`, `TIE_MARGIN`, the 50-residue `no_alignment` floor);
  - the round-trip tables from Step 2;
  - a "Constants" line saying that changing any of them needs a superseding record.

  Run `.venv/bin/python -m pytest tests/python/test_decisions.py -q`. Expected: pass.

- [ ] **Step 4: Commit:** `test(germline): V gapping round trips, tuned constants and decision 0037`.

---

### Task 5: Python API and `abstar germline gap`

**Files:**
- Modify: `crates/py/src/lib.rs` (add `germline_gap`, register it, and expose `DEFAULT_TEMPLATES` as module attribute `GAP_DEFAULT_TEMPLATES`), `python/abstar/germline.py` (add `gap`), `python/abstar/cli.py` (add the `gap` command to `germline_app`), `python/abstar/errors.py` (add `NO_ALLELES_GAPPED`), `python/abstar/manual/topics/examples.md` (a gap example under "Germline databases")
- Test: `tests/python/test_germline_gap.py`

**Interfaces:**
- Consumes: `abstar_germline::gap::{gap_v, GapResult, DEFAULT_TEMPLATES}`.
- Produces:
  - `_core.germline_gap(queries: list[tuple[str, str]], template_paths: list[str]) -> list[dict]`. Each dict holds the `GapResult` fields; `warnings` is a list of str.
  - `abstar.germline.gap(sequences, templates: list[str] | None = None) -> tuple[list[tuple[str, str]], pl.DataFrame]`. `sequences` is a FASTA path (str/Path, gz allowed) or a list of `(header, sequence)`. The returned FASTA records keep each input's full header line verbatim, gapped records only, in input order. Report columns, in order: `name, status, template, template_species, identity, cdr1_length, cdr2_length, warnings, failure_reason`. `status` is `gapped` or `failed`; `warnings` is joined with `;`. The name passed to Rust is the allele name parsed from the header the way `sources::fasta` parses it (field 1 of a `|`-separated IMGT header, otherwise the first word). Raises `AbstarError("NO_ALLELES_GAPPED", ...)` when no record was gapped. `templates=None` means `GAP_DEFAULT_TEMPLATES`, each resolved with `resolve()`.
  - CLI `abstar germline gap V_FASTA [--templates DB]... [--output PATH] [--report PATH] [--json]`:
    - FASTA goes to `--output`, or to stdout.
    - The report TSV goes to `--report`.
    - The summary (gapped/failed counts, then one line per failed allele with its reason) goes to stdout when `--output` is given, otherwise to stderr. `--json` makes the summary a JSON object `{n_input, n_gapped, n_failed, templates, failed: [{name, failure_reason}]}`.

- [ ] **Step 1: Write the failing tests** in `tests/python/test_germline_gap.py`:
  - `test_gap_regaps_human_ighv`: take every human IGHV allele from `alleles.tsv`, strip it, and call `gap(..., templates=["human"])`. Every output equals the shipped gapped string, and the report has 0 failures. This is the idempotence check from Review Focus. Self-templates are allowed here because it checks reproduction, not accuracy.
  - `test_gap_report_columns_and_order`: 3 records (two with the same header, plus one IGHJ sequence). The report has 3 rows in input order, with the exact columns above. The J row is `failed`/`no_alignment`. Both duplicates are gapped.
  - `test_gap_none_gapped_raises`: only J sequences raise `AbstarError` with code `NO_ALLELES_GAPPED`.
  - `test_cli_gap_stdout_and_report(tmp_path)`: run the CLI via `typer.testing.CliRunner` on a 2-record FASTA with `--report r.tsv`. Stdout is FASTA with the original headers and gapped sequences, and `r.tsv` has 2 rows.
  - `test_cli_gap_json_summary(tmp_path)`: with `--output o.fasta --json`, stdout parses as JSON with `n_input == 2`.
  - `test_gap_templates_replace_default`: `templates=["c57bl6"]` on a human IGKV allele reports `template_species == "mouse"`.

- [ ] **Step 2: Run to verify failure:** rebuild with maturin, then `.venv/bin/python -m pytest tests/python/test_germline_gap.py -q`. Expected: failures (no `gap`).

- [ ] **Step 3: Implement** the binding (release the GIL around `gap_v`, as `germline_build` does with `detach`), `germline.gap`, the CLI command, and the error code. Add the gap example to `examples.md`. `test_manual.py` checks that every new option is documented, and the commands topic is generated from Typer.

- [ ] **Step 4: Run to verify pass:** rebuild, then `.venv/bin/python -m pytest tests/python -q --ignore tests/python/test_accuracy.py`. Expected: all pass.

- [ ] **Step 5: Commit:** `feat: abstar germline gap command and abstar.germline.gap()`.

---

### Task 6: Recipe option `gap_v` in `abstar germline build`

**Files:**
- Modify: `crates/germline/src/recipe.rs` (Source gets `gap_v`, `gap_templates`), `crates/germline/src/build.rs`, `python/abstar/germline.py` (`build` resolves template paths), `python/abstar/manual/topics/examples.md` (recipe example)
- Test: `crates/germline/tests/build.rs` (add cases), `tests/python/test_germline_api.py` (add case)

**Interfaces:**
- Consumes: `gap_v`, `DEFAULT_TEMPLATES`.
- Produces:
  - Recipe keys per `[[source]]`: `gap_v = true` (default false) and `gap_templates = ["name", ...]` (optional).
  - `recipe::read` rejects `gap_v` on a `database` source, and rejects `gap_templates` without `gap_v`, each with `RECIPE_INVALID` naming the source number.
  - `build()` takes the template databases from the existing `database_paths` map: Python's `germline.build` adds each source's `gap_templates`, or `GAP_DEFAULT_TEMPLATES` when `gap_v` is set without it, resolved with `resolve()`.

For a `gap_v` source, every V allele from its files (FASTA and OGRDB) is gapped by `gap_v` before `anchors::assign`:
- A failed allele is dropped and gets a `BuildWarning { allele: id, message: "not gapped: <reason>" }`.
- A gapped allele's warnings become one `BuildWarning` with message `"gapped: <code>, <code>"`.
- The source's `Provenance` gets one `files` entry `gap-template:<db name>` per template, with that database's `alleles_sha256` in `sha256`.

- [ ] **Step 1: Write the failing tests:**
  - Rust `build.rs`, `gap_v_source_is_gapped`: a temp recipe with `species = "human"`, `gap_v = true`, `gap_templates = ["human"]`, and a V FASTA of 3 stripped human IGHV alleles (pass the human path in `database_paths`). The built database's V alleles equal the shipped gapped strings, all have anchors, and provenance contains `gap-template:human`.
  - Rust, `gap_v_drops_failures_with_warning`: the same with one J sequence in the V FASTA. 3 alleles are built, and there's a warning `not gapped: no_alignment` for the J record.
  - Rust, `gap_v_recipe_validation`: `gap_v` on a database source, and `gap_templates` without `gap_v`, each give `RECIPE_INVALID`.
  - Python `test_build_gap_v_default_templates(tmp_path)`: a recipe with `gap_v = true` and no `gap_templates` builds, and `info().provenance` lists `gap-template:human`, `gap-template:c57bl6` and `gap-template:balbc`.

- [ ] **Step 2: Run to verify failure:** `cargo test -p abstar_germline --release --test build`, then the Python test after a rebuild. Expected: failures.

- [ ] **Step 3: Implement** the recipe fields, validation, build hook and Python path resolution. Add a `gap_v` recipe example to `examples.md`.

- [ ] **Step 4: Run to verify pass:** the full verification block from Global Constraints. Expected: all pass.

- [ ] **Step 5: Commit:** `feat(germline): gap_v recipe option gaps V sources during build`.

---

### Task 7: Ship macaque again (#137)

**Files:**
- Modify: `germline-sources/recipes/macaque.toml` (add `gap_v = true`; replace the description's "does not follow IMGT numbering" sentence with "V genes IMGT-gapped by abstar germline gap against human, c57bl6 and balbc (decision 0037)"), `scripts/build_builtin_germlines.py` (stop skipping macaque, and keep macaque out of any template list), `tests/python/test_builtin_databases.py` (`DATABASES` gains `macaque`; `LEGACY_SOURCES` maps it to the legacy macaque directory; segment counts; legacy-match dropouts; replace `test_macaque_is_not_builtin`), `python/abstar/decisions/0019-builtin-databases.md` and `0036-imgt-gapped-databases.md` (dated amendments that quote the old text), `README.md` (built-in list), `AGENTS.md` (if it lists built-ins)
- Create: `python/abstar/germline_dbs/macaque/` (built by the script)

**Interfaces:**
- Consumes: Task 6 (`gap_v` recipes).

- [ ] **Step 1: Build:** `.venv/bin/python scripts/build_builtin_germlines.py macaque`. Record the number of V alleles built, the failures with reasons, and the warning counts by code.

- [ ] **Step 2: Write the tests** in `test_builtin_databases.py`:
  - `test_macaque_v_anchored`: among macaque V alleles whose ungapped sequence has an open reading frame (no stop in at least one frame), at least 99% have an anchor.
  - `test_sequences_match_legacy[macaque-*]`: the existing ungapped comparison, with a module-level `MACAQUE_DROPOUTS = {name: reason}` set from Step 1 that is allowed to be missing from `built` (an empty dict if nothing failed).
  - `test_segment_counts[macaque]`: the measured counts.
  - `test_builtins_load[macaque]`: passes database validation (record 0036).
  - `test_macaque_never_a_default_template`: `"macaque" not in _core.GAP_DEFAULT_TEMPLATES`.

- [ ] **Step 3: Run:** the full verification block. Expected: all pass. If the 99% bar fails, inspect the failed and `no_cys104` alleles and report them instead of loosening the bar.

- [ ] **Step 4: Docs:**
  - 0019 amendment: macaque is built-in again, its V genes gapped by decision 0037.
  - 0036 amendment: ungapped databases can now be made valid with `abstar germline gap` or `gap_v`.
  - README built-in list.
  - Run `.venv/bin/python -m pytest tests/python/test_decisions.py tests/python/test_manual.py -q`. Expected: pass.

- [ ] **Step 5: Commit:** `feat: ship macaque again, V genes IMGT-gapped by abstar (#137)`.
