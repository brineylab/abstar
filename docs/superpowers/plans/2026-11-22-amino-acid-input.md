# Amino-Acid Input Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `abstar run --sequence-type aa` annotates protein sequences (V/J calls, IMGT regions, CDR3/junction, amino-acid mutations, indels, masks, productivity) into the existing AIRR table, with nucleotide-only fields null.

**Architecture:** A second `Annotate` implementation, `AaEngine` (new module `crates/annotate/src/aa/`), built from the existing germline database by translating V (IMGT-gapped, codon-wise) and J (anchor frame) once per run, plus a new scalar BLOSUM62 aligner (`align/protein.rs`). The pipeline picks `Engine` or `AaEngine` from `RunOptions.sequence_type`; everything downstream (writers, schema, `run.json`, Python, CLI) is shared. The nucleotide path is not modified except for the schema (`sequence` nullable, four new columns).

**Tech Stack:** Rust (abstar_annotate, PyO3 in abstar_py), Python (abstar API/CLI, typer, polars), pytest, cargo test.

**Spec:** `docs/superpowers/specs/2026-11-22-amino-acid-input-design.md`. Read it before any task.

## Global Constraints

- New code only in `crates/` and `python/abstar/`; never touch `legacy/`.
- Nucleotide-mode values must not change: existing cargo tests and the existing nucleotide accuracy floors (`tests/python/data/accuracy_floors.toml`) pass with no floor rewrite. The only allowed nucleotide-output changes are `sequence` becoming nullable and four new (null) columns.
- No new dependencies. BLOSUM62 is embedded as a constant (block-aligner is only a dependency on SIMD targets).
- Fields are single-sourced in `crates/annotate/fields.toml`; every new `abstar`-origin field has `decision = "0039"`.
- Default `sequence_type` is `"nt"`; values are exactly `"nt"` and `"aa"`.
- Gap costs: open −11, extend −1. Amino-acid alphabet accepted on input: `ACDEFGHIKLMNPQRSTVWY` + `XBZJUO*`, any case.
- Build/test: `maturin develop --release --extras test`, `cargo test --workspace --release`, `python -m pytest -q`.
- Every commit message ends with a blank line and `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.

## Review Focus

1. **A trailing `*`** (FASTA of translations): stop after FWR4 must not make the chain nonproductive; `stop_codon` counts only `*` inside `sequence_vdj_aa`. Test in Task 6.
2. **A leader/signal peptide or tag before FWR1**: V starts at FWR1; `sequence_vdj_aa` and the masks exclude the leader. Test in Task 6.
3. **Fab heavy chain with CH1** (structure data): trailing constant residues are ignored; regions and junction identical to the bare VH. Test in Task 6.
4. **scFv (VH–linker–VL) in one sequence**: one row, `annotated`, no panic. Test in Task 6.
5. **A nucleotide sequence passed in `aa` mode** (ACGT are valid residues): `unassigned`/`no_v_hit`, not a bogus call or a crash. Test in Task 6.

---

### Task 1: Schema, record type and decision record

**Files:**
- Modify: `crates/annotate/fields.toml` (`sequence` entry; four new entries)
- Modify: `crates/annotate/src/record.rs:28,128`, `crates/annotate/src/pipeline.rs:101,183,195`, `crates/annotate/src/stub.rs`, any other `sequence:` constructor the compiler flags
- Modify: `crates/annotate/tests/writers.rs:38-39` and any test comparing `r.sequence` as `String`
- Create: `python/abstar/decisions/0039-amino-acid-input.md`
- Modify: `AGENTS.md` (Specs line: add the amino-acid spec and plan)

**Interfaces:**
- Produces: `Record.sequence: Option<String>`; new `Record` fields `v_insertions_aa`, `v_deletions_aa`, `j_insertions_aa`, `j_deletions_aa: Option<Vec<String>>`.

- [ ] **Step 1: Update the schema test first.** In `schema_matches_toml`, replace the `sequence` non-nullable assertion and the count:

```rust
assert!(s.field_with_name("sequence").unwrap().is_nullable());
assert_eq!(s.fields().iter().filter(|f| !f.is_nullable()).count(), 2);
for f in ["v_insertions_aa", "v_deletions_aa", "j_insertions_aa", "j_deletions_aa"] {
    assert!(matches!(s.field_with_name(f).unwrap().data_type(), DataType::List(_)), "{f}");
}
```

- [ ] **Step 2: Run it to see it fail.** `cargo test --release -p abstar_annotate --test writers schema_matches_toml` → FAIL (`sequence` not nullable / no field `v_insertions_aa`).

- [ ] **Step 3: Edit `fields.toml`.** `sequence`: `nullable = true`, description `"The query nucleotide sequence. Null for amino-acid input (sequence_type aa), whose input is in sequence_aa."`, `decision = "0039"`. Add `v_insertions_aa`/`v_deletions_aa` right after `v_deletions`, and `j_insertions_aa`/`j_deletions_aa` right after `j_deletions`: `type = "list<string>"`, `nullable = true`, `origin = "abstar"`, `coordinates = "germline: IMGT residue number (V), J residue number from the Trp/Phe118 frame (J)"`, `decision = "0039"`. Descriptions: `"Insertions in V relative to the germline protein, as \"IMGT{pos}:ins:RESIDUES\", POS the germline residue the insertion follows. Amino-acid input only; null for nucleotide input."` (deletions: `"IMGT{pos}:del:RESIDUES"`, POS the first deleted germline residue; J versions use `{pos}` without `IMGT`). Also append `" For amino-acid input: <meaning>."` to the descriptions of `sequence_aa` (the input, uppercased), `productive` (rule in spec §Fields), `stop_codon` (`*` in `sequence_vdj_aa`), `vj_in_frame` (null), `germline_alignment_aa` and `gene_segment_mask_aa` (no D residues).

- [ ] **Step 4: Fix compile errors** (`sequence: seq.into()` → `Some(seq.into())`; `tsv_safe(&mut r.sequence)` → `if let Some(s) = r.sequence.as_mut() { tsv_safe(s) }`; tests comparing `r.sequence.as_str()` → `r.sequence.as_deref()`).

- [ ] **Step 5: Write `0039-amino-acid-input.md`** in the format of `python/abstar/decisions/README.md` (Status `Accepted (2026-11-22)`). Context/Decision/Alternatives/Consequences come from the spec's Intent, Design (mode selection, germline translation, fields table, null list, indel fields, no D/C), Approaches, and the schema note (nucleotide Parquet `sequence` becomes optional; output is not `airr`-valid in aa mode).

- [ ] **Step 6: Verify.** `cargo test --workspace --release` → PASS. `maturin develop --release --extras test && python -m pytest -q -m "not slow"` → PASS (includes `test_decisions.py`).

- [ ] **Step 7: Commit.** `git commit -am "feat: nullable sequence and amino-acid indel fields (record 0039)"` (add the new decision file first).

---

### Task 2: Protein aligner

**Files:**
- Create: `crates/annotate/src/align/protein.rs`
- Modify: `crates/annotate/src/align/mod.rs` (add `mod protein; pub use protein::{align_protein, blosum62, ProteinScoring};`)
- Test: `crates/annotate/tests/align_protein.rs`

**Interfaces:**
- Produces:
  - `pub struct ProteinScoring { pub gap_open: i32, pub gap_extend: i32 }` with `Default` = `{ -11, -1 }`.
  - `pub fn blosum62(a: u8, b: u8) -> i32`: NCBI BLOSUM62 over `ARNDCQEGHILKMFPSTWYVBZX*`, case-insensitive; `J`, `O`, `U` and any other byte score as `X`.
  - `pub fn align_protein(query: &[u8], germline: &[u8], s: &ProteinScoring) -> Alignment` (existing `align::Alignment`/`Op`; coordinates in residues).

- [ ] **Step 1: Write failing tests** in `tests/align_protein.rs`:

```rust
use abstar_annotate::align::{align_protein, blosum62, Op, ProteinScoring};
const S: ProteinScoring = ProteinScoring { gap_open: -11, gap_extend: -1 };

#[test] fn blosum62_values() {
    assert_eq!((blosum62(b'W', b'W'), blosum62(b'A', b'A'), blosum62(b'C', b'C')), (11, 4, 9));
    assert_eq!((blosum62(b'a', b'R'), blosum62(b'*', b'*'), blosum62(b'J', b'A')), (-1, 1, blosum62(b'X', b'A')));
}
#[test] fn identical_is_one_match_run() {
    let a = align_protein(b"QVQLVQSGAEV", b"QVQLVQSGAEV", &S);
    assert_eq!((a.q_start, a.q_end, a.g_start, a.g_end, a.ops.clone()), (0, 11, 0, 11, vec![(Op::Match, 11)]));
}
#[test] fn free_ends_on_both_sides() {
    // leader before the germline, germline running past the query's end
    let a = align_protein(b"MGWSQVQLVQSGAEV", b"QVQLVQSGAEVKKPG", &S);
    assert_eq!((a.q_start, a.q_end, a.g_start, a.g_end), (4, 15, 0, 11));
}
#[test] fn deletion_in_a_repeat_goes_leftmost() {
    let a = align_protein(b"MKVLLAAWGQGTLV", b"MKVLLLAAWGQGTLV", &S);
    assert_eq!(a.ops, vec![(Op::Match, 3), (Op::Del, 1), (Op::Match, 11)]);
}
#[test] fn insertion_in_a_repeat_goes_leftmost() {
    let a = align_protein(b"MKVLLLAAWGQGTLV", b"MKVLLAAWGQGTLV", &S);
    assert_eq!(a.ops, vec![(Op::Match, 3), (Op::Ins, 1), (Op::Match, 11)]);
}
#[test] fn x_never_matches() {
    let a = align_protein(b"QVQXVQSGAEV", b"QVQLVQSGAEV", &S);
    assert_eq!(a.ops, vec![(Op::Match, 3), (Op::Mismatch, 1), (Op::Match, 7)]);
}
#[test] fn unrelated_is_empty() {
    assert_eq!(align_protein(b"WWWW", b"PPPP", &S).ops, vec![]);
}
```

- [ ] **Step 2: Run** `cargo test --release -p abstar_annotate --test align_protein` → FAIL (unresolved imports).

- [ ] **Step 3: Implement.** Scalar Gotoh overlap alignment (query and germline ends free), BLOSUM62 substitution, affine gaps. `score` is the DP optimum; terminal gaps are not reported in `ops` (same contract as `align_semiglobal`). Optimum ≤ 0 → empty alignment (all zero, no ops). `Op::Match` iff residues equal case-insensitively and neither is `X`. Traceback prefers diagonal over gap on equal scores, so gaps land leftmost. Add `// ponytail: scalar O(n·m) per allele; AA k-mer prefilter + block-aligner AAMatrix if protein volume grows`.

- [ ] **Step 4: Run** the same command → PASS.

- [ ] **Step 5: Commit** `feat: scalar BLOSUM62 protein aligner`.

---

### Task 3: Amino-acid germlines

**Files:**
- Create: `crates/annotate/src/aa/mod.rs` (module root; only `pub mod germline;` for now), `crates/annotate/src/aa/germline.rs`
- Modify: `crates/annotate/src/lib.rs` (`pub mod aa;`)
- Test: `crates/annotate/tests/aa_germline.rs`

**Interfaces:**
- Consumes: `abstar_germline::{Allele, Locus, Segment}`, `derive::translate`.
- Produces:

```rust
pub struct AaGroup {
    /// Database indices of every allele with this exact translation, ascending.
    pub alleles: Vec<usize>,
    pub locus: Locus,
    /// Ungapped residues.
    pub seq: Vec<u8>,
    /// V: ungapped residue index -> IMGT gapped residue index (0-based); J: identity.
    pub imgt: Vec<usize>,
    /// V: ungapped index of Cys104 (IMGT 104); J: index of Trp/Phe118. None if the allele has no anchor.
    pub anchor: Option<usize>,
}
pub struct AaSegment { pub groups: Vec<AaGroup> }
impl AaSegment {
    pub fn new(alleles: &[Allele], segment: Segment, loci: Option<&[Locus]>, nonfunctional: bool) -> AaSegment;
}
pub fn translate_v(gapped_nt: &str) -> Option<Vec<u8>>; // IMGT-gapped residues with '.'; None if no '.' or a codon mixes '.' and bases
pub fn translate_j(ungapped_nt: &str, anchor: usize) -> (Vec<u8>, usize); // residues from anchor % 3, anchor residue index anchor / 3
```

Groups are keyed by (locus, gapped residues), in ascending order of their first allele's database index. Allele filtering matches `SegmentIndex::for_loci` exactly. V alleles where `translate_v` is `None`, and J alleles without `anchor`, are excluded. V `anchor` = number of non-`.` residues before gapped index 103, when the nucleotide allele has an anchor.

- [ ] **Step 1: Failing tests** (`tests/aa_germline.rs`, using `common::human_db()`):
  - `translate_v_keeps_imgt_numbering`: for IGHV1-2*02, `translate_v(gapped)` has `C` at index 103, `.` where the nt gapped codon is `...`, and its non-`.` residues start `QVQLVQSGAEVKKPGASVKVSCKAS`.
  - `translate_v_rejects_ungapped`: `translate_v("CAGGTGCAG")` is `None`.
  - `translate_j_from_anchor`: IGHJ4*02 → residues `YFDYWGQGTLVTVSS`, anchor index 4 (the `W`).
  - `synonymous_alleles_share_a_group`: in `AaSegment::new(&db.alleles, Segment::V, Some(&[Locus::IGH]), false)`, the group containing IGHV1-2*02 also holds every other IGHV1-2 allele with an identical translation (assert by comparing translations), and no group mixes loci.
  - `nonfunctional_excluded_by_default`: IGKV2-29*01 (NF) is in no group unless `nonfunctional = true`.

- [ ] **Step 2: Run** `cargo test --release -p abstar_annotate --test aa_germline` → FAIL.
- [ ] **Step 3: Implement** `germline.rs` to the interface above (non-ACGT codons → `X` via `translate`).
- [ ] **Step 4: Run** → PASS.
- [ ] **Step 5: Commit** `feat: translate V and J germlines for amino-acid annotation`.

---

### Task 4: Amino-acid V/J assignment

**Files:**
- Create: `crates/annotate/src/aa/assign.rs`
- Modify: `crates/annotate/src/aa/mod.rs` (`pub mod assign;`). The nucleotide `assign::v::trim` is not reused (its column scale is ±`matched`); write a private `trim_protein` here.
- Test: `crates/annotate/tests/aa_assign.rs`

**Interfaces:**
- Consumes: `align_protein`, `blosum62`, `ProteinScoring` (Task 2); `AaSegment`, `AaGroup` (Task 3); `assign::primary_by_tiebreak(tied: &[usize], alleles: &[Allele]) -> usize`.
- Produces:

```rust
pub const MIN_V_SCORE_AA: i32 = 100; // starting value; Task 8 calibrates
pub const MIN_J_SCORE_AA: i32 = 25;  // starting value; Task 8 calibrates
pub struct AaHit {
    pub primary: usize,      // database index
    pub ties: Vec<usize>,    // database indices, primary first, then record 0008 order
    pub group: usize,        // index into AaSegment.groups of the primary's group
    pub aln: Alignment,      // query vs group.seq, trimmed to the segment-derived residues
}
pub struct AaCalls { pub v: AaHit, pub j: Option<AaHit>, pub v_end: usize, pub j_start: usize }
pub fn assign_aa(q: &[u8], v: &AaSegment, j: &AaSegment, alleles: &[Allele], s: &ProteinScoring) -> Option<AaCalls>; // None = no_v_hit
```

Rules (spec §Assignment):
- **V:** align against every group; best score must reach `MIN_V_SCORE_AA`. Trim the best alignment to its maximum-scoring contiguous stretch of columns (BLOSUM62 per aligned column, affine cost per gap run; equal scores: longest, then leftmost). Ties: groups of the best's locus with the fewest differences (mismatched + gapped columns of their own alignment) inside the best's trimmed query range; tied alleles = union of those groups' `alleles`; primary by `primary_by_tiebreak`; `aln` = the primary group's alignment, trimmed.
- **J:** search `q[from..]`, `from = min(v.aln.q_end, cys)` where `cys` is the query residue aligned to the V group's `anchor` (if aligned); only groups of the V locus; best ≥ `MIN_J_SCORE_AA`; trim the start only (stretch must end at the last column; equal scores: shortest). Ties/primary as V. `aln` coordinates are on the full query.
- **Boundary:** `v_end = v.aln.q_end`, `j_start = j.aln.q_start`; if `j_start < v_end`, cut at `c` in `j_start..=v_end` maximizing (V column scores with query < c) + (J column scores with query ≥ c), earliest `c` on ties, and set both to `c`. Without J: `v_end = v.aln.q_end`, `j_start = q.len()`.

- [ ] **Step 1: Failing tests** (`tests/aa_assign.rs`; `IGH_AA` = `"QVQLVQSGAEVKKPGASVKVSCKASGYTFTGYYMHWVRQAPGQGLGWMGWINPNSGGTNYAQKFQGRVTMTRDTSISTAYMELSRLRSDDTAVYYCARDRITIFGVVIIPDYFDYWGQGTLVTVSS"`):
  - `igh_calls`: V primary gene `IGHV1-2`, `IGHV1-2*02` in ties; J gene `IGHJ4`, `IGHJ4*02` in ties; `(v_end, j_start) == (98, 111)`.
  - `synonymous_alleles_tie`: every allele in the primary V group's `alleles` is in `ties`.
  - `no_j`: `IGH_AA[..96]` → `j` is `None`, `j_start == 96`.
  - `unrelated_protein_no_v`: ubiquitin `"MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"` → `None`.
  - `overlap_cut_earliest`: a read whose V germline tail and J head both match the same residues gets `v_end == j_start` at the earliest maximal cut (build from IGHV and IGHJ translations with a 2-residue overlap; assert the cut).

- [ ] **Step 2: Run** `cargo test --release -p abstar_annotate --test aa_assign` → FAIL.
- [ ] **Step 3: Implement** `assign.rs`.
- [ ] **Step 4: Run** → PASS.
- [ ] **Step 5: Commit** `feat: amino-acid V and J assignment`.

---

### Task 5: Amino-acid regions, mutations, masks and alignment strings

**Files:**
- Create: `crates/annotate/src/aa/derive.rs`
- Modify: `crates/annotate/src/aa/mod.rs` (`pub mod derive;`), `crates/annotate/src/derive/mod.rs` (`fn columns` → `pub(crate) fn columns`), `crates/annotate/src/derive/mutations.rs` (`fn segment` → `pub(crate) fn segment`; re-export both via `derive` as `pub(crate)`)
- Test: `crates/annotate/tests/aa_derive.rs`

**Interfaces:**
- Consumes: `AaCalls`, `AaHit` (Task 4); `AaSegment` (Task 3); `derive::columns(aln: &Alignment, span: Range<usize>) -> Vec<Col>`; `derive::segment(q, cols, g, label, None) -> Option<SegmentMutations>` (its `nt`, `insertions`, `deletions`, `count`, `frequency` are the residue-level values here).
- Produces:

```rust
pub struct AaRegions {
    pub fwr1: Option<(usize, usize)>, pub cdr1: Option<(usize, usize)>, pub fwr2: Option<(usize, usize)>,
    pub cdr2: Option<(usize, usize)>, pub fwr3: Option<(usize, usize)>, pub cdr3: Option<(usize, usize)>,
    pub fwr4: Option<(usize, usize)>,
    pub junction: Option<(usize, usize)>, // Cys104..=Trp/Phe118 as half-open
    pub vdj: (usize, usize),              // first V residue .. J end (query end without J)
    pub labels: Vec<Option<usize>>,       // region index 0..7 per query residue
}
pub fn aa_regions(q: &[u8], c: &AaCalls, v: &AaSegment, j: &AaSegment) -> AaRegions;
pub fn aa_mutations(q: &[u8], c: &AaCalls, v: &AaSegment, j: &AaSegment) -> [Option<SegmentMutations>; 2]; // V, J
pub struct AaMasks { pub cdr: Option<String>, pub segment: String, pub nongermline: String }
pub fn aa_masks(q: &[u8], c: &AaCalls, v: &AaSegment, j: &AaSegment, r: &AaRegions) -> AaMasks;
pub fn aa_alignments(q: &[u8], c: &AaCalls, v: &AaSegment, j: &AaSegment) -> (String, String); // sequence_alignment_aa, germline_alignment_aa
```

Rules (spec §Regions and junction, §Fields): region ends on IMGT gapped residue index (exclusive, 0-based) `[26, 38, 55, 65, 104]`; V columns = `columns(&v.aln, v.aln.q_start..v_end)`, J columns = `columns(&j.aln, j_start..j.aln.q_end)`; Cys104 projection ≤ 16 residues, accepted only if `C` and before `j_start`; Trp/Phe118 projection along J's first diagonal accepted if the residue equals the J germline's anchor residue, else look back up to 5 residues from the last-diagonal projection; CDR3 null unless both anchors; inserted residue takes the previous residue's region. Mutation labels: V `IMGT{imgt[u]+1}`, J `{u+1}`. Masks over `vdj`: segment `V`/`N`/`J`; nongermline `1` for N, inserted or mismatched residues; cdr mask `0`/`1`/`2`/`3`, null when CDR3 is null. Alignment strings: V with the group's IMGT gaps (`.` in both), N residues against `-`, J as aligned, insert → `-` in germline, delete → `-` in sequence.

- [ ] **Step 1: Failing tests** (`tests/aa_derive.rs`, `IGH_AA` from Task 4):
  - `igh_regions`: `fwr1..fwr4` slices equal `QVQLVQSGAEVKKPGASVKVSCKAS`, `GYTFTGYY`, `MHWVRQAPGQGLGWMGW`, `INPNSGGT`, `NYAQKFQGRVTMTRDTSISTAYMELSRLRSDDTAVYYC`, `ARDRITIFGVVIIPDYFDY`, `WGQGTLVTVSS`; junction slice `CARDRITIFGVVIIPDYFDYW`.
  - `igh_mutations`: V `nt` (residue list) `["IMGT51:E>G"]`, count 1, insertions/deletions empty; J empty, count 0.
  - `igh_masks`: segment `"V"*98 + "N"*13 + "J"*15`; cdr `"0"*25+"1"*8+"0"*17+"2"*8+"0"*38+"3"*19+"0"*11`; nongermline has exactly 14 `1`s, 13 of them at 98..111.
  - `v_insertion_and_deletion`: IGH_AA with `"GS"` inserted after IMGT 61's residue and one FWR3 residue deleted → V `insertions == ["IMGT{p}:ins:GS"]`, `deletions == ["IMGT{d}:del:{R}"]` with `p`, `d`, `R` computed from the IGHV1-2*02 gapped translation in the test.
  - `projected_cys`: IGH_AA with the V's last 6 residues before `C` mutated so V's trimmed alignment ends early → `cdr3` still starts after the `C`.
  - `no_cys_no_cdr3`: IGH_AA with residues `90..98` (`TAVYYCAR`) removed: V ends before IMGT 104, no `C` at the projected position → `cdr3 == None`, `junction == None`, `aa_masks(..).cdr == None`.
  - `alignment_strings`: `sequence_alignment_aa.replace('.', "") == IGH_AA`; `germline_alignment_aa` has `-` exactly at the 13 N residues.

- [ ] **Step 2: Run** `cargo test --release -p abstar_annotate --test aa_derive` → FAIL.
- [ ] **Step 3: Implement** `aa/derive.rs`; make `columns`/`segment` `pub(crate)` (no behavior change).
- [ ] **Step 4: Run** → PASS; also `cargo test --release -p abstar_annotate --test derive --test mutations` → PASS (unchanged).
- [ ] **Step 5: Commit** `feat: amino-acid regions, mutations, indels and masks`.

---

### Task 6: `AaEngine` and pipeline wiring

**Files:**
- Modify: `crates/annotate/src/aa/mod.rs` (`AaEngine`, record assembly)
- Modify: `crates/annotate/src/engine.rs` (`Annotate::failed` default method)
- Modify: `crates/annotate/src/pipeline.rs` (`SequenceType`, `RunOptions.sequence_type`, engine factory, failure records via `failed`, `tsv_safe` on `sequence_aa`, `run.json` parameter)
- Test: `crates/annotate/tests/aa_engine.rs`

**Interfaces:**
- Consumes: Tasks 3–5.
- Produces:

```rust
// engine.rs
pub trait Annotate: Sync {
    fn annotate(&self, id: &str, seq: &str) -> Record;
    /// A record with only identity and status; the raw input goes where this engine reports it.
    fn failed(&self, id: &str, seq: &str, status: &str, reason: &str) -> Record { crate::record::failed(id, seq, status, reason) }
}
// pipeline.rs
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum SequenceType { #[default] Nt, Aa }
impl std::str::FromStr for SequenceType { /* "nt" | "aa"; else Err("sequence_type must be 'nt' or 'aa', got {s:?}") */ }
impl std::fmt::Display for SequenceType { /* "nt" | "aa" */ }
pub struct RunOptions { /* existing */ pub sequence_type: SequenceType }
// aa/mod.rs
pub struct AaEngine<'a> { /* alleles, v: AaSegment, j: AaSegment, scoring: ProteinScoring */ }
impl<'a> AaEngine<'a> { pub fn new(db: &'a Database, opts: &EngineOptions) -> Result<AaEngine<'a>, RunError>; }
impl Annotate for AaEngine<'_> { /* annotate + failed (input into sequence_aa, sequence None) */ }
```

`AaEngine::annotate`: uppercase input; invalid characters/empty → `failed(.., "invalid_input", "invalid_characters"|"empty_sequence")`; `assign_aa` None → `no_v_hit`; fill every field per spec §Fields (scores: BLOSUM62 sum over reported columns with affine gaps; identity: matches / aligned columns); status `no_j_hit` / `junction_undefined` / `annotated` as `record::assemble` does. `productive` = junction present ∧ no `*` in `sequence_vdj_aa` ∧ `junction_aa` starts with `C`, ends with `W` or `F` ∧ no `X` in it. `vj_in_frame` = None. In the pipeline, `run_to_files`, `run_to_writer` and `run_in_memory` build `Box<dyn Annotate + '_>` from `opts.sequence_type`; `annotate_guarded` and the `invalid_encoding` branch call `ann.failed(..)`; `run.json` `parameters.sequence_type` = `"nt"`/`"aa"`.

- [ ] **Step 1: Failing tests** (`tests/aa_engine.rs`; engine = `AaEngine::new(human_db, &EngineOptions::default())`):
  - `aa_end_to_end_igh`: on `IGH_AA`, status `annotated`; `sequence == None`; `sequence_aa == IGH_AA`; `locus IGH`; `v_gene IGHV1-2`, `j_gene IGHJ4`; regions `_aa` and `_length` as Task 5 with `_start/_end` all None; `junction_aa`; `v_mutations_aa == ["IMGT51:E>G"]`, `v_mutation_count_aa == 1`; `v_insertions_aa == Some(vec![])`; `productive == Some(true)`, `stop_codon == Some(false)`, `vj_in_frame == None`; `sequence_vdj_aa == IGH_AA`.
  - `aa_null_contract`: for the same record, every field in the spec's null list is None (list them explicitly: `sequence, rev_comp, vj_in_frame, v_cigar, d_cigar, j_cigar, v_sequence_start, …, junction, junction_length, np1, np1_length, np2, np2_length, sequence_vdj, sequence_vdjc, sequence_vdjc_aa, sequence_alignment, germline_alignment, fwr1 … fwr4 / cdr1 … cdr3 (nt) and their _start/_end, cdr_mask, gene_segment_mask, nongermline_mask, v_mutations, v_insertions, v_deletions, v_mutation_count, v_mutation_frequency, j_* equivalents, all d_* and c_*`).
  - `aa_statuses`: `""` → invalid_input/empty_sequence; `"QVQ-LV"` → invalid_input/invalid_characters, `sequence_aa == "QVQ-LV"`; ubiquitin → unassigned/no_v_hit; `IGH_AA[..96]` → unassigned/no_j_hit; IGH_AA with residues `90..98` removed (as in Task 5's `no_cys_no_cdr3`) → unassigned/junction_undefined.
  - `x_in_junction_is_nonproductive`: CDR3 residue replaced by `X` → `productive == Some(false)`.
  - Review Focus tests: `trailing_stop_is_productive` (`IGH_AA + "*"` → productive, stop_codon false); `leader_is_excluded` (`"MDWTWRILFLVAAATGAHS" + IGH_AA` → same regions, `sequence_vdj_aa == IGH_AA`); `fab_with_ch1` (`IGH_AA + "ASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSC"` → same regions/junction as IGH_AA, productive); `scfv_one_row` (`IGH_AA + "GGGGSGGGGSGGGGS" + "DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYSTPLTFGGGTKVEIK"` → status annotated, no panic); `nucleotide_in_aa_mode` (the nucleotide `igh_read()` sequence from `tests/engine.rs`, rebuilt here → unassigned/no_v_hit). If this fails at the starting thresholds, leave it failing with `#[ignore = "calibrated in Task 8"]` and report it; Task 8 adds it to the negative controls and removes the ignore.
  - `pipeline_aa_mode`: `run_in_memory` with `RunOptions { sequence_type: SequenceType::Aa, ..Default::default() }` over `Input::Records(vec![("a".into(), IGH_AA.into()), ("b".into(), "QV\tQ".into())])` → one batch, row a annotated with null `sequence`; row b invalid_input with `sequence_aa` tab-free.

- [ ] **Step 2: Run** `cargo test --release -p abstar_annotate --test aa_engine` → FAIL.
- [ ] **Step 3: Implement** `AaEngine`, `failed`, `SequenceType`, the factory and `run.json` field.
- [ ] **Step 4: Run** → PASS, then `cargo test --workspace --release` → PASS (nucleotide unchanged).
- [ ] **Step 5: Commit** `feat: AaEngine and sequence_type run option`.

---

### Task 7: Python API, bindings and CLI

**Files:**
- Modify: `crates/py/src/lib.rs` (`options(..., sequence_type: String)`; `run_to_files`, `run_to_stdout`, `run_in_memory` gain a trailing `sequence_type: String` parameter)
- Modify: `python/abstar/_core.pyi`, `python/abstar/api.py` (`run(..., sequence_type: str = "nt")`, `_classify(item, aa: bool)`, `_normalize(sequences, aa: bool = False)`, `_AA` alphabet set), `python/abstar/cli.py` (`--sequence-type`)
- Test: `tests/python/test_run_api.py`, `tests/python/test_cli.py`

**Interfaces:**
- Consumes: `SequenceType::from_str` (Task 6).
- Produces: `abstar.run(sequences, output_dir=None, *, ..., sequence_type="nt")`; CLI `abstar run --sequence-type {nt,aa}`.

- [ ] **Step 1: Failing tests:**
  - `test_run_api.py::test_aa_string_is_a_sequence`: `abstar.run("QVQLVQSGAEVKKPGASVKVSCKASGYTFTGYYMHWVRQAPGQGLGWMGWINPNSGGTNYAQKFQGRVTMTRDTSISTAYMELSRLRSDDTAVYYCARDRITIFGVVIIPDYFDYWGQGTLVTVSS", sequence_type="aa")` → one row, `v_gene == "IGHV1-2"`, `sequence` null, `cdr3_aa == "ARDRITIFGVVIIPDYFDY"`.
  - `test_aa_string_rejected_in_nt_mode`: the same string without `sequence_type` → `AbstarError` code `INPUT_NOT_FOUND` (unchanged behavior).
  - `test_bad_sequence_type`: `sequence_type="protein"` → `AbstarError` `INVALID_ARGUMENT`.
  - `test_aa_files_run_json`: FASTA of the IGH protein to `tmp_path` → `run.json["parameters"]["sequence_type"] == "aa"`.
  - `test_cli.py::test_run_sequence_type_aa`: `abstar run --seq <IGH_AA> --sequence-type aa` streams TSV whose `v_gene` column is `IGHV1-2`.

- [ ] **Step 2: Run** `maturin develop --release --extras test && python -m pytest tests/python/test_run_api.py tests/python/test_cli.py -q` → the new tests FAIL.
- [ ] **Step 3: Implement.** Validate `sequence_type` in `api.run` (`INVALID_ARGUMENT`, message `sequence_type must be 'nt' or 'aa', got {x!r}`) before calling `_core`; Rust `options()` also rejects with the same code. In `aa` mode `_classify` accepts strings that are a subset of `_AA` (both cases); error text `not an existing file and not an amino-acid sequence`. Update the `run` docstring.
- [ ] **Step 4: Run** the same command → PASS; `python -m pytest -q -m "not slow"` → PASS.
- [ ] **Step 5: Commit** `feat: sequence_type in abstar.run and --sequence-type in the CLI`.

---

### Task 8: Amino-acid accuracy gate and threshold calibration

**Files:**
- Modify: `tests/python/test_accuracy.py`
- Create: `tests/python/data/accuracy_floors_aa.toml`
- Modify (calibration only): `crates/annotate/src/aa/assign.rs` (`MIN_V_SCORE_AA`, `MIN_J_SCORE_AA`), possibly `ProteinScoring::default` gap costs

**Interfaces:**
- Consumes: `abstar.generate`, `abstar.run(..., sequence_type="aa")`.
- Produces: `test_accuracy_floors_aa`, floors file.

Gate (spec §Accuracy gate): a module fixture `truths` holds `abstar.generate(N, locus, seed=20261121, **opts)` per `SETS` entry, shared by `_measure` and the new test (no second generation). Per set: keep rows with `vj_in_frame == True`; run `abstar.run(list(zip(ids, truth["sequence_aa"])), sequence_type="aa")`; grade with a local `_grade_aa(truth, calls)` returning `{field: (accuracy, n)}` under these rules:
- `v_call`/`j_call`: truth allele ∈ call `*_call_ties`.
- `v_gene`/`j_gene`: truth gene ∈ genes of the call's tie set (gene = allele name before `*`).
- `locus`, `productive`, `stop_codon`: equality.
- `cdr3_aa`, `junction_aa`, `fwr1_aa`…`fwr3_aa`, `cdr1_aa`, `cdr2_aa` and their `_length`: equality, only on rows whose `v_call` is correct by the rule above (`fwr4_aa`/`fwr4_length`: rows whose `j_call` is correct; `cdr3`/`junction`: both).
- `v_mutations_aa`: entries with IMGT position 1–100 equal truth's entries with IMGT 1–100, on rows with correct `v_call`.
Floors: same writer/checker logic as `test_accuracy_floors` against `accuracy_floors_aa.toml` (`ABSTAR_WRITE_FLOORS=1` writes it), same `TARGET`, `MIN_N`, measured − 0.5 pp.

- [ ] **Step 1: Write the test** (fixture refactor, `_grade_aa`, `test_accuracy_floors_aa`), plus `test_aa_negative_controls`: ubiquitin, GFP (`MSKGEELFTGVVPVLVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK`), a nucleotide read (`truths["igh_memory"]["sequence"][0]`), hen lysozyme (`KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL`) and 20 seeded uniform-random 120-residue proteins → all `failure_reason == "no_v_hit"`.
- [ ] **Step 2: Calibrate.** For `MIN_V_SCORE_AA` in {50, 75, 100, 150, 200} and `MIN_J_SCORE_AA` in {15, 25, 35}, rebuild (`maturin develop --release`) and run `ABSTAR_FULL_TABLE=1 python -m pytest tests/python/test_accuracy.py -k "aa" -s` with floors temporarily ignored; record `v_gene`, `j_gene`, `junction_aa` accuracy per set and whether the negative controls pass. Pick the most permissive pair (lowest values) that maximizes true accuracy and still passes the negative controls; null rates are informational. Keep gap costs −11/−1 unless accuracy clearly improves with another BLAST protein preset (−10/−1, −9/−2). Record the table in the PR description. If Task 6 left `nucleotide_in_aa_mode` ignored, remove the ignore once the chosen thresholds pass it.
- [ ] **Step 3: Measure throughput**: time `abstar.run` on the 2,000-read `igh_memory` amino-acid set with default threads; put sequences/s in the PR description and in `0039`'s Consequences.
- [ ] **Step 4: Write floors** `ABSTAR_WRITE_FLOORS=1 python -m pytest tests/python/test_accuracy.py -k aa -q`, then run without the flag → PASS. Confirm `tests/python/data/accuracy_floors.toml` is unchanged (`git diff --stat` shows only the aa file).
- [ ] **Step 5: Full suite.** `cargo test --workspace --release` and `python -m pytest -q` → PASS.
- [ ] **Step 6: Commit** `test: amino-acid accuracy gate; calibrate amino-acid thresholds`.
