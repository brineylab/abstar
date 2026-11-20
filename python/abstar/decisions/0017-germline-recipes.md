# 0017: Germline builds from TOML recipes; segment from the recipe key; OGRDB J anchors win

## Status

Accepted (2026-11-20)

## Context

Germline databases (0005) are built from OGRDB germline-set JSON, IMGT-style gapped FASTA, and existing databases (donor builds, 0009), in any mix of species (0007). A build has to say which files go into which database, which segment each FASTA holds, which species each source is, and how mixed-species IDs are prefixed. The build must also set V and J anchors (0016: 0-based offsets into the ungapped allele, first base of the V Cys104 or J Trp118/Phe118 codon), which later drive junction and CDR3 identification.

## Decision

**Recipe format.** A build is driven by one TOML recipe; paths are relative to the recipe file. Top-level keys: `name`, `description`, optional `keep` (a donor keep-list file) and `[prefixes]` (species → short prefix, required for every species iff the database has more than one species). Each `[[source]]` is either a file source (`species` plus any of `ogrdb_json`, `v_fasta`, `d_fasta`, `j_fasta`, `c_fasta`, and an optional free-text `description`) or a `database` source (a name looked up in the map Python supplies, else a directory path; it takes no species or files). Unknown keys are rejected. Each source yields one `Provenance` entry with file basenames and their sha256; OGRDB sources append each germline set's name, release version and release date to the description. An allele's `source` is its file basename, or `database:<name>`.

**Segment from the recipe key.** A FASTA's segment comes from the key it is listed under (`v_fasta` → V, and so on), never from the allele name: the IMGT constant `IGHD*01` would otherwise read as a D gene. Locus is the first three characters of the name and must be a known locus. IMGT `|`-delimited headers give the name from field 2 and functionality from field 4 (`F`/`ORF`/`P`, else `unknown`); other headers give the first token and `unknown`. OGRDB alleles take name, segment and locus from `label`, `sequence_type` and `locus`, functionality `F` if `functional` is true and `NF` otherwise, the V sequence from the IMGT delineation's `aligned_sequence`, and the D/J sequence from `coding_sequence`.

**Anchors.** V: the anchor is found only when IMGT positions 310-312 of the gapped sequence are `TGT`/`TGC`; it is the count of non-gap characters before position 310. J: an OGRDB `j_cdr3_end` gives the anchor `j_cdr3_end - 1` (in OGRDB's coding-sequence coordinates this is the first base of Trp/Phe118; verified on IGHJ6*03 = 29 and IGHJ4*02 = 14); otherwise the first match in any frame of `(TGG|TTT|TTC)GG[ACGT]{4}GG[ACGT]`. D and C have no anchor.

**Disagreement policy.** When OGRDB supplies `j_cdr3_end` and the motif search finds a different position (or none), the OGRDB value is kept and the build warns `motif anchor X disagrees with OGRDB j_cdr3_end Y`. Curated coordinates outrank a heuristic, and the warning makes every disagreement visible. A V or J without an anchor is reported by name with one warning per allele: the specific anchor-search message, which replaces the generic validation warning for that allele.

**IUPAC ambiguity codes.** Germline sequences accept the same alphabet as queries: `ACGTN` plus the IUPAC codes `RYSWKMBDHV`, case-insensitive, plus `.` as the IMGT gap. Anything else (`X`, `-`, `*`, `U`, digits) is rejected. Real OGRDB sets contain ambiguity codes: C57BL/6 `IGHV0-KAOG*00` has an `S` at aligned position 130, and legacy kept it verbatim. Anchor detection stays strict: an ambiguity code never satisfies the V `TGT`/`TGC` check or a `[ACGT]` motif position.

## Alternatives considered

- Parse the segment from the allele name: fails on constant-region names (`IGHD*01`, `IGHM*01`) and on non-IMGT names.
- A Python-side build script per database: the recipe keeps builds declarative and reproducible from the CLI, and the same Rust code builds both built-in and user databases.
- Always use the J motif, ignoring `j_cdr3_end`: discards OGRDB curation, and the motif can match in the wrong frame on unusual alleles.
- Fail the build on an OGRDB/motif disagreement: would block whole databases over one questionable allele; a named warning is enough for review.
- Accept V anchors by searching for a Cys codon near the 3' end of ungapped sequences: plausible but can pick the wrong codon; ungapped V genes wait for alignment-based gapping (spec build mode 2).

## Consequences

- One recipe file fully describes a database build, including provenance and checksums of every input.
- V genes whose gapping is not IMGT-numbered get no anchor and a named warning rather than a silently wrong anchor.
- OGRDB J anchors are authoritative; every place the motif disagrees is listed in the manifest warnings.
- Germline sequences may contain IUPAC ambiguity codes. The annotator must treat them as ambiguous bases (neither a match nor a mutation by default), not as literal characters.
- A recipe `name` must be a single plain path component: not empty, `.` or `..`, and free of `/`, `\` and control characters (`+`, as in `human+c57bl6`, is fine). The name becomes the database directory, so `name = "../victim"` with `overwrite` would otherwise replace and delete an unrelated directory, and a nested name would create a database that `list()` and `resolve()` cannot find.
