# Known edge cases and how abstar handles them

This catalog describes handling implemented in the active MMseqs assignment, annotation,
preprocessing, and output paths in this checkout. It groups cases by biological
mechanism or by the shared procedure that handles them. It includes established indel,
truncation, ambiguity, and preprocessing behavior alongside the more recent junction and
region-boundary recoveries.

“Handled” has several meanings here: recovering a supported annotation, preserving an
ambiguous call, emitting an explicit unassigned record, retaining a null value when
evidence is absent, or rejecting an unsupported annotation with an inspectable failure.
Each entry explains the outcome and its limits. An alignment-supported annotation is not
independent proof of biological correctness. This is a source-and-test-backed catalog,
not a claim of exhaustive coverage of all possible rearrangements.

Examples are hypothetical or schematic unless explicitly identified as tested regression
values. Short fragments illustrate a local mechanism and are not necessarily annotatable
standalone reads. `...` omits sequence; `-` is an alignment gap; `.` is IMGT numbering
padding. NT means nucleotide and AA means amino acid. NP1/NP2 denote the intersegment
sequence attributed to neither retained flanking germline segment; a null D call does
not establish the biological absence of D-derived bases.

Unless stated otherwise, query/reference intervals are **zero-based, half-open**
(`[start, end)`, suitable for Python slicing). Query coordinates address the
**oriented** input. IMGT positions and reported mutation/indel positions are one-based;
low-level region alignment endpoints and MMseqs coordinates have their own explicitly
noted conventions. AIRR TSV converts public intervals to one-based closed coordinates
once. The aligned rows used as evidence are called the *retained alignment* below.

Source links identify the implementation; named tests identify representative contracts
and exact regression examples. Source-only details are identified where a cited test
does not directly exercise that particular branch. Legacy commented implementations and
bypassed helper-only behavior are not presented as active pipeline recovery.

## Contents

- [Assignment ambiguity, receptor identity, and short search windows](#assignment-ambiguity-receptor-identity-and-short-search-windows)
- [Germline evidence, alignment coordinates, and IMGT numbering](#germline-evidence-alignment-coordinates-and-imgt-numbering)
- [Mutation and indel representation](#mutation-and-indel-representation)
- [Region extraction and continuous sequence assembly](#region-extraction-and-continuous-sequence-assembly)
- [Junction anchors, missing boundaries, and repeats](#junction-anchors-missing-boundaries-and-repeats)
- [Productivity and mask consistency](#productivity-and-mask-consistency)
- [Input collections and record identity](#input-collections-and-record-identity)
- [Preprocessing variation](#preprocessing-variation)
- [Explicit outcomes and recoverability](#explicit-outcomes-and-recoverability)
- [Strand, alignment, and serialization](#strand-alignment-and-serialization)
- [Custom-reference integrity and interrupted database builds](#custom-reference-integrity-and-interrupted-database-builds)

## Assignment ambiguity, receptor identity, and short search windows

### Equal allele evidence and competing high-identity hits

**Description.** A read may cover only sequence shared by several alleles, or a short,
nearly perfect hit may compete with a longer, better-supported hit. Taking the first
search result or the highest identity alone can produce an unstable or misleading call.

**Example.** These hypothetical V hits illustrate the two decisions:

| Hit | Bit score | E-value | Identity | Outcome |
| --- | ---: | ---: | ---: | --- |
| Short allele A match | 180 | 1e-8 | 0.99 | Loses to the stronger alignment |
| Longer allele B match | 220 | 1e-10 | 0.97 | Wins despite lower identity |
| Alleles B and C with identical values for **all** ranking metrics | 220 each | 1e-10 each | 0.97 each | Both calls retained |

**Handling.** MMseqs hits are ranked by descending bit score, ascending E-value, then
descending fractional identity, target coverage, query coverage, and alignment length.
Exact ties across those metrics become a sorted, unique, comma-delimited call, such as
`IGHV3-15*01,IGHV3-15*07`. Representative alignment details come from the alphabetically
first tied call. Equal-evidence rows for that same call are ordered by query start,
query end, query sequence, and remaining detail columns, making the result independent
of MMseqs row order. This retains ties among returned candidates; it is not a claim that
every possible allele was returned by the search.

**Evidence.** [mmseqs.py](../abstar/assigners/mmseqs.py), `select_best_hits`;
[test_mmseqs.py](../abstar/tests/test_mmseqs.py),
`test_select_best_hits_uses_alignment_evidence_not_output_order`,
`test_select_best_hits_retains_exact_allele_ties_deterministically`, and
`test_select_best_hits_same_call_ties_choose_details_independent_of_row_order`.

### Cross-locus matches and receptor-specific D-segment eligibility

**Description.** Short J or D sequences can resemble genes from another chain. The same
pipeline must also distinguish VJ rearrangements from VDJ rearrangements and BCR from
TCR reference databases.

**Example.** A read assigned `IGKV1-5*01` has candidate J matches to `IGHJ4*02` and
`IGKJ2*01`. Only the IGK candidate is compatible. In a TCR run, a TRB read can receive a
`TRBD` call, while a TRA read must not acquire a D call simply because its junction
resembles `TRBD1`.

**Handling.** The requested receptor selects the database tree and is carried through
assignment and realignment. J, D, and C candidate lists are filtered to the three-letter
locus of the selected V call before choosing their best hits. D-search queries are
generated only for IGH, TRB, and TRD; secondary D realignment uses the matching locus
and receptor too. IGK, IGL, TRA, and TRG skip D assignment. This is a compatibility
filter relative to the selected V locus, not an independent classifier for chimeric or
cross-locus rearrangements.

**Evidence.** [mmseqs.py](../abstar/assigners/mmseqs.py), `filter_compatible_locus`,
`build_dquery_fasta`; [germline.py](../abstar/annotation/germline.py), `reassign_dgene`;
[test_mmseqs.py](../abstar/tests/test_mmseqs.py),
`test_filter_compatible_locus_removes_cross_locus_hits`;
[test_germline.py](../abstar/tests/test_germline.py),
`test_reassign_dgene_uses_tcr_locus_database` and
`test_reassign_dgene_skips_loci_without_d_genes`;
[test_tcr_e2e.py](../abstar/tests/test_tcr_e2e.py), `test_tcr_goldens`.

### Very short post-V sequence and short or undetectable D remnants

**Description.** Read truncation can leave too little sequence for a J search.
Exonuclease trimming during recombination can leave only a few D-derived bases, or no
distinguishable D evidence at all.

**Example.** A V alignment ends at base 100 of a 112-nt read: only 12 nt remain.
Separately, a V-to-J interval might be `ACGTA` (five bases) or `ACGT` (four bases), even
in a D-bearing chain. A null D call alone cannot distinguish complete D trimming from
insufficient evidence.

**Handling.** Post-V J queries shorter than 14 nt are excluded because those queries can
break MMseqs with the configured J parameters. D queries require at least five extracted
bases; the D search uses a five-base minimum alignment, three-base k-mer, and permissive
E-value threshold (`1e6`) to allow short matches. Empty query sets are handled without
invoking a search on them. Missing J evidence follows the unassigned-record path.
Missing D evidence does not by itself prevent annotation: the annotator has a
locus-specific local-alignment fallback, and otherwise represents the interval without a
D call. These thresholds are search safeguards, not biological minimum lengths.

**Evidence.** [mmseqs.py](../abstar/assigners/mmseqs.py), `assign_germlines`,
`build_jquery_fasta`, `build_dquery_fasta`;
[annotator.py](../abstar/annotation/annotator.py), `annotate_single_sequence` (D-gene
block); [test_mmseqs.py](../abstar/tests/test_mmseqs.py),
`test_build_jquery_fasta_filters_short_sequences` and
`test_build_dquery_fasta_uses_d_bearing_loci_and_query_length`.

### Reverse-strand reads and inclusive search coordinates

**Description.** A reverse-complement read has decreasing MMseqs coordinates. The
sequence after V in biological orientation lies before V in the original input string.
An off-by-one extraction can also include the last V base in the next search or omit the
first available junction base.

**Example.** MMseqs coordinates are one-based inclusive:

```text
forward V: input bases 1..5       next query = input[5:]
reverse V: input bases 30..16     next query = input[:15]
           biological direction <--------------------
```

**Handling.** J/C query extraction branches on start versus end: forward post-segment
queries use `seq[end:]`; reverse ones use `seq[:end - 1]`. D extraction takes the pre-J
part of the already restricted J query, using `seq[:start - 1]` forward or `seq[start:]`
reverse. V orientation sets `rev_comp`, and annotation operates on `sequence_oriented`.
Conversion to public output coordinate conventions occurs separately at serialization.

**Evidence.** [mmseqs.py](../abstar/assigners/mmseqs.py), `build_jquery_fasta`,
`build_dquery_fasta`, `build_cquery_fasta`, `_write_assignment_outputs`;
[test_mmseqs.py](../abstar/tests/test_mmseqs.py),
`test_build_jquery_fasta_respects_inclusive_mmseqs_coordinates`,
`test_build_cquery_fasta_respects_inclusive_mmseqs_coordinates`, and
`test_build_dquery_fasta_uses_d_bearing_loci_and_query_length`;
[test_tcr_e2e.py](../abstar/tests/test_tcr_e2e.py), `test_tcr_goldens`.

### Missing constant sequence and deliberately absent optional database segments

**Description.** An amplicon can end within J, and a custom reference may deliberately
contain only V and J genes. Neither necessarily means the input or database is corrupt.
Conversely, a partially copied C or D index must not silently look like an intentionally
omitted segment.

**Example.** `V--junction--J|read ends` has no post-J sequence to search for C. A
light-chain database with complete V/J FASTAs and indexes but no D/C files is also
intentional. A database with only `mmseqs/c.lookup` remaining is incomplete.

**Handling.** C assignment is skipped when its database is absent or there is no post-J
query, leaving C fields null. Runtime database validation permits complete omission of
D/C, but once any component is present, requires the corresponding source FASTAs and
core search components together. Tests exercise VJ-only databases with IGK, IGL, TRA,
and TRG. This does not establish D-bearing chain support for a database that omits D.

**Evidence.** [mmseqs.py](../abstar/assigners/mmseqs.py), `assign_germlines`;
[abstar.py](../abstar/core/abstar.py), `_validate_assignment_database`;
[test_failure_contracts.py](../abstar/tests/test_failure_contracts.py),
`test_partial_optional_database_segment_is_invalid_input` and
`test_vj_only_database_without_optional_components_annotates_vj_chains`.

### Ambiguous reference names: species suffixes, substring matches, and IgD

**Description.** A gene name can be a substring of another name; combined-species
databases can contain species-qualified versions of the same allele name; and `IGHD`
names can describe diversity genes or the IgD constant region.

**Example.** A convenience lookup for `IGHV1-2` also matches `IGHV1-24`.
`IGHV1-2*02__homo_sapiens` retains information absent from the bare allele label. An
`IGHD` constant-region query must be directed to `c.fasta`, rather than the default
`d.fasta` inferred from its fourth character.

**Handling.** `get_germline` provides `exact_match`, `truncate_species`, and
`force_constant` controls. Appending `*` to a gene-level substring query separates
`IGHV1-2*` from `IGHV1-24*`; full allele lookups can use exact matching. Annotation
realignment uses exact matching with species truncation disabled, and constant-region
realignment forces the C reference file. Species stripping for display must therefore
not be confused with reference selection.

**Evidence.** [germline.py](../abstar/annotation/germline.py), `get_germline` and
`realign_germline`; [annotator.py](../abstar/annotation/annotator.py),
`annotate_single_sequence` (constant-region block). Basic exact, missing, and
multiple-match contracts are tested in
[test_germline.py](../abstar/tests/test_germline.py), `test_get_single_germline`,
`test_get_multiple_germlines`, and `test_get_single_germline_nonunique`; the lookup
controls above are also documented directly in the implementation.

## Germline evidence, alignment coordinates, and IMGT numbering

### Competing local alignments, terminal clipping, and repeated J sequence

**Description.** Semiglobal clipping can hide a better local endpoint, and repeated
J-like sequence can yield two equal scoring matches.

**Example.**

```text
oriented query:  V ... NP ... [J-like copy 1] ... [J-like copy 2]
local scores:                      65                  65
retained J:                   earlier endpoint
```

**Handling and limits.** Current V/J annotation obtains boundary evidence by local
alignment of the entire supplied query against the entire ungapped named reference,
independently of semiglobal clipping (`local_full_query=True`). Its scoring policy is
match `2`, mismatch `-3`, gap open `-12`, gap extension `-2`. V's query origin is the
oriented input origin; J's origin is the beginning of the supplied downstream query. The
retained local alignment supplies the endpoints, with alignment-library inclusive ends
converted once to half-open ends. The tested Parasail tie policy selects the earliest
query endpoint, then target endpoint. Semiglobal alignments remain available for other
annotation mappings. This is a deterministic scoring policy, not proof that the earliest
repeated copy is biologically correct. Older callable branches of
`process_vgene_alignment` and `process_jgene_alignment` extend to semiglobal flanks; the
main pipeline's `local_full_query=True` bypasses them. Current V/J annotation does not
always restore terminal mismatches to full germline length.

**Evidence.** [germline.py](../abstar/annotation/germline.py), `VJ_BOUNDARY_PARAMS`,
`realign_germline`, `process_vgene_alignment`, `process_jgene_alignment`;
[annotator.py](../abstar/annotation/annotator.py), V/J realignment calls;
[test_germline.py](../abstar/tests/test_germline.py),
`test_real_bcr_full_query_boundaries_match_authenticated_traces` (36 stored cases ×
V/J), `test_real_bcr_equal_score_j_repeat_selects_earlier_query_endpoint` (equal-score J
intervals `[452,487)` and `[493,528)`, retaining `[452,487)`).

### Indels make query and reference spans different

**Description.** A retained alignment containing an insertion or deletion cannot use one
sequence's length as the other's coordinate span.

**Example.** A D alignment may retain query `[11,14)` (three bases) but reference
`[2,7)` (five bases). Separately, an insertion-only interval can consume one query base
and zero reference bases.

**Handling and limits.** D query endpoints use query alignment coordinates; reference
endpoints use target coordinates, independently. C sequence endpoints add the
downstream-of-J origin and semiglobal/local offsets exactly once, avoiding duplicated
offsets on truncated constant reads. `alignment_columns_for_span` walks the paired trace
and records `(query_position, reference_position)` at every column boundary. Both
coordinate spaces must identify a represented boundary; missing boundaries, reversed
spans, unequal row lengths, and double-gap columns raise `ValueError`. A trace cannot be
silently combined with coordinates from another alignment. Coordinates here are
zero-based, half-open; the indel descriptions below instead use IMGT numbering.

**Evidence.** [germline.py](../abstar/annotation/germline.py),
`process_dgene_alignment`, `process_cgene_alignment`;
[positions.py](../abstar/annotation/positions.py), `alignment_columns_for_span`;
[test_germline.py](../abstar/tests/test_germline.py),
`test_d_germline_end_uses_target_span_when_alignment_has_indel`,
`test_c_sequence_end_does_not_double_count_query_offset`;
[test_positions.py](../abstar/tests/test_positions.py),
`test_retained_span_maps_both_coordinate_spaces_to_alignment_columns`,
`test_retained_span_rejects_a_boundary_absent_from_its_trace`.

### IMGT dots and biological insertions are different kinds of gaps

**Description.** A query can contain inserted bases adjacent to IMGT padding or at the
retained reference end. Advancing an IMGT template cursor for every query base would
lose or shift residues.

**Example.**

```text
query alignment:       ACAAAGT
reference alignment:   AC---GT
full IMGT template:    ..A.C..GT.
gapped query output:   ..A.CAAA..GT
```

**Handling and limits.** `get_gapped_sequence` advances through the template only when
the alignment consumes a reference residue. It emits intervening dots before their
following reference residue; insertion columns append the query character without
consuming the reference. It retains deletion dashes, leading IMGT dots only when the
retained reference begins at offset zero, and insertions after the last retained
reference residue. Trailing template dots without a following consumed residue are not
emitted. Nonzero starts begin at the first retained reference residue. The helper
rejects unequal alignment lengths, starts outside the reference, or alignment residues
inconsistent with the template. This is projection of an existing alignment, not de novo
detection of an insertion.

**Evidence.** [positions.py](../abstar/annotation/positions.py), `get_gapped_sequence`;
[test_positions.py](../abstar/tests/test_positions.py),
`test_gapped_sequence_preserves_insertions_offsets_and_reference_end`;
[test_real_bcr.py](../abstar/tests/test_real_bcr.py),
`test_full_constant_insertion_preserves_terminal_alignment_columns` (constant reference
starts 0, 1, and 2).

### A read starting partway through a codon still uses full-reference AA numbering

**Description.** Starting a retained V or C reference slice at nucleotide 1 or 2 must
not shift all subsequent amino-acid mutation positions by one.

**Example.** A V slice beginning at reference nucleotide offset `2` uses one-based slice
frame `2`; its first complete codon begins at full-reference nucleotide `3`, hence AA
offset `1`, not `0`.

**Handling and limits.** V and C frames derive from the retained ungapped reference
start: `(3 - start % 3) % 3 + 1`. `translated_reference_start(start, frame)` converts to
`(start + frame - 1) // 3`, accounting for skipped partial leading codons before AA
mutation/region mapping. No incomplete leading codon is invented. The exact real BCR
regression retains a V at oriented-query start `138`, reference start `2`, frame `2`,
and mutation `29:T>I` rather than `28:T>I`.

**Evidence.** [germline.py](../abstar/annotation/germline.py),
`translated_reference_start`, `process_vgene_alignment`, `process_cgene_alignment`;
[test_germline.py](../abstar/tests/test_germline.py),
`test_translated_reference_origin_accounts_for_partial_first_codon`;
[test_real_bcr.py](../abstar/tests/test_real_bcr.py),
`test_partial_v_codon_keeps_full_imgt_amino_acid_positions`.

### Identity includes indels, and short D evidence can lack AA identity

**Description.** An insertion/deletion must reduce sequence identity, while a retained
short D span may carry useful nucleotide evidence without containing a complete codon in
its frame.

**Examples.** Query `AC-GTA` versus reference `ACTG-A` has four matches in six columns:
identity `4/6`, not 100%. D query/reference `CCGG` in frame 3 retains NT identity `1.0`
but AA identity `None`; query `AC` versus reference `AT` in frame 1 gives `(0.5, None)`.

**Handling and limits.** `calculate_alignment_identity` counts exact nongap matches over
alignment columns, including single-gap columns in the denominator. Double-gap columns
are ignored defensively; empty evidence returns `None`; unequal alignment lengths raise
`ValueError`. `_segment_identities` uses retained NT alignment evidence when supplied
(or a global alignment for standalone calls). It separately translates the retained
query and reference in the requested slice frame and globally aligns those AA strings.
If either translation is empty, it retains NT identity and returns no AA identity. AA
identity therefore derives from translated segment alignment, not division of NT
identity or proof of complete D translation. Exact-match counting does not assign
partial credit to ambiguous nucleotide codes.

**Evidence.** [annotator.py](../abstar/annotation/annotator.py),
`calculate_alignment_identity`, `_segment_identities`;
[test_annotator.py](../abstar/tests/test_annotator.py),
`test_alignment_identity_counts_substitutions_and_indels`,
`test_alignment_identity_ignores_double_gap_columns`,
`test_alignment_identity_rejects_different_alignment_lengths`;
[test_real_bcr.py](../abstar/tests/test_real_bcr.py),
`test_pilot_loss_empty_d_translation_preserves_nt_identity`.

## Mutation and indel representation

### Substitutions do not absorb adjacent insertions or deletions

**Description.** An alignment can contain a substitution near a gap and a reference with
IMGT padding. Counting alignment columns as reference positions would misnumber the
substitution or call the gap a substitution.

**Example.**

```text
query:       ATGCTTCC
reference:   ATGCATGC
mutations:   5:A>T and 7:G>C (ungapped template)
```

Adding an IMGT dot before the fifth reference residue shifts its reported IMGT position
to 6.

**Handling and limits.** `annotate_mutations` walks paired alignment columns. Matches
advance the raw reference cursor; insertions (reference `-`) do not; deletions (query
`-`) advance it without creating a mutation. Other differences are recorded as
`IMGT_position:reference>query`. Conversion counts IMGT dots in the full numbering
template. V/C wrappers store `|`-joined events and counts for NT and AA inputs. The
algorithm records mismatches, including ambiguous characters; it does not establish that
a mismatch is a true somatic mutation rather than sequencing error. The low-level helper
assumes correctly paired inputs and uses `zip`, so it is not itself an unequal-length
alignment validator.

**Evidence.** [mutations.py](../abstar/annotation/mutations.py), `annotate_mutations`,
`annotate_v_mutations`, `annotate_c_mutations`;
[test_mutations.py](../abstar/tests/test_mutations.py),
`test_annotate_mutations_single_mutation_single_germline_gap`,
`test_annotate_mutations_multi_mutation_germline_gaps`;
[test_annotator.py](../abstar/tests/test_annotator.py), `test_v_mutations_annotated`.

### Multiple indel runs, IMGT-spanning deletions, and leading insertions

**Description.** A deletion can cross an IMGT dot; an insertion can precede the first
reference base; several indel runs can occur in one retained alignment.

**Examples:**

```text
query:       A--D
reference:   ABCD
template:    AB.CD
deletion:    2-4:2>BC!   (2 bases deleted, spanning IMGT positions 2 to 4)

query:       AA
reference:   -A
insertion:   0:1>A!     (before the first reference base)
```

Letters in the first schematic are positional placeholders. An in-frame nucleotide
insertion might be `5:3>GGG`; a separate two-base run receives its own event with `!`.

**Handling and limits.** Maximal runs of `-` in the reference identify insertions; runs
in the query identify deletions. Prefix counts exclude alignment gaps and add the
retained reference start before conversion to IMGT coordinates. Insertion positions
denote the preceding reference residue, including the explicit position-zero boundary.
Deletion start and inclusive end are converted independently, so intervening IMGT dots
do not corrupt the end position. `|` separates events; payload length is the number of
inserted/deleted bases, not the numeric width of the IMGT interval. `!` means that this
event's length is not divisible by three. No indels produce an empty string. These
helpers serialize supplied alignment evidence; they neither resolve alternative gap
placements nor correct the query. V calls use the IMGT V template; constant-region calls
currently pass the complete ungapped constant reference as their numbering template.

**Evidence.** [indels.py](../abstar/annotation/indels.py), `annotate_insertions`,
`annotate_deletions`; [positions.py](../abstar/annotation/positions.py),
`get_gapped_position_from_raw`; [test_indels.py](../abstar/tests/test_indels.py),
`test_deletion_end_accounts_for_intervening_imgt_gap`,
`test_leading_insertion_precedes_first_germline_base`,
`test_annotate_multiple_scattered_deletions`, `test_annotate_insertions_complex_case`;
[test_properties.py](../abstar/tests/test_properties.py),
`test_indel_descriptions_reconstruct_aligned_pair`, `test_imgt_gapped_round_trip`.

### Compensating frameshift-sized V indels

**Description.** A one-base insertion and later one-base deletion each receive `!` but
have zero net frame offset.

**Example.** The tested derived read emits `105:1>A!` and `114:1>T!`; these two events
alone do not add `out-of-frame indel(s)` to its productivity issues.

**Handling and limits.** The V annotation step sums insertion lengths and deletion
lengths from the retained alignment. It adds the V frameshift issue and sets
`v_frameshift=True` only when `(total_inserted - total_deleted) % 3 != 0`. Per-event
exclamation marks remain. Compensating indels can still alter intervening amino acids,
introduce a stop codon, or cause other productivity failures. Net compensation is not an
assertion that the read is functional. This specific cumulative rule is applied to V
indels; constant-region indel serialization does not run this V productivity rule.

**Evidence.** [annotator.py](../abstar/annotation/annotator.py), cumulative V indel
processing; [test_real_bcr.py](../abstar/tests/test_real_bcr.py),
`test_compensating_v_indels_use_one_retained_alignment`.

## Region extraction and continuous sequence assembly

### A read begins inside a region or after an entire region

**Example.** A 5′-truncated read starts halfway through FWR1, or begins inside FWR3 with
FWR1–CDR2 entirely absent. Missing sequence must not be replaced with germline bases.

**Handling.** `get_region_sequence()` returns the aligned, ungapped portion for a
partial region and `RegionSequence(None, None, "")` for a region entirely upstream of
the read. Annotation gives absent regions empty nucleotide/amino-acid strings and null
query intervals. A reviewed FWR3 example remains valid after removing 320 or 370 input
bases: its remaining FWR3 intervals are `[0, 94)` and `[0, 44)`, respectively, with the
original junction preserved. This is distinct from a read with upstream V regions
present but no support for locating FWR3.

**Evidence.** [`get_region_sequence()`](../abstar/annotation/regions.py);
[`test_get_region_sequence_truncated`,
`test_get_region_sequence_missing_region`](../abstar/tests/test_regions.py);
[`test_read_beginning_inside_fwr3_remains_a_valid_partial_region`](../abstar/tests/test_missing_fwr3.py).

### Inserted bases fall exactly between two regions

**Description.** A reference gap at a framework/CDR boundary can cause both regions to
omit the inserted query bases: the upstream region stops at its last reference base,
while the downstream start skips reference gaps.

**Example.** Suppose the nominal cut follows reference `ACG`:

```text
query:       ACGAAA|TTC
reference:   ACG---|TTC
columns:     012345 678
```

**Handling.** `get_aligned_position_from_ungapped(..., is_end_position=True)` extends an
inclusive region endpoint through the immediately following reference-gap run. In this
example, the last upstream reference base maps to column 2, but the region endpoint
becomes column 5; the next reference base maps to column 6. The inserted `AAA` belongs
to the preceding region and is not lost. This nucleotide ownership rule is separate from
the downstream ownership of a codon that crosses an already established nucleotide
boundary during continuous AA projection. The later FWR3 synchronization gives the
established junction anchor priority at the FWR3/CDR3 cut.

**Evidence.** [positions.py](../abstar/annotation/positions.py),
`get_aligned_position_from_ungapped`; [regions.py](../abstar/annotation/regions.py),
`get_region_sequence`. This specific boundary policy is documented in the active
implementation; the illustrative helper outputs are 5 and 6 for the two endpoints above.

### A complete-codon deletion straddles a framework/CDR boundary

**Example.** In an ungapped-reference coordinate system, FWR1 nominally ends at 78 and
CDR1 starts at 78. Deleting reference bases `[76, 79)` splits one three-base deletion
across that boundary:

```text
reference interval:    ... [76 77 | 78] 79 80 | 81 ...
query alignment:       ... [ -  - |  -] 79 80 | 81 ...
nominal boundary:               78
normalized boundary:                         81 (alignment column)
```

**Handling.** `get_region_sequence()` recognizes a contiguous deletion whose full length
is divisible by three even though its portion at the region boundary is not. It moves
the shared alignment boundary forward by three columns until the residual gap is
resolved. In this example, FWR1 uses inclusive alignment columns `0–80`, and CDR1
`81–113`; ungapped lengths are 78 and 33. The two regions preserve query order without
duplicating residues. Start normalization also runs when an FWR3 endpoint is missing and
will subsequently be recovered.

**Evidence.** [`get_region_sequence()`,
`_region_start_after_split_deletion()`](../abstar/annotation/regions.py);
[`test_complete_codon_deletion_spanning_regions_preserves_order`](../abstar/tests/test_regions.py);
[`test_recovered_fwr3_preserves_codon_deletion_boundary`](../abstar/tests/test_missing_fwr3.py).

### Separate one-base insertion/deletion events compensate across a region

**Example.** FWR1 contains a one-base insertion before reference position 70 and a
one-base deletion at reference position 77. Its net length change is zero, but its
aligned final column is a gap.

**Handling.** The complete-codon boundary rebalance above applies only when the entire
contiguous deletion has length divisible by three. It does not move a boundary simply
because its last one or two alignment columns are gaps. In the tested example, FWR1
retains inclusive columns `0–78`, CDR1 starts at column 79, and concatenation reproduces
the original ungapped query exactly. The terminal single-base deletion therefore does
not cause the next region's residues to be copied into FWR1.

**Evidence.** [`get_region_sequence()`](../abstar/annotation/regions.py);
[`test_compensating_single_base_indels_at_region_end_do_not_duplicate_next_region`](../abstar/tests/test_regions.py).

### Retained V evidence ends inside an otherwise present FWR3

**Example.** The retained V hit reaches FWR3 but stops before its final anchor codon,
while the oriented read continues through a supported junction and J region:

```text
query:         ... CDR2 | FWR3................TGC | CDR3 ... J
retained V:    ===================]
region FWR3:            [........................]
```

**Handling.** When ordinary extraction maps the FWR3 start but not its end, annotation
uses that retained start and `junction_start + 3` as the end. It requires
`v_sequence_start <= fwr3_start < v_sequence_end <= junction_start + 3 <=
j_sequence_end`. The sequence comes directly from that oriented-query interval. V calls,
retained V alignment, mutation evidence, and the established junction remain unchanged.
A conflicting interval raises a record annotation error. If neither FWR3 boundary is
mapped despite upstream regions being present, unsupported recovery also raises an error
unless the separately constrained joint-boundary recovery succeeds.

**Evidence.** [`annotate_single_sequence()`](../abstar/annotation/annotator.py);
[`test_missing_fwr3_recovered_from_oriented_query`,
`test_unsupported_fwr3_is_an_explicit_record_error`,
`test_unsupported_fwr3_diagnostics_preserve_other_records`](../abstar/tests/test_missing_fwr3.py).

### FWR3 extraction and the junction disagree about ownership of anchor-adjacent bases

**Example.** A repeated anchor or an insertion immediately after the established anchor
leads ordinary region extraction to end FWR3 at 315, while the established junction
begins at 309. The consistent FWR3/CDR3 cut is 312, immediately after the three-base
anchor.

**Handling.** Annotation synchronizes `fwr3_end` to `junction_start + 3`, which is also
`cdr3_start`, and reslices FWR3 from the oriented query. It validates `v_sequence_start
<= fwr3_start < junction_start + 3 <= j_sequence_end`. It does not relocate the junction
or rewrite retained V evidence to make the region fit. Invalid intervals fail before
slicing. In AIRR's one-based closed coordinates, adjacent FWR3/CDR3 fields consequently
satisfy `fwr3_end + 1 == cdr3_start`.

**Evidence.** [`annotate_single_sequence()`](../abstar/annotation/annotator.py);
[`test_real_fwr3_cdr3_boundary_preserves_junction_and_v_evidence`,
`test_inconsistent_boundaries_are_rejected_before_slicing`,
`test_shared_boundary_survives_public_outputs`](../abstar/tests/test_fwr3_cdr3_boundaries.py).

### A nucleotide region boundary cuts through a translated codon

**Example.** With coding origin zero, `ATGAAATTT` translates to `MKF`. Two adjacent
nucleotide intervals `[0, 4)` and `[4, 9)` contain `ATGA` and `AATTT`. Translating those
strings independently would lose the codon spanning the cut.

**Handling.** All FWR/CDR proteins are slices of the continuous assembled V(D)J
translation. For query interval `[s, e)` and coding origin `v_sequence_start + frame -
1`, both translated boundaries use floor division by three, clamped at zero. Here the
protein slices are `[0, 1)` (`M`) and `[1, 3)` (`KF`): the crossing `AAA` codon belongs
to the downstream region. Leading bases before the coding origin and incomplete terminal
codons do not translate. `cdr3_length` is updated to the resulting `len(cdr3_aa)`, and
amino-acid CDR masks use the same partition. This also applies to out-of-frame
sequences; locally translated `junction_aa` and CDR3 gene subdivisions retain separate
semantics and can differ from continuous-query `cdr3_aa`.

**Evidence.** [`annotate_single_sequence()`](../abstar/annotation/annotator.py);
[`test_real_aa_regions_partition_continuous_query`,
`test_partial_leading_codon_does_not_enter_aa_regions`,
`test_aa_region_partition_public_outputs`,
`test_no_project_python_return_uses_continuous_regions`](../abstar/tests/test_aa_region_boundaries.py).


## Junction anchors, missing boundaries, and repeats

### A gap near the FWR3 anchor would be lost by independently realigning only a short region

**Example.** The full V alignment supports an insertion just before the final FWR3
codon, but a short FWR3-only alignment could prefer mismatches over reopening the gap:

```text
query:      ... GGA ATC TGT
reference:  ... GGA --- TGT
```

**Handling.** Ordinary junction-start mapping includes the full IMGT codon 104, rather
than stopping immediately before it. It first aligns the FWR3 reference slice
`gapped_v[195:312]` without IMGT dots to the already aligned semiglobal V reference,
using gap-open −2 to recover existing gap placement. It then aligns the already aligned
V query to this reconstructed FWR3 reference and converts its endpoint back to ungapped
query coordinates. Including the anchor also supplies context for deletions immediately
before it. An unmappable endpoint activates the conservative raw-query recovery below;
this is not a claim that every nearby indel admits a unique biological interpretation.

**Evidence.** Active FWR3/junction block in
[`annotate_single_sequence()`](../abstar/annotation/annotator.py); the exact biological
and retained-evidence expectations in
[`test_real_fwr3_cdr3_boundary_preserves_junction_and_v_evidence`](../abstar/tests/test_fwr3_cdr3_boundaries.py)
and [`test_lc_anchor_exact_saved_assignment`](../abstar/tests/test_junction_anchor.py).

### The normal FWR3 endpoint is unmappable, but raw query bases uniquely support an anchor

**Example.** Let `U` be a uniquely matching upstream FWR3 sequence. A query window is `U
+ TGG + AAATGT...`, whereas the reference ends `U + TGT`. The homologous observed anchor
is `TGG`, despite the later cysteine codon.

**Handling.** `recover_fwr3_anchor()` requires a complete reference IMGT anchor and a
retained mapping of the upstream FWR3 boundary. It fits the raw query from that boundary
up to, but excluding, `j_sequence_start` against ungapped FWR3 through codon 104.
Leading ends are fixed; the unused query suffix is free. Defaults are match +3, mismatch
−2, gap-open −35, and extension −1. Recovery requires a positive score and the same
complete, contiguous three-base anchor projection across all optimal alignments. Any
already retained anchor-base mapping must agree. Coordinates are oriented-query,
zero-based half-open, and upstream insertions/deletions of one, two, or three bases
shift them without forcing frame preservation.

**Evidence.** [`recover_fwr3_anchor()`, `recover_junction_anchor()`,
`_recover_query_positions()`](../abstar/annotation/junction.py);
[`test_lc_anchor_exact_saved_assignment`,
`test_query_origin_preserves_coordinate_space`,
`test_upstream_insertions_shift_anchor_without_imposing_frame`,
`test_upstream_deletions_shift_anchor_without_imposing_frame`,
`test_recovery_does_not_override_conflicting_retained_anchor`](../abstar/tests/test_junction_anchor.py).

### A recovered anchor is mutated, ambiguous, or inconsistent with productivity

**Example.** The aligned homologous anchor can be `TGA` (stop), `TGG` (tryptophan), or
`TGN` (ambiguous), even when a downstream `TGT` could make a more plausible
cysteine-starting junction.

**Handling.** Recovery chooses coordinates from alignment evidence, with no
cysteine-motif or productivity preference, and retains the observed bases. Ordinary
downstream productivity assessment records the consequences. An exact tested fallback
keeps interval `[365, 398)` for all three examples: `TGA` produces `*AYATDGTLDF` and
stop/conserved-C issues; `TGG` produces `WAYATDGTLDF` and a conserved-C issue; `TGN`
produces `XAYATDGTLDF` and ambiguity/conserved-C issues. Recovery also preserves a
junction frameshift instead of changing its length to a multiple of three. A missing
biological anchor therefore does not license scanning forward to the next cysteine;
depending on alignment support, the homologous non-C bases are retained or recovery
fails.

**Evidence.** [`recover_junction_anchor()`](../abstar/annotation/junction.py);
[`test_real_fallback_preserves_nonproductive_anchor`,
`test_deleted_anchor_does_not_select_a_downstream_cysteine`,
`test_fallback_preserves_junction_frameshift`](../abstar/tests/test_junction_anchor.py).

### Multiple optimal alignments imply different anchor coordinates

**Example.** A query ending in six consecutive `T` bases is fitted to a reference ending
in seven. Several equally scoring placements of the deleted base may place the three
anchor bases differently. In contrast, an upstream seven-versus-eight-base `A` run can
admit several gap placements while leaving the downstream anchor unchanged.

**Handling.** Recovery propagates the requested coordinate projection across optimal
paths, retaining up to two distinct projections to establish ambiguity. Different anchor
projections raise `JunctionAnchorError`; multiple tracebacks with the same anchor
projection are equivalent and can succeed. Deleted, interrupted, or incomplete anchor
codons and nonpositive alignment support also fail explicitly. These errors use the
ordinary record-failure accounting, rather than fabricating a completed annotation.

**Evidence.** [`_best_projection()`, `_recover_query_positions()`,
`recover_junction_anchor()`](../abstar/annotation/junction.py);
[`test_competing_optimal_anchor_projections_are_rejected`,
`test_equivalent_upstream_gap_placements_do_not_make_anchor_ambiguous`,
`test_missing_anchor_or_unsupported_upstream_fails_explicitly`](../abstar/tests/test_junction_anchor.py).

### A short retained V hit does not reach FWR3, or contradicts the proposed FWR3 endpoint

**Example.** The retained V alignment ends inside FWR2. Raw query bases still cover
CDR2, FWR3, and the anchor, but recovering the anchor alone would leave intermediate
region boundaries unsupported.

**Handling.** `recover_v_region_boundaries()` is narrowly triggered when a hit that
begins no later than FWR3 ends at/before the FWR3 start, or when the retained reference
stops short of the anchor while its query end lies past the proposed anchor end. It
starts at the last region boundary that maps to a retained query base, fits through
codon 104 before J, and fixes every retained reference-to-query mapping in that window,
including retained deletions. Every requested downstream region start and all three
anchor bases must have one joint projection across optimal alignments. Intervals must
increase strictly; missing upstream support, incomplete reference anchors, deleted
boundaries, interrupted anchors, and overlapping reference boundaries fail explicitly.
Earlier retained regions and assignment/mutation evidence are preserved.

**Example with exact offsets.** Suppose the reference has no IMGT dots and the query
matches it without indels through the anchor, after a prefix of length 11. With
retained reference bases `[0, 150)`, the last mapped region start is FWR2 at reference
offset 114. The recovered intervals are FWR2 `[125, 176)`, CDR2 `[176, 206)`, FWR3
`[206, 323)`, with anchor `[320, 323)`.

**Evidence.** [`recover_v_region_boundaries()`](../abstar/annotation/junction.py);
[`test_extended_regions_do_not_prefer_productive_anchor_codons`,
`test_recovery_requires_unique_boundaries_and_preserves_assignment`,
`test_missing_mapping_requires_valid_reference_and_upstream_support`,
`test_deleted_boundary_in_retained_evidence_is_not_fabricated`,
`test_extended_mapping_rejects_interrupted_anchor_codon`](../abstar/tests/test_region_boundary_recovery.py).

### The anchor is unique, but intervening region boundaries remain ambiguous

**Example.** An upstream homopolymer allows two equally supported CDR2/FWR3 boundaries
while both paths agree on the anchor. A reviewed case has competing joint projections
`(244, 265, 379, 380, 381)` and `(247, 268, 379, 380, 381)`.

**Handling.** Joint recovery rejects this record even though anchor bases 379–381 agree.
It does not independently select a favorable boundary for each region. Fixed retained
mappings can resolve an otherwise ambiguous repeat placement, but unresolved competing
joint projections are logged as a record failure. Public API/CLI regressions preserve
the other two recoverable records in the three-record fixture and the failed record's
diagnostic; strict mode aborts with the recorded failure.

**Evidence.** [`_recover_query_positions()`,
`recover_v_region_boundaries()`](../abstar/annotation/junction.py);
[`test_region_tie_is_not_hidden_by_a_unique_anchor`,
`test_retained_mapping_constrains_competing_repeat_placements`,
`test_public_recovery_preserves_ambiguous_failure_and_other_records`,
`test_strict_recovery_still_aborts_on_ambiguous_regions`](../abstar/tests/test_region_boundary_recovery.py).

### Anchor recovery must use the correct receptor reference and strand

**Example.** A reverse-complement TRA read, or a TRB/TRD/TRG read, needs its own V
reference and oriented-query coordinates; a BCR reference or original-input offset would
misplace the anchor.

**Handling.** Annotation orients the input once, infers receptor/locus from the assigned
V call, and propagates receptor identity through reference lookup and realignment. Both
recovery mechanisms use the supplied receptor's actual IMGT-gapped V reference. TCR
regressions cover TRA, TRB, TRD, and TRG; the basic fallback is deliberately forced in
those controls, and the broader recovery is exercised with artificially shortened
retained evidence. These establish recovery behavior for the controls, not evidence that
every natural TCR failure is recoverable. Reverse-complement/prefix regressions check
that adding a five-base oriented prefix adds five to junction coordinates without
changing junction sequence.

**Evidence.** [`annotate_single_sequence()`](../abstar/annotation/annotator.py),
[`recover_fwr3_anchor()`,
`recover_v_region_boundaries()`](../abstar/annotation/junction.py);
[`test_fallback_projects_supported_tcr_anchors`,
`test_recovered_anchor_reverse_complement_and_prefix`](../abstar/tests/test_junction_anchor.py);
[`test_broader_recovery_uses_tcr_reference_boundaries`](../abstar/tests/test_region_boundary_recovery.py).

### A downstream J-like repeat wins a local J alignment

**Example.** In a reviewed IGHJ5 case, an unrestricted downstream local hit occupies
`[487, 528)`, wholly beyond the primary junction end at 458. The primary retained J
interval is `[455, 487)`.

**Handling.** Ordinary FWR4 mapping preserves gap context from the semiglobal J
reference and searches only the query interval from the established junction start to
the current J end. If the subsequently considered local J hit starts at or after
`junction_end`, annotation bounds a new local search at `junction_end - 3 +
len(germ_fr4_sequence)` and reprocesses the primary J boundaries. FWR4 reference length
is 34 bases for IGH/TRA/TRD and 31 for the other supported loci. In the exact example,
the final reference interval is `[17, 49)` and the retained query/reference both equal
`TGGGGCCAGGGAACCCTGGTCACCGTCTCCTC`.

**Evidence.** [`annotate_single_sequence()`](../abstar/annotation/annotator.py);
[`test_primary_j_boundary_rejects_downstream_j5_repeat`](../abstar/tests/test_real_bcr.py).

### Final J boundary selection or read truncation shortens FWR4

**Example.** FWR4 was located using a larger semiglobal alignment, but the retained J
endpoint is subsequently shortened by primary-J selection or by a 3′-truncated read.

**Handling.** Final region assembly sets `fwr4_start = junction_end - 3` and `fwr4_end =
j_sequence_end`, then takes that exact oriented-query slice. The FWR4 locator
establishes the anchor; it does not supply a stale final endpoint. In a tested truncated
kappa read, final J is `[393, 407)`, FWR4 is `TTCGGCCAA` (`FGQ`), and the junction is
unchanged. Tests check native and AIRR serialization, forward/reverse inputs, retained
gene evidence, and region/mask assembly.

**Evidence.** [`annotate_single_sequence()`](../abstar/annotation/annotator.py);
[`test_fwr4_uses_final_j_endpoint_without_changing_gene_evidence`,
`test_truncated_j_preserves_retained_endpoint`,
`test_final_fwr4_endpoint_survives_public_serialization`](../abstar/tests/test_fwr4_endpoints.py).

### CDR3 or D-like motifs occur more than once in the assembled sequence

**Example.** An explanatory assembled sequence `CCCXXXAAABBBCCCJJJ` contains `CCC` at
offsets 0 and 12. Given `v_sequence_start = 100`, `junction_start = 109`, and
`junction_end = 118`, the CDR3 is `[12, 15)`, the second `CCC`. With
`v_sequence_end = j_sequence_start = 112`, that entire CDR3 lies in retained J.

**Handling.** `identify_cdr3_regions()` computes assembled-sequence coordinates by
subtracting `v_sequence_start` from oriented-query coordinates; it does not call a
substring search. It validates that the resulting CDR3 slice lies inside the assembly
and exactly reproduces `ab.cdr3`. D subdivisions likewise derive from retained D
coordinates, so repeated sequence elsewhere in V cannot steal their positions. The
example assigns `CCC` to CDR3 J and leaves CDR3 V/N1 empty. The letters `X`, `B`, and
`J` in this toy illustrate positions and are not proposed biological input.

**Evidence.** [`identify_cdr3_regions()`](../abstar/annotation/regions.py);
[`test_identify_cdr3_uses_alignment_coordinates_with_repeated_motifs`](../abstar/tests/test_regions.py).

### A V, D, or J segment contributes no complete codon to CDR3

**Example.** Retained V ends at the conserved junction-start codon, J begins only at the
final junction anchor, or a two-base D segment falls between them. The gene evidence can
exist even when its CDR3 amino-acid contribution is empty.

**Handling.** CDR3 subdivisions use the local CDR3 frame. V contribution is trimmed to
complete codons at its 3′ end and clamped to empty if it ends before CDR3 starts. J's
beginning advances to the next CDR3-frame boundary and its contribution is clamped to
empty when necessary. D begins/ends at retained coordinates and trims partial codons at
both ends; surrounding bases enter the N1/N2 subdivisions. These local subdivisions must
not be used to reconstruct continuous-query region proteins on out-of-frame reads.
Independently, a nonempty short D retains nucleotide evidence and nucleotide identity
even if translation is empty; an actually empty D alignment clears D-specific evidence
and assigns the entire V–J gap to NP1.

**Evidence.** [`identify_cdr3_regions()`](../abstar/annotation/regions.py);
[`_segment_identities()`, `_clear_empty_d_alignment()`,
`annotate_single_sequence()`](../abstar/annotation/annotator.py);
[`test_identify_cdr3_uses_alignment_coordinates_with_repeated_motifs`](../abstar/tests/test_regions.py);
[`test_pilot_loss_empty_d_clears_stale_evidence`,
`test_pilot_loss_empty_d_translation_preserves_nt_identity`,
`test_pilot_loss_empty_d_translation_retains_nucleotide_evidence`](../abstar/tests/test_real_bcr.py).
The V/J codon trimming rules are source-established; the named repeated-motif test
directly exercises the empty-V case.

## Productivity and mask consistency

### Stops, ambiguity, missing evidence, and locus-specific junction motifs

**Description.** An apparently annotated record may contain a stop codon, a non-ACGT
base, incompatible V/J calls, or a missing/mutated junction anchor.

**Examples.** `CAF` is accepted by the motif check for each supported TCR locus; `CAW`
is rejected there with `junction does not end with conserved F`. Heavy-chain `CAW` has
the expected terminal residue. A query containing `R`, `Y`, or `N` fails the
unambiguous-DNA check even when translation lacks `*`. The standalone productivity
helper also rejects `U` and dashes; public nucleotide input validation rejects those
earlier.

**Handling and limits.** `assess_productivity` accumulates distinct issues, preserves
earlier issues, and serializes them with `|`. It flags missing AA sequence, translated
`*` (also sets `stop_codon=True`), absent sequence or any non-ACGT nucleotide after
uppercasing, missing V/J calls, different three-letter V/J loci, and missing/truncated
junction AA. The junction must begin with C and end with W for IGH, or F for IGK, IGL,
TRA, TRB, TRD, and TRG; unsupported motif loci are issues. When junction NT is supplied,
fewer than six bases, nonmultiple-of-three length, and ambiguity are explicitly
recorded. Any issue yields `productive=False`. These are software productivity criteria,
not an expression or binding assay. The standalone helper skips NT-specific checks if
`junction is None` and frame-specific checks if `frame is None`; partially populated
objects therefore do not receive every check.

**Evidence.** [productivity.py](../abstar/annotation/productivity.py),
`assess_productivity`, `JUNCTION_MOTIFS`;
[test_productivity.py](../abstar/tests/test_productivity.py),
`test_stop_codon_antibody`, `test_locus_mismatch_antibody`,
`test_all_non_acgt_bases_are_ambiguous`,
`test_empty_or_truncated_junction_is_nonproductive`,
`test_tcr_junction_requires_terminal_phenylalanine`,
`test_tcr_junction_rejects_terminal_tryptophan`.

### Junction length and phase are separate checks; a 5′ prefix must not change phase

**Description.** A nine-base junction is a multiple of three but can still begin out of
phase with the translated V region. Conversely, adding upstream sequence before the V
should not change productivity.

**Example.** With V start `137`, junction start `437`, and slice frame `1`, the
difference is `300` and the junction is in frame. Moving only the junction start to
`438` makes it out of frame. A junction `TGTAAAT` has length seven and fails regardless
of its start.

**Handling and limits.** Valid frames are 1, 2, and 3. Phase uses `(junction_start -
v_sequence_start - (frame - 1)) % 3 == 0` in oriented-query coordinates. Invalid frames,
out-of-phase starts, or invalid junction lengths add explicit issues and set
`vj_in_frame=False`. Motif and stop-codon failures can coexist with `vj_in_frame=True`;
the frame flag is not synonymous with `productive`.

**Evidence.** [productivity.py](../abstar/annotation/productivity.py),
`junction_is_in_frame`, `assess_productivity`;
[test_productivity.py](../abstar/tests/test_productivity.py),
`test_junction_length_must_be_a_multiple_of_three`,
`test_junction_frame_uses_v_region_origin`,
`test_productivity_uses_v_region_frame_origin`,
`test_invalid_frame_remains_nonproductive_with_v_origin`.

### Repeated sequence motifs and absent D segments cannot shift gene masks

**Description.** Searching for a CDR3 string can find an earlier repeated copy; using
CDR3-only subdivisions can omit bases outside the CDRs. Light chains and some TCR loci
have no D segment.

**Example.**

```text
assembled NT:  AAAA CC GGGG T CCCC
segment mask:  VVVV NN DDDD N JJJJ
```

For V=`AAA`, no N/D, J=`TTT`, the mask is `VVVJJJ`.

**Handling and limits.** Gene masks concatenate the actual assembled V, NP1, optional D,
optional NP2, and J span lengths, in that order. They never search for a CDR3 substring.
Missing D/NP2 contribute zero length. The total must equal ungapped assembled sequence
length or a `ValueError` is raised. CDR masks instead repeat labels for the established
region lengths (`0` for all FWRs, including FWR4; `1`, `2`, `3` for CDR1–3). Mask
generation labels existing intervals; it does not validate their biological placement.
The gene-mask check validates length, not equality of segment concatenation to sequence
content.

**Evidence.** [mask.py](../abstar/annotation/mask.py), `_generate_gene_segment_mask_nt`,
`_generate_cdr_mask_nt`; [test_mask.py](../abstar/tests/test_mask.py),
`test_gene_segment_mask_does_not_search_for_repeated_cdr3`,
`test_gene_segment_mask_uses_complete_assembled_spans`,
`test_gene_segment_mask_rejects_inconsistent_assembly`, `test_generate_cdr_mask_nt`.

### Codons spanning segment boundaries and untranslated partial codons

**Description.** An AA can receive bases from more than one gene segment, even when
there are no nucleotide N additions.

**Example.** NT segment labels `VVV VNN DDD DNJ JJJ` in frame 1 produce AA labels
`VNDNJ`. A V/D boundary codon `VDD` is labeled `N` even if every nucleotide matches its
assigned germline.

**Handling and limits.** The AA gene mask walks complete three-base windows starting at
`frame - 1` in the assembled NT mask. A codon receives V/D/J only if all three labels
agree; any mixed codon is `N`. Leading bases before the frame and terminal incomplete
codons do not produce labels. Result length must equal ungapped assembled AA length.
Thus `N` in an AA segment mask means a mixed or N-containing codon, not necessarily an
entire amino acid encoded by newly added nucleotides.

**Evidence.** [mask.py](../abstar/annotation/mask.py), `_generate_gene_segment_mask_aa`;
[test_mask.py](../abstar/tests/test_mask.py), `test_generate_gene_segment_mask_aa`
(frames 1, 2, and 3), `test_gene_segment_mask_uses_complete_assembled_spans`.

### Query deletions and N additions in non-germline masks

**Description.** Alignment rows contain deletion columns that have no query base to
label. N regions must remain non-germline even if a placeholder aligned reference
happens to match.

**Example.** Query `A-`, reference `AT`, segment mask `V` yields non-germline mask `0`;
no extra character is added for the deleted T. With identical query/reference
`AAAACCGGGGTCCCC` and segment mask `VVVVNNDDDDNJJJJ`, the result remains
`000011000010000`.

**Handling and limits.** `generate_nongermline_mask` first requires equal alignment
lengths and a gene mask matching the ungapped query length. It skips query-gap columns
before indexing the gene mask; N labels and query/reference mismatches (including
insertion columns) become `1`, and remaining matches become `0`. Skipping terminal
deletions avoids indexing past an already consumed mask. The same algorithm works on AA
alignments. This mask reports alignment disagreement or N/mixed origin; it is not a
probability of somatic mutation.

**Evidence.** [mask.py](../abstar/annotation/mask.py), `generate_nongermline_mask`;
[test_mask.py](../abstar/tests/test_mask.py),
`test_nongermline_mask_handles_terminal_deletion_after_mask_is_consumed`,
`test_nongermline_mask_rejects_short_segment_mask`, `test_generate_nongermline_mask_nt`,
`test_generate_nongermline_mask_aa`.

## Input collections and record identity

### Opaque identifiers, duplicate identifiers, and asynchronous completion

**Description.** An identifier can look numeric (`00123`) or scientific (`10E8`), and
multiple biological records can share the same identifier. Running chunks in parallel
must not make identifiers join keys or reorder results by completion time.

**Example.** Four reads arrive in this order:

```text
input ordinal    sequence_id    internal identity
0                10E8           abstar_0_0
1                00123          abstar_0_1
2                duplicate      abstar_0_2
3                duplicate      abstar_0_3

worker completion: chunk 2, chunk 0, chunk 3, chunk 1 (one record per chunk)
published order:   original record ordinals 0, 1, 2, 3
```

**Handling.** Assignment gives each record an immutable sample/record ordinal key and
carries the external ID as a string. Annotation results are restored to numeric
sample/record order; duplicate internal keys are errors. Public objects, AIRR, and
Parquet omit `row_id`, but failure diagnostics retain it to distinguish duplicates. This
handles duplicate labels, not a deduplication or clonotyping request.

**Evidence.** [`MMseqs.prepare_input_files`](../abstar/assigners/mmseqs.py),
[`_sort_annotation_workframe` and `run`](../abstar/core/abstar.py);
[`test_run_conserves_duplicate_opaque_ids_and_order_across_workers`](../abstar/tests/test_pipeline.py)
asserts exact IDs/statuses for `(n_processes, chunksize)` values `(1,1), (1,3), (2,1),
(2,3)`;
[`test_api_input_and_worker_matrix_matches_authenticated_serial_baseline`](../abstar/tests/test_public_e2e.py)
compares annotations and order across public input/worker forms.

### One-shot iterators and data-dependent output counts

**Description.** A generator may be consumable only once, and a multi-record input may
contain only one successful annotation or no successful annotations.

**Example.** A generator yields `[valid antibody, N]`. The second record has no V
assignment. The API returns a two-element list, including the explicit unassigned
record. If one of two records instead raises a recoverable annotation exception, its
omission does not turn the other result into a scalar.

**Handling.** `_process_inputs` deliberately materializes an iterable once and accepts
an iterable of `abutils.Sequence` objects. Return shape follows the original record
count: one retained record from a single-record input gives a `Sequence`; multi-record
inputs give a list, even when fewer records survive. `as_dataframe=True` gives a Polars
DataFrame; a project-backed run writes outputs and returns `None`. All-failed
recoverable annotation gives an empty list/DataFrame accompanied by a warning and
persistent diagnostics. Empty input itself is rejected.

**Limits.** This iterable branch requires `Sequence` elements, rather than arbitrary
file paths or raw strings. Materialization is deliberate and does not imply
bounded-memory streaming.

**Evidence.** [`_process_inputs` and final return logic in
`run`](../abstar/core/abstar.py); [`test_api_arbitrary_iterable_is_consumed_once`,
`test_api_one_record_return_shape`,
`test_api_mixed_outcomes_keep_multiple_input_return_shape`](../abstar/tests/test_public_e2e.py);
[`test_no_project_all_failed_retains_diagnostics`](../abstar/tests/test_recoverable_annotation.py).

### Repeated filenames in nested sample directories

**Description.** Different samples can both be named `reads.fasta`; suffix differences
can also produce the same apparent stem.

**Example.** `batch/sample1/reads.fasta` and `batch/sample2/reads.fasta` become
`sample1__reads` and `sample2__reads` when the common root is `batch`. A further
collision is resolved with `__2`, `__3`, etc.

**Handling.** Directory discovery is recursive and naturally sorted. Unique stems remain
simple; colliding stems use paths relative to the common root, with separators replaced
by `__`, followed by deterministic numeric disambiguation if needed. With
`copy_inputs_to_project=True`, source-relative subdirectories are preserved instead of
flattening copies.

**Limits.** `run` directory discovery explicitly lists `fasta`, `fa`, `fastq`, `fq`,
`fasta.gz`, and `fastq.gz`; other compressed suffix variants are not all included in
directory discovery. Output names depend on the whole discovered collection/common root.

**Evidence.** [`_process_inputs`, `_get_sample_names`,
`_copy_inputs_to_project`](../abstar/core/abstar.py);
[`test_sample_names_are_unique_for_nested_and_extension_collisions`,
`test_copy_inputs_preserves_relative_directories`](../abstar/tests/test_pipeline.py);
[`test_api_project_returns_none_and_preserves_nested_input_paths`](../abstar/tests/test_public_e2e.py).

### Empty files versus malformed or truncated sequence records

**Description.** An empty sample is different from a nonempty file that cannot be parsed
or contains invalid sequence characters. FASTQ sequences and qualities may span multiple
lines.

**Example.** A directory contains `blank.fastq` holding only whitespace and
`valid.fastq` holding an antibody read. Only the valid sample is annotated. Conversely:

```text
@00123
ACGT
+
III          # three qualities for four nucleotides: invalid
```

The example's explanatory comment is not literal FASTQ data. A valid multiline
equivalent can split `ACGT` as `AC`/`GT` and `IIII` as `II`/`II`.

**Handling.** Whitespace-only files are filtered out before processing; a wholly empty
collection raises `ValueError`. Nonempty malformed files are retained for input
validation. The MMseqs input preparation validates all FASTQ records through Biopython's
strict `FastqGeneralIterator` before normal parsing, so a malformed later record is not
silently dropped. Non-IUPAC sequence characters raise an input error. The public
controller reports these as `preprocess/invalid_input`, writes a diagnostic, and does
not publish successful sample output.

**Limits.** This does not repair malformed FASTQ or infer missing qualities; corruption
is rejected. Record conservation below begins with parsed annotation input (after
optional read merging), not all raw paired-end reads.

**Evidence.** [`_has_input_content`, `_process_inputs`,
`run`](../abstar/core/abstar.py),
[`MMseqs.prepare_input_files`](../abstar/assigners/mmseqs.py);
[`test_api_mixed_empty_and_valid_files_keep_the_valid_annotation`,
`test_api_nonempty_malformed_files_keep_structured_failure_artifacts`](../abstar/tests/test_public_e2e.py);
[`test_invalid_sequence_content_retains_distinct_input_diagnostic`,
`test_public_valid_multiline_fastq_retains_annotation`](../abstar/tests/test_failure_contracts.py).

## Preprocessing variation

### End-specific UMIs with anchors, offsets, and reverse orientation

**Description.** A UMI can be specified by a fixed terminal length, an adjacent
conserved sequence, or flanking conserved sequences. Extra bases can shift an anchor
away from the literal read end.

**Example.** `ATGC AAAA` with pattern `ATGC[UMI]`, length `4` returns `AAAA`. `ACGTAC
TTGGCC` with `[UMI]TTGGCC` and no length infers `ACGTAC`. With input `GGGGTTTTGCAT`,
pattern `ATGC[UMI]`, length `-4`, reverse complementation exposes `ATGCAAAACCCC` and
gives `AAAA`.

**Handling.** Without a pattern, positive/negative lengths slice the start/end of the
original uppercased read. With a pattern, a negative length selects the
reverse-complement view unless `ignore_strand=True`. Single-flank patterns use the
absolute requested length when supplied. Two-flank patterns instead extract the entire
interval between the aligned flanks: six intervening bases yield a six-base UMI even
if the requested length is four. That length still affects the search window.
Conserved flanks use local alignment and mismatch accounting against all conserved
pattern positions. Search is restricted to an end window containing
conserved-flank length, requested UMI length (if any), and `extra_length_for_alignment`
(default 25). A trailing anchor can determine a UMI with unspecified length; a pattern
ending in `[UMI]` requires a length. Multiple accepted UMI components are joined with
`+`.

**Limits.** Patternless negative length means a suffix slice, not reverse
complementation. This is terminal UMI extraction, not an unrestricted internal-motif
search. An anchor at offset 60 is outside the default window in the regression and is
found when the extra window is expanded to 60. Invalid/missing pattern/length
specifications are rejected; no inference of an unbounded terminal UMI.

**Evidence.** [`UMI.__init__`, `process_sequence`, `get_umi`, `get_mismatches`,
`parse_umis`](../abstar/annotation/umi.py);
[`test_negative_pattern_length_uses_absolute_slice_after_reverse_complement`,
`test_pattern_without_length_infers_umi_before_trailing_anchor`,
`test_pattern_ending_in_umi_requires_length`,
`test_pattern_search_is_limited_to_sequence_end`,
`test_end_search_window_can_be_extended`,
`test_parse_umis_builtin_pattern_defaults`](../abstar/tests/test_umi.py). The two-flank
length behavior above is established by the active `get_umi` implementation.

### Reads lacking an acceptable UMI

**Description.** A UMI may be absent, incomplete, or fail the conserved-pattern mismatch
threshold even while the immunoglobulin sequence remains annotatable.

**Example.** For `ATGC[UMI]`, length `4`, mismatch allowance `0`, inputs `ATGCAAAA` and
`TTTTCCCC` yield UMI values `AAAA` and `None`. Both records remain present.

**Handling.** Single-sequence `parse_umis` returns `None` when nothing is detected. The
iterable helper retains each record and adds `umi=None` for a miss; the file helper
retains its unchanged identifier when no UMI is found. `run(..., umi_pattern=...,
umi_length=...)` retains the antibody annotation with a null UMI. Successful standalone
file extraction appends the UMI to the output identifier; the main annotation pipeline
preserves the original identifier and writes the UMI field. Neither path here removes
the UMI bases from the source sequence. File extraction defaults to a sibling
`.umis.<format>` output and rejects an explicitly identical input/output path.

**Evidence.** [`_parse_umis_from_single_sequence`, `_parse_umis_from_sequences`,
`_parse_umis_from_file`](../abstar/annotation/umi.py),
[`annotate_single_sequence`](../abstar/annotation/annotator.py);
[`test_iterable_retains_record_without_detected_umi`,
`test_file_retains_record_without_detected_umi`,
`test_file_rejects_explicit_in_place_output`,
`test_annotation_pipeline_conserves_umi_present_and_absent_records`](../abstar/tests/test_umi.py).
The last regression asserts exact IDs, `["ACGT", None]`, both annotated statuses, and
unchanged V/D/J/C calls versus baseline.

### Paired-end reads split across lanes or interleaved in one file

**Description.** A sample can have R1/R2 files for multiple sequencing lanes, index-read
files in the same directory, or interleaved paired records.

**Example.** `sample_S1_L001_R1_001.fastq`, its matching R2, and a second L002 pair form
one sample. `sample_S1_L001_I1_001.fastq` is not a mate. A lane with two R1 files and
one R2 is ambiguous.

**Handling.** Illumina and Element naming parsers identify samples/reads/lanes. Index
reads are filtered out; every represented lane must have exactly one R1 and one R2,
otherwise pairing raises `ValueError`. Lanes are naturally sorted, merged independently
through fastp, then concatenated into one sample result. Interleaved inputs are
normalized to temporary one-line sequence/quality records before fastp and temporary
input is removed afterward. The paired and interleaved regression fixtures each yield
exactly their two expected merged sequences, IDs, and order.

**Limits.** Public `merge_fastqs` supports fastp; an old standalone vsearch function is
not evidence that the active public selector supports it. Mergeability, overlap mismatch
thresholds, adapter trimming, and quality filtering remain fastp-controlled; raw pairs
are not guaranteed to produce annotation records, and unmerged reads are not guaranteed
to be retained. Filename validation does not establish biological mate correctness.

**Evidence.** [`group_paired_fastqs`, `MergeGroup._validate_pairs`, `MergeGroup.merge`,
`merge_fastqs`](../abstar/preprocess/merging.py);
[`test_group_paired_fastqs_rejects_missing_or_duplicate_reads`,
`test_fastp_paired_merge_conserves_records_and_order`,
`test_fastp_interleaved_merge_conserves_records_and_cleans_temporary_input`,
`test_interleaved_merge_normalizes_text_and_cleans_temporary_file`](../abstar/tests/test_merging.py).

## Explicit outcomes and recoverability

### Nonassignment, record annotation errors, and infrastructure failures

**Description.** A short/ambiguous read may have no assignable V gene; a V-only read may
lack J; an assigned record can fail during detailed annotation; an entire worker/tool
can fail. These are different outcomes.

**Example.** In a three-record sample: read A annotates, read B is `N`, and read C
triggers an annotation exception. Default output contains A and an explicit unassigned
B, while C is counted in `failures.tsv`. A following sample still runs.

```text
input record
  ├─ annotated ───────────────────────────> output row
  ├─ no compatible V/J ──────────────────> unassigned output row
  └─ record annotation exception ────────> diagnostic + indexed failure

input count = annotated + unassigned + indexed record failures
```

**Handling.** Biological nonassignment is represented by
`annotation_status="unassigned"`, a specific `failure_reason`, and null unsupported
annotation/productivity values. A V-only record retains its V call and reports `no
compatible J gene assignment`. Recoverable annotation exceptions omit incomplete rows,
persist diagnostics, warn API callers, report counts, and continue across
chunks/samples. `strict=True`/`--strict` raises `AnnotationRunError` for record errors
before final output for the failing sample; default CLI recovery exits successfully
while reporting failures. Record conservation is explicitly checked per sample and per
run.

**Limits.** Recovery describes failures caught around individual annotation operations,
not biological correction. `OSError` and `MemoryError` bypass record recovery; worker
crashes, external-tool errors, and output/diagnostic-storage failures are fatal. A
worker crash is represented with failure entries for its input records and surviving
partial work is exposed, rather than presented as a successful sample. An all-failed
sample can produce an empty table only with explicit failure accounting and persistent
diagnostics.

**Evidence.** [`annotate`](../abstar/annotation/annotator.py), [`run`,
`_assert_record_conservation`, `_chunk_failure_result`](../abstar/core/abstar.py);
[`test_single_unassignable_record_returns_explicit_sequence_outcome`,
`test_v_assigned_j_unassigned_record_skips_empty_downstream_search`](../abstar/tests/test_pipeline.py);
[`test_failed_records_do_not_discard_chunk_or_later_samples`,
`test_no_project_all_failed_retains_diagnostics`,
`test_strict_mode_still_raises_for_record_error`,
`test_cli_strict_mode_keeps_nonzero_exit_and_diagnostics`,
`test_infrastructure_errors_are_not_record_failures`](../abstar/tests/test_recoverable_annotation.py);
[`test_public_mixed_worker_failures_account_for_every_row`](../abstar/tests/test_failure_contracts.py).

### Failed records with duplicate, path-like, or very long IDs, including reruns

**Description.** The same failed ID can occur twice; an ID can contain `/`, `%`, or
hundreds of characters. Diagnosing a rerun must not erase earlier evidence.

**Example.** `duplicate` appears twice in `sample.fasta`, while another ID is
`../../escape`. The two duplicates receive distinct `...__abstar_0_<ordinal>.failed`
paths, and the path-like ID is encoded as one filename component.

**Handling.** Diagnostics encode path syntax and bound long components with a hash
suffix. Each filename combines the encoded external ID and internal row key.
`logs/failures.tsv` retains exact external IDs, input paths, categories, exception
details, and diagnostic locations. `logs/run.json` records parameters, package versions,
and source hashes. Reruns allocate fresh sample log namespaces and archive prior failure
indexes/run metadata. No-project API cleanup retains failure logs; warnings point to the
persistent failure index, even when all records fail. Promotion of the warning to an
exception does not erase those files.

**Limits.** Safe diagnostic filenames do not change the original identifier in the
index/record. This retention policy is for diagnostics; historical successful output
versions are not all archived.

**Evidence.** [`safe_component`, `failure_path`, `sample_directory`,
`initialize_diagnostics`, `index_failures`](../abstar/core/diagnostics.py),
[`_project_workspace`](../abstar/core/abstar.py);
[`test_record_diagnostics_encode_ids_and_preserve_duplicate_records`,
`test_diagnostic_reruns_preserve_prior_index_and_record`,
`test_promoted_warning_retains_no_project_diagnostics`](../abstar/tests/test_recoverable_annotation.py).

### Failed external tools and failed output publication

**Description.** A merge/MMseqs executable can fail, or output storage can become
unwritable after earlier samples completed.

**Example.** Sample 1 successfully writes AIRR/Parquet; sample 2 fails to serialize. The
raised error lists surviving sample 1 files and current diagnostics as partial results.
A failed fastp rerun leaves a preexisting merged destination intact.

**Handling.** External tools use checked argument lists, preserving
command/return-code/stdout/stderr evidence. fastp writes into staging and replaces the
destination only after success. Final AIRR/Parquet outputs are both staged before
publication; surviving promoted files are listed if a later promotion fails.
`AnnotationRunError` carries structured failures and surviving `partial_output_paths`.
If project diagnostics cannot be written, a separate temporary diagnostic location is
attempted; diagnostic retention failures do not replace the original error.

**Limits.** Two final formats are not a single atomic transaction: publication is one
file at a time, and a failure after the first promotion can leave that file as partial
output. The whole multilane merge and entire multi-sample run are not globally atomic.
There is no successful recovery when required storage itself fails.

**Evidence.** [`merge_fastqs_fastp`](../abstar/preprocess/merging.py),
[`_write_sample_outputs`, `_raise_pipeline_failure`,
`_project_workspace`](../abstar/core/abstar.py);
[`test_public_merge_failure_preserves_existing_destination_and_cleans_staging`,
`test_run_reports_fastp_failure_as_structured_preprocess_error`](../abstar/tests/test_merging.py);
[`test_public_output_failure_retains_diagnostics_without_final_success`,
`test_later_sample_failure_reports_earlier_files_as_partial`,
`test_initial_project_storage_failure_is_structured_and_preserves_caller_file`,
`test_fallback_diagnostic_failure_preserves_original_output_error`](../abstar/tests/test_failure_contracts.py).

## Strand, alignment, and serialization

### Reverse-strand reads, leading context, and file/API sequence meanings

**Description.** An input read can be reverse-complemented relative to V(D)J, and the
read can contain bases before/after the retained rearrangement.

**Example.** Input `AACG` has oriented query `CGTT` when `rev_comp=True`. Internal
oriented interval `[1,4)` addresses `GTT`; AIRR writes interval `[2,4]` while keeping
`sequence="AACG"`. This is a coordinate illustration, not an annotatable four-base
antibody.

**Handling.** Final AIRR and Parquet `sequence` preserve the full original input;
annotation coordinates and alignments address `sequence_oriented`. AIRR performs a
one-time zero-based half-open to one-based closed conversion; native Parquet retains
zero-based half-open intervals and native booleans/nulls. No-project API
objects/DataFrames deliberately retain legacy assembled V(D)J meanings for `sequence`
and the legacy amino-acid fields. Output normalization tests compare every public field
under the documented file representation conversions.

**Limits.** AIRR intervals cannot be directly sliced against the original reverse-strand
`sequence`, and Parquet intervals must not be converted twice. External AIRR readers may
normalize coordinates on loading (the tests' AIRR reader does), so distinguish bytes on
disk from reader-returned values. API/file sequence meanings are an intentional
compatibility difference, not universal identical serialization.

**Evidence.** [`to_airr_interval`, `to_airr_row`](../abstar/annotation/airr.py),
[`_write_sample_outputs`](../abstar/core/abstar.py);
[`test_row_conversion_is_pure_and_converts_every_coordinate`,
`test_airr_region_coordinates_slice_the_oriented_query`,
`test_normalized_tsv_and_parquet_agree_every_public_field`](../abstar/tests/test_airr.py);
[`test_no_project_api_keeps_internal_sequence_meanings_and_return_shapes`](../abstar/tests/test_output_parity.py);
[`test_unassigned_reverse_strand_preserves_oriented_sequence`](../abstar/tests/test_failure_contracts.py).

### Null/partial coordinates, quoted IDs, and forbidden TSV delimiters

**Description.** Some regions are absent; identifiers may contain quotation marks; a
tab/newline inside a field would corrupt a TSV record.

**Example.** Internal `[0,1)` becomes AIRR `[1,1]`. `(None,9)` becomes a null pair
rather than a fabricated interval. ID `clone"A` round-trips through ordinary CSV-style
quote escaping; an ID containing a tab is rejected.

**Handling.** Missing either interval endpoint yields two nulls. Non-null intervals must
be integer, nonnegative, and nonempty with end greater than start. AIRR emits `T`/`F`,
empty nulls, LF line endings, and properly escaped literal quotes. Every row is
converted and delimiter-validated before the destination is opened, preventing a bad
later field from truncating an existing TSV.

**Limits.** This is validation, not sanitization of IDs. Embedded TSV delimiters cannot
be exported in a field. Null productivity represents unavailable evidence, distinct from
`False`.

**Evidence.** [`to_airr_interval`, `write_airr_tsv`](../abstar/annotation/airr.py);
[`test_half_open_to_closed`, `test_invalid_intervals_raise`,
`test_writer_boolean_null_lf_and_reference_validation`,
`test_literal_quotes_round_trip_all_fields_and_rows`,
`test_delimiter_prevalidation_does_not_create_or_truncate_file`](../abstar/tests/test_airr.py).

### CIGARs for retained indels and absent germline bases in NP regions

**Description.** Insertions/deletions and non-templated junction bases require explicit
alignment evidence; a mismatch is not an insertion.

**Example.** Retained query `AC-GT`, germline `ACCGT`, query start `2`, germline start
`1` produces `2S1N2M1D2M`. Here `S` skips unaligned query prefix and `N` skips the
reference prefix. A non-templated `AAA` between retained segments contributes `AAA` to
query alignment and `---` to germline alignment.

**Handling.** `build_cigar` compresses retained columns into `M`, `I`, or `D`, with
prefix query clips/reference skips. `M` includes both matches and mismatches. Unequal
alignment lengths, negative/noninteger origins, or a double-gap column raise errors. The
assembled V(D)J alignment pads NP bases with germline gaps rather than inventing
templated residues. Independent test replay reconstructs exact retained segment
sequences, intervals, identities, and assembled NP columns.

**Evidence.** [`build_cigar`](../abstar/annotation/airr.py),
[`annotate_single_sequence`](../abstar/annotation/annotator.py);
[`test_cigar_run_lengths`, `test_invalid_cigar_evidence_raises`,
`test_cigar_replay_and_np_columns_preserve_retained_traces`](../abstar/tests/test_airr.py).

### Amino-acid alignment columns spanning gaps, NP boundaries, or incomplete codons

**Description.** Codon-sized gaps, non-triplet gaps, NP segments beginning within a
codon, and partial terminal codons cannot be handled by independently ungapping and
translating both alignment rows.

**Example.** These are exact serializer regression values:

```text
query NT:       ATG AAA CCC GGG
germline NT:    ATG A-- --C GGG
query AA:        M   K   P   G
germline AA:     M   X   X   G

query NT:       ATG AAA CCC
germline NT:    ATG --- CCC
query AA:        M   K   P
germline AA:     M   -   P
```

**Handling.** Final AIRR/Parquet aligned amino-acid fields use the same triplets of
alignment columns in both rows after skipping the leading partial query codon. An
all-gap codon is `-`, an all-dot IMGT spacer is `.`, and mixed-gap or ambiguous triplets
are `X`. A terminal incomplete triplet is dropped from both aligned rows. Missing
context or no complete codon gives null paired AA values; invalid frame, unequal row
length, or double gaps are rejected.

**Limits.** These are translations of retained aligned rows, not independent
translations of ungapped sequences. They are also distinct from continuous-query FWR/CDR
amino-acid partitioning and local junction-frame subdivisions described in the
biological sections.

**Evidence.** [`translate_airr_alignment`](../abstar/annotation/airr.py),
[`_write_sample_outputs`](../abstar/core/abstar.py);
[`test_official_paired_amino_acids_follow_shared_codon_columns`,
`test_official_paired_translation_rejects_invalid_evidence`,
`test_public_official_amino_acids_map_to_source_nt_and_preserve_internals`](../abstar/tests/test_airr.py).

### Full-input amino-acid phase with flanking sequence or a truncated first codon

**Description.** The V alignment can start inside a larger query, or the read can begin
midway through a codon. Translating the original read from position zero can select the
wrong frame.

**Example.** Oriented input `CATGAAACCCGGGTAA`, V query start `4`, and V frame `1` gives
full-query `sequence_aa="MKPG*"`. The same result follows with V start `2`, frame `3`. A
read `AATGAAACC` with start `0`, frame `2` yields `MK`, omitting leading/terminal
incomplete codons.

**Handling.** Official file `sequence_aa` translates the entire oriented query using
phase `(v_sequence_start + frame - 1) % 3`, preferring `v_frame` when available. Thus
full input AA can include flanking translated context while `sequence` remains original
input. Missing query/origin/frame yields null rather than guessing a phase, and
unassigned official AA fields are null.

**Evidence.** [`translate_airr_query`, `to_airr_row`](../abstar/annotation/airr.py);
[`test_official_query_translation_uses_input_phase_and_orientation`,
`test_unknown_query_translation_context_emits_null`,
`test_unassigned_official_amino_acids_are_null`](../abstar/tests/test_airr.py);
[`test_final_parquet_official_sequences_match_independent_source_evidence`](../abstar/tests/test_output_parity.py).

## Custom-reference integrity and interrupted database builds

### Mixed gapped/ungapped references and biologically incompatible inputs

**Description.** A custom database can mix new ungapped V alleles, existing IMGT-gapped
V alleles, and D/J genes. Adding gaps indiscriminately would alter already-numbered
genes or introduce gaps into segments that do not use V-style IMGT numbering. Duplicate
labels or mismatched receptor loci can make later lookup ambiguous.

**Example.** Suppose the input contains a gapped V fragment `ACG...TTC`, an ungapped V
allele, and a J gene `TTCGGT`. The first V must keep its existing dots and the J must
remain ungapped. Two entries both called `IGKV1-5*01`, or an `IGKD...` diversity gene,
are invalid inputs.

**Handling.** The builder validates unique whitespace-free identifiers, nucleotide
alphabets, receptor-compatible loci, D-bearing loci, and the presence of V and J for
each supplied VDJ locus. IMGT dots are allowed in V segments; D/J dots are rejected.
Already-gapped V sequences pass through. Ungapped V sequences are semiglobally aligned
to gapped references from the same locus and receptor with an IMGT-specific matrix.
Query gaps become dots; trailing dots are removed, while leading dots preserve
numbering. D/J genes pass through unchanged. Ungapped database FASTAs are derived by
removing dots. Staged validation checks matching unique IDs and exact gapped-to-ungapped
sequence equivalence. Reference-based gapping supplies a numbering procedure, not
independent biological validation of every novel allele.

**Evidence.** [core/germline.py](../abstar/core/germline.py), `validate_germlines`,
`add_imgt_gaps`, `validate_staged_database`;
[test_custom_germline.py](../abstar/tests/test_custom_germline.py),
`test_validate_germlines_rejects_unsafe_inputs` and
`test_add_imgt_gaps_handles_mixed_inputs_with_locus_reference`;
[test_database_integrity.py](../abstar/tests/test_database_integrity.py),
`test_packaged_germline_database_integrity`.

### Partial builds, replacement failures, and competing builders

**Description.** Reference gapping or MMseqs indexing can fail halfway through a build.
Overwriting an existing database in place would leave a mixture of old and new files.
Two processes might also build the same destination concurrently.

**Example.** A replacement has completed V and J indexes when D indexing exits nonzero.
The installed reference should remain the old complete database, rather than an
apparently usable directory with new V genes and stale D genes.

**Handling.** Builds take place in a private staging directory. Source FASTAs, IDs,
sequences, nonempty regular-file components, and indexes are validated before
publication. Publication uses a per-database lock and checks that the destination's
filesystem identity has not changed during the build. Existing databases are moved to a
backup before the staged directory is renamed into place. A failed publication restores
the backup when possible; if rollback also fails, a structured error retains both errors
and the backup recovery path. Backup-cleanup or secondary logging failures do not turn a
completed publication into a reported build failure. Staging cleanup and stale-file
removal are covered separately. This procedure handles the tested build and publication
errors; it is not a claim of power-loss durability or uninterrupted visibility to
readers between the two replacement renames.

**Evidence.** [core/germline.py](../abstar/core/germline.py), `build_germline_database`,
`publish_database`, `database_build_lock`;
[test_custom_germline.py](../abstar/tests/test_custom_germline.py),
`test_mmseqs_failure_at_each_segment_and_command_preserves_destination`,
`test_incomplete_staged_index_is_rejected_before_publication`,
`test_publication_failure_rolls_back_existing_database`,
`test_approved_destination_identity_change_is_not_touched`, and
`test_publish_and_rollback_failure_retains_both_errors_and_recovery_path`.
