# Adjudicated original BCR contigs

These 36 original human BCR contigs are from **“Functional antibodies exhibit
light chain coherence,” Figshare record `19617633`, version `3`**, available at
<https://figshare.com/articles/preprint/Functional_antibodies_exhibit_light_chain_coherence/19617633/3>.
The source license is **MIT**. Copied fixtures are redistributed under the
repository MIT license. Original dataset and contig IDs are retained exactly.
The small copied subset contains 21,559 nucleotides: 20 IGH, eight IGK and
eight IGL contigs. No external corpus path or dependency on that corpus is
required to load these fixtures.

`sequences.fasta` and `cases.json` have identical record order. A FASTA header
is the original contig ID; the corresponding JSON row supplies its original
dataset. The pair identifies the source record, including when different
datasets reuse an external contig ID. Each case pins the SHA-256 of its exact
ASCII sequence, original Cell Ranger fields, source-file checksums, germline
checksums, alignment evidence and curation rationale. The loader returns fresh,
recursively immutable cases. `case.as_sequence()` creates an independent mutable
`abutils.Sequence` when a test needs one.

## Selection

The first eight cases are all eight original missing-output records from dataset
`1279068`. Both junction disagreements from that complete 1,935-record pilot
cohort are included. The other candidates come from the bounded 200-record-per-
dataset discovery sweep, using seed `abstar-real-bcr-v1`, algorithm version 1.
Rank is SHA-256 of UTF-8 `1\0abstar-real-bcr-v1\0dataset\0sequence_id`.
After the mandatory records, each bucket uses the lowest rank passing review:

| Bucket | Accepted |
| --- | ---: |
| Original pilot losses | 8 |
| Original pilot junction disagreements | 2 |
| Concordant V/J, junction and productivity controls | 2 per locus |
| Productivity disagreements | 2 per locus |
| Tied V/J calls (including within-gene allele ties) | 2 per locus |
| V insertions | 2 |
| V deletions | 2 |
| IGH with no defensible D assignment | 2 |
| Shortest and longest concordant junction in the sweep | 1 each |

There are no duplicated records across buckets. Full ranks and one-based bucket
positions are in `source.selection`. The first deletion candidate,
`1287189/CATGGCGCAGCGTTCG-1_contig_2`, was rejected: permissive nucleotide
alignment split the loss into one- and two-base gaps, whereas stronger penalties
placed a contiguous three-base loss beside substitutions. The retained deletion
representatives are ranks two and three. Repeat-related shifts of an intact
indel remain inspectable in the evidence; their coordinates use the leftmost
optimal placement under the documented stronger-penalty alignment.

The reviewed historical discovery report was schema version 1, SHA-256
`ec01dfd40aad781019d5eff61d64e076804a6d24f63f48f3be05acb4d3a8a752`.
The complete pilot diagnostic and the ordered accepted/rejected decision report
were working evidence outside the repository. Their full source data are not
needed for fixture integrity tests.

## Biological expectations and boundaries

These are curated expectations, not regenerated abstar output. Current abstar
and Cell Ranger were comparison evidence. Each original sequence and CSV row
was re-read; every accepted V/J call was checked against the packaged human BCR
ungapped germlines, and the V anchor was mapped through the corresponding
IMGT-gapped sequence. Direct exhaustive locus-compatible nucleotide alignment
used match 2, mismatch -3, gap-open 12 and gap-extension 2; comparison at gap-open
5, extension 1 exposed unstable gap structures. Full traces and top candidate
scores remain in each case. Both query strands were compared. All originals
support the forward orientation.

Gene calls are sorted **gene-level allowed sets**, even for a singleton. None
claims an exact allele. Equal V-gene alternatives are retained, including
`IGKV1-13`/`IGKV1D-13`; reviewed near alternatives differing by one substitution
are retained where appropriate. Scores are evidence of sequence compatibility,
not proof of donor genotype. C-gene calls and unreviewed IGH D-gene calls are
omitted from expectations. Light chains have null D. In the two no-D IGH cases,
only four bases remain between retained V and J alignments: null D means no
reliable assignment, not proof of a biologically D-free rearrangement.

Expected segment and junction coordinates are **zero-based, half-open in the
oriented full original input**. Segment boundaries describe the retained local
matches, allowing genuine terminal mismatch clipping; they are not claims about
exact recombination breakpoints. If V and J local matches overlap, V owns the
shared bases and J starts immediately afterward. The raw alignment coordinates
remain in evidence. A consumer writing AIRR must convert starts by adding one;
ends retain their numeric value. See the [AIRR coordinate and junction
specification](https://docs.airr-community.org/en/stable/datarep/rearrangements.html).

The V junction anchor is the homologous IMGT position 104, located by counting
ungapped bases before nucleotide offset 309 in the gapped V reference. The J
anchor is the W/F of its homologous W/F-G-X-G motif. The [IMGT numbering
reference](https://www.imgt.org/IMGTScientificChart/Numbering/IMGTIGVLsuperfamily.html)
identifies those conserved positions. A mutation does not move the boundary to
a more convenient downstream residue. Junction includes the anchor codons;
CDR3 excludes them. Where a mutated J anchor lies just outside the local match,
its codon is projected through the contiguous immediately adjacent J sequence,
with the projection stated in evidence.

Productivity here tests an unambiguous, stop-free V-through-J ORF, an in-frame
junction and the package's locus-appropriate conserved-anchor rule. It does not
establish antibody expression, antigen binding, regulatory integrity or function.
Thirty-three cases satisfy that rule; three have an altered C/W/F anchor. The
six deliberately selected productivity disagreements have an intact ORF and
correct anchors. Their coding-origin arithmetic is explicit, so an erroneous
comparison of a full-query junction coordinate to a trimmed-query frame cannot
become an expected failure.

Special pilot adjudications are documented per case:

- `ATCATCTTCAGCAACT-1_contig_1`: the homologous V cysteine is GGT (glycine).
  Its junction begins `GARDE...`; a downstream C does not replace the V anchor.
- `GTTACAGCACATAACC-1_contig_2`: four optimal insertion placements all map
  the conserved C to query 365, supporting `CSSYCNSYTSSSTLYVF`.
- `CTAAGACAGCAATCTC-1_contig_2`: the J F anchor is ATC (isoleucine).
- `CTGTTTACAGGTGCCT-1_contig_1`: the J W anchor is GTT (valine).
- `CCATGTCCAGTCTTCC-1_contig_1`: two in-frame J-like tracts occur in the
  original sequence. The primary junction ends at the first complete J motif
  (W at 455); the later repeat is kept in evidence. Both IGHJ4 and IGHJ5 are
  compatible with the primary tract. The duplication's molecular origin remains
  unresolved; the fixture does not claim to resolve it.

Indel expectations are structured lists with `query_start`, `query_end` and
`sequence`. Insertions span observed query bases; deletions have zero query
width and name the removed representative-germline bases. These are explicit
fixture coordinates, not abstar's IMGT-numbered indel strings. Equivalent
repeat rotations require normalization before a consumer compares them.

## Integrity check

Run `python -m pytest abstar/tests/test_real_bcr.py -m "not e2e" -q`.
Integrity tests check record conservation, ordering, immutable state, provenance,
hashes, expectation types and mandatory membership. They must never launch
MMseqs. Loading the fixtures does not establish that current annotation passes
these biological expectations; the original missing outputs remain diagnostic
regression targets.
