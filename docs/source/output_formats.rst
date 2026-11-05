.. _output-formats:

Output Formats
==============

abstar targets the `AIRR Data Standards 2.0 Rearrangement schema`_. The test
suite validates TSV output with the official AIRR Python library 2.0.0 and
checks biological coordinate and alignment semantics. Two file formats
are supported:

- ``airr``: Tab-delimited TSV file with header row
- ``parquet``: Columnar binary format, more space-efficient for large datasets

.. _AIRR Data Standards 2.0 Rearrangement schema: https://docs.airr-community.org/en/stable/datarep/rearrangements.html

.. warning::

   **Breaking change in 0.9.0.** Every public output (AIRR TSV, Parquet, and
   the dataframe or ``Sequence`` objects returned by ``abstar.run()``) now has
   the same fields, in the same order, with the same meanings. Compared with
   0.8.x:

   - **Coordinates are 1-based closed intervals everywhere.** Every
     ``*_start`` value is one greater than in 0.8.x; every ``*_end`` value is
     unchanged.
   - **``sequence`` is the original input**, and ``sequence_aa`` translates the
     full oriented query (see `Official Amino-acid Fields`_). In 0.8.x both
     held the assembled V(D)J. The V(D)J is now in the new
     ``sequence_vdj`` and ``sequence_vdj_aa`` fields, which the masks index.
     A returned ``Sequence`` object's ``.sequence`` is therefore the input.
   - ``sequence_alignment_aa`` and ``germline_alignment_aa`` translate the
     paired nucleotide alignment columns instead of aligning the two
     translations.
   - TSV booleans are ``T``/``F``, and the AIRR-required fields lead the field
     order in every output.

   For example, record ``AAACCTGAGACCACGA-1_contig_2`` (535 nt input, 343 nt
   V(D)J):

   .. list-table::
      :header-rows: 1

      * - Field
        - 0.8.x
        - 0.9.0
      * - ``sequence``
        - 343 nt V(D)J (``CAGGTGCAGCTG...``)
        - 535 nt input (``GGGAGCATCACC...``)
      * - ``sequence_vdj``
        - (absent)
        - 343 nt V(D)J (``CAGGTGCAGCTG...``)
      * - ``v_sequence_start`` / ``v_sequence_end``
        - ``121`` / ``417``
        - ``122`` / ``417``
      * - ``cdr3_start`` / ``cdr3_end``
        - ``409`` / ``430``
        - ``410`` / ``430``

   To compare 0.8.x output with 0.9.0, add one to every 0.8.x ``*_start`` and
   use ``sequence_vdj`` wherever 0.8.x code used ``sequence``. To slice a 0.9.0
   coordinate pair in Python, use ``seq[start - 1:end]``. See
   :ref:`changelog`.


Coordinates
-----------

Every public output uses **1-based closed intervals**, as AIRR requires.
Annotation internals (``Antibody`` objects) use 0-based half-open intervals and
are converted once, when each record's public row is built. For example, an
internal interval ``[137, 439)`` is published as start ``138``, end ``439``;
``[0, 1)`` is published as ``1, 1``. Missing intervals are null (empty TSV
cells). This applies to V/D/J/C query and germline coordinates and to
``fwr1/cdr1/fwr2/cdr2/fwr3/cdr3/fwr4_start/end``.

Query coordinates address ``sequence_oriented``. Germline coordinates address
the ungapped reference of the corresponding gene. In every output,
``sequence`` is always the original, unmodified ``sequence_input``. When
``rev_comp`` is true, all
alignments and query coordinates refer to its reverse complement. Opaque
identifiers, including duplicate IDs and leading zeroes, retain their input
order. The private internal ``row_id`` is never serialized.

The first columns of every output, in order, are ``sequence_id``, ``sequence``, ``rev_comp``,
``productive``, ``v_call``, ``d_call``, ``j_call``, ``sequence_alignment``,
``germline_alignment``, ``junction``, ``junction_aa``, ``v_cigar``, ``d_cigar``,
and ``j_cigar``. All remaining public fields follow in the declared output
schema order, without duplicates. Boolean values are ``T`` or ``F``, nulls
are empty cells, and lines end with LF. Values containing a tab, newline or
carriage return raise ``ValueError`` before the destination is opened. Literal
quote characters use standard CSV quoting with doubled quotes, so IDs such as
``"quoted"`` and ``"unterminated`` round-trip through the official AIRR reader
without losing characters or merging records.

``sequence_alignment`` and ``germline_alignment`` retain the same V/D/J
alignment columns used as annotation evidence. Non-templated NP query bases
align to germline gaps. Segment CIGARs use ``M`` for aligned residue pairs
(including substitutions), ``I`` for query insertions, and ``D`` for query
deletions. Leading ``S`` operations skip query bases and leading ``N``
operations skip reference bases. For example, ``AC-GT`` aligned to ``ACCGT``
at query offset 2 and reference offset 1 produces ``2S1N2M1D2M``.
The C CIGAR describes the separately retained constant-region alignment.
``junction`` includes the conserved endpoint codons; ``cdr3`` excludes them.


Official Amino-acid Fields
~~~~~~~~~~~~~~~~~~~~~~~~~~

In every output, ``sequence_aa`` translates the full
``sequence_oriented`` query.
The retained V frame is one-based within ``v_sequence``. Its first complete
codon starts at the 0-based oriented-query offset
``(v_sequence_start - 1) + (v_frame - 1)``; the full-query translation therefore
starts at that offset modulo 3. The ``frame`` alias is used when
``v_frame`` is absent. Leading and terminal partial codons are omitted.
A reverse-complement input and its forward counterpart have the same
``sequence_aa``. Translation includes input flanks; productivity is still
assessed over the retained V(D)J coding sequence.

``sequence_alignment_aa`` and ``germline_alignment_aa`` translate the paired
nucleotide alignment over one shared coding window. After skipping the leading
partial query codon, both rows use the same triplets of alignment columns.
Complete gap codons become ``-``, complete IMGT-dot spacers become ``.``, and
mixed-gap or unresolved codons become ``X``. Thus ``ATGAAACCC`` aligned to
``ATG---CCC`` becomes ``MKP`` aligned to ``M-P``; NP bases are never copied into
the inferred germline amino acids. An indel within a codon can produce ``X``
without an invented amino-acid correction. These fields translate the aligned
columns; independently translating the ungapped rows may give different
results for frameshifts. A terminal incomplete column triplet is omitted from
both rows together.

``sequence_aa`` is null if the oriented query, retained V origin or coding
frame is unavailable, or no complete codon remains. The paired AA fields are
both null when either alignment row or its frame is unavailable, or the shared
window has no complete triplet. Invalid known frames or inconsistent alignment
rows raise an error rather than producing inferred residues.

Unassigned records retain ``sequence_input``, ``sequence_oriented``,
``annotation_status=unassigned``, and an inspectable ``failure_reason``.
Their ``productive`` value and unavailable annotation fields remain null.
Sequence annotation exceptions have no normal output row. They are recorded in
``logs/failures.tsv`` and individual diagnostic files while other records continue.
An all-failed sample has an empty output with explicit failure counts and logs.
``strict=True`` / ``--strict`` instead raises a structured run error. Input, worker,
external-tool, and output/storage errors always remain fatal.

V(D)J and assembly fields
~~~~~~~~~~~~~~~~~~~~~~~~~

``sequence_vdj`` is the assembled V(D)J nucleotide sequence (retained V, NP1,
optional D and NP2, and J; no constant region) and ``sequence_vdj_aa`` is its
translation in the V coding frame. Both are null for unassigned records. They
are abstar extensions with the same names for every locus, including light
chains and TCRs without a D gene. The ``germline``, ``*_gapped``, ``*_vdjc``
and ``*_lvdjc`` assembly fields are also abstar extensions, including their
legacy NP content; they must not be used as substitutes for the paired AIRR
alignment fields. ``cdr3_length`` is also an abstar extension and counts
amino acids.

When comparing TSV with Parquet or Python returns, decode raw TSV numeric cells
and ``T``/``F`` and reconcile empty strings with nulls; every field then compares
directly. The official ``airr.read_rearrangement(..., validate=True)`` reader
converts ``*_start`` values to 0-based offsets, so add one back before comparing
its rows with Parquet. ``abstar.annotation.airr.to_public_row()`` is the single
pure mapper from an internal annotation row to the public row.


Specifying Output Format
------------------------

**Command Line:**

.. code-block:: bash

    # Default: AIRR TSV
    abstar run sequences.fasta output/

    # Parquet format
    abstar run sequences.fasta output/ -o parquet

    # Both formats
    abstar run sequences.fasta output/ -o airr -o parquet

**Python:**

.. code-block:: python

    import abstar

    # Write to files
    abstar.run("sequences.fasta", "output/", output_format=["airr", "parquet"])

    # Return as DataFrame (no file output)
    df = abstar.run("sequences.fasta", as_dataframe=True)


Output Directory Structure
--------------------------

.. code-block:: text

    output/
    ├── airr/                 # AIRR TSV files
    │   └── sequences.tsv
    ├── parquet/              # Parquet files
    │   └── sequences.parquet
    └── logs/                 # Log files
        └── abstar.log


Output Fields
-------------

Core Identification
~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Field
     - Type
     - Description
   * - ``sequence_id``
     - String
     - Unique sequence identifier
   * - ``sequence_input``
     - String
     - Original input sequence
   * - ``sequence_oriented``
     - String
     - Sequence in V->J orientation
   * - ``rev_comp``
     - Boolean
     - True if sequence was reverse-complemented
   * - ``quality``
     - String
     - Quality scores (if FASTQ input)
   * - ``umi``
     - String
     - Unique molecular identifier (if parsed)
   * - ``locus``
     - String
     - Locus (e.g., IGH, IGK, IGL, TRA, TRB)
   * - ``species``
     - String
     - Species from germline database
   * - ``germline_database``
     - String
     - Name of germline database used


Gene Calls
~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Field
     - Type
     - Description
   * - ``v_call``
     - String
     - V gene assignment with allele (e.g., IGHV1-2*02)
   * - ``d_call``
     - String
     - D gene assignment
   * - ``j_call``
     - String
     - J gene assignment
   * - ``c_call``
     - String
     - C gene (isotype) assignment
   * - ``v_gene``
     - String
     - V gene without allele (e.g., IGHV1-2)
   * - ``d_gene``
     - String
     - D gene without allele
   * - ``j_gene``
     - String
     - J gene without allele
   * - ``c_gene``
     - String
     - C gene without allele
   * - ``v_support``
     - Float
     - V gene assignment E-value
   * - ``d_support``
     - Float
     - D gene assignment E-value
   * - ``j_support``
     - Float
     - J gene assignment E-value
   * - ``c_support``
     - Float
     - C gene assignment E-value


Regions
~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Field
     - Type
     - Description
   * - ``fwr1``, ``fwr1_aa``
     - String
     - Framework region 1 (nucleotide, amino acid)
   * - ``cdr1``, ``cdr1_aa``
     - String
     - CDR1
   * - ``fwr2``, ``fwr2_aa``
     - String
     - Framework region 2
   * - ``cdr2``, ``cdr2_aa``
     - String
     - CDR2
   * - ``fwr3``, ``fwr3_aa``
     - String
     - Framework region 3
   * - ``cdr3``, ``cdr3_aa``
     - String
     - CDR3
   * - ``fwr4``, ``fwr4_aa``
     - String
     - Framework region 4
   * - ``junction``, ``junction_aa``
     - String
     - Junction region (conserved C to conserved W/F)
   * - ``cdr3_length``
     - Integer
     - CDR3 length in amino acids
   * - ``np1``, ``np2``
     - String
     - N-nucleotide regions (non-templated)
   * - ``np1_length``, ``np2_length``
     - Integer
     - Length of N-regions


Junction Components
~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Field
     - Type
     - Description
   * - ``cdr3_v``, ``cdr3_v_aa``
     - String
     - V gene contribution to CDR3
   * - ``cdr3_n1``, ``cdr3_n1_aa``
     - String
     - N1 region (V-D junction)
   * - ``cdr3_d``, ``cdr3_d_aa``
     - String
     - D gene contribution to CDR3
   * - ``cdr3_n2``, ``cdr3_n2_aa``
     - String
     - N2 region (D-J junction)
   * - ``cdr3_j``, ``cdr3_j_aa``
     - String
     - J gene contribution to CDR3


Quality Metrics
~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Field
     - Type
     - Description
   * - ``productive``
     - Boolean
     - True if sequence is productive
   * - ``vj_in_frame``
     - Boolean
     - True if the V-to-J junction is in frame and has a codon-complete length
   * - ``productivity_issues``
     - String
     - List of productivity issues (if any)
   * - ``stop_codon``
     - Boolean
     - True if stop codon present
   * - ``complete_vdj``
     - Boolean
     - True if V, D (heavy only), and J assigned
   * - ``v_identity``
     - Float
     - V gene alignment identity (0-1), including indel columns
   * - ``v_identity_aa``
     - Float
     - V gene amino acid alignment identity, including indel columns
   * - ``d_identity``
     - Float
     - D gene alignment identity, including indel columns
   * - ``d_identity_aa``
     - Float
     - D gene amino acid alignment identity, including indel columns
   * - ``j_identity``
     - Float
     - J gene alignment identity, including indel columns
   * - ``j_identity_aa``
     - Float
     - J gene amino acid alignment identity, including indel columns
   * - ``c_identity``
     - Float
     - Constant-region alignment identity, including indel columns
   * - ``c_identity_aa``
     - Float
     - Constant-region amino acid alignment identity, including indel columns
   * - ``frame``
     - Integer
     - One-based reading frame within the retained V sequence (1, 2, or 3)


D and J amino acid identities are null when either retained alignment row has no
complete codon in its germline reading frame. A short nonempty nucleotide
alignment retains its gene call, nucleotide identity, and alignment evidence.
A D alignment with an empty nucleotide query or germline is discarded; the full
interval between retained V and J sequences is represented as NP1.


Mutations
~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Field
     - Type
     - Description
   * - ``v_mutations``
     - String
     - V gene mutations (format: "pos:ref>alt")
   * - ``v_mutations_aa``
     - String
     - V gene amino acid mutations
   * - ``v_mutation_count``
     - Integer
     - Number of V gene mutations
   * - ``v_mutation_count_aa``
     - Integer
     - Number of V gene AA mutations
   * - ``v_insertions``
     - String
     - Non-templated insertions in V
   * - ``v_deletions``
     - String
     - Non-templated deletions in V
   * - ``v_frameshift``
     - Boolean
     - True if frameshift in V region
   * - ``c_mutations``
     - String
     - C gene mutations
   * - ``c_mutation_count``
     - Integer
     - Number of C gene mutations


Sequences
~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Field
     - Type
     - Description
   * - ``sequence``
     - String
     - Original input query
   * - ``sequence_vdj``
     - String
     - Assembled V(D)J sequence, without the constant region (abstar extension)
   * - ``germline``
     - String
     - Corresponding germline sequence
   * - ``sequence_aa``
     - String
     - Full oriented query translated in the V-derived coding phase
   * - ``sequence_vdj_aa``
     - String
     - Translation of ``sequence_vdj`` (abstar extension)
   * - ``germline_aa``
     - String
     - Germline amino acid sequence
   * - ``sequence_gapped``
     - String
     - IMGT-gapped sequence
   * - ``germline_gapped``
     - String
     - IMGT-gapped germline
   * - ``sequence_alignment``
     - String
     - Aligned sequence (with gaps from alignment)
   * - ``germline_alignment``
     - String
     - Aligned germline


Masks
~~~~~

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Field
     - Type
     - Description
   * - ``cdr_mask``
     - String
     - CDR region mask (0=FWR, 1=CDR1, 2=CDR2, 3=CDR3)
   * - ``gene_segment_mask``
     - String
     - Assembled V(D)J segment mask (V, N, D, J); N denotes NP1/NP2 bases
   * - ``gene_segment_mask_aa``
     - String
     - Codon segment mask; codons spanning different segments are N
   * - ``nongermline_mask``
     - String
     - Mutation position mask


Gene-segment masks contain one label per nucleotide or amino acid in
``sequence_vdj`` or ``sequence_vdj_aa``. The ``sequence`` and ``sequence_aa``
fields include input flanks, so these masks are not indexed over them.
Nucleotide labels follow the
retained V, NP1, optional D/NP2, and J spans. Amino acid labels use complete
codons in the assembled query reading frame; a codon receives V, D, or J only
when all three nucleotides belong to that segment. These masks exclude the
constant region and do not depend on CDR3 or framework subdivisions.

Independent V(D)J consistency checks
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``sequence_alignment.replace("-", "")`` provides the nucleotide V(D)J
query without assembling FWR/CDR sequences. For a reference independent of
both regional and gene-segment assembly, slice the oriented input directly.
For any public row (Parquet, Python return, or decoded TSV):

.. code-block:: python

   from Bio.Seq import Seq

   vdj_nt = row["sequence_oriented"][
       row["v_sequence_start"] - 1:row["j_sequence_end"]
   ]
   coding = vdj_nt[row["frame"] - 1:]
   vdj_aa = str(Seq(coding[:len(coding) // 3 * 3]).translate())
   regions = ("fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4")
   assert row["sequence_vdj"] == vdj_nt
   assert row["sequence_alignment"].replace("-", "") == vdj_nt
   assert "".join(row[r] for r in regions) == vdj_nt
   assert "".join(row[r + "_aa"] for r in regions) == vdj_aa

These are diagnostic assertions for annotated records, not a guarantee that
every existing annotation passes. In particular, independent regional
translations can disagree across frameshifted boundaries. Compare sequence
contents as well as all six mask lengths: equal lengths can conceal duplicated
or mistranslated residues. The direct slice validates consistency with the
annotated V/J coordinates; it does not establish that those coordinates are
biologically correct. Rows from the official AIRR reader have 0-based starts,
so do not apply this slice to them directly.

To audit a directory of 0.9+ Parquet files without rerunning annotation:

.. code-block:: bash

   python scripts/audit_annotation_consistency.py /path/to/parquet \
     --report /tmp/abstar-consistency.json

The report must be new and outside the input and repository trees. It contains
counts and candidate records identified by source Parquet path, zero-based
Parquet row ordinal, and exact sequence ID. Unassigned records are counted
separately. Exit status 1 indicates discovered inconsistencies; status 0 means
all implemented checks passed. The audit does not verify biological truth or
the contents of gene-segment and non-germline masks.


Position Coordinates
~~~~~~~~~~~~~~~~~~~~

V/J boundaries use nucleotide Smith-Waterman alignment with match 2,
mismatch -3, gap-open penalty 12, and gap-extension penalty 2. Boundary
alignment uses the complete supplied query, independently of the semiglobal
alignment used for IMGT and junction mapping. Equal optimal scores select the
earlier query endpoint. A J-like repeat wholly downstream of the mapped
primary junction is excluded from the primary J boundary search. Coordinates
refer to the selected reference; an accepted alternative reference may have
a different supported boundary.

These boundary corrections can change segment identities, mutation counts,
N-region subdivisions, and fallback D calls when the retained V–J interval
changes. ``v_score`` and ``j_score`` describe the boundary alignments;
``v_support`` and ``j_support`` retain search evidence. Productivity evaluates
the junction start relative to the V-region origin and its one-based
``frame``. Internal annotation coordinates are zero-based, half-open; public
outputs are 1-based closed.

If the junction-start alignment endpoint falls in a query gap, annotation
attempts recovery by aligning raw oriented-query FWR3 bases through the complete
IMGT 104 reference codon. The upstream boundary comes from the retained V
alignment. Recovery requires a unique, contiguous three-base anchor projection
across all optimal alignments and agreement with retained anchor evidence.
Competing boundaries, incomplete codons, and conflicting evidence remain
explicit per-sequence failures. The diagnostic log records the recovered
interval and alignment score.

Recovery does not search for a nearby cysteine or prefer a productive reading
frame. Observed mutations, ambiguous bases, and frameshifts remain subject to
the usual productivity checks. Sequences whose original endpoint maps to a
query base continue through the existing junction calculation.

``c_sequence_gapped`` and ``c_germline_gapped`` now contain the retained C
query/reference pair in the same alignment columns. Previously, the latter
could contain the entire constant reference for a partial read. Consumers
needing the full C reference should retrieve it from the germline database.
Their amino acid counterparts use the C-region frame; C identity includes
mismatch and insertion/deletion columns. Insertion columns consume query bases
without advancing the reference template, including at the final reference
base or when the retained reference starts inside a codon.

All retained V nucleotide alignments, mutation and indel events, identities,
and mapped region sequences use the same boundary-alignment trace. Re-aligning
that span with a different gap penalty could previously erase compensating
indels or make annotation fail. Region boundary adjustment applies only to
complete-codon deletions spanning a boundary, preserving each query base once.

If retained V evidence ends before the FWR3 endpoint can be mapped, but the
FWR3 start is mapped, FWR3 is recovered directly from the oriented query through
the existing junction-start codon. This restores a previously empty FWR3 and
its CDR-mask labels without extending V assignment evidence or moving the
junction. Its amino-acid sequence follows the shared query-frame projection
used for all FWR/CDR regions. The diagnostic log records the recovered interval.

If retained V evidence stops before FWR3 begins, or a partial retained anchor
contradicts the proposed junction boundary, a constrained recovery alignment
starts at the last mapped upstream region boundary and fits through IMGT codon
104 before the assigned J region. This is separate anatomical mapping evidence:
all retained reference-to-query base coordinates (including deletions) are
fixed, and gene calls, assignment scores, identities, mutations, and indels
remain based on the retained alignment. Recovery can fill incomplete FWR2,
CDR2, and FWR3 intervals and correct an unsupported junction anchor.

Every recovered region boundary and all three anchor bases must agree across
optimal alignments to the assigned primary V reference. A unique junction anchor alone is insufficient if the
intermediate region boundaries remain ambiguous. Deleted boundaries, interrupted
anchors, and unsupported mappings fail explicitly; no conserved-motif or
productivity preference is used to resolve ties. Diagnostics include examples
of competing coordinate projections. Default mode records these failures and
continues; ``strict=True`` / ``--strict`` aborts as usual. Recovered sequences
then follow ordinary junction construction, AA projection, and productivity
assessment; a corrected anchor can consequently change those derived fields.

FWR3 ends immediately after the established junction-start codon, sharing its
boundary with CDR3. A retained V alignment can place an insertion immediately
after that codon in FWR3 or map a repeated/deleted anchor differently from
junction finding. The annotated region uses the junction boundary in both
cases, preventing duplicated or omitted query bases. Retained V alignments,
mutation/indel evidence, gene calls, and junction selection remain unchanged;
this synchronization is not a new biological adjudication of ambiguous anchors.
Invalid shared intervals produce explicit per-record errors rather than a
silently truncated slice. Diagnostics record each endpoint adjustment.

FWR4 starts at the existing junction-ending codon and ends at the final retained
J endpoint. If J-boundary selection shortens the initial alignment (for example,
to reject a downstream repeat), FWR4 and its CDR-mask labels stop at that same
endpoint. This synchronization does not change the J assignment or junction.

All FWR/CDR amino-acid fields are projected from their established nucleotide
intervals into the continuous V(D)J translation, using ``v_sequence_start`` and
``frame`` as the coding origin. A codon crossing a region boundary belongs to
the downstream region. Leading bases before the coding origin and terminal
incomplete codons are excluded. Consequently, contiguous nucleotide regions
partition the translated sequence exactly once, even in nonproductive reads;
``cdr_mask_aa`` uses that same partition. The amino-acid projection itself does not adjudicate nucleotide region
boundaries; FWR3/CDR3 synchronization happens before that projection.

Migration note: independently aligned V-region AA boundaries and locally
translated CDR3/FWR4 fields could previously duplicate, omit, or mistranslate
residues. Their corrected values may change region labels even when total AA
assembly was previously correct. ``cdr3_length`` follows the corrected
``cdr3_aa``. ``junction_aa`` and the ``cdr3_v_aa``, ``cdr3_n1_aa``,
``cdr3_d_aa``, ``cdr3_n2_aa``, and ``cdr3_j_aa`` subdivisions retain their local
junction-frame semantics for biological interpretation; on out-of-frame reads,
they need not assemble into the continuous-query ``cdr3_aa``. Gene calls,
retained alignment evidence, nucleotide coordinates, and productivity checks
are unchanged. These region semantics apply to Python returns, Parquet, and
AIRR TSV.

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Field
     - Type
     - Description
   * - ``v_sequence_start``
     - Integer
     - V region start in input sequence
   * - ``v_sequence_end``
     - Integer
     - V region end in input sequence
   * - ``v_germline_start``
     - Integer
     - Start position in V germline
   * - ``v_germline_end``
     - Integer
     - End position in V germline
   * - ``j_sequence_start``
     - Integer
     - J region start in input sequence
   * - ``j_sequence_end``
     - Integer
     - J region end in input sequence
