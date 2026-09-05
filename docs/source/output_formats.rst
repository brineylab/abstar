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


AIRR TSV and Python Coordinates
-------------------------------

AIRR TSV uses **1-based closed intervals**. Python annotations, dataframe
returns, and final or temporary Parquet retain **0-based half-open intervals**.
For example, an
internal interval ``[137, 439)`` is written as start ``138``, end ``439``;
``[0, 1)`` is written as ``1, 1``. Missing intervals have empty start/end cells.
This applies to V/D/J/C query and germline coordinates and to
``fwr1/cdr1/fwr2/cdr2/fwr3/cdr3/fwr4_start/end``.

Query coordinates address ``sequence_oriented``. Germline coordinates address
the ungapped reference of the corresponding gene. In both final file formats,
``sequence`` is always the original, unmodified ``sequence_input``. When
``rev_comp`` is true, all
alignments and query coordinates refer to its reverse complement. Opaque
identifiers, including duplicate IDs and leading zeroes, retain their input
order. The private internal ``row_id`` is never serialized.

The first TSV columns, in order, are ``sequence_id``, ``sequence``, ``rev_comp``,
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

In AIRR TSV and final Parquet files, ``sequence_aa`` translates the full
``sequence_oriented`` query.
The retained V frame is one-based within ``v_sequence``. Its first complete
codon starts at the internal oriented-query offset
``v_sequence_start + v_frame - 1``; the full-query translation therefore starts
at ``(v_sequence_start + v_frame - 1) % 3``. The ``frame`` alias is used when
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
Internal annotation errors raise a structured run error instead of producing
successful empty output.

Compatibility note
~~~~~~~~~~~~~~~~~~

Earlier TSV output exposed Python coordinate offsets, Python boolean text,
and the assembled V(D)J sequence as ``sequence``. Consumers must use the AIRR
conventions above for TSV. Final Parquet files now also expose the original
input as ``sequence`` and derive the three official AA fields from their
corresponding nucleotide query/alignment fields. Consumers of older Parquet
files must account for these four changed field meanings. Parquet booleans
and nulls remain native, and its coordinates remain zero-based half-open.

Python annotation objects, dataframe returns, and temporary work Parquets
retain the assembled V(D)J ``sequence``, its ``sequence_aa``, and their existing
AA alignments for compatibility with masks and productivity. The ``germline``,
``*_gapped`` and ``*_vdjc`` assembly fields in both file formats
remain abstar extensions, including their legacy NP content. They must not be
used as substitutes for the paired AIRR alignment fields. ``cdr3_length`` is
also an abstar extension and counts amino acids.

The assembly extensions such as
``sequence_vdjc_aa`` and ``germline_vdjc_aa`` retain their existing meanings;
no AA fields are silently changed inside annotation objects.

When comparing the two final formats, decode raw TSV numeric cells and
``T``/``F``, reconcile empty strings with nulls, and subtract one from raw TSV
starts. All shared fields, including the four official sequence fields, then
compare directly without recomputing translations or changing gene calls,
evidence, identifiers, order, or outcome reasons. The official
``airr.read_rearrangement(..., validate=True)`` reader already normalizes known
coordinate starts to Python offsets; do not subtract again. The pure
``abstar.annotation.airr.to_airr_row()`` mapper converts an internal annotation
row to AIRR sequence, coordinate and AA conventions.


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
     - AIRR TSV/final Parquet: original input query. Python returns: assembled V(D)J sequence.
   * - ``germline``
     - String
     - Corresponding germline sequence
   * - ``sequence_aa``
     - String
     - AIRR TSV/final Parquet: full oriented query in the V-derived coding phase. Python returns: assembled V(D)J translation.
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


Gene-segment masks contain one label per nucleotide or amino acid in the
ungapped assembled V(D)J sequence (the Python-returned ``sequence`` or
``sequence_aa``). Both final file formats' query fields include input flanks,
so these masks are not indexed over their full ``sequence`` or ``sequence_aa``.
Nucleotide labels follow the
retained V, NP1, optional D/NP2, and J spans. Amino acid labels use complete
codons in the assembled query reading frame; a codon receives V, D, or J only
when all three nucleotides belong to that segment. These masks exclude the
constant region and do not depend on CDR3 or framework subdivisions.


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
``frame``. Current internal coordinates remain zero-based, half-open.

``c_sequence_gapped`` and ``c_germline_gapped`` now contain the retained C
query/reference pair in the same alignment columns. Previously, the latter
could contain the entire constant reference for a partial read. Consumers
needing the full C reference should retrieve it from the germline database.
Their amino acid counterparts use the C-region frame; C identity includes
mismatch and insertion/deletion columns. Insertion columns consume query bases
without advancing the reference template, including at the final reference
base or when the retained reference starts inside a codon.

All retained V nucleotide alignments, mutation and indel events, identities,
and region sequences now use the same boundary-alignment trace. Re-aligning
that span with a different gap penalty could previously erase compensating
indels or make annotation fail. Region boundary adjustment applies only to
complete-codon deletions spanning a boundary, preserving each query base once.

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
