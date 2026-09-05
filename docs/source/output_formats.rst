.. _output-formats:

Output Formats
==============

abstar outputs annotations in `AIRR-compatible format`_. Two file formats
are supported:

- ``airr``: Tab-delimited TSV file with header row
- ``parquet``: Columnar binary format, more space-efficient for large datasets

.. _AIRR-compatible format: https://docs.airr-community.org/en/stable/datarep/rearrangements.html


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
     - Reading frame (0, 1, or 2)


V, D, and J amino acid identities are null when either retained segment has no
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
     - V(D)J sequence (no leader, no constant)
   * - ``germline``
     - String
     - Corresponding germline sequence
   * - ``sequence_aa``
     - String
     - Amino acid sequence
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
ungapped assembled ``sequence`` or ``sequence_aa``. Nucleotide labels follow the
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
