.. _umis:

UMI Support
===========

abstar can detect and extract unique molecular identifiers (UMIs) from
sequencing data. UMIs are short random sequences added during library
preparation that enable error correction and duplicate identification.


Basic Usage
-----------

**Command Line:**

.. code-block:: bash

    # Using a built-in pattern
    abstar run sequences.fasta output/ --umi_pattern smartseq-human-bcr

    # Custom pattern with specified length
    abstar run sequences.fasta output/ --umi_pattern "[UMI]TCAGCGGGAAGACATT" --umi_length 12

    # UMI by position only (no pattern matching)
    abstar run sequences.fasta output/ --umi_length 12

**Python:**

.. code-block:: python

    import abstar

    # Using built-in pattern
    sequences = abstar.run(
        "sequences.fasta",
        umi_pattern="smartseq-human-bcr"
    )

    # Custom pattern
    sequences = abstar.run(
        "sequences.fasta",
        umi_pattern="[UMI]TCAGCGGGAAGACATT",
        umi_length=12
    )


Pattern Format
--------------

UMI patterns use ``[UMI]`` as a placeholder to indicate where the UMI
sequence is located, with conserved flanking sequences for anchoring:

**Pattern examples:**

``"[UMI]TCAGCGGGAAGACATT"``
    UMI at 5' end, followed by conserved sequence ``TCAGCGGGAAGACATT``

``"ATGCATGC[UMI]"``
    Conserved sequence ``ATGCATGC`` followed by UMI

``"ATGC[UMI]GCTA"``
    UMI flanked by conserved sequences on both sides

When the pattern has a trailing conserved sequence (like ``[UMI]TCAG...``),
the UMI length can be inferred from the alignment. When the pattern ends
with ``[UMI]`` (no trailing sequence), ``--umi_length`` is required.


UMI Position (5' vs 3')
-----------------------

The sign of ``--umi_length`` determines which end of the sequence to search:

- **Positive length**: UMI is at the 5' end of the sequence
- **Negative length**: UMI is at the 3' end of the sequence

.. code-block:: bash

    # 12bp UMI at 5' end
    abstar run seqs.fasta out/ --umi_length 12

    # 8bp UMI at 3' end (sequence is reverse-complemented before matching)
    abstar run seqs.fasta out/ --umi_length -8

When a negative length is used, the sequence is automatically reverse-complemented
before pattern matching, so you can write patterns in the 5'->3' orientation
of your primers.


Built-in Patterns
-----------------

smartseq-human-bcr
~~~~~~~~~~~~~~~~~~

For `Takara's Smart-Seq Human BCR kit`_ with UMIs. This pattern set includes
primers for heavy chain (IgG, IgM, IgA, IgD, IgE) and light chains (kappa, lambda),
each with a 12bp UMI.

.. code-block:: bash

    abstar run sequences.fasta output/ --umi_pattern smartseq-human-bcr

The patterns automatically handle mixed samples containing heavy, kappa,
and lambda chains by attempting to match each chain-specific primer pattern.

.. _Takara's Smart-Seq Human BCR kit: https://www.takarabio.com/products/next-generation-sequencing/immune-profiling/human-repertoire/smart-seq-human-bcr-with-umis


Multiple UMIs
-------------

If a sequence contains multiple UMIs (e.g., at both ends, or with different
primers), they are concatenated with ``+``:

.. code-block:: text

    # Single UMI
    umi: "ATCGATCGATCG"

    # Multiple UMIs
    umi: "ATCGATCGATCG+GCTAGCTAGCTA"


Mismatch Tolerance
------------------

By default, up to 1 mismatch is allowed when matching the conserved flanking
sequences. The ``smartseq-human-bcr`` pattern allows 2 mismatches.


Output
------

UMIs are stored in the ``umi`` field of the annotation output:

**In Python:**

.. code-block:: python

    sequences = abstar.run("input.fasta", umi_pattern="smartseq-human-bcr")
    for seq in sequences:
        print(seq["umi"])  # e.g., "ATCGATCGATCG"

**In AIRR TSV output:**

The UMI appears in the ``umi`` column.

**In Parquet output:**

The UMI is stored in the ``umi`` field.


Sequences Without UMIs
----------------------

Sequences where the UMI pattern cannot be matched (due to too many mismatches
or the pattern not being found) will have a ``null``/empty UMI field but are
still annotated normally.
