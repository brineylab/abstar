.. _tcr:

TCR Annotation
==============

abstar fully supports T-cell receptor (TCR) sequence annotation using
the ``--receptor tcr`` option.


Basic Usage
-----------

**Command Line:**

.. code-block:: bash

    abstar run tcr_sequences.fasta output/ --receptor tcr

    # Human TCR (default germline database)
    abstar run tcr.fasta output/ --receptor tcr --germline_database human

    # Mouse TCR
    abstar run tcr.fasta output/ --receptor tcr --germline_database mouse

**Python:**

.. code-block:: python

    import abstar

    # Annotate TCR sequences
    sequences = abstar.run("tcr_sequences.fasta", receptor="tcr")

    # Return as DataFrame
    df = abstar.run("tcr_sequences.fasta", receptor="tcr", as_dataframe=True)


TCR Chain Types
---------------

abstar supports all TCR chain types:

- **Alpha chain (TRA)**: V and J gene segments
- **Beta chain (TRB)**: V, D, and J gene segments
- **Gamma chain (TRG)**: V and J gene segments
- **Delta chain (TRD)**: V, D, and J gene segments

D gene assignment applies only to beta and delta chains.


Available Germline Databases
----------------------------

TCR germline databases are available for:

- ``human``: Human TCR (alpha, beta, gamma, delta)
- ``mouse``: Mouse TCR (alpha, beta, gamma, delta)


Output Fields
-------------

TCR annotations include the same fields as BCR annotations:

- Gene calls: ``v_call``, ``d_call`` (beta/delta only), ``j_call``
- Regions: CDR1, CDR2, CDR3, FWR1-4
- Mutations and quality metrics
- Junction analysis

The ``locus`` field indicates the chain type (TRA, TRB, TRG, TRD).


Custom TCR Databases
--------------------

Build custom TCR databases using ``build_germline_database``:

**Command Line:**

.. code-block:: bash

    abstar build_germline_database my_tcr_db \
        -f tcr_v_genes.fasta \
        -f tcr_d_genes.fasta \
        -f tcr_j_genes.fasta \
        --receptor tcr

**Python:**

.. code-block:: python

    import abstar

    abstar.tl.build_germline_database(
        name="my_tcr_db",
        fastas=["tcr_v.fasta", "tcr_d.fasta", "tcr_j.fasta"],
        receptor="tcr"
    )

Then use the custom database:

.. code-block:: bash

    abstar run tcr.fasta output/ --receptor tcr --germline_database my_tcr_db


Differences from BCR Annotation
-------------------------------

- **No isotype assignment**: TCRs do not have isotypes, so ``c_call`` is
  not assigned for TCR sequences
- **D gene**: Only beta and delta chains have D gene segments
- **Constant regions**: TCR constant regions are not annotated
