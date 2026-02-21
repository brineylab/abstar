.. _quickstart:

Quick Start
===========

This guide will get you annotating antibody and TCR sequences in minutes.


Basic BCR Annotation
--------------------

**Command Line:**

.. code-block:: bash

    # Annotate a FASTA file (human BCR, AIRR TSV output)
    abstar run sequences.fasta output_dir/

    # Annotate a directory of FASTA/FASTQ files
    abstar run fastq_directory/ output_dir/

**Python:**

.. code-block:: python

    import abstar

    # Return annotated Sequence objects
    sequences = abstar.run("sequences.fasta")

    # Access annotation fields
    for seq in sequences:
        print(seq.id, seq["v_call"], seq["cdr3_aa"])


Different Species
-----------------

**Command Line:**

.. code-block:: bash

    # C57BL/6 mouse sequences
    abstar run sequences.fasta output_dir/ --germline_database c57bl6

    # Macaque sequences
    abstar run sequences.fasta output_dir/ --germline_database macaque

**Python:**

.. code-block:: python

    sequences = abstar.run("sequences.fasta", germline_database="c57bl6")


TCR Annotation
--------------

**Command Line:**

.. code-block:: bash

    abstar run tcr_sequences.fasta output_dir/ --receptor tcr

**Python:**

.. code-block:: python

    sequences = abstar.run("tcr_sequences.fasta", receptor="tcr")


Output Formats
--------------

**Command Line:**

.. code-block:: bash

    # Parquet format (efficient for large datasets)
    abstar run sequences.fasta output_dir/ -o parquet

    # Both AIRR TSV and Parquet
    abstar run sequences.fasta output_dir/ -o airr -o parquet

**Python:**

.. code-block:: python

    # Return as polars DataFrame
    df = abstar.run("sequences.fasta", as_dataframe=True)

    # Write to files
    abstar.run("sequences.fasta", "output_dir/", output_format=["airr", "parquet"])


Paired-End Read Merging
-----------------------

For paired-end FASTQ files from Illumina or Element sequencers:

**Command Line:**

.. code-block:: bash

    abstar run paired_reads_directory/ output_dir/ --merge

**Python:**

.. code-block:: python

    sequences = abstar.run("paired_reads_directory/", merge=True)


Common Options
--------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Option
     - Description
   * - ``--germline_database``
     - Species database (BCR): ``human``, ``macaque``, ``c57bl6``, ``balbc``, ``human+c57bl6``
   * - ``--receptor``
     - ``bcr`` (default) or ``tcr``
   * - ``-o``, ``--output_format``
     - ``airr`` (TSV) or ``parquet``
   * - ``-m``, ``--merge``
     - Merge paired-end FASTQ files before annotation
   * - ``--n_processes``
     - Number of parallel annotation workers
   * - ``--umi_pattern``
     - Pattern for UMI extraction
   * - ``--debug``
     - Retain temp files and enable detailed logging


Output Directory Structure
--------------------------

After running abstar, your output directory will contain:

.. code-block:: text

    output_dir/
    ├── airr/           # AIRR-formatted TSV files
    │   └── sequences.tsv
    ├── parquet/        # Parquet files (if requested)
    │   └── sequences.parquet
    ├── logs/           # Log files
    │   └── abstar.log
    └── tmp/            # Temporary files (removed unless --debug)


Next Steps
----------

- :doc:`cli` - Complete CLI reference
- :doc:`python_api` - Python API documentation
- :doc:`output_formats` - Understanding output fields
- :doc:`germline_dbs` - Custom germline databases
- :doc:`umis` - UMI extraction
