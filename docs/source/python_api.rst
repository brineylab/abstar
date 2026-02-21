.. _python-api:

Python API
==========

abstar provides a Python API for integrating sequence annotation into
custom analysis pipelines.


The abstar.run() Function
-------------------------

The main entry point for annotation:

.. code-block:: python

    import abstar

    # Basic usage - returns annotated Sequence objects
    sequences = abstar.run("sequences.fasta")

    # Return as polars DataFrame
    df = abstar.run("sequences.fasta", as_dataframe=True)

    # Write to project directory
    abstar.run("sequences.fasta", "project/", output_format=["airr", "parquet"])

    # TCR annotation
    sequences = abstar.run("tcr.fasta", receptor="tcr")

    # C57BL/6 mouse sequences with built-in germline database
    sequences = abstar.run("sequences.fasta", germline_database="c57bl6")


Parameters
~~~~~~~~~~

``sequences``
    Input sequences. Can be:

    - Path to a FASTA/FASTQ file
    - Path to a directory of FASTA/FASTQ files
    - A single ``abutils.Sequence`` object
    - An iterable of ``Sequence`` objects

``project_path`` (optional)
    Directory for output files. If provided, results are written to disk
    and the function returns ``None``. If not provided, annotated sequences
    are returned.

``germline_database``
    Germline database name. Default: ``"human"``

    Built-in options for ``receptor="bcr"``:
    ``human``, ``macaque``, ``c57bl6``, ``balbc``, ``human+c57bl6``

    Built-in option for ``receptor="tcr"``: ``human``

``receptor``
    Receptor type: ``"bcr"`` (default) or ``"tcr"``

``output_format``
    Output format(s): ``"airr"`` (TSV), ``"parquet"``, or a list of both.
    Default: ``"airr"``

``as_dataframe``
    If ``True``, return a polars DataFrame instead of Sequence objects.
    Default: ``False``

``umi_pattern``
    Pattern for UMI extraction. See :doc:`umis` for details.

``umi_length``
    UMI length. Positive for 5' end, negative for 3' end.

``merge``
    Merge paired-end FASTQ files before annotation. Default: ``False``

``merge_kwargs``
    Additional arguments for the merge function as a dict.

``chunksize``
    Sequences per annotation batch. Default: ``500``

``mmseqs_chunksize``
    Sequences per MMseqs2 batch. Default: ``1000000``

``mmseqs_threads``
    Threads for MMseqs2. Default: auto-detected

``n_processes``
    Parallel annotation workers. Default: CPU count

``verbose``
    Print progress information. Default: ``False``

``debug``
    Retain temp files and enable detailed logging. Default: ``False``


Return Types
~~~~~~~~~~~~

**When project_path is None (default):**

Returns annotated ``abutils.Sequence`` objects:

.. code-block:: python

    sequences = abstar.run("input.fasta")
    for seq in sequences:
        print(seq.id)
        print(seq["v_call"])      # V gene assignment
        print(seq["cdr3_aa"])     # CDR3 amino acid sequence
        print(seq["productive"])  # Productivity status

**When as_dataframe=True:**

Returns a polars DataFrame:

.. code-block:: python

    import polars as pl

    df = abstar.run("input.fasta", as_dataframe=True)

    # Filter productive sequences
    productive = df.filter(pl.col("productive") == True)

    # Group by V gene
    v_gene_counts = df.group_by("v_gene").len()

**When project_path is provided:**

Returns ``None``; writes files to project directory:

.. code-block:: python

    abstar.run("input.fasta", "project/")

    # Output files:
    #   project/airr/input.tsv
    #   project/logs/abstar.log


Module Namespaces
-----------------

abstar.gl - Germline Functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Access germline sequences and database paths.

**Get database path:**

.. code-block:: python

    import abstar

    # Get path to built-in human BCR database
    path = abstar.gl.get_germline_database_path("human", receptor="bcr")

    # Get path to custom database
    path = abstar.gl.get_germline_database_path("my_custom_db")

**Get germline sequences:**

.. code-block:: python

    # Get a specific allele
    vgene = abstar.gl.get_germline("IGHV1-2*02", "human")
    print(vgene.sequence)

    # Get all alleles of a gene (returns list)
    alleles = abstar.gl.get_germline("IGHV1-2", "human")

    # Get IMGT-gapped sequence
    vgene_gapped = abstar.gl.get_germline("IGHV1-2*02", "human", imgt_gapped=True)


abstar.pp - Preprocessing
~~~~~~~~~~~~~~~~~~~~~~~~~

Paired-end read merging.

.. code-block:: python

    import abstar

    # Merge paired FASTQ files in a directory
    merged_files = abstar.pp.merge_fastqs(
        "fastq_directory/",
        "merged_output/",
        schema="illumina"  # or "element"
    )

    # With quality trimming options
    merged_files = abstar.pp.merge_fastqs(
        "fastq_directory/",
        "merged_output/",
        minimum_overlap=30,
        quality_cutoff=20,
        trim_adapters=True
    )

**Parameters:**

- ``schema``: Filename schema (``"illumina"`` or ``"element"``)
- ``minimum_overlap``: Minimum overlap for merging (default: 30)
- ``allowed_mismatches``: Allowed mismatches in overlap (default: 5)
- ``trim_adapters``: Trim adapters (default: True)
- ``quality_trim``: Quality trim (default: True)
- ``quality_cutoff``: Quality threshold (default: 20)


abstar.tl - Tools
~~~~~~~~~~~~~~~~~

Utility functions for database building and UMI parsing.

**Build custom germline database:**

.. code-block:: python

    import abstar

    abstar.tl.build_germline_database(
        name="my_database",
        fastas=["v_genes.fasta", "d_genes.fasta", "j_genes.fasta"],
        constants=["c_genes.fasta"],
        receptor="bcr"
    )

**Parse UMIs from sequences:**

.. code-block:: python

    # Parse UMIs and return annotated sequences
    umi = abstar.tl.parse_umis(
        "sequence_string_or_file",
        pattern="[UMI]TCAGCGGGAAGACATT",
        length=12
    )

See :doc:`umis` for detailed UMI documentation.


Examples
--------

**Basic annotation pipeline:**

.. code-block:: python

    import abstar

    # Annotate sequences
    sequences = abstar.run("sequences.fasta")

    # Filter productive sequences
    productive = [s for s in sequences if s["productive"]]

    # Extract CDR3 sequences
    cdr3_sequences = [s["cdr3_aa"] for s in productive]

**Large-scale processing:**

.. code-block:: python

    import abstar

    # Process with multiple output formats
    abstar.run(
        "large_dataset/",
        "output/",
        output_format=["airr", "parquet"],
        n_processes=16,
        mmseqs_threads=8
    )

**DataFrame analysis:**

.. code-block:: python

    import abstar
    import polars as pl

    df = abstar.run("sequences.fasta", as_dataframe=True)

    # Analyze V gene usage
    v_usage = (
        df.filter(pl.col("productive") == True)
        .group_by("v_gene")
        .len()
        .sort("len", descending=True)
    )
    print(v_usage)
