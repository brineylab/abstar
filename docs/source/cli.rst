.. _cli:

CLI Reference
=============

abstar provides two main commands: ``run`` for annotating sequences and
``build_germline_database`` for creating custom germline databases.


abstar run
----------

Annotate antibody or TCR sequences.

**Usage:**

.. code-block:: bash

    abstar run INPUT_PATH PROJECT_PATH [OPTIONS]

**Arguments:**

``INPUT_PATH``
    Path to a FASTA/Q file or a directory of FASTA/Q files.
    Gzip-compressed files (``.gz``) are supported.

``PROJECT_PATH``
    Directory for output files, logs, and temporary data.
    Created if it does not exist.

**Basic Examples:**

.. code-block:: bash

    # Single file
    abstar run sequences.fasta output_dir/

    # Directory of files
    abstar run fastq_directory/ output_dir/

    # Mouse sequences
    abstar run sequences.fasta output_dir/ --germline_database mouse

    # TCR sequences
    abstar run tcr.fasta output_dir/ --receptor tcr


Germline and Receptor Options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``--germline_database TEXT``
    Name of the germline database for assignment.

    Built-in options: ``human`` (default), ``mouse``, ``macaque``, ``humouse``

    Custom databases in ``~/.abstar/germline_dbs/`` are also available.

    .. code-block:: bash

        abstar run seqs.fasta out/ --germline_database mouse

``--receptor [bcr|tcr]``
    Receptor type. Default: ``bcr``

    .. code-block:: bash

        abstar run tcr_seqs.fasta out/ --receptor tcr


Output Options
~~~~~~~~~~~~~~

``-o, --output_format [airr|parquet]``
    Output format. Can be specified multiple times for multiple formats.
    Default: ``airr``

    .. code-block:: bash

        # Parquet only
        abstar run seqs.fasta out/ -o parquet

        # Both formats
        abstar run seqs.fasta out/ -o airr -o parquet

``--copy-inputs/--no-copy-inputs``
    Copy input files to project directory. Default: ``--copy-inputs``


UMI Options
~~~~~~~~~~~

``--umi_pattern TEXT``
    Pattern for UMI extraction. Use ``[UMI]`` as placeholder for the UMI sequence,
    with surrounding conserved sequences for anchoring.

    Built-in patterns: ``smartseq-human-bcr``

    .. code-block:: bash

        # Built-in pattern
        abstar run seqs.fasta out/ --umi_pattern smartseq-human-bcr

        # Custom pattern
        abstar run seqs.fasta out/ --umi_pattern "[UMI]TCAGCGGGAAGACATT" --umi_length 12

``--umi_length INT``
    Length of the UMI sequence to extract.

    - Positive value: UMI at 5' end of sequence
    - Negative value: UMI at 3' end of sequence

    .. code-block:: bash

        # 12bp UMI at 5' end
        abstar run seqs.fasta out/ --umi_length 12

        # 8bp UMI at 3' end
        abstar run seqs.fasta out/ --umi_length -8


Read Merging Options
~~~~~~~~~~~~~~~~~~~~

``-m, --merge``
    Merge paired-end FASTQ files before annotation using fastp.

    .. code-block:: bash

        abstar run paired_reads/ out/ --merge

``--merge_kwargs TEXT``
    Additional arguments for the merge function.
    Format: ``key1=val1,key2=val2``

    .. code-block:: bash

        abstar run reads/ out/ --merge --merge_kwargs "minimum_overlap=20,quality_cutoff=25"

``--interleaved_fastq``
    Input FASTQ files are interleaved (R1 and R2 alternating in single file).
    Implies ``--merge``.


Performance Options
~~~~~~~~~~~~~~~~~~~

``-c, --chunksize INT``
    Number of sequences per annotation batch. Default: ``500``

``--mmseqs_chunksize INT``
    Number of sequences per MMseqs2 search batch. Default: ``1000000``

``--mmseqs_threads INT``
    Number of threads for MMseqs2 searches. Default: auto-detected

``--n_processes INT``
    Number of parallel annotation workers. Default: number of CPU cores


Logging Options
~~~~~~~~~~~~~~~

``--verbose/--quiet``
    Enable or disable verbose output. Default: ``--verbose``

``--debug``
    Retain temporary files and enable detailed logging.
    Useful for troubleshooting.


abstar build_germline_database
------------------------------

Build a custom germline database from FASTA or JSON files.

**Usage:**

.. code-block:: bash

    abstar build_germline_database NAME [OPTIONS]

**Arguments:**

``NAME``
    Name for the new germline database.

**Options:**

``-f, --fasta PATH``
    FASTA file containing gapped germline gene sequences.
    Can be specified multiple times.

``-j, --json PATH``
    JSON file containing germline gene sequences.
    Can be specified multiple times.

``-c, --constant PATH``
    FASTA file containing constant region sequences.
    Can be specified multiple times.

``-m, --manifest PATH``
    Plain text file with database metadata (origin, date, etc.).

``-r, --receptor [bcr|tcr]``
    Receptor type. Default: ``bcr``

``-l, --location PATH``
    Custom location for the database.
    Default: ``~/.abstar/germline_dbs/``

``--reference TEXT``
    Reference species for adding IMGT gaps to ungapped sequences.
    Default: ``human``

``--include_species_in_name/--exclude_species_from_name``
    Include species in sequence names (e.g., ``IGHV1-2*02__homo_sapiens``).
    Default: ``--include_species_in_name``

**Example:**

.. code-block:: bash

    abstar build_germline_database my_custom_db \
        -f v_genes.fasta \
        -f d_genes.fasta \
        -f j_genes.fasta \
        -c constant_genes.fasta \
        --receptor bcr


Related Topics
--------------

.. toctree::
    :maxdepth: 1

    output_formats
    read_merging
    umis
    germline_dbs
