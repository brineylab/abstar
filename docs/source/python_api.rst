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

    # C57BL/6 mouse sequences
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

    BCR options: ``human``, ``macaque``, ``c57bl6``, ``balbc``, and
    ``human+c57bl6``. TCR currently provides ``human``.

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

``copy_inputs_to_project``
    Copy source files into ``project_path/input/``. Default: ``False`` in
    Python (the CLI defaults to copying). For directory inputs, preserve paths
    relative to the original input directory, including nested directories.

``verbose``
    Print progress information. Default: ``False``

``strict``
    Abort on individual sequence annotation exceptions. Default: ``False``

``debug``
    Retain temp files and enable detailed logging. Default: ``False``


Return Types
~~~~~~~~~~~~

**When project_path is None (default):**

One input record that annotates or is unassigned returns an ``abutils.Sequence``; multiple input records return
a list of ``Sequence`` objects in input order. Lists, iterators, and generators
are supported, and arbitrary iterables are consumed once. Directory inputs are
discovered recursively in natural path order, with record order retained within
each file. This order is stable across worker counts and annotation chunk sizes.

Every returned or written row receives an ``annotation_status``. An ordinary
biological non-assignment returns ``"unassigned"`` with a ``failure_reason``;
it remains in the result and does not change the return shape. Sequence annotation
exceptions have no normal result row: they produce persistent diagnostics and a
warning, while remaining records continue. An all-failed input returns an empty
list/DataFrame. With ``strict=True``, sequence exceptions instead raise
``abstar.AnnotationRunError``. Worker, input, external-tool, and output/storage
failures always raise and appear in the exception's structured ``failures``.
``partial_output_paths`` contains only diagnostic or partial artifacts that
could be retained; it may be empty when storage fails before an artifact can be
preserved.

Visible identifiers are data, not join keys. Duplicate identifiers, leading
zeroes, and values such as ``10E8`` are preserved exactly. abstar uses a unique
private ``row_id`` while work is split and joined, then removes it at the public
boundary. Results remain in deterministic input order across process counts and
chunk sizes.

The required ``abutils.tl.translate`` capability is checked before starting
workers or creating a project; an incompatible installation raises
``preprocess/internal_error`` with installation guidance in the failure's
``message``. Malformed nonempty FASTA/FASTQ and non-IUPAC bases raise
``preprocess/invalid_input``. The accepted nucleotide alphabet is
``ACGTRYSWKMBDHVN`` (case-insensitive). FASTQ validation uses Biopython and supports
multiline sequence and quality data. Missing required components in the selected
germline database raise ``assignment/invalid_input``, naming the database,
receptor, and component. Content-empty inputs still raise ``ValueError`` before
creating a caller project.

Final writer failures raise ``output/internal_error``. AIRR and Parquet files
are staged together before publication; failed staging files are removed, while
annotation work and diagnostics remain inspectable. If publication itself fails
after a file was promoted, that file is listed in ``partial_output_paths`` and
must be treated as partial run output.
Files published for earlier samples are also listed as partial if a later
sample fails.

Without a caller project, ordinary temporary workspaces are cleaned after
success with no record errors. Recovered record failures retain their diagnostics
in the API workspace, with the index path reported in the warning. On a fatal
run error, failed work and logs are transferred to an
``abstar-failed-*`` directory in the system temporary directory, and the
surviving partial files and failure logs are listed in ``partial_output_paths``.
If that transfer fails, abstar keeps the surviving owned workspace files and
preserves the original structured error and its ``failures``. The error's
``retention_diagnostics`` tuple and message describe the secondary storage
failure. Only surviving files are listed in ``partial_output_paths``.
Caller projects retain their diagnostic and partial work files; current-run
logs are listed without including historical failures from earlier runs.
``debug=True`` retains the complete API workspace, including on success, and
lists surviving failure logs when the run raises.

Initial project, log, temporary, and output directory failures use
``output/internal_error``. If the project cannot store its diagnostic, abstar
retains a fallback log in an ``abstar-failed-*`` temporary directory and lists
its path in ``partial_output_paths``. If fallback storage also fails, the
original structured failure and exception cause remain available; its message
also describes the diagnostic storage errors, and no incomplete fallback file
is reported.

For a file containing multiple records:

.. code-block:: python

    sequences = abstar.run("input.fasta")
    for seq in sequences:
        print(seq.id)
        print(seq["v_call"])      # V gene assignment
        print(seq["cdr3_aa"])     # CDR3 amino acid sequence
        print(seq["productive"])  # Productivity status

**When as_dataframe=True:**

Returns a polars DataFrame, including for one input record:

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

Empty input iterables, empty or whitespace-only raw strings, and directories
without supported FASTA/FASTQ files raise ``ValueError`` before creating the
requested project directory. Files containing only whitespace or no content
(including gzip-compressed files) are skipped; if no nonempty files remain, the
run also raises before project creation. This preflight does not parse sequence
records: malformed nonempty files continue through the parser and retain their
structured failure diagnostics.

Process and chunk counts must be positive integers (not Booleans);
``n_processes=None`` selects the CPU count. Unsupported or empty output formats
also raise before project creation.


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


Coordinate conventions
----------------------

Python annotations, dataframe returns, and Parquet use zero-based half-open
query/reference coordinates. Region coordinates and V/D/J/C sequence coordinates
address ``sequence_oriented``. AIRR TSV converts these to one-based closed
intervals. Both final file formats write the original input as ``sequence``. See
:doc:`output_formats` for sequence, CIGAR, null, and migration details.
The official TSV/final Parquet ``sequence_aa`` uses the full oriented-query coding phase;
its paired AA alignment fields share the nucleotide alignment columns.
Python annotation objects and dataframe returns retain their assembled
``sequence`` and translations for compatibility. Project mode returns ``None``;
the file mapping does not change no-project API return shapes or values.
Use ``abstar.annotation.airr.to_airr_row()`` for the explicit serialization mapping.


Migration examples
------------------

A biological non-assignment could previously disappear from the result. Given
a two-record ``mixed.fasta`` with an assignable antibody named ``assigned`` and
a query ``N`` named ``unassigned``, the legacy result contained only the first
row:

.. code-block:: python

    # Before
    result = abstar.run("mixed.fasta")
    assert result.id == "assigned"  # one object despite two input records

The conserved result now contains both records in input order and makes the
biological outcome explicit:

.. code-block:: python

    import abstar

    # After
    result = abstar.run("mixed.fasta", n_processes=1, mmseqs_threads=1)
    assert len(result) == 2
    assert [row.id for row in result] == ["assigned", "unassigned"]
    assert [row["annotation_status"] for row in result] == ["annotated", "unassigned"]
    assert result[1]["sequence_input"] == "N"
    assert result[1]["failure_reason"] == "no compatible V gene assignment"
    assert result[1]["v_call"] is result[1]["j_call"] is None

If both input records use the identifier ``duplicate``, both returned IDs stay
``duplicate`` and remain in the same order; the private row identity is absent.

Sequence annotation exceptions are recoverable by default. Failed records have no
normal output row; the remaining records and samples continue. A ``RuntimeWarning``
reports the failure count and persistent ``logs/failures.tsv`` location. Without a
project directory, the temporary workspace containing diagnostics is retained;
ordinary scratch Parquets are cleaned unless ``debug=True``. An entirely failed
input returns an empty list or a schema-bearing empty DataFrame, with the warning
and diagnostics distinguishing it from biological non-assignment. Multiple-input
object returns remain lists even when only one record survives.

Migration: callers relying on annotation exceptions must now pass ``strict=True``.
Worker, input, external-tool, and output/storage failures remain fatal in either mode.

.. code-block:: python

    import abstar

    # Default: keep successful/unassigned records and log sequence errors.
    result = abstar.run("input.fasta")

    # Strict: retain the previous exception policy.
    try:
        abstar.run("input.fasta", "project/", strict=True)
    except abstar.AnnotationRunError as error:
        print(error.failures, error.partial_output_paths)
        raise

``failures`` contains immutable ``RecordFailure`` values. Stages are
``preprocess``, ``assignment``, ``annotation``, or ``output``; categories are
``invalid_input``, ``unassigned``, ``external_tool``, or ``internal_error``.
Ordinary unassigned rows are returned as data and do not raise. Paths in
``partial_output_paths`` are diagnostic or already-promoted artifacts, not a
successful complete result. The tuple is empty when no artifact could be
retained.
