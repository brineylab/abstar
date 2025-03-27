.. _output-formats:

output formats
==============

``abstar`` can output annotations in two different file formats, each of which contains identical data 
in `AIRR-compatible format <https://docs.airr-community.org/en/stable/datarep/rearrangements.html>`. 
The two formats are:

- ``airr``: a tab-delimited file with a header row
- ``parquet``: a columnar, binary file format that is more space-efficient for large datasets and 
    faster to read into memory

By default, ``abstar`` will output annotations in the ``airr`` format. The format is specified using 
the ``--output_format`` (or ``-o``) option:

.. code-block:: bash

    abstar run path/to/sequences.fasta path/to/output/ --output_format parquet

Multiple output formats can be specified by passing the ``--output_format`` (``-o``) option multiple times:

.. code-block:: bash

    abstar run path/to/sequences.fasta path/to/output/ -o airr -o parquet

In this case, ``abstar`` will produce annotations in both ``airr`` and ``parquet`` formats.


