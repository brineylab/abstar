abstar: scalable AIRR annotation
=================================

abstar is a tool for VDJ germline gene assignment and antibody/TCR sequence
annotation. It performs germline gene assignment using MMseqs2 and detailed
sequence annotation including mutations, indels, regions (CDR/FWR), and
productivity assessment. Scalable from single sequences to billions.


Key Features
------------

- **Fast**: MMseqs2-powered germline assignment scales to billions of sequences
- **AIRR-compliant**: Full compatibility with `AIRR data standards`_
- **BCR and TCR**: Support for both B-cell and T-cell receptor sequences
- **Flexible output**: AIRR TSV or Parquet formats
- **UMI support**: Built-in UMI extraction and parsing
- **Read merging**: Automatic paired-end read merging with fastp
- **Custom germline databases**: Build databases from OGRDB, IgDiscover, Digger, or FASTA

.. _AIRR data standards: https://docs.airr-community.org/en/latest/


Quick Example
-------------

**Command Line:**

.. code-block:: bash

    # Install
    pip install abstar

    # Annotate human BCR sequences
    abstar run sequences.fasta output_dir/

    # TCR sequences
    abstar run tcr.fasta output_dir/ --receptor tcr

**Python:**

.. code-block:: python

    import abstar

    # Return annotated Sequence objects
    sequences = abstar.run("sequences.fasta")

    # Return as polars DataFrame
    df = abstar.run("sequences.fasta", as_dataframe=True)


Input and Output
----------------

**Input**: FASTA or FASTQ files (gzip-compressed supported)

**Output**: AIRR-compliant TSV or Parquet files containing:

- V(D)J gene assignments
- CDR/FWR region sequences
- Mutation and indel annotations
- Productivity assessment
- Junction/CDR3 analysis


Germline Databases
------------------

Built-in databases: ``human``, ``mouse``, ``macaque``, ``humouse``

Human and mouse databases are based on the `OGRDB`_ germline reference sets.
Custom databases can be built from FASTA/JSON files.

.. _OGRDB: https://ogrdb.airr-community.org/


.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Getting Started

   installation
   quickstart


.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: User Guide

   cli
   python_api
   tcr


.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: API Reference

   api


.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: About

   license


.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Related Projects

   abutils <https://github.com/brineylab/abutils>
   scab <https://github.com/brineylab/scab>


Index
-----

* :ref:`genindex`
* :ref:`modindex`
