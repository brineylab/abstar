 
.. _germline-dbs:

germline databases
=========================

``abstar`` comes pre-packaged with built-in germline databases for human, macaque, and mouse. 
The default germline database is human, but a different germline database can be specified with 
``--germline_database``:

.. code-block:: bash

    abstar run --germline_database mouse path/to/sequences.fasta path/to/output/


``abstar`` can also create custom germline databases, either for a species that is not included in 
the set of built-in germline databases or for donor-specific databases created using tools like 
`IgDiscover <https://www.nature.com/articles/ncomms13642>`_ or `Digger <https://academic.oup.com/bioinformatics/article/40/3/btae144/7628126>`_.
The ``build_germline_database`` command can be used to create germline databases using one of two
different types of input files (or even a mix of the two):

* FASTA-formatted files containing IMGT-gapped germline gene sequences
* JSON-formatted files containing germline gene sequences in `AIRR <https://docs.airr-community.org/en/latest/>`_ 
  format, such as those from `OGRDB <https://ogrdb.airr-community.org/>`_.

FASTA-formatted files can be supplied using ``--fasta`` (or ``-f``), which can be used 
multiple times to specify multiple files. JSON-formatted files can be supplied using  ``--json`` 
(or ``-j``), which can also be used multiple times to specify multiple files. The files can contain
a mix of V, D or J gene sequences. Constant region sequences can be supplied as FASTA-formatted file(s) 
using ``--constants`` (or ``-c``), which can also be used multiple times to specify multiple files. 
An example command for creating a database named ``my_germline_db`` might look like this:

.. code-block:: bash

    abstar build_germline_database my_germline_db -f germlines.fasta -j more_germlines.json -c constants.fasta

|

**Germline database location**

By default, germline databases are deposited in ``~/.abstar/germline_dbs/``. This can be changed
using the ``-l`` (or ``--location``) option, which can be used to specify a alternative location. 
``abstar`` will only look for custom germline databases in ``~/.abstar/germline_dbs/``, so 
the option to specify a custom location is provided primarily for testing purposes.

.. warning::
    When running ``abstar``, databases in ``~/.abstar/germline_dbs/`` will have priority over 
    the built-in databases. This means that if a custom database named ``human`` exists in 
    ``~/.abstar/germline_dbs/``, it will be used instead of the built-in human database.

|

**Receptor type**

The ``-r`` (or ``--receptor``) option can be used to specify the receptor type for the germline database. 
The default receptor type is ``bcr``, but ``tcr`` can also be specified.

|

**Manifest files**

An optional manifest file can be supplied using the ``--manifest`` (or ``-m``) option. A manifest file 
is a text file (of any format) that contains supplementary information about the germline database. For example,
the the manifest file could contain information about the source of the germline sequences, the download 
date of the germline sequences, or any other relevant information. An example using the ``-m`` option might look 
like this:

.. code-block:: bash

    abstar build_germline_database my_germline_db -f germlines.fasta -m manifest.txt

|

**Species names**

The ``--include_species_in_name`` option can be used to include the species name in the 
name of each sequence in the germline database. This option is provided primarily to simplify the creation of multi-species databases that 
may result in duplicate germline gene names. This is useful when analyzing data from, for example, transgenic 
mouse models that contain one or more human sequences in addition to the mouse sequences. 

.. note::
    The ``--include_species_in_name`` option is only applicable when using JSON-formatted files as input.

The resulting 
germline database will have unique sequence names like so: ``IGHV1-2*02__homo_sapiens``. When processing 
data with a multi-species database, ``abstar`` will automatically remove the species when populating the 
germline call fields, and the species name will be included in the ``species`` field. For example, ``IGHV1-2*02__homo_sapiens`` 
will be truncated to ``IGHV1-2*02`` when populating the ``v_call`` field, and the ``species`` field will be populated 
with ``homo_sapiens``. For example:

.. code-block:: bash

    abstar build_germline_database my_germline_db -j human.json -j mouse.json --include_species_in_name













