cli
===============

Running ``abstar`` from the command line is designed to be relatively straightforward, even for users with
minimal experience with command-line applications. In the most
basic case, with a single input file of human antibody sequences, ``abstar`` can be run with:

.. code-block:: bash

    abstar run path/to/sequences.fasta path/to/output/

``abstar`` will process all sequences contained in ``sequences.fasta`` and the
results, including annotations and logs, will be deposited into ``path/to/output/``. If
 ``path/to/output/`` does not exist, it will be created.

If you have a directory of FASTA/Q-formatted files for ``abstar`` to process, you 
can pass a directory instead of a single file, and all files in the directory will be processed:

.. code-block:: bash

    abstar run path/to/input/ path/to/output/



.. toctree::
    :hidden:
    
    outputs <output_formats>
    read merging <read_merging>
    umis <umis>
    germline databases <germline_dbs>

