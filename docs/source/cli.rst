cli
===============

Running ``abstar`` from the command line is designed to be relatively straightforward, even for users with
minimal experience with command-line applications. In the most
basic case, with a single input file of human antibody sequences, ``abstar`` can be run with:

.. code-block:: bash

    abstar run path/to/sequences.fasta path/to/output/

``abstar`` will process all sequences contained in ``sequences.fasta`` and the
results, including annotations and logs, will be deposited into ``path/to/output/``. 

.. note::
    If ``path/to/output/`` does not exist, it will be created.

If you have a directory of FASTA/Q-formatted files for ``abstar`` to process, you 
can pass a directory instead of a single file, and all files in the directory will be processed:

.. code-block:: bash

    abstar run path/to/input/ path/to/output/


`abstar` has many other options that can be used to customize the analysis or to handle more complex cases. 
This includes:

- :ref:`merging paired-end reads <read_merging>` prior to annotation
- :ref:`UMI detection <umis>` and incorporation into the annotated output
- creation and use of :ref:`custom germline databases <germline_dbs>`


.. toctree::
    :hidden:
    
    outputs <output_formats>
    read merging <read_merging>
    umis <umis>
    germline databases <germline_dbs>

