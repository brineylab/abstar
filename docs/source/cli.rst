cli
===============

Running ``abstar`` from the command line is designed to be relatively straightforward, even for users with
minimal experience with command-line applications. In the most
basic case, with a single input file of human antibody sequences, ``abstar`` can be run with::

    $ abstar run path/to/sequences.fasta path/to/output/

``abstar`` will process all sequences contained in ``sequences.fasta`` and the
results, including annotations and logs, will be deposited into ``path/to/output/``. If
 ``path/to/output/`` does not exist, it will be created.

If you have a directory of FASTA/Q-formatted files for ``abstar`` to process, you 
can pass a directory instead of a single file, and all files in the directory will be processed::

    $ abstar run path/to/input/ path/to/output/

|

**Merging paired-end reads**

For input directories that contain paired FASTQ files that need to be merged
prior to processing, passing the ``--merge`` (``-m``) flag instructs ``abstar`` to merge
paired files with `fastp <https://github.com/OpenGene/fastp>`_::

    $ abstar run --merge path/to/input/ path/to/output/

Reads will be merged prior to annotation, and the merged reads will be deposited into a ``merged`` 
sub-directory located within the output directory. ``abstar`` can also accomodate interleaved 
FASTQ files, in which paired-end reads are stored in a single file. This is common when using 
data downloaded from the `SRA <https://www.ncbi.nlm.nih.gov/sra>`_. To process interleaved 
FASTQ files, pass the ``--interleaved_fastq`` flag::

    $ abstar run --interleaved_fastq path/to/input/ path/to/output/

The ``--interleaved_fastq`` flag will automatically trigger read merging, so the ``--merge`` 
flag is not necessary.

|

**Germline database creation**

``abstar`` can be used to create germline databases for use with ``abstar``. The ``build_germline_database`` 
command can be used to create germline databases from a variety of input formats.


