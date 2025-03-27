.. _read-merging:

merging paired-end reads
=========================

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
|
  .. note::
    The ``--interleaved_fastq`` flag will automatically trigger read merging, so the ``--merge`` 
    flag is not necessary.