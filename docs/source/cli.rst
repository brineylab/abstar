Commandline Use
===============

Running AbStar from the command line is simple. In the most
basic case, with a single input file of human antibody sequences::

    $ abstar -i /path/to/mydata.fasta -t /path/to/temp/ -o /path/to/output/

AbStar will process all sequences contained in ``mydata.fasta`` and the
results will be written to ``/path/to/output/mydata.json``. If either (or both)
of ``/path/to/temp/`` or ``/path/to/output/`` don't exist, they will be created.

If you have a directory of FASTA/Q-formatted files for AbStar to process, you 
can pass a directory via ``-i`` and all files in the directory will be processed::

    $ abstar -i /path/to/input/ -t /path/to/temp/ -o /path/to/output/

For input directories that contain paired FASTQ files that need to be merged
prior to processing, passing the ``-m`` flag instructs AbStar to merge
paired files with PANDAseq::

    $ abstar -i /path/to/input/ -t /path/to/temp/ -o /path/to/output/ -m

By default, AbStar will use PANDAseq's ``simple_bayesian`` merging algorithm, and 
alternate merging algorithms can be selected with ``pandaseq-algo``. 

For data generated with Illumina sequencers, AbStar can directly interface with
BaseSpace to download raw sequencing data. First, you need BaseSpace access token. 
The easiest way to do this is to set up a BaseSpace developer account following
`these instructions <https://support.basespace.illumina.com/knowledgebase/articles/403618-python-run-downloader>`_. 
Once you have your credentials, you can generate a BaseSpace credentials file by::

    $ make_basespace_credfile

and following the instructions. The simplest way to use AbStar with BaseSpace is
to pass a project directory using ``-p``. AbStar will then create all of the
necessary directories within the project directory. Running AbStar with the ``-b``
option indicates that input data should be downloaded from BaseSpace::

    $ abstar -p /path/to/project_dir/ -b

A list of available BaseSpace projects will be displayed and you can select the 
appropriate project. If downloading data from BaseSpace, ``-m`` is assumed and
paired-end reads will be merged. AbStar uses a human germline database by default, 
but germline databases are also provided for macaque, mouse and rabbit. To process 
macaque antibody sequences (from BaseSpace)::

    $ abstar -p /path/to/project_dir/ -b -s macaque
