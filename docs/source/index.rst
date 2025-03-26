abstar: scalable AIRR annotation
===================================================================


Continuous improvements in the throughput of next-generation sequencing platforms 
have made adaptive immune receptor repertoire (AIRR) sequencing an increasingly 
important tool for detailed characterization of the immune response to infection 
and immunization. Accordingly, there is a need for open, scalable software for 
the genetic analysis of repertoire-scale antibody sequence data.

abstar is a core component of the ab[x] toolkit for antibody sequence analysis. 
abstar performs V(D)J germline gene assignment and antibody sequence annotation, and can readily scale
from a single sequence to billions of sequences. abstar is fully compliant with `AIRR <https://docs.airr-community.org/en/latest/>`_
data standards and produces annotated sequence data in `AIRR <https://docs.airr-community.org/en/latest/>`_ format.


usage
--------

abstar can be used both as a command-line tool and as a Python API. The command-line interface 
provides a straightforward way to process files containing antibody sequences, while the Python API 
allows for integration of abstar's annotation capabilities into custom analysis pipelines. For large 
datasets, abstar offers distributed processing capabilities to accelerate annotation of many sequences 
in parallel. Detailed usage instructions are available in the CLI and API documentation sections.


file formats
---------------

**Input Formats**: abstar accepts antibody sequences in FASTA and FASTQ formats, which are standard 
formats for storing nucleotide sequences. These can be raw sequencing output files (for example, paired-end reads
from Illumina or Element sequencing platforms) or pre-processed sequence data.

**Output Formats**: by default, abstar generates annotation results in AIRR-compliant TSV format. We 
also offer the ability to generate output in Parquet format, which is a columnar storage format that 
is more space-efficient for large datasets and can be faster for certain types of analysis. In either
case (TSV or Parquet), all output adheres to the standardized AIRR schema, ensuring interoperability 
with other tools in the immunoinformatics ecosystem.


getting started
---------------

.. toctree::
   :maxdepth: 1
   :caption: getting started

   installation


usage
-----

.. toctree::
   :maxdepth: 1
   :caption: usage

   cli
   api
   

about
-----

.. toctree::
   :maxdepth: 1
   :caption: about

   license
   news


related projects
----------------

.. toctree::
   :maxdepth: 1
   :caption: related projects

   abutils <https://github.com/brineylab/abutils>
   scab <https://github.com/brineylab/scab>


Index
-----

* :ref:`modindex`
* :ref:`search`


