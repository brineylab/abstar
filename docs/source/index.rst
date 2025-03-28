abstar: scalable AIRR annotation
===================================================================


Continuous improvements in the throughput of next-generation sequencing platforms 
have made adaptive immune receptor repertoire (AIRR) sequencing an increasingly 
important tool for detailed characterization of the immune response to infection 
and immunization. Accordingly, there is a need for open, scalable software for 
the genetic analysis of repertoire-scale antibody sequence data.

``abstar`` is a core component of the ab[x] toolkit for antibody sequence analysis. 
``abstar`` performs V(D)J germline gene assignment and antibody sequence annotation, and can readily scale
from a single sequence to billions of sequences. ``abstar`` is fully compliant with `AIRR <https://docs.airr-community.org/en/latest/>`_
data standards and produces annotated sequence data in `AIRR <https://docs.airr-community.org/en/latest/>`_ format.


usage
--------

``abstar`` can be used both as a command-line tool and as a Python API. The command-line interface 
provides a straightforward way to process files containing antibody sequences, while the Python API 
allows for integration of ``abstar``'s annotation capabilities into custom analysis pipelines. For large 
datasets, ``abstar`` offers distributed processing capabilities to accelerate annotation of many sequences 
in parallel. Detailed usage instructions are available in the CLI and API documentation sections.


file formats
---------------

**Input Formats**: ``abstar`` accepts antibody sequences in FASTA and FASTQ formats, which are standard 
formats for storing nucleotide sequences. These can be raw sequencing output files (for example, paired-end reads
from Illumina or Element sequencing platforms) or pre-processed sequence data.

**Output Formats**: by default, ``abstar`` generates annotation results in AIRR-compliant TSV format. We 
also offer the ability to generate output in Parquet format, which is a columnar storage format that 
is more space-efficient for large datasets and can be faster for certain types of analysis. In either
case (TSV or Parquet), all output adheres to the standardized AIRR schema, ensuring interoperability 
with other tools in the immunoinformatics ecosystem.


germline databases
-------------------

``abstar`` comes pre-packaged with built-in germline databases for human, macaque, and mouse. The human and mouse databases
are based on the `Open Germline Receptor Database  (OGRDB) <https://ogrdb.airr-community.org/>`_ germline 
reference sets. ``abstar`` also supports the use and creation of custom germline databases, which can 
be used to annotate sequences from species that are not included in the built-in databases or to use
donor-specific databases created using tools like `IgDiscover <https://www.nature.com/articles/ncomms13642>`_ 
or `Digger <https://academic.oup.com/bioinformatics/article/40/3/btae144/7628126>`_.



.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: getting started

   installation


.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: usage

   cli
   api
   

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: about

   license


.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: related projects

   abutils <https://github.com/brineylab/abutils>
   scab <https://github.com/brineylab/scab>


index
-----

* :ref:`genindex`
* :ref:`modindex`


