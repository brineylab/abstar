Overview
========

With technical breakthroughs in the throughput and read-length of 
next-generation sequencing platforms, antibody repertoire sequencing 
is becoming an increasingly important tool for detailed characterization 
of the immune response to infection and immunization. Accordingly, 
there is a need for open, scalable software for the genetic analysis of 
repertoire-scale antibody sequence data.

We built AbStar to be a modular, scalable component of these analyses. 
AbStar is engineered to be highly scalable, capable of processing a single 
sequence or billions of sequences and scaling from individual laptops and
large clusters of cloud computing instances.

Workflows
---------

In addition to V(D)J germline gene assignment and primary antibody
sequence annotation, AbStar contains additional utilities for
sequence acquisition, pre-processing, and database import. To ease
integration of AbStar into currently existing antibody analysis
pipelines based on IMGT, AbStar can optionally produce output
that mimics the IMGT-HighV/Quest Summary output file. AbStar also
exposes a public API to many of the core functions, which allows
users to easily construct custom analysis workflows using multiple
AbStar utilities as well as other third-party tools.

File formats
------------

AbTools accepts standard FASTA or FASTQ files and produces, by default,
JSON-formatted output. This format allows us to build output using
data structures that match the way we process data programatically: 
equivalents to Python dictionaries, lists, etc. JSON is also easily 
importable into NoSQL databases like MongoDB. We have found NoSQL databases
to be very advantageous when performing secondary analyses, as the
flexible schema allows for easy updating of sequence records with
additional annotation.
