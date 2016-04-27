Overview
========

With technical breakthroughs in the throughput and read-length of 
next-generation sequencing platforms, antibody repertoire sequencing 
is becoming an increasingly important tool for detailed characterization 
of the immune response to infection and immunization. Accordingly, 
there is a need for open, scalable software for the genetic analysis of 
repertoire-scale antibody sequence data.

We built AbStar to be a modular component of these analyses. 
AbStar is engineered to be highly flexible, capable of processing a single 
sequence or billions of sequences and scaling from individual laptops to
large clusters of cloud computing instances.

Workflows
---------

In addition to V(D)J germline gene assignment and primary antibody
sequence annotation, AbStar contains utilities for
sequence acquisition, pre-processing, and database import. AbStar also
exposes a high-level public API to many of the core functions, which allows
users to easily construct :ref:`custom analysis workflows <APIexamples>`
using multiple AbStar utilities as well as other third-party tools. To 
ease integration of AbStar into currently existing antibody analysis
pipelines based on IMGT, AbStar can optionally produce output
that mimics the IMGT-HighV/Quest Summary output file.

File formats
------------

AbTools accepts standard FASTA or FASTQ files and produces, by default,
JSON-formatted output. This output format allows us to build the output using
data structures that match the way we process data programatically.
JSON is also easily importable into NoSQL databases like MongoDB. We have 
found NoSQL databases to be very well suited for performing downstream 
analyses of antibody repertoire data, as the flexible schema allows for easy 
updating of sequence records with additional annotation information. Although 
additional data can be added to relational databases, querying this data 
often involves joining tables, which can require significant optimization 
for very large datasets.

Scalability
-----------

Cloud computing has dramatically changed the landscape of high-performance
computing (HPC), and has allowed small academic labs to 'rent' access
to computational resources that would have been previously far outside their 
budget. AbStar is tightly integrated with AbCloud_, which provides tools
for launching, configuring and managing clusters of compute instances on
Amazon's Elastic Compute Cloud (EC2). Using the Celery distributed task queue,
jobs are distributed to worker nodes and processed in parallel.

In order to maintain compatability with alternate cloud computing platforms
with minimal effort, an `AbStar Docker image`_ is also provided.

.. _AbCloud: https://github.com/briney/abcloud
.. _AbStar Docker Image: https://hub.docker.com/r/briney/abstar/
