![](https://img.shields.io/pypi/v/abstar.svg?colorB=blue)
[![Build Status](https://travis-ci.com/briney/abstar.svg?branch=master)](https://travis-ci.com/briney/abstar)
[![Documentation Status](https://readthedocs.org/projects/abstar/badge/?version=latest)](https://abstar.readthedocs.io/en/latest/?badge=latest)
![](https://img.shields.io/pypi/pyversions/abstar.svg)
![](https://img.shields.io/badge/license-MIT-blue.svg)

# abstar  
  
VDJ assignment and antibody sequence annotation. Scalable from a single sequence to billions of sequences.  
  
  - Source code: [github.com/briney/abstar](https://github.com/briney/abstar)  
  - Documentation: [abstar.readthedocs.org](http://abstar.readthedocs.org)  
  - Download: [pypi.python.org/pypi/abstar](https://pypi.python.org/pypi/abstar)  
  - Docker: [hub.docker.com/r/briney/abstar/](https://hub.docker.com/r/briney/abstar/)  
  
### install  
`pip install abstar`  
  
### use  

To run abstar on a single FASTA or FASTQ file:  
`abstar -i <input-file> -o <output-directory> -t <temp-directory>`

To iteratively run abstar on all files in an input directory:  
`abstar -i <input-directory> -o <output-directory> -t <temp-directory>`
  
To run abstar using the included test data as input:  
`abstar -o <output-directory> -t <temp-directory> --use-test-data`  
  
When using the abstar test data, note that although the test data file contains 1,000 sequences, one of the test sequences is not a valid antibody recombination. Only 999 sequences should be processed successfully.  

When using BaseSpace as the input data source, you can optionally provide all of the required directories:  
`abstar -i <input-directory> -o <output-directory> -t <temp-directory> -b`  
  
Or you can simply provide a single project directory, and all required directories will be created in the project directory:  
`abstar -p <project_directory> -b`  
  
  
### additional options  
`-l LOG_LOCATION, --log LOG_LOCATION` Change the log directory location. Default is the parent directory of `<output_directory>`.  
  
`-m, --merge` Input directory should contain paired FASTQ (or gzipped FASTQ) files. Paired files will be merged with PANDAseq prior to processing with abstar. Note that when using the BaseSpace option (`-b, --basespace`), this option is implied.  
  
`-b, --basespace` Download a sequencing run from BaseSpace, which is Illumina's cloud storage environment. Since Illumina sequencers produce paired-end reads, `--merge` is implied.  
  
`-u N, --uaid N` Sequences contain a unique antibody ID (UAID, or molecular barcode) of length N. The uaid will be parsed from the beginning of each input sequence and added to the JSON output. Negative values result in the UAID being parsed from the end of the sequence.  
  
`-s SPECIES, --species SPECIES` Select the species from which the input sequences are derived. Supported options are 'human', 'mouse', and 'macaque'. Default is 'human'.  
   
`-c, --cluster` Runs abstar in distributed mode on a Celery cluster.  
  
`-h, --help` Prints detailed information about all runtime options.
  
`-D --debug` Much more verbose logging.  
  

### api  
Most core abstar functions are available through a public API, making it easier to run abstar as a component of integrated analysis pipelines. See the abstar [documentation](http://abstar.readthedocs.org) for more detail about the API.  
  
  
### helper scripts  
A few helper scripts are included with abstar:  
`batch_mongoimport` automates the import of multiple JSON output files into a MongoDB database.  
`build_abstar_germline_db` creates abstar germline databases from IMGT-gapped FASTA files of V, D and J gene segments.  
`make_basespace_credfile` makes a credentials file for BaseSpace, which is required if downloading sequences from BaseSpace with abstar. Developer credentials are required, and the process for obtaining them is explained [here](https://support.basespace.illumina.com/knowledgebase/articles/403618-python-run-downloader)  
  

### testing  
To run the test suite, clone or download the repository and run `pytest ./` from the top-level directory.  
  
### requirements  
Python 3.8+   
abutils  
biopython  
celery  
nwalign3  
pymongo  
pytest  
scikit-bio  

All of the above dependencies can be installed with pip, and will be installed automatically when installing abstar with pip.  
If you're new to Python, a great way to get started is to install the [Anaconda Python distribution](https://www.continuum.io/downloads), which includes pip as well as a ton of useful scientific Python packages.
  
sequence merging requires [PANDAseq](https://github.com/neufeld/pandaseq)  
batch_mongoimport requires [MongoDB](http://www.mongodb.org/)  
BaseSpace downloading requires the [BaseSpace Python SDK](https://github.com/basespace/basespace-python-sdk)  
