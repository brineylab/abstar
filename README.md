![](https://img.shields.io/pypi/v/abstar.svg?colorB=blue)
[![tests](https://github.com/brineylab/abstar/actions/workflows/pytest.yml/badge.svg)](https://github.com/brineylab/abstar/actions/workflows/pytest.yml)
[![Documentation Status](https://readthedocs.org/projects/abstar/badge/?version=latest)](https://abstar.readthedocs.io/en/latest/?badge=latest)
![](https://img.shields.io/pypi/pyversions/abstar.svg)
![](https://img.shields.io/badge/license-MIT-blue.svg)

# abstar  
  
VDJ assignment and antibody sequence annotation. Scalable from a single sequence to billions of sequences.  
  
  - Source code: [github.com/brineylab/abstar](https://github.com/brineylab/abstar)  
  - Documentation: [abstar.readthedocs.org](http://abstar.readthedocs.org)  
  - Download: [pypi.python.org/pypi/abstar](https://pypi.python.org/pypi/abstar)  
  - Docker: [hub.docker.com/r/brineylab/datascience/](https://hub.docker.com/r/brineylab/datascience/)  
  
## install  
``` bash
pip install abstar
```  
  
## use  

To run abstar on a single FASTA or FASTQ file, you need to supply the input file and the project directory (into which output and logs will be written):  
``` bash
abstar path/to/sequences.fasta path/to/project_directory
```

To iteratively run abstar on all files in an input directory, pass a directory containing FASTA or FASTQ files instead of the path to a single file:  
``` bash
abstar path/to/input_directory path/to/project_directory
```
  
  
## additional options  
`abstar` contains a number of additional options and tools, including merging paired-end reads, parsing unique molecular identifiers (UMIs), and building/using custom germline databases. These are described in the `abstar` [documentation](http://abstar.readthedocs.org).  
  

## api  
Most core `abstar` functions are available through a Python API, making it easier to run `abstar` as a component of integrated analysis pipelines or to run `abstar` interactively (e.g. in a Jupyter notebook). See the `abstar` [documentation](http://abstar.readthedocs.org) for more detail about the API.  
  

## testing  
You can run the complete `abstar` test suite by first installing `pytest`:
``` bash
pip install pytest
```

followed by:

``` bash
git clone https://github.com/brineylab/abstar
cd abstar
pytest
```

This test suite is automatically run after every commit, and is tested against all supported versions of Python.

  
## requirements  
**Python 3.10+**   
  
abutils  
click  
matplotlib  
numpy  
pandas  
parasail  
polars  
pyarrow  
pytest  
