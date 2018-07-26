Install
=======

The easiest way to install abstar locally (on OSX or Linux) is to use pip::

    $ pip install abstar

If you don't have pip, the Anaconda_ Python distribution contains pip along 
with a ton of useful scientific Python packages and is a great way to get 
started with Python.

abstar does not run natively on Windows, but Windows users can run abstar with Docker_::

    $ docker pull briney/abstar
    $ docker run -it briney/abstar

Stable_ and development_ versions of abstar can also be downloaded from Github. 
You can manually install the latest development version of abstar with::

    $ git clone https://github.com/briney/abstar
    $ cd abstar/
    $ python setup.py install

.. note::

    If installing manually via setup.py and you don't already have scikit-bio installed, 
    you may get an error when setuptools attempts to install scikit-bio. This can be fixed 
    by first installing scikit-bio with pip::

        $ pip install scikit-bio

    and then retrying the manual install of abstar. Starting with version 0.5, scikit-bio 
    dropped support for Python 2.7, so install scikit-bio on Python 2.7 with::

        $ pip install scikit-bio<=0.4.2


Requirements
------------

* Python 2.7 or 3.5+
* abutils_
* biopython_
* celery_
* pymongo_
* pytest_
* `scikit bio`_


Optional dependencies
---------------------

Several optional abstar components have additional dependencies:

* ``abstar.preprocessing`` requires FASTQC_, cutadapt_ and sickle_
* sequence merging requires PANDAseq_
* downloading data from BaseSpace requires the `BaseSpace Python SDK`_
* ``batch_mongoimport`` requires MongoDB_

If using Docker, all of the the optional dependencies are included.


.. _Docker: https://www.docker.com/
.. _Anaconda: https://www.continuum.io/downloads
.. _stable: https://github.com/briney/abstar/releases
.. _development: https://github.com/briney/abstar
.. _abutils: https://github.com/briney/abutils
.. _biopython: http://biopython.org/
.. _scikit bio: http://scikit-bio.org/
.. _pymongo: https://api.mongodb.org/python/current/
.. _celery: http://www.celeryproject.org/
.. _PANDAseq: https://github.com/neufeld/pandaseq
.. _FASTQC: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
.. _cutadapt: https://github.com/marcelm/cutadapt/
.. _sickle: https://github.com/najoshi/sickle
.. _BaseSpace Python SDK: https://github.com/basespace/basespace-python-sdk
.. _MongoDB: https://www.mongodb.org/
.. _pytest: https://docs.pytest.org/en/latest/
