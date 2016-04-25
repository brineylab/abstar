Install
=======

The easiest way to install AbStar locally (on OSX or Linux) is to use pip::

    $ pip install abstar

If you don't have pip, the Anaconda_ Python distribution contains pip along 
with a ton of useful scientific Python packages and is a great way to get 
started with Python.

AbStar does not run natively on Windows, but Windows users can run AbStar with Docker_::

    $ docker pull briney/abstar
    $ docker run -it briney/abstar

Stable_ and development_ versions of AbStar can also be downloaded from Github. 
You can manually install the latest development version of AbStar with::

    $ git clone https://github.com/briney/abstar
    $ cd abstar/
    $ python setup.py install

.. note:

    If installing manually via setup.py and you don't already have scikit-bio installed, 
    you may get an error when setuptools attempts to install scikit-bio. This can be fixed 
    by first installing scikit-bio with pip::

        $ pip install scikit-bio

    and then retrying the manual install of AbStar.


Requirements
------------

* Python 2.7.x (Python 3 compatability is in the works)
* abtools_
* biopython_
* `scikit bio`_
* pymongo_
* celery_


Optional dependencies
---------------------

Several optional AbStar components have additional dependencies:

* ``abstar.preprocessing`` requires FASTQC_, cutadapt_ and sickle_
* sequence merging requires PANDAseq_
* downloading data from BaseSpace requires the `BaseSpace Python SDK`_
* ``batch_mongoimport`` requires MongoDB_

If using Docker, all of the the optional dependencies are included.


.. _Docker: https://www.docker.com/
.. _Anaconda: https://www.continuum.io/downloads
.. _stable: https://github.com/briney/abstar/releases
.. _development: https://github.com/briney/abstar
.. _abtools: https://github.com/briney/abtools
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
