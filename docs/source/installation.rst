install
=======

The easiest way to install ``abstar`` locally (on macOS or Linux) is to use pip::

    $ pip install abstar

``abstar`` does not run natively on Windows, but Windows users can use Docker_ 
(the brineylab datascience_ Docker container contains the entire ab[x] toolkit,
which includes ``abstar``)::

    $ docker pull brineylab/datascience
    $ docker run -it brineylab/datascience

Stable_ and development_ versions of ``abstar`` can also be downloaded from Github. 
You can manually install the latest development version of ``abstar`` with::

    $ git clone https://github.com/briney/abstar
    $ cd abstar/
    $ python setup.py install


.. _Docker: https://www.docker.com/
.. _Anaconda: https://www.continuum.io/downloads
.. _stable: https://github.com/briney/abstar/releases
.. _development: https://github.com/briney/abstar
