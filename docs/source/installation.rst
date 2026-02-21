Installation
============

Requirements
------------

abstar requires Python 3.10 or later.

MMseqs2_ is used internally for germline gene assignment, but it is bundled with
abutils_ (a dependency of abstar), so no separate installation is required.


Install abstar
--------------

The easiest way to install abstar is via pip:

.. code-block:: bash

    pip install abstar


Docker
------

abstar is included in the brineylab datascience_ Docker container:

.. code-block:: bash

    docker pull brineylab/datascience
    docker run -it brineylab/datascience

This container includes abstar and all dependencies pre-configured.


Development Installation
------------------------

To install from source for development:

.. code-block:: bash

    git clone https://github.com/briney/abstar
    cd abstar/
    pip install -e .


Verify Installation
-------------------

To verify that abstar is installed correctly:

.. code-block:: bash

    python -c "import abstar; print(abstar.__version__)"


.. _MMseqs2: https://github.com/soedinglab/MMseqs2
.. _abutils: https://github.com/briney/abutils
.. _datascience: https://hub.docker.com/repository/docker/brineylab/datascience/general
