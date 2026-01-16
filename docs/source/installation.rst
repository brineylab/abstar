Installation
============

Requirements
------------

abstar requires:

- Python 3.9 or later
- MMseqs2_ (external binary for germline gene assignment)

Install abstar
--------------

The easiest way to install abstar is via pip:

.. code-block:: bash

    pip install abstar


Install MMseqs2
---------------

MMseqs2 must be installed separately and available in your PATH.

**Using conda (recommended):**

.. code-block:: bash

    conda install -c conda-forge -c bioconda mmseqs2

**On Ubuntu/Debian:**

.. code-block:: bash

    sudo apt install mmseqs2

**On macOS with Homebrew:**

.. code-block:: bash

    brew install mmseqs2

**From source:**

See the `MMseqs2 installation guide`_ for additional options.


Docker
------

abstar is included in the brineylab datascience_ Docker container:

.. code-block:: bash

    docker pull brineylab/datascience
    docker run -it brineylab/datascience

This container includes MMseqs2 and all other dependencies pre-configured.


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

    # Check abstar version
    python -c "import abstar; print(abstar.__version__)"

    # Check MMseqs2 is available
    mmseqs --help


.. _MMseqs2: https://github.com/soedinglab/MMseqs2
.. _MMseqs2 installation guide: https://github.com/soedinglab/MMseqs2#installation
.. _datascience: https://hub.docker.com/repository/docker/brineylab/datascience/general
