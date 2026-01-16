Installation
============

Requirements
------------

abstar requires Python 3.9 or later.

MMseqs2_ is used internally for germline gene assignment, but it is bundled with
abutils_ (a dependency of abstar), so no separate installation is required.


Install abstar
--------------

The easiest way to install abstar is via pip:

.. code-block:: bash

    pip install abstar


Using a Custom MMseqs2 Binary
-----------------------------

If you prefer to use a specific version of MMseqs2 rather than the bundled binary,
you can specify the path to your custom binary using the ``mmseqs_binary`` parameter:

.. code-block:: python

    import abstar

    sequences = abstar.run("sequences.fasta", mmseqs_binary="/path/to/mmseqs")

See the `MMseqs2 installation guide`_ for information on installing MMseqs2.


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
.. _MMseqs2 installation guide: https://github.com/soedinglab/MMseqs2#installation
.. _datascience: https://hub.docker.com/repository/docker/brineylab/datascience/general
