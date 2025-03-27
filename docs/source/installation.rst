install
=======

The easiest way to install ``abstar`` locally (on macOS or Linux) is to use pip:

.. code-block:: console

    pip install abstar

``abstar`` does not run natively on Windows, but Windows users can use Docker_ 
(the brineylab datascience_ Docker container contains the entire ab[x] toolkit,
which includes ``abstar``):

.. code-block:: console

    docker pull brineylab/datascience
    docker run -it brineylab/datascience

Stable and development versions of ``abstar`` can also be downloaded from GitHub_. 
You can manually install the latest version of ``abstar`` with:

.. code-block:: console

    git clone https://github.com/briney/abstar
    cd abstar/
    python setup.py install


.. _Docker: https://www.docker.com/
.. _datascience: https://hub.docker.com/repository/docker/brineylab/datascience/general
.. _GitHub: https://github.com/briney/abstar
