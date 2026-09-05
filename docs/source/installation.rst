Installation
============

Requirements
------------

abstar supports Python 3.10 through 3.13. The release matrix tests both the
lowest declared dependency profile and the newest versions allowed by these
bounds:

- ``abutils>=0.6,<0.7``
- ``polars>=1.5,<2``
- ``pyarrow>=16.1,<26``
- ``biopython>=1.80``

MMseqs2_ is used internally for germline gene assignment, but it is bundled with
abutils_ (a dependency of abstar), so no separate installation is required.


Install abstar
--------------

The easiest way to install abstar is via pip:

.. code-block:: bash

    pip install abstar


External executables
--------------------

The normal installation obtains the MMseqs2 and fastp executables through
``abutils``. abstar resolves MMseqs2 with the public ``abutils.bin.get_path``
accessor and executes it with checked argument lists. A failed external command
raises a structured ``AnnotationRunError`` that includes its stage and category;
retained diagnostics include the exit status, stdout, and stderr.


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

    git clone https://github.com/brineylab/abstar
    cd abstar/
    python -m pip install -r requirements-test.txt


Development verification
------------------------

The marker scopes are intentionally separate. The bulk published corpus is not
needed for any of these commands:

.. code-block:: bash

    # Fast unit and property tests
    python -m pytest -m "not integration and not e2e and not slow" -q

    # Real component boundaries, public entry points, and reserved slow checks
    python -m pytest -m "integration" -q
    python -m pytest -m "e2e" -q
    python -m pytest -m "slow" -q

    # Complete suite
    python -m pytest -q

    # AIRR and packaged-database gates
    python -m pytest abstar/tests/test_airr.py -q
    python -m pytest abstar/tests/test_database_integrity.py -q

    # Statement/branch package floor and critical-module floors
    python -m pytest --cov=abstar --cov-branch --cov-report=term-missing \
      --cov-report=json:/tmp/abstar-coverage.json -q
    python scripts/check_coverage.py /tmp/abstar-coverage.json coverage-floors.json

    # Documentation warnings are errors
    python -m sphinx -W --keep-going -b html docs/source docs/_build/html

Optional corpus discovery requires explicit read-only source paths and an output
outside the source, input, and repository trees. It never runs in ordinary push
or pull-request CI:

.. code-block:: bash

    python scripts/discover_bcr_cases.py \
      --fasta-dir /path/to/bcr_fastas \
      --manifest /path/to/sample_manifest.csv \
      --cellranger-root /path/to/cellranger \
      --output /tmp/abstar-bcr-candidates.jsonl \
      --per-dataset 25 --n-processes 2

The scheduled nightly workflow is separate. It downloads an explicitly
provisioned artifact with ``bcr_fastas/``, ``sample_manifest.csv``, and
``cellranger/``, limits the per-dataset cohort, and uploads its outcome report.


Verify Installation
-------------------

To verify that abstar is installed correctly:

.. code-block:: bash

    python -c "import abstar; print(abstar.__version__)"


.. _MMseqs2: https://github.com/soedinglab/MMseqs2
.. _abutils: https://github.com/briney/abutils
.. _datascience: https://hub.docker.com/repository/docker/brineylab/datascience/general
