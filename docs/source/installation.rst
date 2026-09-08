Installation
============

Requirements
------------

abstar supports Python 3.10 through 3.13. The release matrix tests both the
lowest declared dependency profile and the newest versions allowed by these
bounds:

- ``abutils>=0.6,<0.7``
- ``polars>=1.6,<2``
- ``pyarrow>=16.1,<26``
- ``biopython>=1.80``

MMseqs2_ is used internally for germline gene assignment, but it is bundled with
abutils_ (a dependency of abstar), so no separate installation is required.


Install abstar
--------------

The easiest way to install abstar is via pip:

.. code-block:: bash

    pip install abstar


macOS source-build prerequisites
--------------------------------

On Apple Silicon, ``parasail==1.3.4`` is installed from source because it does
not provide a compatible wheel. Install the Xcode Command Line Tools
(``xcode-select --install``) and Homebrew_, then prepare its build tools before
installing abstar:

.. code-block:: bash

    brew install autoconf automake libtool m4
    export M4="$(brew --prefix m4)/bin/m4"
    python -m pip install abstar

Selecting ``M4`` explicitly keeps the build on Homebrew's GNU M4 even when
Parasail prepends ``/usr/bin`` to ``PATH``. The macOS wheel-install CI job
uses these prerequisites and checks a native Parasail alignment after installation.


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

    # Install documentation dependencies, then treat warnings as errors
    python -m pip install -r docs/doc_requirements.txt
    python -m sphinx -W --keep-going -b html docs/source docs/_build/html

The dedicated committed-corpus workflow runs a fixed, enriched BCR subset on
pushes and pull requests. It runs once outside the pytest version matrix and
coverage job. Its input and baseline are committed under
``test_data/bcr_corpus/`` and are excluded from installed distributions.
Reproduce it with Python 3.12 and a new output directory outside the checkout:

.. code-block:: bash

    python -m pip install -r requirements-corpus.txt
    POLARS_MAX_THREADS=2 OMP_NUM_THREADS=2 python scripts/run_corpus.py \
      --output /tmp/abstar-corpus-check

The output includes a comparison report and persistent diagnostics. Unexpected
biological changes and new record failures fail the job. Follow
``abstar/tests/README.md`` and ``test_data/bcr_corpus/README.md`` before changing
baseline expectations. The five-minute runtime target requires measurement on
the hosted runner; local benchmarks alone do not establish it.

Optional full-corpus discovery requires explicit read-only source paths and an output
outside the source, input, and repository trees. It never runs in ordinary push
or pull-request CI:

.. code-block:: bash

    python scripts/discover_bcr_cases.py \
      --fasta-dir /path/to/bcr_fastas \
      --manifest /path/to/sample_manifest.csv \
      --cellranger-root /path/to/cellranger \
      --output /tmp/abstar-bcr-candidates.jsonl \
      --per-dataset 25 --n-processes 2

Full external corpus discovery is a manual workflow; there is no scheduled corpus
job. The scheduled documentation linkcheck runs independently.


Verify Installation
-------------------

To verify that abstar is installed correctly:

.. code-block:: bash

    python -c "import abstar; print(abstar.__version__)"


.. _MMseqs2: https://github.com/soedinglab/MMseqs2
.. _Homebrew: https://brew.sh/
.. _abutils: https://github.com/briney/abutils
.. _datascience: https://hub.docker.com/repository/docker/brineylab/datascience/general
