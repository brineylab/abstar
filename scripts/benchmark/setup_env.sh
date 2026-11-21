#!/usr/bin/env bash
# Benchmark-only env: IgBLAST, MMseqs2 and legacy abstar (installed from legacy/, read-only).
# Not a CI or package dependency. Usage: scripts/benchmark/setup_env.sh
set -euo pipefail
MM=${MICROMAMBA:-$HOME/.local/bin/micromamba}
ENV=${BENCH_ENV:-abstar-bench}
ROOT=$(cd "$(dirname "$0")/../.." && pwd)

"$MM" create -y -n "$ENV" -c conda-forge -c bioconda python=3.12 igblast blast mmseqs2
# pip builds in-tree, so install from a temporary copy: legacy/ is never written to.
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
cp -r "$ROOT/legacy" "$TMP/legacy"
"$MM" run -n "$ENV" python -m pip install --no-cache-dir -r "$TMP/legacy/requirements.txt"
"$MM" run -n "$ENV" python -m pip install --no-cache-dir --no-deps "$TMP/legacy"
"$MM" run -n "$ENV" sh -c 'igblastn -version; mmseqs version; python -c "import abstar; print(\"legacy abstar\", abstar.__version__)"'
