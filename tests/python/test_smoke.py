import subprocess, sys
import abstar

def test_version_matches_core():
    from abstar import _core
    assert abstar.__version__ == "1.0.0.dev0"
    assert _core.version() == "1.0.0-dev.0"

def test_import_is_light():
    code = "import sys, abstar; print(any(m in sys.modules for m in ('polars', 'typer', 'abstar._core')))"
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True)
    assert out.stdout.strip() == "False"
