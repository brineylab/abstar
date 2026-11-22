import subprocess, sys
import abstar

def test_version_matches_core():
    from abstar import _core
    # Cargo's 1.0.0-rc.1 is Python's 1.0.0rc1: equal once separators are dropped
    strip = lambda v: v.replace(".", "").replace("-", "")
    assert strip(abstar.__version__) == strip(_core.version())

def test_import_is_light():
    code = "import sys, abstar; print(any(m in sys.modules for m in ('polars', 'typer', 'abstar._core')))"
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True)
    assert out.stdout.strip() == "False"
