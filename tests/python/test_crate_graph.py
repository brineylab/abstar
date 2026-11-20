import tomllib
from pathlib import Path

CRATES = Path(__file__).resolve().parents[2] / "crates"


def _deps(crate):
    d = tomllib.loads((CRATES / crate / "Cargo.toml").read_text())
    return set(d.get("dependencies", {})) | set(d.get("build-dependencies", {}))


def test_generate_never_depends_on_annotate():
    crates = ["table"] + (["generate"] if (CRATES / "generate").exists() else [])
    for c in crates:
        assert "abstar_annotate" not in _deps(c), c
    assert not [d for d in _deps("table") if d.startswith("abstar_")]
