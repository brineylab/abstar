import tomllib
from pathlib import Path

CRATES = Path(__file__).resolve().parents[2] / "crates"
TABLES = ("dependencies", "build-dependencies", "dev-dependencies")


def _abstar_deps(crate):
    """Real package names of every abstar_* dependency, including target-specific tables and renames."""
    d = tomllib.loads((CRATES / crate / "Cargo.toml").read_text())
    tables = [d.get(t, {}) for t in TABLES]
    tables += [tgt.get(t, {}) for tgt in d.get("target", {}).values() for t in TABLES]
    names = set()
    for t in tables:
        for key, spec in t.items():
            names.add(spec.get("package", key) if isinstance(spec, dict) else key)
    return {n for n in names if n.startswith("abstar_")}


def test_generate_never_depends_on_annotate():
    if (CRATES / "generate").exists():
        extra = _abstar_deps("generate") - {"abstar_germline", "abstar_table"}
        assert not extra, extra
    assert not _abstar_deps("table")


def _tables(crate):
    """Every dependency table of `crate` (plain and target-specific) as (table name, table)."""
    d = tomllib.loads((CRATES / crate / "Cargo.toml").read_text())
    out = [(t, d.get(t, {})) for t in TABLES]
    out += [(t, tgt.get(t, {})) for tgt in d.get("target", {}).values() for t in TABLES]
    return out


def test_annotate_uses_generate_only_in_tests():
    """abstar_generate is a dev-dependency of abstar_annotate only, never a library dependency."""
    found = set()
    for table, deps in _tables("annotate"):
        for key, spec in deps.items():
            if (spec.get("package", key) if isinstance(spec, dict) else key) == "abstar_generate":
                found.add(table)
    assert found == {"dev-dependencies"}, found
    assert "abstar_annotate" not in _abstar_deps("generate")
    used = [p for p in (CRATES / "annotate" / "src").rglob("*.rs") if "abstar_generate" in p.read_text()]
    assert not used, used
