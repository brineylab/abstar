"""Germline databases: resolve, build, list, and inspect."""
from __future__ import annotations

import builtins
import gzip
import os
import re
import tomllib
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Literal

from .errors import AbstarError, from_core

_BUILTIN_DIR = Path(__file__).parent / "germline_dbs"
# Leftovers of an interrupted build (crates/germline io.rs `sibling`); never databases.
_LEFTOVER = re.compile(r"\.(staging|old)-\d+$")

if TYPE_CHECKING:
    import polars as pl

Kind = Literal["builtin", "user", "path"]


@dataclass(frozen=True)
class DatabaseInfo:
    name: str
    kind: Kind
    path: Path
    species: list[str]
    loci: list[str]
    n_alleles: int
    counts: dict[str, int]
    description: str
    built: str = ""
    abstar_version: str = ""
    alleles_sha256: str = ""
    prefixes: dict[str, str] = field(default_factory=dict)
    provenance: builtins.list[dict] = field(default_factory=builtins.list)
    warnings: builtins.list[dict] = field(default_factory=builtins.list)
    error: str | None = None  # set by list() for a database that failed to load


@dataclass(frozen=True)
class BuildReport:
    name: str
    path: Path
    n_alleles: int
    counts: dict[str, int]
    restricted: list[str]
    passed_through: list[str]
    warnings: list[dict]


def _core():
    from . import _core

    return _core


def abstar_home() -> Path:
    """`$ABSTAR_HOME`, or `~/.abstar`."""
    return Path(os.environ.get("ABSTAR_HOME") or Path.home() / ".abstar")


def _names(d: Path) -> builtins.list[str]:
    if not d.is_dir():
        return []
    return sorted(p.name for p in d.iterdir() if (p / "manifest.toml").is_file() and not _LEFTOVER.search(p.name))


def _user_dir() -> Path:
    return abstar_home() / "germline"


def builtin_names() -> builtins.list[str]:
    """Names of the packaged databases."""
    return _names(_BUILTIN_DIR)


def resolve(name_or_path: str | Path) -> tuple[Path, Kind]:
    """A database directory, then a built-in name, then a user name."""
    p = Path(name_or_path)
    if (p / "manifest.toml").is_file():
        return p, "path"
    name = str(name_or_path)
    for kind, d in (("builtin", _BUILTIN_DIR), ("user", _user_dir())):
        if name in _names(d):
            return d / name, kind
    available = ", ".join(builtin_names() + _names(_user_dir())) or "none"
    raise AbstarError("DATABASE_NOT_FOUND", f"no germline database {name!r} (available: {available})")


def build(
    recipe: str | Path, *, output_dir: str | Path | None = None, overwrite: bool = False
) -> BuildReport:
    """Build a database from a recipe; by default into `abstar_home()/germline/<name>`."""
    recipe = Path(recipe)
    try:
        spec = tomllib.loads(recipe.read_text())
    except (OSError, UnicodeDecodeError, tomllib.TOMLDecodeError) as e:
        raise AbstarError("RECIPE_INVALID", f"{recipe}: {e}") from None
    name = spec.get("name")
    if not isinstance(name, str) or not name:
        raise AbstarError("RECIPE_INVALID", f"{recipe}: missing name")
    if output_dir is None:
        if name in builtin_names():
            raise AbstarError(
                "BUILTIN_NAME_CONFLICT",
                f"{name!r} is a built-in database name; rename the recipe or pass output_dir",
            )
        output_dir = _user_dir() / name
    output_dir = Path(output_dir)
    sources = spec.get("source", [])
    databases = [s["database"] for s in sources if isinstance(s, dict) and isinstance(s.get("database"), str)]
    paths = {db: str(resolve(db)[0]) for db in databases}
    core = _core()
    try:
        r = core.germline_build(str(recipe), str(output_dir), overwrite, paths)
    except core.CoreError as e:
        raise from_core(e) from None
    return BuildReport(**{**r, "path": Path(r["path"])})


def _info(path: Path, kind: Kind) -> DatabaseInfo:
    core = _core()
    try:
        s = core.germline_summary(str(path))
    except core.CoreError as e:
        raise from_core(e) from None
    return DatabaseInfo(
        name=s["name"],
        kind=kind,
        path=path,
        species=s["species"],
        loci=s["loci"],
        n_alleles=s["n_alleles"],
        counts=s["counts"],
        description=s["description"],
        built=s["built"],
        abstar_version=s["abstar_version"],
        alleles_sha256=s["alleles_sha256"],
        prefixes=s.get("prefixes", {}),
        provenance=s["provenance"],
        warnings=s["warnings"],
    )


def _listed(path: Path, kind: Kind) -> DatabaseInfo:
    """`_info`, but a database that fails to load is listed with its error instead of raising."""
    try:
        return _info(path, kind)
    except AbstarError as e:
        return DatabaseInfo(path.name, kind, path, [], [], 0, {}, "", error=str(e))


def info(name_or_path: str | Path) -> DatabaseInfo:
    return _info(*resolve(name_or_path))


def list() -> builtins.list[DatabaseInfo]:
    """Built-in databases, then user databases; each sorted by name. A database that
    fails to load is still listed, with `n_alleles` 0 and the reason in `error`."""
    return [_listed(_BUILTIN_DIR / n, "builtin") for n in builtin_names()] + [
        _listed(_user_dir() / n, "user") for n in _names(_user_dir())
    ]


def _read_fasta(path: Path) -> builtins.list[tuple[str, str]]:
    try:
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt") as f:
            lines = f.read().splitlines()
    except (OSError, UnicodeDecodeError) as e:
        raise AbstarError("INPUT_UNREADABLE", f"{path}: {e}") from None
    out: builtins.list[tuple[str, str]] = []
    for line in lines:
        if line.startswith(">"):
            out.append((line[1:], ""))
        elif out:
            out[-1] = (out[-1][0], out[-1][1] + line.strip())
    return out


def _allele_name(header: str) -> str:
    """The allele name the way crates/germline `sources::fasta` takes it."""
    if "|" in header:
        f = header.split("|")
        return f[1].strip() if len(f) > 1 else ""
    return header.split()[0] if header.split() else ""


def gap(
    sequences: "str | Path | builtins.list[tuple[str, str]]", templates: builtins.list[str] | None = None
) -> "tuple[builtins.list[tuple[str, str]], pl.DataFrame]":
    """IMGT-gap V sequences by transferring positions from the best template allele.

    Returns the gapped `(header, sequence)` records (full input headers, input order) and a
    per-record report. Raises `NO_ALLELES_GAPPED` when nothing could be gapped."""
    import polars as pl

    core = _core()
    records = _read_fasta(Path(sequences)) if isinstance(sequences, (str, Path)) else builtins.list(sequences)
    queries = [(_allele_name(h), s) for h, s in records]
    paths = [str(resolve(t)[0]) for t in (templates if templates is not None else core.GAP_DEFAULT_TEMPLATES)]
    try:
        results = core.germline_gap(queries, paths)
    except core.CoreError as e:
        raise from_core(e) from None
    out = [(h, r["gapped"]) for (h, _), r in zip(records, results) if r["gapped"] is not None]
    if not out:
        raise AbstarError("NO_ALLELES_GAPPED", f"none of {len(records)} sequences could be gapped")
    report = pl.DataFrame(
        {
            "name": [r["name"] for r in results],
            "status": ["failed" if r["failure"] else "gapped" for r in results],
            "template": [r["template"] for r in results],
            "template_species": [r["template_species"] for r in results],
            "identity": [r["identity"] for r in results],
            "cdr1_length": [r["cdr1_length"] for r in results],
            "cdr2_length": [r["cdr2_length"] for r in results],
            "warnings": [";".join(r["warnings"]) for r in results],
            "failure_reason": [r["failure"] for r in results],
        },
        schema={
            "name": pl.String, "status": pl.String, "template": pl.String, "template_species": pl.String,
            "identity": pl.Float64, "cdr1_length": pl.Int64, "cdr2_length": pl.Int64,
            "warnings": pl.String, "failure_reason": pl.String,
        },
    )
    return out, report
