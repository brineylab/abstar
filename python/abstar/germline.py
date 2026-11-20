"""Germline databases: resolve, build, list, and inspect."""
from __future__ import annotations

import builtins
import os
import tomllib
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from .errors import AbstarError, from_core

_BUILTIN_DIR = Path(__file__).parent / "germline_dbs"

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
    return sorted(p.name for p in d.iterdir() if (p / "manifest.toml").is_file())


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
    output_dir.parent.mkdir(parents=True, exist_ok=True)
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
    )


def info(name_or_path: str | Path) -> DatabaseInfo:
    return _info(*resolve(name_or_path))


def list() -> builtins.list[DatabaseInfo]:
    """Built-in databases, then user databases; each sorted by name."""
    return [_info(_BUILTIN_DIR / n, "builtin") for n in builtin_names()] + [
        _info(_user_dir() / n, "user") for n in _names(_user_dir())
    ]
