"""Recombination models for `abstar.generate()`: resolve, list, and inspect."""
from __future__ import annotations

import builtins
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

from .errors import AbstarError, from_core
from .germline import abstar_home, _is_path

_BUILTIN_DIR = Path(__file__).parent / "models"

Kind = Literal["builtin", "user", "path"]


@dataclass(frozen=True)
class ModelInfo:
    name: str
    kind: Kind
    path: Path
    species: str
    locus: str
    description: str
    sources: builtins.list[dict] = field(default_factory=builtins.list)
    error: str | None = None  # set by list() for a model that failed to load


def _core():
    from . import _core

    return _core


def _names(d: Path) -> builtins.list[str]:
    if not d.is_dir():
        return []
    return sorted(p.name for p in d.iterdir() if (p / "model.toml").is_file())


def _user_dir() -> Path:
    return abstar_home() / "models"


def resolve(name_or_path: str | Path | None, species: str, locus: str) -> tuple[Path, Kind]:
    """The model file for `name_or_path` (a built-in name, a user name, then a model.toml or a
    directory holding one), or the built-in `f"{species}-{locus.lower()}"` when it is None. A path
    (`germline._is_path`) skips the names, so a local `human-igh/` never replaces the built-in."""
    if name_or_path is None:
        name_or_path = f"{species}-{locus.lower()}"
    name = str(name_or_path)
    builtins_, users = _names(_BUILTIN_DIR), _names(_user_dir())
    if not _is_path(name_or_path):
        if name in builtins_:
            if name in users:
                raise AbstarError(
                    "BUILTIN_NAME_CONFLICT",
                    f"user model {_user_dir() / name} has the name of a built-in model; rename it",
                )
            return _BUILTIN_DIR / name / "model.toml", "builtin"
        if name in users:
            return _user_dir() / name / "model.toml", "user"
    p = Path(name_or_path)
    if p.is_file():
        return p, "path"
    if (p / "model.toml").is_file():
        return p / "model.toml", "path"
    raise AbstarError(
        "MODEL_NOT_FOUND", f"no recombination model {name!r} (available: {', '.join(builtins_ + users) or 'none'})"
    )


def _info(path: Path, kind: Kind) -> ModelInfo:
    core = _core()
    try:
        s = core.model_summary(str(path))
    except core.CoreError as e:
        raise from_core(e) from None
    if kind != "path":
        s["name"] = path.parent.name  # the name it is resolved by, whatever the file says
    return ModelInfo(kind=kind, path=path, **s)


def _listed(path: Path, kind: Kind) -> ModelInfo:
    """`_info`, but a model that fails to load is listed with its error instead of raising."""
    try:
        return _info(path, kind)
    except AbstarError as e:
        return ModelInfo(path.parent.name, kind, path, "", "", "", error=str(e))


def info(name_or_path: str | Path) -> ModelInfo:
    return _info(*resolve(name_or_path, "", ""))


def list() -> builtins.list[ModelInfo]:
    """Built-in models, then user models (minus any that shadow a built-in name); each sorted
    by name. A model that fails to load is still listed, with the reason in `error`."""
    builtins_ = _names(_BUILTIN_DIR)
    return [_listed(_BUILTIN_DIR / n / "model.toml", "builtin") for n in builtins_] + [
        _listed(_user_dir() / n / "model.toml", "user") for n in _names(_user_dir()) if n not in builtins_
    ]
