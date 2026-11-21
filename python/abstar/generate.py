"""`abstar.generate()`: synthetic V(D)J reads with truth, and the recombination models."""
from __future__ import annotations

import builtins
import secrets
import sys
import tomllib
import types
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Literal

from .errors import AbstarError, from_core
from .germline import abstar_home

if TYPE_CHECKING:
    import polars as pl

_BUILTIN_DIR = Path(__file__).parent / "models"
_FORMATS = {"airr", "parquet"}

Kind = Literal["builtin", "user", "path"]


@dataclass(frozen=True)
class GenerateSummary:
    n: int
    productive: int
    attempts: int
    seed: int
    fasta: Path | None
    truth: dict[str, Path]
    generate_json: Path | None
    warnings: list[str]


@dataclass(frozen=True)
class ModelInfo:
    name: str
    kind: Kind
    path: Path
    species: str
    locus: str
    description: str
    sources: list[dict]


def _core():
    from . import _core

    return _core


def _names(d: Path) -> builtins.list[str]:
    if not d.is_dir():
        return []
    return sorted(p.name for p in d.iterdir() if (p / "model.toml").is_file())


def _user_dir() -> Path:
    return abstar_home() / "models"


def resolve_model(name_or_path: str | Path | None, species: str, locus: str) -> tuple[Path, Kind]:
    """The model file for `name_or_path` (a model.toml, a directory holding one, a built-in
    name, then a user name), or the built-in `f"{species}-{locus.lower()}"` when it is None."""
    if name_or_path is None:
        name_or_path = f"{species}-{locus.lower()}"
    p = Path(name_or_path)
    if p.is_file():
        return p, "path"
    if (p / "model.toml").is_file():
        return p / "model.toml", "path"
    name = str(name_or_path)
    builtins_, users = _names(_BUILTIN_DIR), _names(_user_dir())
    if name in builtins_:
        if name in users:
            raise AbstarError(
                "BUILTIN_NAME_CONFLICT",
                f"user model {_user_dir() / name} has the name of a built-in model; rename it",
            )
        return _BUILTIN_DIR / name / "model.toml", "builtin"
    if name in users:
        return _user_dir() / name / "model.toml", "user"
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


def model_info(name_or_path: str | Path) -> ModelInfo:
    return _info(*resolve_model(name_or_path, "", ""))


def list_models() -> builtins.list[ModelInfo]:
    """Built-in models, then user models (minus any that shadow a built-in name); each sorted by name."""
    builtins_ = _names(_BUILTIN_DIR)
    return [_info(_BUILTIN_DIR / n / "model.toml", "builtin") for n in builtins_] + [
        _info(_user_dir() / n / "model.toml", "user") for n in _names(_user_dir()) if n not in builtins_
    ]


def _species(db_path: Path) -> builtins.list[str]:
    try:
        return tomllib.loads((db_path / "manifest.toml").read_text())["species"]
    except (OSError, tomllib.TOMLDecodeError, KeyError) as e:
        raise AbstarError("DATABASE_INVALID", f"{db_path}: {e}") from None


def _prepare(n, locus, germline_db, model, preset, seed, threads, options) -> tuple[str, str, dict]:
    if isinstance(n, bool) or not isinstance(n, int) or n < 0:
        raise AbstarError("INVALID_ARGUMENT", f"n must be an integer >= 0, got {n!r}")
    if threads is not None and (isinstance(threads, bool) or not isinstance(threads, int) or threads < 1):
        raise AbstarError("INVALID_ARGUMENT", f"threads must be None or an integer >= 1, got {threads!r}")
    from . import germline

    db_path = germline.resolve(germline_db)[0]
    species = ""
    if model is None:
        sp = _species(db_path)
        if len(sp) != 1:
            raise AbstarError(
                "INVALID_ARGUMENT", f"germline database {germline_db!r} has species {', '.join(sp)}; pass a model"
            )
        species = sp[0]
    model_path = resolve_model(model, species, locus)[0]
    if seed is None:
        seed = secrets.randbits(63)
    params = {**options, "n": n, "locus": locus, "preset": preset, "seed": seed, "threads": threads}
    return str(db_path), str(model_path), params


def generate(
    n: int,
    locus: str,
    *,
    germline_db: str = "human",
    model: str | Path | None = None,
    preset: str | None = None,
    seed: int | None = None,
    output_dir: str | Path | None = None,
    name: str = "generated",
    formats: str | Iterable[str] = ("airr",),
    threads: int | None = None,
    **options,
) -> "pl.DataFrame | GenerateSummary":
    """Simulate `n` V(D)J reads of `locus` with their truth (`abstar manual generate`).

    `model` defaults to the built-in model for the database's species and `locus`
    (a database with more than one species needs `model`). `preset` (naive,
    memory, bnab, tcr) sets the starting parameters and `options` override them:
    `mutation`, `long_cdr3_weight`, `indels`, `frameshift_indels`,
    `indel_duplication`, `truncate_5p`, `truncate_3p`, `constant` (nt of constant
    region, or None), `isotype_weights`, `productive_only`, `shm_table` (a
    directory of S5F CSVs) and `shm_model`. `mutation` and the truncations take a
    number or `fixed:X`, `uniform:LO-HI`, `beta:A,B`. Without `seed` one is
    drawn (and recorded in generate.json).

    With `output_dir=None` the truth table is returned as a Polars DataFrame.
    Otherwise `<name>.fasta`, `<name>.truth.tsv` / `.parquet` (`formats`) and
    `generate.json` are written and a `GenerateSummary` is returned.
    """
    fmts = [formats] if isinstance(formats, str) else list(formats)
    if not set(fmts) <= _FORMATS:
        raise AbstarError("INVALID_ARGUMENT", f"formats must be a subset of {sorted(_FORMATS)}, got {fmts!r}")
    db_path, model_path, params = _prepare(n, locus, germline_db, model, preset, seed, threads, options)
    core = _core()
    try:
        if output_dir is None:
            import polars as pl

            return pl.DataFrame(core.generate_in_memory(db_path, model_path, params))
        d = core.generate_to_files(db_path, model_path, params, str(output_dir), name, fmts)
    except core.CoreError as e:
        raise from_core(e) from None
    return _summary(d)


def _summary(d: dict) -> GenerateSummary:
    path = lambda p: None if p is None else Path(p)  # noqa: E731
    return GenerateSummary(
        n=d["n"],
        productive=d["productive"],
        attempts=d["attempts"],
        seed=d["seed"],
        fasta=path(d["fasta"]),
        truth={k: Path(v) for k, v in d["truth"].items()},
        generate_json=path(d["generate_json"]),
        warnings=d["warnings"],
    )


class _CallableModule(types.ModuleType):
    def __call__(self, *args, **kwargs):
        return generate(*args, **kwargs)


# Importing this submodule sets `abstar.generate` to the module, shadowing the function;
# a callable module keeps both `abstar.generate(...)` and `abstar.generate.list_models()`.
sys.modules[__name__].__class__ = _CallableModule
