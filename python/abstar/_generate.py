"""`abstar.generate()`: synthetic V(D)J reads with truth. Models are in `abstar.model`."""
from __future__ import annotations

import secrets
import tomllib
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

from .errors import AbstarError, from_core

if TYPE_CHECKING:
    import polars as pl

_FORMATS = {"airr", "parquet"}


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


def _core():
    from . import _core

    return _core


def _species(db_path: Path) -> list[str]:
    try:
        return tomllib.loads((db_path / "manifest.toml").read_text())["species"]
    except (OSError, tomllib.TOMLDecodeError, KeyError) as e:
        raise AbstarError("DATABASE_INVALID", f"{db_path}: {e}") from None


def _prepare(n, locus, germline_db, model, preset, seed, threads, options) -> tuple[str, str, dict]:
    if isinstance(n, bool) or not isinstance(n, int) or n < 0:
        raise AbstarError("INVALID_ARGUMENT", f"n must be an integer >= 0, got {n!r}")
    if threads is not None and (isinstance(threads, bool) or not isinstance(threads, int) or threads < 1):
        raise AbstarError("INVALID_ARGUMENT", f"threads must be None or an integer >= 1, got {threads!r}")
    from . import germline, model as _model

    db_path = germline.resolve(germline_db)[0]
    species = ""
    if model is None:
        sp = _species(db_path)
        if len(sp) != 1:
            raise AbstarError(
                "INVALID_ARGUMENT", f"germline database {germline_db!r} has species {', '.join(sp)}; pass a model"
            )
        species = sp[0]
    model_path = _model.resolve(model, species, locus)[0]
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
) -> pl.DataFrame | GenerateSummary:
    """Simulate `n` V(D)J reads of `locus` with their truth (`abstar manual generate`).

    `model` defaults to the built-in model for the database's species and `locus`
    (a database with more than one species needs `model`). `preset` (naive,
    memory, bnab, tcr) sets the starting parameters and `options` override them:
    `mutation`, `long_cdr3_weight`, `indels`, `frameshift_indels`,
    `indel_duplication`, `truncate_5p`, `truncate_3p`, `constant` (nt of constant
    region, or None), `isotype_weights`, `productive_only`, `shm_table` (a
    directory of S5F CSVs) and `shm_model`. `mutation` takes a number or `fixed:X`,
    `uniform:LO-HI`, `beta:A,B`; the truncations (nt) a number, `fixed:N` or
    `uniform:LO-HI`. Without `seed` one is
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

