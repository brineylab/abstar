"""`abstar.run()`: the single public entry point."""
from __future__ import annotations

import errno
import os
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

from .errors import AbstarError, from_core

if TYPE_CHECKING:
    import polars as pl

_EXTS = (".fasta", ".fa", ".fna", ".fastq", ".fq")
_NT = frozenset("ACGTNRYSWKMBDHVacgtnryswkmbdhv")
_AA = frozenset("ACDEFGHIKLMNPQRSTVWYXBZUOJ*acdefghiklmnpqrstvwyxbzuoj")
# Errors that mean "this string cannot be a path" (a long sequence; an invalid name on Windows).
_NOT_A_PATH = (errno.ENAMETOOLONG, errno.EINVAL)
_FORMATS = {"airr", "parquet"}
_STATUSES = ("annotated", "unassigned", "invalid_input", "error")


@dataclass(frozen=True)
class InputSummary:
    source: str
    outputs: dict[str, Path]
    counts: dict[str, int]


@dataclass(frozen=True)
class RunSummary:
    inputs: list[InputSummary]
    counts: dict[str, int]
    run_json: Path


def _is_seq_file(name: str) -> bool:
    """A FASTA/FASTQ name, optionally compressed. Only gzip is readable; the others are
    included so the core reports them (INPUT_UNREADABLE) instead of being skipped."""
    name = name.lower()
    for z in (".gz", ".bz2", ".xz", ".zst"):
        name = name.removesuffix(z)
    return name.endswith(_EXTS)


def _exists(p: Path) -> bool:
    try:
        os.stat(p)
        return True
    except (FileNotFoundError, NotADirectoryError, ValueError):  # ValueError: embedded NUL
        return False
    except OSError as e:
        if e.errno in _NOT_A_PATH:
            return False
        raise AbstarError("INPUT_UNREADABLE", f"{p}: {e.strerror or e}") from None


def _classify(item, aa: bool = False) -> tuple[str, object]:
    """('path', str) | ('seq', str) | ('rec', (id, seq))."""
    if isinstance(item, tuple):
        if len(item) == 2 and all(isinstance(x, str) for x in item):
            return "rec", item
        raise AbstarError("INVALID_ARGUMENT", f"records must be (id, sequence) string pairs, got {item!r}")
    if item == "-":  # stdin (the CLI's `-`); never a valid sequence
        return "path", Path("-")
    if isinstance(item, Path):
        if not _exists(item):
            raise AbstarError("INPUT_NOT_FOUND", f"{item} does not exist")
        return "path", item
    if not isinstance(item, str):
        raise AbstarError("INVALID_ARGUMENT", f"unsupported input of type {type(item).__name__}")
    if item and _exists(Path(item)):
        return "path", Path(item)
    if not (_AA if aa else _NT).issuperset(item):
        hint = " (protein? use sequence_type='aa')" if not aa and _AA.issuperset(item) else ""
        raise AbstarError("INPUT_NOT_FOUND", f"{item}: not an existing file and not {'an amino-acid' if aa else 'a nucleotide'} sequence{hint}")
    return "seq", item


def _expand(p: Path) -> list[str]:
    if p.is_dir():
        files = sorted(
            (f for f in p.rglob("*") if f.is_file() and _is_seq_file(f.name)),
            key=lambda f: f.relative_to(p).parts,
        )
        return [str(f) for f in files]
    return [str(p)]


def _normalize(sequences, aa: bool = False) -> tuple[list[str], list[tuple[str, str]]]:
    items = [sequences] if isinstance(sequences, (str, Path)) else list(sequences)  # consumed once
    kinds = [_classify(i, aa) for i in items]
    paths = [v for k, v in kinds if k == "path"]
    if paths and len(paths) != len(kinds):
        raise AbstarError("INVALID_ARGUMENT", "inputs mix file paths with sequences or records")
    if paths:
        return [f for p in paths for f in _expand(p)], []
    return [], [v if k == "rec" else (str(n), v) for n, (k, v) in enumerate(kinds, 1)]


def run(
    sequences,
    output_dir: str | Path | None = None,
    *,
    germline_db: str = "human",
    locus: str | list[str] | None = None,
    include_nonfunctional: bool = False,
    formats: str | Iterable[str] = ("airr",),
    threads: int | None = None,
    strict: bool = False,
    sequence_type: str = "nt",
) -> "pl.DataFrame | RunSummary":
    """Annotate antibody/TCR sequences.

    `sequences` is a FASTA/FASTQ path (optionally gzipped), a directory of them
    (searched recursively, sorted by relative path), a single sequence string, an
    iterable of any one of those kinds, or an iterable of `(id, sequence)` tuples.
    A string that exists on disk is a path; otherwise, if it consists only of
    nucleotide characters (`ACGTN` and IUPAC `RYSWKMBDHV`, any case) it is a
    sequence, with ID `"1"`..`"N"` by position, and anything else raises
    `INPUT_NOT_FOUND`. With `sequence_type="aa"` the input is protein: a string
    is a sequence if it consists only of amino-acid letters (20 standard, `XBZUOJ`
    and `*`, any case). A path that cannot be checked (for example, permission
    denied) raises `INPUT_UNREADABLE`. `"-"` is stdin. Mixing paths and
    sequences raises `INVALID_ARGUMENT`. Iterables are consumed once.

    With `output_dir=None` a Polars DataFrame is returned and `formats` is
    ignored. Otherwise files are written (`formats` from `airr`, `parquet`) and a
    `RunSummary` is returned. `threads` is None or >= 1.

    Alleles the database marks non-functional (pseudogene `P`, OGRDB `NF`) are
    not searched unless `include_nonfunctional` is true; `F`, `ORF` and
    `unknown` alleles always are.

    `sequence_type` is `"nt"` (default) or `"aa"`, any case; anything else raises
    `INVALID_ARGUMENT`.
    """
    fmts = [formats] if isinstance(formats, str) else list(formats)
    if not fmts or not set(fmts) <= _FORMATS:
        raise AbstarError("INVALID_ARGUMENT", f"formats must be a non-empty subset of {sorted(_FORMATS)}, got {fmts!r}")
    if threads is not None and (isinstance(threads, bool) or not isinstance(threads, int) or threads < 1):
        raise AbstarError("INVALID_ARGUMENT", f"threads must be None or an integer >= 1, got {threads!r}")
    sequence_type = sequence_type.lower() if isinstance(sequence_type, str) else sequence_type
    if sequence_type not in ("nt", "aa"):
        raise AbstarError("INVALID_ARGUMENT", f"sequence_type must be 'nt' or 'aa', got {sequence_type!r}")
    loci = [locus] if isinstance(locus, str) else (list(locus) if locus is not None else None)

    from . import _core, germline

    inputs, records = _normalize(sequences, sequence_type == "aa")
    db_path = str(germline.resolve(germline_db)[0])
    try:
        if output_dir is None:
            import polars as pl

            return pl.DataFrame(_core.run_in_memory(inputs, records, db_path, threads, strict, loci, bool(include_nonfunctional), sequence_type))
        d = _core.run_to_files(inputs, records, str(output_dir), db_path, fmts, threads, strict, loci, bool(include_nonfunctional), sequence_type)
    except _core.CoreError as e:
        raise from_core(e) from None
    return RunSummary(
        inputs=[
            InputSummary(i["source"], {k: Path(v) for k, v in i["outputs"].items()}, dict(i["counts"]))
            for i in d["inputs"]
        ],
        counts={k: d["counts"][k] for k in _STATUSES},
        run_json=Path(d["run_json"]),
    )
