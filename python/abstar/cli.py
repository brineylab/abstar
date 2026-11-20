"""The `abstar` command line. stdout carries results only; stderr carries everything else."""
import json
import os
import sys
from dataclasses import asdict
from pathlib import Path
from typing import Annotated, Optional

import typer

from .errors import AbstarError, from_core

app = typer.Typer(
    help="Annotate antibody and TCR sequences.",
    no_args_is_help=True,
    add_completion=False,
    pretty_exceptions_enable=False,
)
germline_app = typer.Typer(help="Build, list and inspect germline databases.", no_args_is_help=True)
app.add_typer(germline_app, name="germline")


def _version(value: bool) -> None:
    if value:
        from . import __version__

        print(__version__)
        raise typer.Exit()


@app.callback()
def _root(
    version: Annotated[bool, typer.Option("--version", callback=_version, is_eager=True, help="Show the version.")] = False,
) -> None:
    pass


def _json(obj) -> str:
    return json.dumps(asdict(obj) if hasattr(obj, "__dataclass_fields__") else obj, default=str)


def _err(*a) -> None:
    print(*a, file=sys.stderr)


def _counts_line(counts: dict) -> str:
    return " ".join(f"{k}={v}" for k, v in counts.items())


@app.command()
def run(
    sequences: Annotated[Optional[str], typer.Argument(help="FASTA/FASTQ file, directory, or '-' for stdin.")] = None,
    output_dir: Annotated[Optional[Path], typer.Argument(help="Write files here; omit to stream TSV to stdout.")] = None,
    seq: Annotated[Optional[str], typer.Option("--seq", help="A single sequence (ID '1') instead of SEQUENCES.")] = None,
    germline_db: Annotated[str, typer.Option("--germline-db", "-g", help="Database name or path.")] = "human",
    locus: Annotated[Optional[list[str]], typer.Option("--locus", help="Restrict to a locus (repeatable).")] = None,
    format: Annotated[Optional[list[str]], typer.Option("--format", help="airr or parquet (repeatable).")] = None,
    threads: Annotated[Optional[int], typer.Option("--threads", help="Worker threads (>= 1).")] = None,
    strict: Annotated[bool, typer.Option("--strict", help="Fail (exit 3) on any record error.")] = False,
    json_out: Annotated[bool, typer.Option("--json", help="Machine-readable output.")] = False,
    quiet: Annotated[bool, typer.Option("--quiet", "-q", help="Suppress the summary.")] = False,
    verbose: Annotated[bool, typer.Option("--verbose", "-v", help="Print the database and inputs to stderr.")] = False,
) -> None:
    """Annotate sequences. Without OUTPUT_DIR, TSV streams to stdout."""
    if seq is not None and sequences is not None:
        if output_dir is not None or (os.path.exists(sequences) and not os.path.isdir(sequences)):
            raise AbstarError("INVALID_ARGUMENT", "SEQUENCES and --seq are mutually exclusive")
        sequences, output_dir = None, Path(sequences)  # `--seq ACGT out`: the lone positional is OUTPUT_DIR
    if (sequences is None) == (seq is None):
        raise AbstarError("INVALID_ARGUMENT", "give exactly one of SEQUENCES or --seq")
    fmts = format or ["airr"]
    if threads is not None and threads < 1:
        raise AbstarError("INVALID_ARGUMENT", f"threads must be >= 1, got {threads}")
    if output_dir is None and "parquet" in fmts:
        raise AbstarError("INVALID_ARGUMENT", "--format parquet requires OUTPUT_DIR")

    from . import api, germline

    # SEQUENCES is always a path (or `-`), never a sequence: that is what --seq is for.
    items = [("1", seq)] if seq is not None else sequences if sequences == "-" else Path(sequences)
    if output_dir is not None:
        if verbose:
            _err(f"database: {germline.resolve(germline_db)[0]}\ninputs: {items}")
        s = api.run(items, output_dir, germline_db=germline_db, locus=locus, formats=fmts, threads=threads, strict=strict)
        if json_out:
            print(_json(s))
        elif not quiet:
            for i in s.inputs:
                print(f"{i.source}: {_counts_line(i.counts)}")
                for fmt, p in i.outputs.items():
                    print(f"  {fmt}: {p}")
            print(f"total: {_counts_line(s.counts)}")
        return

    from . import _core

    inputs, records = api._normalize(items)
    db_path = str(germline.resolve(germline_db)[0])
    if verbose:
        _err(f"database: {db_path}\ninputs: {inputs + [r[0] for r in records]}")
    sys.stdout.flush()
    sys.stderr.flush()
    try:
        counts = _core.run_to_stdout(inputs, records, db_path, threads, strict, locus)
    except _core.CoreError as e:
        raise from_core(e) from None
    finally:
        # Rust wrote to fd 1 directly; a closed pipe must not raise at interpreter shutdown.
        os.dup2(os.open(os.devnull, os.O_WRONLY), 1)
    if json_out:
        _err(json.dumps({"counts": counts}))
    elif not quiet:
        _err(_counts_line(counts))


@germline_app.command("build")
def germline_build(
    recipe: Annotated[Path, typer.Argument(help="Recipe TOML file.")],
    output: Annotated[Optional[Path], typer.Option("--output", help="Output directory (default: user germline dir).")] = None,
    overwrite: Annotated[bool, typer.Option("--overwrite", help="Replace an existing database.")] = False,
    json_out: Annotated[bool, typer.Option("--json", help="Machine-readable output.")] = False,
) -> None:
    """Build a germline database from a recipe."""
    from . import germline

    r = germline.build(recipe, output_dir=output, overwrite=overwrite)
    if json_out:
        print(_json(r))
        return
    print(f"built {r.name}: {r.n_alleles} alleles -> {r.path}")
    for w in r.warnings:
        _err(f"warning: {w}")


@germline_app.command("list")
def germline_list(json_out: Annotated[bool, typer.Option("--json", help="Machine-readable output.")] = False) -> None:
    """List built-in and user databases."""
    from . import germline

    dbs = germline.list()
    if json_out:
        print(json.dumps([asdict(d) for d in dbs], default=str))
    for d in [] if json_out else dbs:
        print(f"{d.name}\t{d.kind}\t{d.n_alleles}\t{d.path}")


@germline_app.command("info")
def germline_info(
    name: Annotated[str, typer.Argument(help="Database name or path.")],
    json_out: Annotated[bool, typer.Option("--json", help="Machine-readable output.")] = False,
) -> None:
    """Show details of one database."""
    from . import germline

    d = germline.info(name)
    if json_out:
        print(_json(d))
        return
    print(f"{d.name} ({d.kind}): {d.path}\n{d.description}\nspecies: {', '.join(d.species)}\nloci: {', '.join(d.loci)}")
    for k, v in d.counts.items():
        print(f"  {k}\t{v}")


@app.command()
def manual(
    topic: Annotated[Optional[str], typer.Argument(help="Topic name; omit for the whole manual.")] = None,
    subtopic: Annotated[Optional[str], typer.Argument(help="Field name (fields) or decision ID (decisions).")] = None,
) -> None:
    """Print the manual as Markdown: abstar manual [TOPIC] [SUBTOPIC]. Topics: commands, run, germline, manual, fields, conventions, ambiguity, errors, decisions."""
    from . import manual as _manual

    print(_manual.render(topic, subtopic))


def main() -> None:
    try:
        app()
    except AbstarError as e:
        if "--json" in sys.argv[1:]:
            _err(json.dumps({"error": {"code": e.code, "message": e.message}}))
        else:
            _err(f"error[{e.code}]: {e.message}")
        sys.exit(e.exit_code)
