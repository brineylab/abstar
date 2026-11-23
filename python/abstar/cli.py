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
model_app = typer.Typer(help="List and inspect recombination models for abstar generate.", no_args_is_help=True)
app.add_typer(model_app, name="model")


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
    include_nonfunctional: Annotated[bool, typer.Option("--include-nonfunctional", help="Also search pseudogene (P) and non-functional (NF) alleles.")] = False,
    sequence_type: Annotated[str, typer.Option("--sequence-type", help="nt (default) or aa (protein input).")] = "nt",
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
        s = api.run(items, output_dir, germline_db=germline_db, locus=locus, include_nonfunctional=include_nonfunctional, formats=fmts, threads=threads, strict=strict, sequence_type=sequence_type)
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

    inputs, records = api._normalize(items, sequence_type == "aa")
    db_path = str(germline.resolve(germline_db)[0])
    if verbose:
        _err(f"database: {db_path}\ninputs: {inputs + [r[0] for r in records]}")
    sys.stdout.flush()
    sys.stderr.flush()
    try:
        counts = _core.run_to_stdout(inputs, records, db_path, threads, strict, locus, include_nonfunctional, sequence_type)
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
    if r.restricted or r.passed_through:  # a keep-list was applied
        print(f"restricted: {', '.join(r.restricted) or '-'}")
        print(f"passed through: {', '.join(r.passed_through) or '-'}")
    for w in r.warnings:
        _err(f"warning: {w['allele']}: {w['message']}")


@germline_app.command("gap")
def germline_gap(
    v_fasta: Annotated[Path, typer.Argument(help="FASTA of ungapped V alleles (gzip allowed).")],
    templates: Annotated[Optional[list[str]], typer.Option("--templates", help="Template database name or path; repeatable (default: human, c57bl6, balbc).")] = None,
    output: Annotated[Optional[Path], typer.Option("--output", help="Gapped FASTA file (default: stdout).")] = None,
    report: Annotated[Optional[Path], typer.Option("--report", help="Per-allele report TSV.")] = None,
    json_out: Annotated[bool, typer.Option("--json", help="Machine-readable summary.")] = False,
) -> None:
    """IMGT-gap V alleles using the best template allele."""
    from . import _core, germline

    recs = germline._read_fasta(v_fasta)
    try:
        out, rep = germline.gap(recs, templates)
    except AbstarError as e:
        if report and getattr(e, "report", None) is not None:
            e.report.write_csv(report, separator="\t")
        raise
    text = "".join(f">{h}\n{s}\n" for h, s in out)
    if output:
        output.write_text(text)
    else:
        sys.stdout.write(text)
    if report:
        rep.write_csv(report, separator="\t")
    failed = rep.filter(rep["status"] == "failed")
    used = []
    for t in templates if templates is not None else list(_core.GAP_DEFAULT_TEMPLATES):
        i = germline.info(t)
        used.append({"name": i.name, "path": str(i.path), "alleles_sha256": i.alleles_sha256})
    dest = sys.stdout if output else sys.stderr
    if json_out:
        failures = [{"name": n, "failure_reason": r} for n, r in zip(failed["name"], failed["failure_reason"])]
        print(json.dumps({"n_input": len(rep), "n_gapped": len(out), "n_failed": len(failed), "templates": used, "failed": failures}), file=dest)
        return
    print(f"gapped {len(out)}, failed {len(failed)}", file=dest)
    print("templates: " + ", ".join(f"{t['name']} ({t['alleles_sha256'][:8]})" for t in used), file=dest)
    for n, r in zip(failed["name"], failed["failure_reason"]):
        print(f"failed {n}: {r}", file=dest)


@germline_app.command("list")
def germline_list(json_out: Annotated[bool, typer.Option("--json", help="Machine-readable output.")] = False) -> None:
    """List built-in and user databases."""
    from . import germline

    dbs = germline.list()
    if json_out:
        print(json.dumps([asdict(d) for d in dbs], default=str))
    for d in [] if json_out else dbs:
        print(f"{d.name}\t{d.kind}\t{d.n_alleles}\t{d.path}" + (f"\tinvalid: {' '.join(d.error.split())}" if d.error else ""))


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
    print(f"built: {d.built} (abstar {d.abstar_version})\nalleles_sha256: {d.alleles_sha256}")
    if d.prefixes:
        print("prefixes: " + ", ".join(f"{s}={p}" for s, p in d.prefixes.items()))
    print("counts:")
    for k, v in d.counts.items():
        print(f"  {k}\t{v}")
    print("provenance:")
    for p in d.provenance:
        print(f"  {p['species']}: {p['description']}")
        for f, h in zip(p["files"], p["sha256"]):
            print(f"    {f}  sha256 {h}")
    print(f"warnings: {len(d.warnings)}")
    for w in d.warnings:
        print(f"  {w['allele']}: {w['message']}")


def _gen_line(s) -> str:
    return f"n={s.n} productive={s.productive} attempts={s.attempts} seed={s.seed}"


@app.command()
def generate(
    n: Annotated[int, typer.Option("--n", help="Number of reads.")],
    locus: Annotated[str, typer.Option("--locus", help="IGH, IGK, IGL, TRA or TRB with the built-in models.")],
    output_dir: Annotated[Optional[Path], typer.Argument(help="Write files here; omit to stream FASTA to stdout.")] = None,
    germline_db: Annotated[str, typer.Option("--germline-db", "-g", help="Database name or path.")] = "human",
    model: Annotated[Optional[str], typer.Option("--model", help="Model name or path (default: <species>-<locus>).")] = None,
    preset: Annotated[Optional[str], typer.Option("--preset", help="naive, memory, bnab or tcr; options override it.")] = None,
    seed: Annotated[Optional[int], typer.Option("--seed", help="Random seed (default: drawn and recorded).")] = None,
    mutation: Annotated[Optional[str], typer.Option("--mutation", help="V mutation frequency: X, uniform:LO-HI or beta:A,B.")] = None,
    long_cdr3: Annotated[Optional[float], typer.Option("--long-cdr3", help="Per-junction weight of the long-tail N lengths.")] = None,
    indels: Annotated[Optional[float], typer.Option("--indels", help="Per-read probability of one indel.")] = None,
    frameshift_indels: Annotated[bool, typer.Option("--frameshift-indels", help="Allow indels that shift the frame.")] = False,
    truncate_5p: Annotated[Optional[str], typer.Option("--truncate-5p", help="nt removed from the 5' end: N or uniform:LO-HI.")] = None,
    truncate_3p: Annotated[Optional[str], typer.Option("--truncate-3p", help="nt removed from the 3' end: N or uniform:LO-HI.")] = None,
    constant: Annotated[Optional[int], typer.Option("--constant", help="Append this many nt of constant region.")] = None,
    include_nonfunctional: Annotated[bool, typer.Option("--include-nonfunctional", help="Also draw pseudogene (P) and non-functional (NF) alleles.")] = False,
    productive_only: Annotated[Optional[bool], typer.Option("--productive-only/--no-productive-only", help="Keep only productive reads (default: the preset's; bnab turns it on).", show_default=False)] = None,
    shm_table: Annotated[Optional[Path], typer.Option("--shm-table", help="Directory of S5F CSVs to mutate with.")] = None,
    name: Annotated[str, typer.Option("--name", help="Output file name stem.")] = "generated",
    format: Annotated[Optional[list[str]], typer.Option("--format", help="Truth table: airr or parquet (repeatable; needs OUTPUT_DIR).")] = None,
    threads: Annotated[Optional[int], typer.Option("--threads", help="Worker threads (>= 1).")] = None,
    json_out: Annotated[bool, typer.Option("--json", help="Machine-readable output.")] = False,
    quiet: Annotated[bool, typer.Option("--quiet", "-q", help="Suppress the summary.")] = False,
) -> None:
    """Simulate V(D)J reads with truth. Without OUTPUT_DIR, FASTA streams to stdout."""
    if output_dir is None and format:
        raise AbstarError("INVALID_ARGUMENT", "--format requires OUTPUT_DIR (without it only FASTA is written)")
    if threads is not None and threads < 1:
        raise AbstarError("INVALID_ARGUMENT", f"--threads must be >= 1, got {threads}")
    options = {
        "mutation": mutation, "long_cdr3_weight": long_cdr3, "indels": indels, "truncate_5p": truncate_5p,
        "truncate_3p": truncate_3p, "constant": constant, "shm_table": None if shm_table is None else str(shm_table),
        # flags only override a preset when given
        "frameshift_indels": frameshift_indels or None, "productive_only": productive_only,
        "include_nonfunctional": include_nonfunctional or None,
    }
    from . import _generate as gen

    kw = dict(germline_db=germline_db, model=model, preset=preset, seed=seed, threads=threads, **options)
    if output_dir is not None:
        s = gen.generate(n, locus, output_dir=output_dir, name=name, formats=format or ["airr"], **kw)
        if json_out:
            print(_json(s))
        elif not quiet:
            print(_gen_line(s))
            for p in [s.fasta, *s.truth.values(), s.generate_json]:
                print(f"  {p}")
        for w in [] if quiet or json_out else s.warnings:
            _err(f"warning: {w}")
        return

    from . import _core

    db_path, model_path, params = gen._prepare(n, locus, germline_db, model, preset, seed, threads, options)
    if seed is None and not json_out:
        _err(f"seed: {params['seed']}")  # even with -q: a streamed run is reproducible only by its seed
    sys.stdout.flush()
    sys.stderr.flush()
    try:
        s = gen._summary(_core.generate_fasta_to_stdout(db_path, model_path, params))
    except _core.CoreError as e:
        raise from_core(e) from None
    finally:
        # Rust wrote to fd 1 directly; a closed pipe must not raise at interpreter shutdown.
        os.dup2(os.open(os.devnull, os.O_WRONLY), 1)
    if json_out:
        _err(_json(s))
    elif not quiet:
        _err("\n".join([_gen_line(s)] + [f"warning: {w}" for w in s.warnings]))


@model_app.command("list")
def model_list_cmd(json_out: Annotated[bool, typer.Option("--json", help="Machine-readable output.")] = False) -> None:
    """List built-in and user recombination models."""
    from . import model

    ms = model.list()
    if json_out:
        print(json.dumps([asdict(m) for m in ms], default=str))
    for m in [] if json_out else ms:
        print(f"{m.name}\t{m.kind}\t{m.species}\t{m.locus}\t{m.path}" + (f"\tinvalid: {' '.join(m.error.split())}" if m.error else ""))


@model_app.command("info")
def model_info_cmd(
    name: Annotated[str, typer.Argument(help="Model name or path.")],
    json_out: Annotated[bool, typer.Option("--json", help="Machine-readable output.")] = False,
) -> None:
    """Show details of one recombination model."""
    from . import model

    m = model.info(name)
    if json_out:
        print(_json(m))
        return
    print(f"{m.name} ({m.kind}): {m.path}\n{m.description}\nspecies: {m.species}\nlocus: {m.locus}\nsources:")
    for src in m.sources:
        print(f"  {src['used_for']}\n    {src['citation']}")


@app.command()
def manual(
    topic: Annotated[Optional[str], typer.Argument(help="Topic name; omit for the whole manual.")] = None,
    subtopic: Annotated[Optional[str], typer.Argument(help="Field name (fields, truth-fields) or decision ID (decisions).")] = None,
) -> None:
    """Print the manual as Markdown: abstar manual [TOPIC] [SUBTOPIC]. Topics: commands, run, germline, generate, manual, examples, fields, truth-fields, conventions, ambiguity, errors, decisions, migration."""
    from . import manual as _manual

    print(_manual.render(topic, subtopic))


def main() -> None:
    code: int | str | None = 0
    try:
        try:
            app()
        except AbstarError as e:
            if "--json" in sys.argv[1:]:
                _err(json.dumps({"error": {"code": e.code, "message": e.message}}))
            else:
                _err(f"error[{e.code}]: {e.message}")
            code = e.exit_code
        except SystemExit as e:
            # Typer turns EPIPE inside a command into exit 1 after wrapping stdout in
            # (its vendored) click's PacifyFlushWrapper; matched by name, not imported.
            if type(sys.stdout).__name__ == "PacifyFlushWrapper":
                raise BrokenPipeError from None
            code = e.code
        sys.stdout.flush()  # surface a closed pipe here, not at interpreter shutdown
    except BrokenPipeError:
        # The reader went away (`| head`): that ends the output, it is not an error. An error's
        # exit code, already set, stands (#142).
        os.dup2(os.open(os.devnull, os.O_WRONLY), sys.stdout.fileno())
    sys.exit(code)
