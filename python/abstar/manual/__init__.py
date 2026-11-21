"""`abstar manual`: Markdown reference generated from the real sources."""
from __future__ import annotations

import tomllib
from pathlib import Path

from ..errors import ERROR_CODES, EXIT_CODES, AbstarError

TOPICS: dict[str, str] = {
    "commands": "Every command and option",
    "run": "The run command",
    "germline": "The germline commands",
    "generate": "Synthetic reads with truth: the generate and model commands, presets, grading, limits",
    "manual": "The manual command",
    "examples": "Common commands and Python API calls",
    "fields": "Output fields (fields NAME adds the decision record)",
    "truth-fields": "abstar generate truth table fields (truth-fields NAME adds the decision record)",
    "conventions": "Coordinates, TSV encoding, record statuses, sequence IDs",
    "ambiguity": "The primary call and ties",
    "errors": "Error codes and exit codes",
    "decisions": "Decision records (decisions ID prints one)",
    "migration": "Old (0.x) options and what replaces them",
}
_DECISIONS = Path(__file__).parent.parent / "decisions"
_TOPIC_DIR = Path(__file__).parent / "topics"


def _bad(kind: str, got: str, valid) -> AbstarError:
    return AbstarError("INVALID_ARGUMENT", f"unknown {kind} {got!r}; valid: {', '.join(valid)}")


def _tree():
    import typer.main

    from .. import cli

    return typer.main.get_command(cli.app)


def _cmd_md(cmd, path: str, level: int) -> str:
    out = [f"{'#' * level} `{path}`", ""]
    args = [p for p in cmd.params if p.param_type_name == "argument"]
    usage = path + (" [OPTIONS]" if any(p.param_type_name != "argument" for p in cmd.params) else "")
    usage += "".join(f" [{p.name.upper()}]" if not p.required else f" {p.name.upper()}" for p in args)
    if hasattr(cmd, "commands"):
        usage += " COMMAND [ARGS]..."
    out += [f"Usage: `{usage}`", ""]
    if cmd.help:
        out += [cmd.help.strip(), ""]
    for p in cmd.params:
        ann = p.type_descriptor.annotation if hasattr(p, "type_descriptor") else str
        typ = "flag" if getattr(p, "is_flag", False) else getattr(ann, "__name__", None) if isinstance(ann, type) else str(ann)
        names = ", ".join(f"`{o}`" for o in p.opts + p.secondary_opts)
        if p.param_type_name == "argument":
            names = f"`{p.name.upper()}`"
        default = "" if p.default is None or p.is_flag else f", default `{p.default}`"
        out.append(f"- {names} ({typ}{default}): {p.help or ''}")
    out.append("")
    for name, sub in getattr(cmd, "commands", {}).items():
        out.append(_cmd_md(sub, f"{path} {name}", level + 1))
    return "\n".join(out)


def _commands(only: str | None) -> str:
    root = _tree()
    if only is None:
        return _cmd_md(root, "abstar", 2)
    return _cmd_md(root.commands[only], f"abstar {only}", 2)


def _decision_files() -> dict[str, Path]:
    return {p.name[:4]: p for p in sorted(_DECISIONS.glob("[0-9][0-9][0-9][0-9]-*.md"))}


def _title(p: Path) -> str:
    return p.read_text().splitlines()[0].split(": ", 1)[1]


def _status(p: Path) -> str:
    lines = p.read_text().splitlines()
    return next(l for l in lines[lines.index("## Status") + 1:] if l.strip())


def _decisions(sub: str | None) -> str:
    files = _decision_files()
    if sub is None:
        rows = ["| ID | Title | Status |", "|---|---|---|"]
        rows += [f"| {i} | {_title(p)} | {_status(p)} |" for i, p in files.items()]
        return "\n".join(rows) + "\n"
    key = sub.zfill(4)
    if key not in files:
        raise _bad("decision ID", sub, files)
    return files[key].read_text()


def _fields(toml: str, sub: str | None) -> str:
    doc = tomllib.loads(toml)
    fields = {f["name"]: f for f in doc["field"]}
    files = _decision_files()
    if sub is not None and sub not in fields:
        raise _bad("field", sub, fields)

    def one(f: dict) -> str:
        d = f["decision"]
        dec = f"{d}: {_title(files[d])}" if d else "-"
        return "\n".join([
            f"- `{f['name']}`: {f['description']}",
            f"  - type: {f['type']}; nullable: {f['nullable']}; origin: {f['origin']}; "
            f"coordinates: {f['coordinates'] or '-'}",
            f"  - decision: {dec}",
        ])

    if sub is not None:
        f = fields[sub]
        text = f"## Field `{sub}`\n\n{one(f)}\n"
        if f["decision"]:
            text += "\n" + files[f["decision"]].read_text()
        return text
    return f"AIRR schema version {doc['airr_version']}.\n\n" + "\n".join(one(f) for f in fields.values()) + "\n"


def _errors() -> str:
    rows = ["## Exit codes", "", "| Code | Meaning |", "|---|---|"]
    rows += [f"| {n} | {m} |" for n, m in EXIT_CODES.items()]
    rows += ["", "## Error codes", "", "| Code | Exit | Meaning |", "|---|---|---|"]
    rows += [f"| {c} | {n} | {m} |" for c, (n, m) in ERROR_CODES.items()]
    return "\n".join(rows) + "\n"


def _topic(name: str, sub: str | None) -> str:
    if sub is not None and name not in ("fields", "truth-fields", "decisions"):
        raise AbstarError("INVALID_ARGUMENT", f"topic {name!r} takes no subtopic")
    if name == "commands":
        return _commands(None)
    if name in ("run", "germline", "manual"):
        return _commands(name)
    if name in ("fields", "truth-fields"):
        from .. import _core

        return _fields(_core.fields_toml() if name == "fields" else _core.truth_toml(), sub)
    if name == "generate":
        return (_TOPIC_DIR / "generate.md").read_text() + "\n" + _commands("generate") + "\n" + _commands("model")
    if name == "errors":
        return _errors()
    if name == "decisions":
        return _decisions(sub)
    return (_TOPIC_DIR / f"{name}.md").read_text()


def render(topic: str | None = None, subtopic: str | None = None) -> str:
    """Markdown for one topic, or the whole manual when `topic` is None."""
    if topic is None:
        from .. import __version__

        parts = [f"# abstar {__version__} manual\n"]
        parts += [f"# {t}\n\n{_topic(t, None)}" for t in TOPICS]
        return "\n".join(parts)
    if topic not in TOPICS:
        raise _bad("topic", topic, TOPICS)
    return _topic(topic, subtopic)
