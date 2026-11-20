# 0011: Agents as first-class users (manual, stdout/stderr, `--json`, codes)

## Status

Accepted (2026-11-20)

## Context

Everything an agent needs to use abstar correctly must be discoverable from the tool itself. A manual, `--help` or docs that drift from the code are an anti-goal.

## Decision

- `abstar manual` prints the complete reference as one Markdown document; `abstar manual <topic>` narrows it. It is generated from the Typer command definitions and `schema/fields.toml`, never hand-written. Prose that cannot be generated lives in one Markdown source shared with the Sphinx docs.
- Results go to stdout only when writing to stdout; logs and progress go to stderr. Progress bars appear only when stderr is a terminal.
- No interactive prompts, ever.
- `--json` on summaries, `germline list`/`info`, and errors.
- Stable, documented exit codes and stable machine-readable error codes for each failure class.

## Alternatives considered

- An MCP server or bundled agent skill: not built now; both can be thin layers over `abstar manual` later.
- A hand-written manual: rejected, since it drifts.

## Consequences

- The manual cannot drift from the code.
- Scripted and agent use is predictable from streams, codes and JSON.
