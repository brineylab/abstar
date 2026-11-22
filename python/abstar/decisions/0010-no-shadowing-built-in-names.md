# 0010: User databases may not shadow built-in names

## Status

Accepted (2026-11-20)

## Context

User databases live under `~/.abstar/germline/<name>/`. Today a user database with the same name as a built-in one silently overrides it.

## Decision

*Amended 2026-11-21 (#141): A directory in the working directory could still shadow a built-in: `resolve()` tried "a directory holding `manifest.toml`" before names, so `./human/` replaced the built-in `human` for `-g human`. The same held for recombination models, where a local `human-igh/` even replaced the default model when no `--model` was given. Names now come first. A value is taken as a path without trying names only when it is a `Path` object or contains a path separator (`./human`, `dbs/human`). A bare name that matches no built-in or user database still finds a local directory, so `-g mydb` keeps working from the parent directory.*

Building a user database whose name collides with a built-in one is an error. Built-in databases ship inside the package. This was a call made in the design.

## Alternatives considered

- Keep today's silent override: a footgun, since results change with no visible signal.

## Consequences

- The name `human` always means the built-in database; a local `./human` directory must be named as a path.
- Users who want a variant pick a distinct name.
