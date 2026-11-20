# 0010: User databases may not shadow built-in names

## Status

Accepted (2026-11-20)

## Context

User databases live under `~/.abstar/germline/<name>/`. Today a user database with the same name as a built-in one silently overrides it.

## Decision

Building a user database whose name collides with a built-in one is an error. Built-in databases ship inside the package. This was a call made in the design.

## Alternatives considered

- Keep today's silent override: a footgun, since results change with no visible signal.

## Consequences

- The name `human` always means the built-in database.
- Users who want a variant pick a distinct name.
