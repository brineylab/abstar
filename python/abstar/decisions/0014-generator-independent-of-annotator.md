# 0014: The generator may not depend on annotator code

## Status

Accepted (2026-11-20)

## Context

A generator that reuses annotation logic makes shared bugs score as correct. Accuracy benchmarks must not overfit to our own idiosyncrasies.

## Decision

`generate` must not depend on `annotate` or on any alignment or annotation code. The crate dependency graph enforces this, not convention: `generate` depends on `germline` only, and `annotate` depends on `germline`. Code that needs alignment, such as IMGT-gapping new-species germlines during database builds, lives outside `germline` so `generate` stays alignment-free. Exactly where it lives is the builder's call, as long as the rule holds. Evaluation also uses external tools and datasets, not just our generator.

## Alternatives considered

- Keep independence by convention or review: not reliable.

## Consequences

- Shared bugs cannot score as correct through shared code.
- The alignment-dependent build code needs its own home outside `germline`.

Amended by 0022: the generator may also depend on `abstar_table`.
