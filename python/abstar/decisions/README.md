# Decision records

Each file records one consequential design decision, so the docs can explain what every non-AIRR field means and why. They ship inside the package; `abstar manual decisions [ID]` shows them.

## File naming

`NNNN-<slug>.md`: a four-digit number, contiguous from `0001`, then a lowercase slug of letters, digits and hyphens.

## Format

The first line is `# NNNN: <Title>`, followed by these sections, in order:

- `## Status`: `Accepted (YYYY-MM-DD)` or `Superseded by NNNN (YYYY-MM-DD)`.
- `## Context`: the problem and constraints.
- `## Decision`: what was decided.
- `## Alternatives considered`: what else was weighed and why it lost.
- `## Consequences`: what follows from the decision.

## Reversing a decision

Never edit a record away. Write a new record that supersedes it and change the old record's status to `Superseded by NNNN (YYYY-MM-DD)`.

Every `abstar`-origin field in `schema/fields.toml` names its defining record (`decision = "NNNN"`).
