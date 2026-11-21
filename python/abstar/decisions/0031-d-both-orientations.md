# 0031: D searched in both orientations

## Status

Accepted (2026-11-20)

## Context

A D gene can recombine inverted, so the read carries its reverse complement. The generator models this with `d_inversion`, and its truth compares D reverse complements when that rate is above 0 (record 0025). An annotator cannot know the rate that produced a read. The output field `d_inverted` reports the orientation of the D call.

## Decision

- The junction step (record 0030) scores every D allele of the locus both forward and reverse-complemented, at every offset, with the same score.
- An allele joins the D tie set when either orientation reaches the best score.
- `d_inverted` is the orientation of the primary D's best configuration. When the primary reaches the best score in both orientations, the forward one is used.
- `d_germline_start` and `d_germline_end` (`DHit::g_start`, `g_end`) index the reverse-complemented allele when inverted. This matches the generator's convention: `Origin::D(i)` indexes the reverse-complemented D when `d_inverted`. Internally they are 0-based half-open (record 0016).
- `d_inverted` is null without a D call.

## Alternatives considered

- **Forward only:** misses every true inversion, and truth from a model with `d_inversion` > 0 compares the reverse complement. On human-igh reads, which never invert, it would gain the 80 naive reads lost below. Even forward-only, the D call is correct on only about 95% of naive reads, because of the short-D ambiguity in record 0030.
- **An orientation prior** (reverse complement must win by a margin): it would trade true inversions for fewer chance ones. No model of the inversion rate exists yet, so it is left for a later record.

## Consequences

- On models that never invert D, chance reverse-complement placements that outscore the forward D are graded wrong. Truth from such a model compares the forward D only.
- Measured with the built-in human-igh model (`d_inversion` = 0), `d_inverted` is correct on 96.00% of naive reads with a D call. It is 90.05% at 0–30% mutation, 93.95% for memory + 20% indels and 83.37% for bnab.
- Searching both orientations doubles the cost of the D search (record 0030).
