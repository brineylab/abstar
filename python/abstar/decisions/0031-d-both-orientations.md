# 0031: D searched in both orientations, with an inverted-D margin

## Status

Accepted (2026-11-20)

## Context

A D gene can recombine inverted, so the read carries its reverse complement. The generator models this with `d_inversion`. An annotator cannot know the rate that produced a read. The output field `d_inverted` reports the orientation of the D call.

Truth compares every D allele in both orientations, whatever the model's `d_inversion` (record 0025, amended 2026-11-21 by ruling A24). The valid set is what the read supports.

Inversion is rare in real repertoires. A short D has many chance reverse-complement matches, which outscore the forward D more often than true inversions occur (ruling A25).

## Decision

- The junction step (record 0030) scores every D allele of the locus both forward and reverse-complemented, at every offset, with the same score.
- **Inverted-D margin of 2.** The D is called inverted only when the best reverse-complement score beats the best forward score by at least 2. Otherwise it is called forward, so forward wins ties.
- The D tie set holds the alleles reaching the best score in the called orientation. `d_inverted` is that orientation.
- `d_germline_start` and `d_germline_end` (`DHit::g_start`, `g_end`) index the reverse-complemented allele when inverted. This matches the generator's convention: `Origin::D(i)` indexes the reverse-complemented D when `d_inverted`. Internally they are 0-based half-open (record 0016).
- `d_inverted` is null without a D call.

**Margin sweep** (ruling A25). The forward sets are human-igh reads (`d_inversion` = 0). The inversion sets are the same model with `d_inversion` = 0.5, seed 3, 2,000 reads.
- D call: the primary is in `d_call_valid`, or it is null and `d_call_null_valid`, over recoverable reads.
- Inverted recall: reads whose true D is inverted with a run of at least 5 nt that are called inverted.

| M | D call, naive | D call, 0–30% | D call, bnab | Inverted recall, naive inv. set | Inverted recall, 0–30% inv. set | D call, naive inv. set | D call, 0–30% inv. set |
|---|---|---|---|---|---|---|---|
| 0 | 96.05% | 85.53% | 70.45% | 96.42% (944/979) | 87.19% (660/757) | 95.35% | 83.74% |
| 2 | 96.15% | 86.09% | 72.35% | 93.16% (912/979) | 80.32% (608/757) | 94.40% | 82.42% |
| 4 | 96.20% | 86.59% | 72.85% | 86.41% (846/979) | 62.75% (475/757) | 91.70% | 75.93% |
| 6 | 96.20% | 86.53% | 72.85% | 78.96% (773/979) | 50.33% (381/757) | 88.20% | 70.83% |

M = 2 takes most of the forward-set gain (bnab +1.9 points over M = 0, against +2.4 at M = 4). It keeps 93% of inverted Ds on naive reads and 80% at 0–30% mutation. M = 4 adds at most 0.5 points on the forward sets but drops inverted recall at 0–30% mutation to 63%.

## Alternatives considered

- **Forward only:** misses every true inversion, which truth now always accepts.
- **No margin (M = 0):** recall is highest, but chance reverse-complement placements cost 1.9 points of D calls on bnab reads (table above).
- **A larger margin (M = 4 or 6):** no further forward-set gain beyond 0.5 points, and inverted recall collapses on mutated reads.
- **Ties across orientations in one tie set** (the first version): it mixed the orientations of one call, so `d_inverted` could disagree with half the tie set.

## Consequences

- `d_inverted` on human-igh reads (no inversion) is correct on 99.05% of naive reads with a D call. It is 98.18% at 0–30% mutation, 99.25% for memory + 20% indels and 96.65% for bnab.
- About 7% of inverted Ds on naive reads, and 20% at 0–30% mutation, are called forward (or as another D).
- Searching both orientations doubles the cost of the D search (record 0030).
- `tests/junction.rs` `inverted_d_recall` asserts recall ≥ 90% on the naive inversion set.
- The margin is a constant of this decision. Changing it needs a superseding record.
