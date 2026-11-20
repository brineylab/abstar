# 0024: Somatic hypermutation model: abstar's own motif default, user-supplied S5F

## Status

Accepted (2026-11-20)

## Context

`abstar generate` places somatic mutations by sequence context. The reference context model is S5F (Yaari et al. 2013, *Front. Immunol.* 4:358): it gives each 5-mer a mutability and gives the centre base a substitution profile. abstar is MIT-licensed and must not bundle data under GPL/AGPL or other restrictive terms. We checked where the S5F tables can be obtained (retrieved 2026-11-20):

- **The paper is CC BY, but its supplement has no tables.** The article (https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2013.00358/full, PMC3828525) has one supplementary file, `64394_Kleinstein_Presentation1.PDF`, which contains Figures S1-S5 as plots.
- **The numeric tables are CC BY-NC-SA 3.0 US.** The paper's model link is http://clip.med.yale.edu/shm/. It serves `Mutability.csv` and `Substitution.csv` under that license (stated at http://clip.med.yale.edu/shm/license.php). The NonCommercial and ShareAlike terms are incompatible with redistribution in an MIT package.
- **SHazaM is AGPL-3.** Its `HH_S5F` and `HKL_S5F` objects come from the R package SHazaM (https://cran.r-project.org/package=shazam), so they are excluded.
- **HKL_S5F (light chain) has no open license.** It comes from Cui et al. 2016, *J. Immunol.* 197:3566 (PMC5161250). That is an author manuscript with no open license; the article is © AAI.

## Decision

- **Default model: `motif`, built in code.** abstar ships its own 5-mer motif model. Its weights are abstar's heuristic choices from the published hotspot and coldspot motif rules. They are **not fit to data** and do not reproduce S5F rates. Mutability starts at 1.0 and is multiplied by:
  - WRCY / RGYW (centre C / G): ×6;
  - otherwise WRC / GYW: ×4;
  - WA / TW (centre A / T): ×2.5;
  - SYC / GRS (centre C / G): ×0.3.

  The strongest matching motif wins. Hot and cold motifs on the same centre cannot both match.
- **Substitution is transition-biased and context-free.**
  - From C or G: transition 0.6, each transversion 0.2.
  - From A or T: transition 0.5, each transversion 0.25.
- **Recombination models use `shm_model = "motif" | "none"`.** S5F names are not valid values.
- **S5F is user-supplied.** `--shm-table DIR` (Python: `shm_table`) points at a directory with `Mutability.csv` and `Substitution.csv` in the Kleinstein-lab format, which the user downloads under its own license.
  - Columns are found by header name: `Fivemer` and `Mutability`; `Fivemer`, `A`, `C`, `G`, `T`.
  - Fields may be whitespace- or comma-separated and quoted. Rows for 5-mers containing N are ignored, because abstar computes its own mean for those.
  - A missing file, malformed value, duplicate or bad 5-mer, substitution row with self ≠ 0 or sum ≠ 1, or fewer than 1024 ACGT 5-mers is rejected as `INVALID_ARGUMENT` naming the file.
  - Nothing S5F-derived is committed, including in test fixtures (they are synthetic).
- **Placement.**
  - A position is picked among mutable positions with probability ∝ the mutability of its 5-mer. The ends are padded with N.
  - A 5-mer containing N or an IUPAC code gets the table's mean mutability, and the mean substitution row for its centre base.
  - Only A/C/G/T centres are mutated, in either case; a lowercase base stays lowercase.
  - Contexts update after each hit, and repeat hits and back-mutations are allowed.
- **Stopping and the safety cap.**
  - Placement stops when the realized V divergence (V positions differing from their original base) reaches `round(target × V length)`.
  - A safety cap of 10 × target hits stops placement early. A capped run returns short of target. The caller, which records realized divergence for the truth record, detects this and adds a warning.
- **Later work.** We plan to fit our own S5F-style model on data we may redistribute, with the repertoire work (board #148).

## Alternatives considered

- **Bundle the Yale S5F tables:** the NC-SA license is incompatible with MIT.
- **Copy SHazaM's tables:** AGPL-3.
- **Download S5F at run time by default:** the generator would no longer be offline, and the NC terms would pass silently to every user. A user-supplied directory makes that choice explicit.
- **Treat the numbers as uncopyrightable facts:** a legal judgement we chose not to rely on.
- **Fit our own model now:** it needs redistributable training data and much more work, and grading annotators needs hotspot bias, not exact S5F rates.
- **Uniform mutation:** loses the hotspot clustering that makes real reads hard to annotate.

## Consequences

- Default mutation spectra show realistic hotspot and coldspot bias and a transition excess, but per-motif rates are not S5F rates. Benchmarks that need S5F-faithful spectra must use `--shm-table`.
- The motif weights are fixed constants. Changing them changes generated sequences for a given seed, so a change needs a superseding record.
- `mutate` has no warning channel: callers must compare realized divergence with the target and report capped sequences.
- When abstar fits its own model (#148), it is added as a new `shm_model` value; `motif` stays available for reproducibility.
