# 0027: Built-in human recombination models from published figures

## Status

Accepted (2026-11-20)

## Context

`abstar generate` ships five default models: `human-igh`, `human-igk`, `human-igl`, `human-tra` and `human-trb`. The spec requires them to come from published summary distributions (source C), with every table cited. Constraints:

- **No GPL or NC data.** abstar is MIT-licensed. IGoR and OLGA model files are GPL, and no GPL, AGPL or NonCommercial data may be bundled.
- **IGoR paper first.** The user asked that the IGoR paper be used wherever it covers a component: Marcou, Mora & Walczak 2018, *Nat Commun* 9:561, CC BY 4.0. The paper's numbers are published only as figures. The journal hosts no supplementary data file of model parameters.
- **Database fit.** Model gene names must bind to the built-in `human` database (record 0023).

## Decision

- **Digitize figures; never use software model files.** Distributions are read from the CC BY figures of peer-reviewed papers. Both PDFs (IGoR and Elhanati et al. 2015) contain vector figures. Bar heights and line vertices are read from the vector paths, so values are essentially exact. Raster figures are used only for the light chains (Mikocziova et al. 2021 box plots, about ±0.003). `germline-sources/models/SOURCES.md` holds:
  - the citations, licenses, URLs and sha256 hashes;
  - the extraction methods;
  - the digitized tables;
  - the extractors.

  `scripts/build_default_models.py` regenerates every `model.toml`. Each model's `[[source]]` entries name the paper and figure for every component and carry the CC BY attribution.
- **Per locus.**
  - **IGH:** all IGoR — V/J usage, P(D|J), the four trims and the VD/DJ insertions (Supplementary Figs. 4 and 13). N composition comes from Elhanati et al. 2015 (Figs. S9-S10, VD sense).
  - **TRB:** all IGoR — usage (Supplementary Fig. 3), P(D2|J) (Supplementary Fig. 14), V/J trims and insertions (Fig. 2), D trims (Supplementary Fig. 7). N composition is borrowed from IGH.
  - **TRA:**
    - IGoR trims and insertions (Fig. 2).
    - Usage for the six most used V and J genes comes from IGoR Supplementary Fig. 5. The rest is a uniform fallback over the database genes.
    - N composition is borrowed from IGH.
  - **IGK and IGL:**
    - Usage from Mikocziova et al. 2021 (naive repertoires, 100 donors).
    - N length is geometric with the Townsend et al. 2016 mean (fallback shape).
    - The V/J trims are the IGH profiles truncated to ≤ 6 nt and renormalized, a fallback surrogate. No published light-chain trim mean/SD or distribution was found; Collins & Watson 2018 say only that removal is small.
    - Truncation keeps the J F118 anchor, which is at ungapped index 7-8 in IGKJ/IGLJ. The full IGH J profile (P(≥ 8) ≈ 0.25) deleted it in about 23 % of IGK and 25 % of IGL reads. Truncated, the no-junction rate is 0.000 for both, against 0.02-0.04 for IGH/TRA/TRB, and a test asserts < 0.08.
    - N composition is borrowed from IGH.
- **Fit to the database at build time.** Source genes absent from the database are dropped when the models are built and the rest renormalized, rather than left for `bind` to warn about on every run. The dropped genes and the source mass they carried are listed in SOURCES.md and in each file's header. Retained source mass is at least 0.92 for every segment (IGH V 0.921, TRB V 0.924), and the build script asserts ≥ 0.90. The built-in models bind without warnings.
- **Generation-model distributions.** IGoR's distributions describe generation (pre-selection), inferred from nonproductive reads. We use them as published and do not reweight them toward productive repertoires.
- **IGH long tail is engineered.** `insertions_long_tail` for VD and DJ is uniform over 20-60 nt. It is labelled as not literature-derived, with Briney et al. 2012 as the rationale. It is sized so the `bnab` preset's `long_cdr3_weight` can reach the CDRH3 targets: measured ≥ 24 aa / ≥ 35 aa fractions (reads with a junction) are 0.57 / 0.17 at weight 0.3. The `bnab` preset uses 0.08 and turns on `productive_only`. Over productive reads (selection, and no SHM repair of unproductive rearrangements; record 0024) that gives 0.281 / 0.024 (`tests/validation.rs` `bnab_targets`). Before the no-repair rule, 0.05 gave 0.290 / 0.024: SHM had rescued many long-CDRH3 rearrangements with a stop in their N bases.
- **No D inversion.** `d_inversion` is 0 for IGH and TRB, as the spec requires.
- **Mutation model.** `shm_model` is `motif` for IG and `none` for TR (record 0024).

## Alternatives considered

- **Copy IGoR's or OLGA's bundled model files:** GPL.
- **Use the Volpe & Kepler 2008 tables for IGH D usage and N length:** the tables are exact, but the source is older and smaller, and it would mix sources. The user asked for IGoR first.
- **Raster-digitize the IGoR main figures from the journal JPEGs:** unnecessary, because the article PDF figures are vector.
- **Truncated-normal trims for light chains:** no published mean/SD was found, and inventing an SD is worse than borrowing a published B-cell profile.
- **Uniform N composition for non-IGH loci:** a measured TdT composition from the same enzyme is closer to reality.
- **Leave dropped source genes in the models:** every run would print the same warnings for the built-in models.

## Consequences

- **Long CDR3s.** With the `naive` preset (no long tail), productive CDRH3 mean is about 19.6 aa and about 21 % are ≥ 24 aa; selected naive repertoires are about 15-16 aa.
  - The cause is IGoR's generative insertion means of about 10.5 nt per junction, from nonproductive (pre-selection) reads.
  - `naive` is therefore pre-selection-like, and at long-tail weight 0 it is already near the `bnab` ≥ 24 aa fraction. Presets that need productive-like lengths must account for it.
- **Uncovered genes.** Database genes that no source covers (e.g. IGHV4-31, IGLJ4/5) are never generated. SOURCES.md lists them.
- **One N chain.** The format has one N Markov chain, so DJ N regions use the VD-strand composition, although Elhanati et al. show that DJ follows it on the opposite strand.
- **Fallbacks.** TRA gene usage, light-chain trims, light-chain N length and the TR/light-chain N composition are fallbacks or surrogates. Replace them when better open data are found, with a superseding record.
