# Sources for the built-in human recombination models

The five built-in models (`python/abstar/models/human-{igh,igk,igl,tra,trb}/model.toml`)
are built from published, peer-reviewed summary distributions (decision 0027). Nothing here
comes from IGoR or OLGA software or model files (GPL), and no GPL, AGPL or
NonCommercial-licensed data is used. Every source below is under a Creative Commons
Attribution license: CC BY 4.0, except Briney et al. 2012 (CC BY, version unstated by
PLoS ONE at the time), which supplies no numbers.

Rebuild the models with `python scripts/build_default_models.py`. It reads the digitized
tables in `digitized/` and the built-in `human` germline database. The figure extractors
that produced `digitized/` are in `extract/`, for reference. All files were retrieved
on 2026-11-20.

## Sources

### [IGoR] Marcou, Mora & Walczak 2018

Marcou Q, Mora T, Walczak AM. High-throughput immune repertoire analysis with IGoR.
*Nat Commun* 9:561 (2018). doi:10.1038/s41467-018-02832-w.

- **License.** CC BY 4.0. The article text says: "Open Access This article is licensed
  under a Creative Commons Attribution 4.0 International License ... © The Author(s) 2018."
  The Supplementary Information is published with the article.
- **Attribution.** Values are digitized from figures of Marcou, Mora & Walczak (2018),
  licensed CC BY 4.0 (http://creativecommons.org/licenses/by/4.0/). Changes made:
  - values renormalized;
  - restricted to genes in the abstar human database;
  - DNA/RNA series averaged where noted.
- **Files.**
  - Article PDF: https://www.nature.com/articles/s41467-018-02832-w.pdf
    - sha256 `2b8550ed19fb6fff6db7df1ceb5d6ba30d492c5e870072baf060d24cd98ed587`
  - Supplementary Information `41467_2018_2832_MOESM1_ESM.pdf`:
    https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-018-02832-w/MediaObjects/41467_2018_2832_MOESM1_ESM.pdf
    - sha256 `0019b74920a040bc606228a370785c40f222341f18b2f86c9271b06e42ae423c`
    - The Europe PMC copy (PMC5805751) is byte-identical.
- **Not used.** The journal publishes no supplementary data file with model parameters,
  only the SI PDF and the peer-review file. The IGoR GitHub repository and its model
  files were not used.

### [Elhanati] Elhanati et al. 2015

Elhanati Y, Sethna Z, Marcou Q, Callan CG Jr, Mora T, Walczak AM. Inferring processes
underlying B-cell repertoire diversity. *Phil Trans R Soc B* 370:20140243 (2015).
doi:10.1098/rstb.2014.0243.

- **License.** CC BY 4.0.
- **Files.** Electronic supplementary material `rstb20140243supp1.pdf`.
  - Retrieved from Europe PMC (PMC4528420):
    https://www.ebi.ac.uk/europepmc/webservices/rest/PMC4528420/supplementaryFiles
  - The Royal Society site refused automated download.
  - sha256 `f6d6f5a1a9ecba258212f7886792b7d34807b2eb0923a77df26966253777f3bb`

### [Mikocziova] Mikocziova et al. 2021

Mikocziova I, Peres A, Gidoni M, Greiff V, Yaari G, Sollid LM. Germline polymorphisms
and alternative splicing of human immunoglobulin light chain genes. *iScience*
24:103192 (2021). doi:10.1016/j.isci.2021.103192.

- **License.** CC BY 4.0.
- **Data.** Naive B-cell repertoires of 100 individuals.
- **Files.** High-resolution figures:
  - Fig. 3: https://ars.els-cdn.com/content/image/1-s2.0-S2589004221011603-gr3_lrg.jpg
    - sha256 `1cf37ae9430626e6c50fc333bd02fe57b9504381697a8f8dddb96d70ef0b646c`
  - Fig. 4: https://ars.els-cdn.com/content/image/1-s2.0-S2589004221011603-gr4_lrg.jpg
    - sha256 `608c9cde9fa3e5264d35deaf2b6854cd7da6ceca125fcf10fafcfb929c32a594`

### [Townsend] Townsend et al. 2016

Townsend CL, Laffy JMJ, Wu YB, et al. Significant differences in physicochemical
properties of human immunoglobulin kappa and lambda CDR3 regions. *Front Immunol* 7:388
(2016). doi:10.3389/fimmu.2016.00388.

- **License.** CC BY 4.0.
- **Values used.** Two numbers read from the Results text ("lambda mean = 3.321 ...
  kappa mean = 3.046"; IMGT HighV-QUEST N additions).
- **Files.** Full-text XML: https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5037968/fullTextXML
  - sha256 `5ad01d197ddec2bd0b8f8a1dc9542f0ef2b9d7f43d0e7c6e0a1964a494dc2199`

### [Briney] Briney, Willis & Crowe 2012

Briney BS, Willis JR, Crowe JE Jr. Human peripheral blood antibodies with long HCDR3s
are established primarily at original recombination using a limited subset of germline
genes. *PLoS ONE* 7:e36750 (2012). doi:10.1371/journal.pone.0036750.

- **License.** CC BY.
- **Use.** Cited only as the rationale for the IGH long tail; no numbers were taken.

## Extraction methods

**Vector digitization** (IGoR and Elhanati figures).

- **How values were read.**
  - The figures in both PDFs are vector graphics.
  - Each page was converted to SVG with `pdftocairo -svg` (`extract/svgpaths.py`).
  - Bar heights and polyline vertices were read in page coordinates, by series colour.
  - They were mapped to data values with the axis tick marks.
- **Interpolation.** matplotlib drops collinear vertices, so polylines are linearly
  interpolated at integer x.
- **Precision.** Precision is effectively exact to the plotted value, about 1e-4
  absolute. As a check, distributions whose full support is plotted sum to 1.000 ± 0.0002.
  - IGoR Fig. 2a-c: every TRB/TRA curve.
  - IGoR Supplementary Fig. 4: every J column.
  - IGoR Supplementary Fig. 7: the D deletion curves.
  - IGoR Supplementary Fig. 3: the TRBJ and TRBD bars.
- **Truncated ranges.** Supplementary Fig. 13 plots only insertions 0-40 and deletions
  −5..15, so those tables cover the plotted range and are renormalized. Mass outside it
  was 1-6 %.

**Raster digitization** (Mikocziova Figs. 3-4; `extract/boxmedians.py`).

- **How medians were read.**
  - Each box-plot median is the thickest dark horizontal run inside its box, at the
    axis tick of its gene.
  - Medians were calibrated against the major gridlines (about 5,200 px per unit for V,
    about 2,000-4,500 px per unit for J).
- **Precision.** About ±0.003, worse for boxes collapsed near zero.
- **Gene order.** Read from the rotated axis labels.

**Read from text.** The Townsend means and the IGoR Supplementary Fig. 5 per-panel labels
("P(V) = 0.08", 2 decimals) were read from the text.

## Per locus

**Fallback** below means a documented fallback, not a published distribution.

### human-igh

- **v_usage.** [IGoR] Supplementary Fig. 13, IGHV panels, IGoR series (60 genes).
  - Dropped, not in the database: IGHV1-38-4, IGHV1-NL1, IGHV3-22, IGHV3-52, IGHV3-71,
    IGHV4-55, IGHV7-34-1.
  - Retained V mass: 0.921.
  - IGHV3-NL1 has the largest share (about 0.08). This reflects IGoR's IMGT-era
    assignment of reads to that gene in the original data, not a curation choice here.
- **j_usage.** [IGoR] Supplementary Fig. 13, IGHJ panel.
- **d_given_j.** [IGoR] Supplementary Fig. 4, P(D|J) for 27 IGHD genes × IGHJ1-6.
  - Dropped: IGHD4-11, which is identical to IGHD4-4 and absent from the database.
  - Retained D mass: 0.998.
- **trim v_3, d_5, d_3, j_5.** [IGoR] Supplementary Fig. 13, IGoR series, −5..15.
  Negative values are P-nucleotides.
- **insertions vd, dj.** [IGoR] Supplementary Fig. 13, IGoR series, 0..40.
- **n_markov.**
  - start: [Elhanati] Fig. S9, "VD sense" mononucleotide frequencies.
  - Rows: Fig. S10, VD-sense dinucleotide conditionals.
  - Both are means of individuals A and B.
- **insertions_long_tail vd, dj.**
  - Engineered, not literature: uniform over 20-60 nt.
  - Rationale: [Briney] Fig. 3A shows that long HCDR3s come mainly from long N
    additions.
  - The `bnab` preset mixes the tail in with probability `long_cdr3_weight` per junction.
  - Measured with this model (20,000 reads, no SHM or indels), fraction of reads with
    CDRH3 ≥ 24 aa / ≥ 35 aa:

    | `long_cdr3_weight` | ≥ 24 aa | ≥ 35 aa |
    |---|---|---|
    | 0 | 0.236 | 0.005 |
    | 0.3 | 0.566 | 0.166 |
    | 0.5 | 0.732 | 0.301 |

    The targets (≥ 0.25 and ≥ 0.02) are reachable at small weights.
- **Length caveat.** The `naive` preset (no long tail, mutation 0) gives a productive
  CDRH3 mean of about 19.6 aa, with about 21 % of productive reads ≥ 24 aa. Measured
  values: 19.66 aa and 20.7 % at n = 20,000, seed 0; an independent review measured
  21.5 %.
  - Selected naive repertoires are about 15-16 aa. The cause is IGoR's generative
    insertion means of about 10.5 nt per junction, inferred from nonproductive reads
    before selection.
  - The `naive` preset is therefore pre-selection-like. At long-tail weight 0 it is
    already near the `bnab` ≥ 24 aa target.

### human-trb

- **v_usage, j_usage.** [IGoR] Supplementary Fig. 3.
  - Allele bars are summed to genes.
  - For each gene, the mean over the 2 DNA individuals and the mean over the 5 mRNA
    individuals are averaged.
  - Dropped: TRBV1, TRBV12-1, TRBV12-2, TRBV17, TRBV23-1, TRBV26, TRBV3-2, TRBV5-3,
    TRBV5-7, TRBV6-7, TRBV7-1, mostly ORF or pseudogenes absent from the database.
  - Retained V mass: 0.924.
- **d_given_j.** [IGoR] Supplementary Fig. 14, IGoR P(D2|J). P(D1|J) = 1 − P(D2|J). This
  keeps the TRBD2-TRBJ1 exclusion.
- **trim v_3, j_5; insertions vd, dj.** [IGoR] Fig. 2a-c, mean of the "TRB DNA" and "TRB
  RNA" curves.
- **trim d_5, d_3.** [IGoR] Supplementary Fig. 7, "true model" curves. This is IGoR's
  data-inferred TRB model, used there as ground truth.
- **n_markov.** Surrogate: the IGH composition from [Elhanati]. No TRB N composition was
  found in these sources, and TdT is the same enzyme.

### human-tra

- **v_usage, j_usage.** Partial literature plus fallback.
  - [IGoR] Supplementary Fig. 5 panel labels give the 6 most used TRAV alleles: TRAV19
    0.08, TRAV13-1 0.07 + 0.06, TRAV27 0.06, TRAV13-2 0.05, TRAV4 0.04.
  - They give the 6 most used TRAJ alleles at 0.04 each: TRAJ45, TRAJ49, TRAJ37, TRAJ53,
    TRAJ42, TRAJ48.
  - **Fallback:** the remaining mass (V 0.64, J 0.76) is uniform over the other
    generatable database genes (40 V, 43 J).
- **trim v_3, j_5; insertions vj.** [IGoR] Fig. 2a-c, "TRA RNA" curves.
- **n_markov.** Surrogate: the IGH composition from [Elhanati].

### human-igk and human-igl

- **v_usage, j_usage.** [Mikocziova] Fig. 4 (IGK) and Fig. 3 (IGL), box-plot medians,
  renormalized.
  - **"E" genes.** The paper labels identical duplicates "E" (e.g. IGKV1E-12). They are
    split equally over the database copies that exist: IGKV1-12/1D-12, 1-13/1D-13,
    6-21/6D-21. Otherwise the one copy present takes the mass: 2-28, 1D-33, 1-39, 2D-40.
  - **Dropped.** IGKV1E-37 (IGKV1-37/1D-37) has no copy in the database. Retained V mass:
    0.996.
  - **Floor.** Boxes drawn at zero get a floor of 0.001: deleted or "Unknown" in most
    donors, e.g. IGKV2-29, IGKV2D-26, IGKV1D-17, IGKV3D-7, IGLV5-52, IGLV3-22, IGLV4-3.
  - **Pseudogenes.** The paper's suspected pseudogenes IGLJ6 and IGLV3-12 are omitted.
- **insertions vj.** Fallback shape. A geometric distribution with the [Townsend] mean N
  addition, kappa 3.046 and lambda 3.321 nt, truncated where p < 1e-4 and renormalized.
- **trim v_3, j_5.** Fallback surrogate: the IGH V 3′ and J 5′ deletion profiles from
  [IGoR] Supplementary Fig. 13, truncated to ≤ 6 nt and renormalized.
  - **Why truncate.** Every IGKJ/IGLJ *01 allele has its F118 anchor at ungapped
    index 7 (IGKJ2 at 8), while IGHJ anchors sit at 18 or later. The untruncated IGH J
    profile has P(trim ≥ 8) ≈ 0.25.
  - **Effect.** It deleted the anchor in about 23 % of IGK and 25 % of IGL reads (no
    junction). After truncation, the no-junction rate is 0.000 for both at mutation 0,
    measured on 20,000 reads. For comparison: IGH 0.036, TRA 0.035, TRB 0.020.
    `crates/generate/tests/builtin_models.rs` asserts < 0.08.
  - **Why the brief's truncated normal was not used.** No published light-chain trim
    mean/SD was found. Collins & Watson 2018 (*Front Immunol* 9:2249, CC BY) state only
    that light-chain exonuclease removal is small.
  - No open, peer-reviewed human light-chain trim distribution was found.
  - Foster et al. 1997 (*J Clin Invest*) and Farner et al. 1999 (*J Immunol*) report
    them, but they could not be retrieved, and Farner et al. is not open.
- **n_markov.** Surrogate: the IGH composition from [Elhanati].

## Database fit

| Model | Genes not in the source but in the database |
|---|---|
| IGH | IGHV1-69-2, IGHV2-70D, IGHV3-25, IGHV3-38, IGHV3-62, IGHV3-64D, IGHV4-30-4, IGHV4-31, IGHV4-NL1, IGHV8-51-1 |
| IGK | IGKV1-NL1, IGKV1D-42, IGKV2-18, IGKV2D-18, IGKV2D-24, IGKV6D-41 |
| IGL | IGLV1-41, IGLV1-50, IGLV11-55, IGLV2-33, IGLV3-12, IGLV3-13, IGLV5-39, IGLV5-48; IGLJ4, IGLJ5, IGLJ6 |
| TRB | TRBV6-1 (IGoR reports it with TRBV6-2/6-3) |
| TRA | none (fallback covers every gene) |

- **Retained source mass.** `scripts/build_default_models.py` asserts that each segment
  keeps ≥ 0.90 of its source mass after genes absent from the database are dropped.
  Current values: IGH V 0.921, IGH D 0.998, TRB V 0.924, IGK V 0.996, all others 1.0.
  IGH V and TRB V are the closest to the threshold.
- Genes dropped from a source are listed above, and in each `model.toml` header comment.
  Every model binds to the built-in database without warnings
  (`crates/generate/tests/builtin_models.rs`).
- The genes in the table are not generated, because the model gives them no mass.
- IGHV and IGHD names follow OGRDB, and IGHV3-NL1 is shared with IGoR's IMGT-era name.
