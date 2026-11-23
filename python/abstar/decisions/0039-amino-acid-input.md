# 0039: Amino-acid input

## Status

Accepted (2026-11-22)

## Context

Many antibody and TCR sequences exist only as protein: structures, patents, therapeutics, published mAbs. Until now the only option was to back-translate them, which invents codons, mutations and D calls. Users with these datasets need V and J calls, regions, CDR3 and amino-acid mutations from the protein itself.

The output must stay an AIRR Rearrangement table with the same columns, shared writers, `run.json`, Python API and CLI. Nucleotide-mode values must not change.

## Decision

- **Mode.** One run is one sequence type: `sequence_type` `nt` (default) or `aa` (`abstar run --sequence-type`, `abstar.run(..., sequence_type=)`), recorded in `run.json` under `parameters`. There is no per-record detection. Amino-acid mode runs a second annotator, `AaEngine`, beside the nucleotide `Engine`; the pipeline, writers and schema are shared, so the nucleotide path is untouched.
- **Input.** Uppercased. Allowed: the 20 standard residues, `X`, `B`, `Z`, `J`, `U`, `O` and `*`. Anything else is `invalid_characters`; empty is `empty_sequence`.
- **Germline translation.** Built once from the existing nucleotide alleles, honoring `loci` and `include_nonfunctional`. V is translated from the IMGT-gapped sequence codon by codon, a `...` codon becoming `.`, so gapped residue `k` is IMGT position `k + 1`; a V allele that is not IMGT-gapped is excluded. J is translated from its anchor's frame, so the Trp/Phe118 residue is `anchor / 3`; a J allele without an anchor is excluded. Non-ACGT bases translate to `X`. The germline database format is unchanged.
- **Fields.** V and J calls, locus, IMGT regions and lengths, `junction_aa`, `cdr3_aa`, `sequence_vdj_aa`, the amino-acid alignments, mutations, counts and frequencies, masks, `stop_codon` and `productive` are filled. `sequence_aa` is the input. `productive` needs a junction, no stop, `junction_aa` starting with `C` and ending with `W` or `F`, and no `X`. Segment masks use `V`, `N` and `J` only.
- **Null in amino-acid mode.** `sequence`, `rev_comp`, `vj_in_frame`, every nucleotide field (cigars, coordinates, regions, `junction`, `np1`/`np2`, `sequence_vdj`, alignments, masks, mutations, indels), and all `d_*` and `c_*` fields. A field whose AIRR meaning is nucleotide coordinates is never reused for residue positions.
- **Indel fields.** Four new `abstar` fields, `v_insertions_aa`, `v_deletions_aa`, `j_insertions_aa` and `j_deletions_aa`, hold amino-acid indels as `IMGT{pos}:ins:RESIDUES` / `IMGT{pos}:del:RESIDUES` (V) and `{pos}:ins:...` / `{pos}:del:...` (J). They are null in nucleotide mode.
- **No D and no C.** A D contributes a few residues, so a call from them would be noise. C calls are deferred.
- **Schema.** `sequence` becomes nullable in every output.

## Alternatives considered

- **Thread an alphabet parameter through `Engine`, `assign` and `derive`.** Nearly every nucleotide function carries codon or strand logic (orientation, the 2-bit k-mer index, SIMD sentinels, frame and Cys projection, D placement). Each would grow a branch, and each branch is a chance to change nucleotide results.
- **Back-translate and run the nucleotide engine.** It invents codons, so nucleotide mutations, D calls, NP regions and allele calls would all be artifacts.
- **Auto-detect protein input.** Short nucleotide and protein strings can be ambiguous, so the mode is always explicit.
- **A per-mode schema with `sequence` required in nucleotide mode.** It would break the single-sourced `fields.toml` (record 0012).

## Consequences

- **Nucleotide Parquet.** `sequence` is now an optional column in nucleotide-mode Parquet and Arrow schemas too. Values are unchanged, but a reader that checks schemas strictly sees the column as no longer required. Four new `*_aa` indel columns, null in nucleotide mode, are the only other change to nucleotide output.
- **Not AIRR-valid.** Amino-acid output has a null `sequence`, so it does not pass `airr` validation. This is documented, not tested.
- **Heavy chains.** With no D residues, the CDR3 middle is non-germline in `germline_alignment_aa` and `gene_segment_mask_aa`, unlike nucleotide mode.
- **Performance.** Alignment is scalar against the distinct V amino-acid sequences, fine for antibody panels. A protein k-mer prefilter and SIMD are the upgrade for large repertoires. Measured: about 560 sequences/s on the 2,000-read `igh_memory` amino-acid set with default threads (32).
- **Thresholds.** `MIN_V_SCORE_AA = 50` and `MIN_J_SCORE_AA = 15`, gap −11/−1, calibrated on the amino-acid accuracy gate by true accuracy, leaning permissive: V 50–100 scored the same, higher V and J lost true calls, and the negative controls (non-immunoglobulin proteins, random proteins, a nucleotide read) stay `no_v_hit`, scoring under 40.
