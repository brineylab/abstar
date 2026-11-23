# 0038: Non-functional alleles are searched and drawn only on request

## Status

Accepted (2026-11-23)

## Context

Until 1.0.0, the annotator searched every allele in the database, and `abstar generate` drew from every allele its model gives usage to. Functionality only broke ties: F > ORF > P > NF > unknown (records 0008, 0029).

A tie-break cannot help when a non-functional allele wins outright. In the 1.0.0rc2 run over all 2,372,992 `jaffe-paired` chains, 7 kappa chains were called IGKV2-29\*01. OGRDB marks that allele `NF`: it reads `TGA` (stop) at Cys104. The reads have `TGT` there, so they mismatch both \*01 and \*02 (`TGC`) at that codon. \*01 won only because it also matches the reads at IMGT 159. Because \*01 has no Cys104 anchor (record 0019), all 7 chains came out `junction_undefined`. Legacy abstar called each of them with a junction.

Functionality labels are incomplete (record 0019). OGRDB alleles are `F` or `NF`, and IMGT `|` headers give `F`, `ORF` or `P`. TR alleles, every macaque allele and every constant region are `unknown`. In `human`, 12 IG V alleles are `NF`, and no D, J or C allele is `NF` or `P`. `balbc` has 3 `NF` V alleles and `c57bl6` has 7.

## Decision

- **Default.** The annotator leaves out alleles whose functionality is `P` or `NF` (`Functionality::is_nonfunctional`), whatever the segment, when it builds each segment's search index (`SegmentIndex::for_loci`). `F`, `ORF` and `unknown` alleles are always searched. Excluded alleles are absent from calls and from `*_call_ties`.
- **Generator.** `abstar generate` leaves the same alleles out of the alleles it draws from and of the candidates its truth tie sets are computed over (`bind`, `bind_alleles`). Truth and annotation therefore cover the same alleles. A model gene whose alleles are all excluded is dropped without a warning, because the run chose to leave it out. IGHV3-19, whose only allele \*01 is `NF`, is one. Its mass still counts against the 90% coverage check (record 0023).
- **Opt-in.** `include_nonfunctional` restores 1.0.0's behavior: `abstar run --include-nonfunctional` / `abstar.run(..., include_nonfunctional=True)`, and the same flag and option for `abstar generate`. It is recorded in `run.json` and `generate.json` under `parameters`.
- Tie-breaking (record 0008) is unchanged. With the opt-in it still prefers F > ORF > P > NF > unknown.

## Alternatives considered

- **Search only `F` alleles.** This would empty the TR, macaque and constant-region searches, which are all `unknown`, and drop ORF alleles. ORF alleles have defects in recombination or splicing signals, can still be expressed, and appear in functional repertoires.
- **Leave non-functional alleles out of the built-in databases.** Then every user would need a second database to include them. A per-run option keeps one database, and a custom database built from a recipe keeps every allele.
- **Anchor IGKV2-29\*01 by IMGT position despite its stop** (#138). This would give these 7 chains a junction, but they would keep a call to an allele that cannot encode a functional chain, and it fixes one allele rather than the class.
- **Exclude only in the annotator.** The generator would keep drawing NF alleles (1.1% of `human-igh` reads, 0.2% of `human-igk`), which the annotator could never call. Those reads would become guaranteed misses on the accuracy gate.

## Consequences

- **Real data.** Of the 2,372,992 Jaffe chains, 7 had an NF primary V call. All were IGKV2-29\*01, and all now go to IGKV2-29\*02 (6) or IGKV2D-29\*02 (1) with legacy's junction, productive. Another 1,049 chains had an NF allele only among their V ties. They lose it from `v_call_ties`, and their primary call and junction are unchanged. No other chain is affected. With #162's fix, 1,349 of the 1,356 chains that rc1 left `junction_undefined` now have legacy's junction.
- **Simulation.** Default `abstar generate` output changes for any model that gives usage to an NF allele. For `human-igh` the sequences drawn are different, so the pinned neutral-output hashes changed. With `include_nonfunctional` the output is 1.0.0's byte for byte, and a test pins this.
- **Accuracy gate.** The gate's reads changed with the generator default. On identical reads, the annotator's default and opt-in results are the same in every set, except for one heavily mutated igh_bnab_c read (69% V identity). There, removing NF alleles changed which genes passed the k-mer prefilter's 90%-of-leader cut (record 0029), and the call moved to another functional gene. The floors were rewritten on the new reads with `ABSTAR_WRITE_FLOORS` (measured − 0.5 pp): 949 went up and 336 down. 48 dropped by more than 1 pp, nearly all in subsets of about 100 reads (`ambiguous_j`, `has_indel`), which is within resampling noise at that size.
- **Graded on the 1.0.0 reads** (which include NF truth), the new default costs about 1% `v_call` in IGH, on reads whose truth is an NF allele. Users annotating genomic-DNA repertoires, which contain nonproductive pseudogene rearrangements, or doing genotyping or novel-allele work should pass `--include-nonfunctional`. Without it, such reads are called to the nearest functional allele. Productivity is still judged on the read, so a stop codon in the read still makes it nonproductive.
