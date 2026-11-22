# 0019: Built-in databases are built from recipes over the legacy sources

## Status

Accepted (2026-11-20)

## Context

The five built-in databases (`human`, `c57bl6`, `balbc`, `macaque`, `human+c57bl6`) must ship in the new format (0005) and contain exactly the sequences the legacy databases contained, so the rewrite changes how germlines are stored and anchored, not which germlines are used. Their sources follow 0006. The legacy `balbc` and `c57bl6` manifests cited the *Homo sapiens* OGRDB URL, which was wrong.

## Decision

*Amended 2026-11-21 (record 0037, #137): `macaque` is built in again. Its V genes are IMGT-gapped by `abstar germline gap` (recipe option `gap_v`) against human, c57bl6 and balbc; 7 of 1,071 V alleles cannot be gapped (framework insertion) and are dropped. The old note above is kept: "`macaque` is no longer built in, because its V genes are not IMGT-gapped and such a database is refused." The built-in list is `human`, `c57bl6`, `balbc`, `macaque`, `human+c57bl6`.*

*Amended 2026-11-21 (record 0036): `macaque` is no longer built in, because its V genes are not IMGT-gapped and such a database is refused. The built-in list is `human`, `c57bl6`, `balbc`, `human+c57bl6`. The recipe and `germline-sources/macaque/` stay for #137. The macaque paragraphs below describe what shipped before.*

Each built-in database is built by the recipe builder (0017) from `germline-sources/recipes/<name>.toml`. Its source files are copied verbatim from `legacy/abstar/germline_dbs/` into `germline-sources/<db>/`. `scripts/build_builtin_germlines.py [NAME ...]` rebuilds them into `python/abstar/germline_dbs/` with `overwrite=True` and prints every warning. The generated `manifest.toml` and `alleles.tsv` are committed and ship in the wheel.

Sources by database and locus:

- `human` (species `human`). IGH, IGK and IGL V, D and J come from OGRDB *Homo sapiens* IGH_VDJ_extended rev 9, IGKappa_VJ_extended rev 4 and IGLambda_VJ_extended rev 3 (all released 2024-10-12). IG constant regions come from IMGT GENE-DB (legacy `bcr/human/imgt_gapped/c.fasta`). TRA, TRB, TRG and TRD V, D, J and C come from IMGT GENE-DB (downloaded 2025-03-27; F, ORF and in-frame P, IMGT-gapped), from legacy `tcr/human/imgt_gapped/{v,d,j,c}.fasta`. The legacy `__homo_sapiens` header suffix is stripped and the sequences are unchanged.
- `c57bl6` (species `mouse`). This uses the OGRDB *Mus musculus* sets: C57BL/6 IGH rev 3 (2022-06-02), C57BL/6J IGKV rev 1 and C57BL/6J IGLV rev 1 (2022-05-03), and IGKJ (all strains) and IGLJ (all strains) rev 1 (2022-05-31). Constant regions come from IMGT GENE-DB *Mus musculus* (downloaded 2024-08-19).
- `balbc` (species `mouse`). This uses the OGRDB *Mus musculus* sets: BALB/c IGH rev 3 (2022-06-02), BALB/c/ByJ IGKV rev 1 (2022-05-03), BALB/c/ByJ IGLV rev 2 (2024-04-20), and the same all-strains IGKJ/IGLJ sets as `c57bl6`. Constant regions are the same IMGT *Mus musculus* file as `c57bl6`.
- `macaque` (species `macaque`). IGH V, D and J come from KIMDB 1.1. IGK and IGL V and J come from IMGT GENE-DB. IGH, IGK and IGL constant regions come from IMGT GENE-DB. All were downloaded 2024-08-19, and the files are legacy `bcr/macaque/imgt_gapped/{v,d,j,c}.fasta`.
- `human+c57bl6`. This combines the `human` IG sources (no TR) with the `c57bl6` sources, using prefixes `human = "hs"` and `mouse = "mm"`.

Recipe `description`s cite these per-locus sources. The mouse databases now cite the *Mus musculus* OGRDB sets, which corrects the legacy provenance error noted in 0006.

Sequence parity with legacy is a test. For each database and segment, the set of (name, ungapped sequence) equals legacy `ungapped/<seg>.fasta`, with `__<species>` suffixes stripped. Legacy `human` is the union of `bcr/human` and `tcr/human`, and `human+c57bl6` compares (species, name, sequence).

**Functionality.** TR, macaque and all IG constant-region alleles have functionality `unknown`, because their source FASTAs carry no functionality field (no IMGT `|` header). OGRDB alleles are `F` or `NF`.

**Macaque V anchors are mostly missing.** The legacy macaque gapped V FASTA does not follow IMGT gap numbering. Only 8 of 1071 macaque V alleles (all IGK) have `TGT`/`TGC` at IMGT positions 310-312: IGH has 0 of 774, IGK 8 of 169 and IGL 0 of 128. The other 1063 get no anchor and a named warning. No heuristic fallback is added and nothing is regapped. The fix is to regap the macaque V genes with the annotator's aligner (spec build mode 2) once that exists.

**Accepted build warnings** (all recorded in each manifest):

- OGRDB/motif J disagreements, where the OGRDB anchor is kept (0017). `c57bl6`, `balbc` and `human+c57bl6` (as `mm|`) each have IGKJ0-5NUZ\*00, IGKJ0-IIHF\*00 and IGLJ0-UWMQ\*00, all `motif anchor none disagrees with OGRDB j_cdr3_end 8`. The two IGKJ alleles have the non-canonical FSDG motif. These come from the shared all-strains sets.
- V alleles without a Cys104 codon at 310-312:
  - human IGKV2-29\*01 (NF; stop codon at 104) and IGLV3-2\*03 (OGRDB functional; the OGRDB gapping is shifted).
  - human TRAV1-2\*02, TRAV13-1\*03, TRAV14/DV4\*04, TRAV26-2\*02, TRAV30\*04, TRAV6\*06 and TRAV8-1\*02, whose 3' ends are truncated before position 310.
  - human TRBV7-3\*02 and \*03, which have CGT (Arg) at 104.
  - c57bl6 IGHV0-KAOG\*00 (NF), IGHV0-XN3Q\*00 (NF; TCT) and IGHV0-OSZZ\*00 (OGRDB functional; the OGRDB gapping is shifted).
  - balbc IGHV0-NUDM\*00 (OGRDB functional; the OGRDB gapping is shifted) and IGHV0-UDF7\*00 (OGRDB functional; the gapped sequence ends at 305).
- J alleles without the Trp/Phe118 motif:
  - human TRAJ16\*01, TRAJ16\*02 (FARG), TRAJ35\*01 (CGSG), TRBJ2-7\*02 (VGPG), TRGJP1\*01 (FAEG) and TRGJP2\*01 (FAKG).
  - macaque IGLJ4\*01 and IGLJ7\*02.

## Alternatives considered

- Copy the legacy `ungapped/` and `imgt_gapped/` outputs instead of rebuilding from the raw sources: this would carry legacy's derived files forward with no recorded provenance or checksums, and the recipe builder would not be exercised on real data.
- Add a Cys-codon search fallback for macaque V, or regap macaque now: a heuristic can choose the wrong codon (0017). Regapping belongs to the aligner-based build mode, not this task.
- Hand-fix the shifted OGRDB gapping (IGLV3-2\*03, IGHV0-OSZZ\*00, IGHV0-NUDM\*00): this edits curated sources by hand. The warnings make the alleles visible until the aligner can anchor them.
- Generate the databases at install time: build-time generation needs the sources in the sdist and makes the wheel depend on build success. Committed outputs are inspectable and diffable.

## Consequences

- Built-ins are reproducible: `python scripts/build_builtin_germlines.py` regenerates them from committed sources; only the manifest `built` date and `abstar_version` can change.
- About 99% of macaque V alleles have no anchor, so junction/CDR3 identification for macaque must not depend on the V anchor until the V genes are regapped.
- V and J alleles without anchors, including a few that OGRDB marks functional, are listed in each manifest's warnings for review.
- Updating a source means replacing the file under `germline-sources/`, editing the recipe's description, rebuilding, and reviewing the parity and count tests deliberately. A legitimate source change changes the legacy-parity expectation, and that change needs its own decision record.
