# 0036: A germline database must be IMGT-gapped

## Status

Accepted (2026-11-21)

## Context

Regions, mutations and masks are numbered by IMGT gapped position (records 0032–0034), and the V anchor is the Cys104 codon at IMGT 310–312 (0017). A V gene set whose gaps are not IMGT's yields plausible, wrong annotations. The built-in macaque V genes (KIMDB) have `.` gaps, but only 8 of 1,071 have TGT/TGC at gapped positions 309..312 (0-based). Record 0035 handled this with a macaque-specific status. Ruling (2026-11-21): a germline database without IMGT gaps is not a germline database. It must arrive IMGT-gapped, or be gapped by an abstar tool (#152, deferred); the annotator makes no accommodation for an incomplete or improperly gapped one.

## Decision

`load` and `write` (so also `build`) refuse a database that has V alleles unless:

- every V allele's `sequence_gapped` contains `.`, and
- at least 95% (`MIN_ANCHORED_V`, `crates/germline/src/validate.rs`) of its V alleles have TGT/TGC at gapped 0-based `309..312` (`find_v_anchor`).

The error names the database, the failing count, up to five example allele names, and says IMGT-gapped V sequences are required. A database with no V alleles skips the check.

The rule is database-level, not per allele: real IMGT-gapped alleles lack Cys104 (a mutated codon, or a partial sequence: human IGKV2-29*01, IGLV3-2*03, TRBV7-3*02/03, IGHV3-20*02). Built-in coverage is human 577/588, c57bl6 193/196, balbc 263/265, against macaque 8/1,071 as the KIMDB/IMGT sources came (before gapping; with `abstar germline gap` it is now 1,053 of 1,064 V alleles anchored, and 1,028 of 1,030 with an ORF).

*Amended 2026-11-21 (record 0037, #137): ungapped databases can now be made valid with `abstar germline gap` or the recipe option `gap_v`, and macaque ships again that way. The text it replaces is kept: "Macaque is no longer shipped. Its recipe and `germline-sources/macaque/` stay for #137, and `scripts/build_builtin_germlines.py` skips it."*

~~Macaque is no longer shipped. Its recipe and `germline-sources/macaque/` stay for #137, and `scripts/build_builtin_germlines.py` skips it.~~ *Superseded by the 2026-11-21 amendment above (record 0037, #137): macaque ships again, V genes gapped by abstar.*

## Alternatives considered

- **Per-allele refusal:** rejects the real human alleles above.
- **A species or status special case (record 0035's):** the annotator would carry accommodations for databases it should refuse.
- **A Cys-codon search fallback to regap:** a heuristic can choose the wrong codon (0017); regapping belongs to #152.

## Consequences

- Every V allele of a loaded database is gapped, so the annotator has no non-IMGT V path; a V allele without an anchor takes the normal path (`junction_undefined`).
- `abstar run -g macaque` fails with `DATABASE_NOT_FOUND` until #137 ships a regapped database (shipped 2026-11-21, see amendment above).
- A custom database with ungapped V sequences fails at `build`, not at annotation.
