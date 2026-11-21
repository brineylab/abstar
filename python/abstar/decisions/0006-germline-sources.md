# 0006: Germline sources: OGRDB first, best curated source otherwise

## Status

Accepted (2026-11-20)

## Context

Germline sets differ in quality and coverage by source and locus. The choice of source should be explicit and follow the best available data over time.

## Decision

OGRDB supplies V/D/J wherever it covers the locus. Otherwise the best curated alternative is used (for example IMGT GENE-DB, KIMDB). IMGT supplies constant regions. The source is recorded per locus in the manifest, and a locus switches to OGRDB when OGRDB adds it. This was a user decision.

Built-in databases: `human` (IG from OGRDB, TR from IMGT), `macaque` (IGH from KIMDB, IGK/IGL from IMGT), `c57bl6`, `balbc`, and `human+c57bl6`. (Amended 2026-11-21: `macaque` is no longer built in, because its V genes are not IMGT-gapped; see 0036 and #137.) Today's `balbc` and `c57bl6` manifests cite the Homo sapiens OGRDB URL; provenance is fixed when they are rebuilt.

## Alternatives considered

- A single source for everything: not possible, as no one source covers every locus and species.

## Consequences

- Provenance is explicit per locus.
- Databases improve as OGRDB coverage grows, by switching source rather than redesigning.
- Rebuilding the mouse databases corrects their recorded provenance.
