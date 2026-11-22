# 0009: Donor keep-lists restrict only the locus/segment pairs they name

## Status

Accepted (2026-11-20)

## Context

Donor-specific databases are built from a base database, a keep-list of alleles (the donor's genotype), and extra FASTA of novel alleles (for example from IgDiscover or TIgGER). A genotype often covers only some loci.

## Decision

The keep-list restricts only the locus/segment combinations it names. A genotype listing IGHV alleles restricts IGHV and leaves IGHD, IGHJ, light chains and C untouched. The build prints what was restricted and what passed through. A donor database should build in seconds. This was a call made in the design.

## Alternatives considered

- Treat a keep-list as the complete database, dropping everything not named: would silently delete loci the genotype never addressed.

## Consequences

- Partial genotypes do not remove unrelated loci.
- The printed restricted/passed-through report makes the effect auditable.
