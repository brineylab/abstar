# 0016: Plan-level conventions

## Status

Accepted (2026-11-20)

## Context

Several conventions cut across the foundation and must be fixed once so that every component agrees.

## Decision

- A record's `sequence_id` is the FASTA/FASTQ header's first whitespace-delimited token, byte-for-byte. It is always a string.
- Valid sequence characters, case-insensitive: `ACGTN` plus IUPAC `RYSWKMBDHV`.
- `sequence` is written exactly as read; it is not uppercased.
- Exit codes: `0` success; `1` fatal run error; `2` invalid usage or input; `3` strict-mode record failure.
- Germline anchor positions are 0-based offsets into the ungapped allele sequence: the first base of the V Cys104 codon, or of the J Trp118/Phe118 codon.
- Functionality values are `F`, `ORF`, `P`, `NF` and `unknown`.

## Alternatives considered

The design leaves exact exit code values to the builder; this record fixes them. No other alternatives are recorded in the design.

## Consequences

- Exact ID preservation and no silent normalization of input.
- Coordinate spaces are explicit and converted at named boundaries.
- Scripts and agents can rely on the exit-code table.
