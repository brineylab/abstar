# 0002: No external binaries (MMseqs, Parasail dropped)

## Status

Accepted (2026-11-20)

## Context

Today's pipeline shells out to MMseqs and Parasail, which brings subprocess startup, temp files and chunking on the hot path. Accuracy comes first and speed a very close second, and both were limited by this structure.

## Decision

The engine has no external binaries and no temp files on the hot path. Search and alignment are purpose-built in Rust. Germline sets are a few hundred alleles, so a k-mer prefilter followed by exhaustive vectorized alignment is feasible. Fast vector instructions (AVX2 on x86_64, NEON on ARM) are selected at runtime with a portable fallback.

## Alternatives considered

- Keep MMseqs and Parasail and wrap them better: does not remove the temp-file and IPC overhead or the structural bug sources.
- Call a Rust binary as a subprocess: the same pattern in different clothes (see 0001).

## Consequences

- Installation is a single wheel with nothing else to install.
- abstar owns the search and alignment code and its accuracy.
- The throughput targets rest on the engine, so a throwaway performance spike confirms them before anything is built on them.
