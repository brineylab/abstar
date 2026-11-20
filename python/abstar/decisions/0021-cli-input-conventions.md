# 0021: CLI input conventions: stdin and --seq

## Status

Accepted (2026-11-20)

## Context

The CLI must read stdin and annotate a single sequence given on the command line, in pipelines where agents cannot rely on prompts. `abstar run` has two optional positionals (SEQUENCES, OUTPUT_DIR), so `--seq` makes a lone positional ambiguous. 0003 and 0020 define how strings are classified as path or sequence.

## Decision

**Stdin.** The string `-` means stdin, in both `abstar run -` and `abstar.run("-")` (as an exception to the path-or-sequence rule in 0003 and 0020). `-` is never a valid sequence, so nothing is lost. The core reads fd 0 directly.

**`--seq`.** `--seq TEXT` is one sequence with ID `"1"` and is mutually exclusive with SEQUENCES. A lone positional given with `--seq` is OUTPUT_DIR (`abstar run --seq ACGT out`). If that positional is an existing non-directory (such as a FASTA file) it is rejected with `INVALID_ARGUMENT`, because the caller almost certainly meant it as SEQUENCES. Two positionals plus `--seq` is also `INVALID_ARGUMENT`. Nothing is written in either case.

## Alternatives considered

- Require `--seq` to be used without any positional, so it can never write files: it would make single sequences impossible to save.
- Silently treat the positional as OUTPUT_DIR always: `abstar run reads.fasta --seq ACGT` would create a directory named `reads.fasta`.
- A separate `--stdin` flag: `-` is the universal convention.

## Consequences

- Misuse fails fast with a stable code and exit status 2.
- `abstar.run("-")` is not a sequence input; callers wanting a literal `-` have no valid use for it.
