Conventions every abstar output follows. Decision records 0004, 0016, 0018 and 0020 hold the reasoning (`abstar manual decisions 0016`).

## Coordinate spaces

- Output coordinates (TSV and Parquet) follow AIRR: 1-based, closed intervals.
- Germline anchor positions (the first base of the V Cys104 codon, or of the J Trp118/Phe118 codon) are 0-based offsets into the ungapped allele sequence. They are stored in the germline database, not in the output.
- Conversion between the two happens at named boundaries only. The `coordinates` attribute of a field in `abstar manual fields` names its space where it matters.

## TSV and Parquet encoding

- TSV: booleans are `T` and `F`, nulls are empty cells, list fields (such as `v_call_ties`) are comma-joined.
- Parquet carries identical content with native types: bool, int, null, and real lists.
- Both writers read the same in-memory record.
- A tab, CR or LF inside a sequence ID or sequence (possible only from in-memory input or malformed files), or in `source_file` (from an input file name), is replaced by a space, so every record stays one line.
- `sequence` is written exactly as read; it is not uppercased. ASCII spaces and tabs inside sequence lines are removed, since whitespace is not sequence content.

## Record statuses

Every input record produces exactly one output row. Records are never dropped, so record conservation is a row-count check. `annotation_status` and `failure_reason` say what happened:

| `annotation_status` | Meaning |
|---|---|
| `annotated` | The record was annotated. |
| `unassigned` | The input is valid but could not be fully annotated. `failure_reason` says why (table below). |
| `invalid_input` | The input is not a usable sequence: non-nucleotide characters (`invalid_characters`), empty (`empty_sequence`), or ID or sequence bytes that are not UTF-8 (`invalid_encoding`). |
| `error` | An internal per-record error. `failure_reason` carries the message, for example `internal_error: <message>`. |

`--strict` turns any `error` row into a fatal run error (exit 3, `STRICT_RECORD_FAILURE`). Fatal run errors (unreadable input, bad database, I/O failure) stop the run with the codes in `abstar manual errors`.

`unassigned` records keep every field the annotator could establish (record 0035):

| `failure_reason` | Meaning | Fields still filled |
|---|---|---|
| `no_v_hit` | No V gene on either strand (random sequence, adapters, a locus the database lacks). | none |
| `no_j_hit` | A V gene but no J gene, for example a read that ends inside CDR3. | the V call, V coordinates, FWR1–FWR3 (and the partial CDR3), V mutations, `sequence_vdj` to the read's end, masks |
| `junction_undefined` | V and J found, but a junction anchor (the Cys104 or Trp/Phe118 codon) is missing or truncated. | every call and every region the read covers; `junction` is null |
| `macaque_not_imgt_gapped` | The V allele's gaps don't follow IMGT numbering (most built-in macaque IGH V genes, until #137; or a custom database without IMGT gaps). Regions, mutations and masks need IMGT numbering. | `rev_comp`, `locus`, the V call and its alignment |

`annotated` records have a junction. Filtering on `annotation_status == "annotated"` keeps only fully annotated rearrangements.

When writing to disk, `run.json` is written last, beside the outputs: its presence marks a finished run. It records parameters, abstar version, database name and content hash, AIRR schema version, and counts by status.

## `sequence_id` rules

- `sequence_id` is the FASTA or FASTQ header's first whitespace-delimited token, byte-for-byte, and is always a string. The one exception: if a record's ID or sequence is not valid UTF-8, the record is `invalid_input` with `failure_reason` `invalid_encoding`, and its `sequence_id` and `sequence` are written with each invalid byte replaced by U+FFFD (output strings must be UTF-8).
- A bare sequence string is one record with ID `"1"`. A list of sequence strings gets IDs `"1"` to `"N"` by position. A list of `(id, sequence)` pairs keeps its IDs. In `abstar.run()`, a string that is neither an existing file nor made only of valid sequence characters raises `INPUT_NOT_FOUND` (0020). On the command line, SEQUENCES is always a file, directory or `-`; use `--seq` for a sequence.
- Valid sequence characters, case-insensitive: `ACGTN` plus the IUPAC codes `RYSWKMBDHV`.
- Functionality values are `F`, `ORF`, `P`, `NF` and `unknown`.

## Input files

- Compression is detected by content: gzip (including multi-member) is read whatever the file is called. bzip2, xz and zstd are not supported; `.bz2`, `.xz` and `.zst` files fail with `INPUT_UNREADABLE`.
- A leading UTF-8 byte-order mark and leading blank lines (LF or CRLF) are skipped. CRLF line endings are accepted.
