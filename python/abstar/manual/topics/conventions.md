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
- For nucleotide input, `sequence` is written exactly as read; it is not uppercased (amino-acid input is uppercased into `sequence_aa`; see Amino-acid input). ASCII spaces and tabs inside sequence lines are removed, since whitespace is not sequence content.

## Record statuses

Every input record produces exactly one output row. Records are never dropped, so record conservation is a row-count check. `annotation_status` and `failure_reason` say what happened:

| `annotation_status` | Meaning |
|---|---|
| `annotated` | The record was annotated. |
| `unassigned` | The input is valid but could not be fully annotated. `failure_reason` says why (table below). |
| `invalid_input` | The input is not a usable sequence: characters outside the valid set, non-nucleotide for nucleotide input and non-amino-acid for amino-acid input (`invalid_characters`), empty (`empty_sequence`), or ID or sequence bytes that are not UTF-8 (`invalid_encoding`). |
| `error` | An internal per-record error. `failure_reason` carries the message, for example `internal_error: <message>`. |

`--strict` turns any `error` row into a fatal run error (exit 3, `STRICT_RECORD_FAILURE`). Fatal run errors (unreadable input, bad database, I/O failure) stop the run with the codes in `abstar manual errors`.

`unassigned` records keep every field the annotator could establish (record 0035):

| `failure_reason` | Meaning | Fields still filled |
|---|---|---|
| `no_v_hit` | No V gene (on either strand, for nucleotide input) (random sequence, adapters, a locus the database lacks). | none |
| `no_j_hit` | A V gene but no J gene, for example a read that ends inside CDR3. | the V call, V coordinates, FWR1–FWR3 (and the partial CDR3), V mutations, `sequence_vdj` to the read's end, `gene_segment_mask` and `nongermline_mask`; `cdr_mask` only when V reaches the Cys104 codon. These regions run to the read's end, so the partial CDR3, `cdr_mask` and `sequence_vdj` can include J or C bases the J search missed. |
| `junction_undefined` | V and J found, but a junction anchor (the Cys104 or Trp/Phe118 codon) is missing or truncated. | every call and every region the read covers; `junction` is null |

`annotated` records have a junction. Filtering on `annotation_status == "annotated"` keeps only fully annotated rearrangements.

When writing to disk, `run.json` is written last, beside the outputs: its presence marks a finished run. It records parameters, abstar version, database name and content hash, AIRR schema version, and counts by status.

## `sequence_id` rules

- `sequence_id` is the FASTA or FASTQ header's first whitespace-delimited token, byte-for-byte, and is always a string. The one exception: if a record's ID or sequence is not valid UTF-8, the record is `invalid_input` with `failure_reason` `invalid_encoding`, and its `sequence_id` and `sequence` (`sequence_aa`, uppercased, for amino-acid input) are written with each invalid byte replaced by U+FFFD (output strings must be UTF-8).
- A bare sequence string is one record with ID `"1"`. A list of sequence strings gets IDs `"1"` to `"N"` by position. A list of `(id, sequence)` pairs keeps its IDs. In `abstar.run()`, a string that is neither an existing file nor made only of valid sequence characters raises `INPUT_NOT_FOUND` (0020). On the command line, SEQUENCES is always a file, directory or `-`; use `--seq` for a sequence.
- Valid nucleotide sequence characters, case-insensitive: `ACGTN` plus the IUPAC codes `RYSWKMBDHV`. Amino-acid input has its own set (Amino-acid input).
- Functionality values are `F`, `ORF`, `P`, `NF` and `unknown`.

## Input files

- Compression is detected by content: gzip (including multi-member) is read whatever the file is called. bzip2, xz and zstd are not supported; `.bz2`, `.xz` and `.zst` files fail with `INPUT_UNREADABLE`.
- A leading UTF-8 byte-order mark and leading blank lines (LF or CRLF) are skipped. CRLF line endings are accepted.

## Amino-acid input

Protein sequences (structures, patents, published antibodies) are annotated without back-translation with `abstar run --sequence-type aa` or `abstar.run(..., sequence_type="aa")` (decision 0039). The value is case-insensitive, and one run is one sequence type: there is no per-record detection.

- Valid characters, any case: the 20 standard residues `ACDEFGHIKLMNPQRSTVWY` plus `XBZJUO` and `*` (stop). Anything else, including `-` and `.`, is `invalid_input` with `invalid_characters`.
- `sequence` is null; the input, uppercased, is in `sequence_aa`, on every record (failures included).
- Filled: V and J calls and ties, `locus`, `v_score`/`v_identity` and `j_score`/`j_identity` (BLOSUM62 over residues), the `_aa` regions and their `_length`s, `junction_aa`, `sequence_vdj_aa`, `sequence_alignment_aa`, `germline_alignment_aa`, the `_aa` masks, mutations, counts, frequencies and indels (`v_insertions_aa` and the like), `stop_codon` and `productive`. `abstar manual fields` says how each differs from nucleotide mode.
- Null: every nucleotide field (coordinates, cigars, regions and their `_start`/`_end`, `junction` and `junction_length`, `np1`/`np2`, `sequence_vdj`, the nucleotide alignments, masks, mutations and indels), `rev_comp` and `vj_in_frame`.
- No D and no C calls: all `d_*` and `c_*` fields are null. For heavy chains the CDR3 middle is non-germline in `germline_alignment_aa` and `gene_segment_mask_aa`.
- Alleles that differ only by synonymous codons translate to the same protein, so they always tie (`v_call_ties`, `j_call_ties`).
- Germline proteins are translated from the database's nucleotide alleles. V alleles that are not IMGT-gapped, or that have a codon mixing gap and bases, and J alleles without an anchor, are not searched. A custom database with ungapped V alleles loses them in amino-acid mode; gap them first (`abstar germline gap`, or `gap_v = true` in a recipe).
- `unassigned` records keep the amino-acid counterparts of the fields listed under Record statuses.
- The output is not valid AIRR: `sequence` is required by the AIRR schema and is null here, so `airr` validation fails.
