# Copyright (c) 2026 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

"""The AIRR 2.0 Rearrangement TSV serialization boundary.

Annotation/Parquet intervals are zero-based and half-open. Conversion here is
one-way; callers must pass annotation rows, never already serialized AIRR rows.
"""

import csv
from collections.abc import Mapping
from itertools import groupby
from pathlib import Path

import abutils
import polars as pl

from .schema import OUTPUT_SCHEMA


AIRR_SCHEMA_VERSION = "2.0"
AIRR_REQUIRED_FIELDS = (
    "sequence_id", "sequence", "rev_comp", "productive", "v_call", "d_call",
    "j_call", "sequence_alignment", "germline_alignment", "junction",
    "junction_aa", "v_cigar", "d_cigar", "j_cigar",
)
AIRR_FIELDS = AIRR_REQUIRED_FIELDS + tuple(
    field for field in OUTPUT_SCHEMA
    if field not in AIRR_REQUIRED_FIELDS and field != "row_id"
)


def to_airr_interval(start: int | None, end: int | None) -> tuple[int | None, int | None]:
    """Convert an internal half-open interval to a 1-based closed interval."""
    if start is None or end is None:
        return None, None
    if type(start) is not int or type(end) is not int or start < 0 or end <= start:
        raise ValueError(f"invalid half-open interval: [{start}, {end})")
    return start + 1, end


def build_cigar(aligned_query: str, aligned_germline: str, *, query_start: int, germline_start: int) -> str:
    """Encode retained columns, with oriented-query clips and reference skips."""
    if len(aligned_query) != len(aligned_germline):
        raise ValueError("aligned query and germline must have equal lengths")
    if any(type(start) is not int or start < 0 for start in (query_start, germline_start)):
        raise ValueError("alignment starts must be nonnegative integer offsets")
    operations = []
    for query_base, germline_base in zip(aligned_query, aligned_germline):
        if query_base == germline_base == "-":
            raise ValueError("an alignment column cannot contain two gaps")
        operations.append("D" if query_base == "-" else "I" if germline_base == "-" else "M")
    prefix = (f"{query_start}S" if query_start else "") + (f"{germline_start}N" if germline_start else "")
    return prefix + "".join(
        f"{sum(1 for _ in group)}{operation}" for operation, group in groupby(operations)
    )


def _validate_frame(frame: int) -> None:
    if type(frame) is not int or frame not in (1, 2, 3):
        raise ValueError("coding frame must be an integer in (1, 2, 3)")


def translate_airr_query(
    oriented_query: str | None, *, query_start: int | None, frame: int | None,
) -> str | None:
    """Translate the full oriented query in the retained V coding phase.

    ``frame`` is one-based within the retained V query. Its first complete
    codon begins at ``query_start + frame - 1`` in the oriented input, so the
    full-query phase is that offset modulo three. Leading/terminal partial
    codons are omitted. Missing frame/origin/query or no complete codon is null.
    """
    if oriented_query is None or query_start is None or frame is None:
        return None
    _validate_frame(frame)
    if type(query_start) is not int or not 0 <= query_start < len(oriented_query):
        raise ValueError("V query start must address the oriented query")
    full_query_frame = (query_start + frame - 1) % 3 + 1
    return abutils.tl.translate(oriented_query, frame=full_query_frame) or None


def translate_airr_alignment(
    aligned_query: str | None, aligned_germline: str | None, *, frame: int | None,
) -> tuple[str | None, str | None]:
    """Translate paired retained NT columns in one shared coding window.

    Skip the leading partial query codon, then group the same three alignment
    columns in both rows. All-gap codons become '-', all-dot numbering spacers
    become '.', and mixed gaps or ambiguous codons become 'X'. The latter
    explicitly marks disrupted codon columns; no missing NP residues are
    inferred. Omit a final incomplete column triplet symmetrically. These are
    translations of the aligned rows, not independent ungapped translations.
    """
    if aligned_query is None or aligned_germline is None or frame is None:
        return None, None
    _validate_frame(frame)
    if len(aligned_query) != len(aligned_germline):
        raise ValueError("aligned query and germline must have equal lengths")
    if any(q == g == "-" for q, g in zip(aligned_query, aligned_germline)):
        raise ValueError("an alignment column cannot contain two gaps")
    start, remaining = 0, frame - 1
    while remaining and start < len(aligned_query):
        remaining -= aligned_query[start] not in "-."
        start += 1
    end = start + (len(aligned_query) - start) // 3 * 3
    if remaining or start == end:
        return None, None
    return tuple(abutils.tl.translate(row[start:end], allow_dots=True)
                 for row in (aligned_query, aligned_germline))


def to_airr_row(row: Mapping[str, object]) -> dict[str, object]:
    """Copy public annotations and map original sequence and AIRR coordinates."""
    output = {field: row.get(field) for field in AIRR_FIELDS}
    output["sequence"] = row.get("sequence_input")
    prefixes = [f"{segment}_{axis}" for segment in ("v", "d", "j", "c")
                for axis in ("sequence", "germline")]
    prefixes += ["fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4"]
    for prefix in prefixes:
        start_field, end_field = f"{prefix}_start", f"{prefix}_end"
        output[start_field], output[end_field] = to_airr_interval(
            row.get(start_field), row.get(end_field)
        )
    frame = row.get("v_frame")
    if frame is None:
        frame = row.get("frame")
    output["sequence_aa"] = translate_airr_query(
        row.get("sequence_oriented"), query_start=row.get("v_sequence_start"), frame=frame,
    )
    output["sequence_alignment_aa"], output["germline_alignment_aa"] = translate_airr_alignment(
        row.get("sequence_alignment"), row.get("germline_alignment"), frame=frame,
    )
    return output


def write_airr_tsv(frame: pl.DataFrame, path: str | Path) -> None:
    """Write T/F booleans, empty nulls and LF lines, rejecting TSV delimiters."""
    rows = []
    for source in frame.iter_rows(named=True):
        encoded = {}
        for field, value in to_airr_row(source).items():
            if isinstance(value, bool):
                value = "T" if value else "F"
            elif value is None:
                value = ""
            elif any(delimiter in str(value) for delimiter in ("\t", "\n", "\r")):
                raise ValueError(f"AIRR field {field!r} contains a forbidden delimiter")
            encoded[field] = value
        rows.append(encoded)
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=AIRR_FIELDS, delimiter="\t",
                                lineterminator="\n", quoting=csv.QUOTE_MINIMAL, doublequote=True)
        writer.writeheader()
        writer.writerows(rows)
