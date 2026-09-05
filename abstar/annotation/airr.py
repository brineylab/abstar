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
                                lineterminator="\n", quoting=csv.QUOTE_NONE, quotechar=None)
        writer.writeheader()
        writer.writerows(rows)
