# Copyright (c) 2025 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

"""Shared helpers for the abstar test suite."""

from collections.abc import Mapping, Sequence
from pathlib import Path


def assert_same_annotations(
    left: Sequence[Mapping[str, object]],
    right: Sequence[Mapping[str, object]],
    fields: Sequence[str],
) -> None:
    """Compare exact values and types in input order; missing fields fail."""
    assert len(left) == len(right), f"record count: {len(left)} != {len(right)}"
    for index, (expected, actual) in enumerate(zip(left, right)):
        for field in fields:
            assert field in expected and field in actual, f"row {index}: missing {field}"
            assert type(expected[field]) is type(actual[field]), (
                f"row {index}, {field}: types differ "
                f"({type(expected[field]).__name__}, {type(actual[field]).__name__})"
            )
            assert expected[field] == actual[field], (
                f"row {index}, {field}: {expected[field]!r} != {actual[field]!r}"
            )


def read_fasta_records(path: Path) -> dict[str, str]:
    """Read FASTA records, rejecting duplicate identifiers."""
    records = {}
    identifier = None
    sequence = []

    with path.open(encoding="utf-8") as fasta:
        for raw_line in fasta:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    records[identifier] = "".join(sequence)
                identifier = line[1:].split()[0]
                if identifier in records:
                    raise ValueError(
                        f"duplicate FASTA identifier '{identifier}' in {path}"
                    )
                sequence = []
            else:
                sequence.append(line)

    if identifier is not None:
        records[identifier] = "".join(sequence)
    return records


def assert_airr_matches_annotations(expected, path, fields):
    """Compare through AIRR's explicit sequence, null and coordinate boundary.

    The reference reader inversely converts known AIRR starts to Python offsets.
    Public extensions are parsed using their declared annotation types.
    """
    import airr
    import polars as pl
    from abstar.annotation.schema import OUTPUT_SCHEMA

    reader = airr.read_rearrangement(str(path), validate=True)
    try:
        actual = list(reader)
    finally:
        reader.close()
    assert all('row_id' not in row for row in actual)
    normalized_expected, normalized_actual = [], []
    for rows, destination, is_airr in ((expected, normalized_expected, False),
                                       (actual, normalized_actual, True)):
        for row in rows:
            normalized = {}
            for field in fields:
                assert field in row
                value = row['sequence_input'] if field == 'sequence' and not is_airr else row[field]
                if value == '' or value is None:
                    value = None
                elif is_airr:
                    dtype = OUTPUT_SCHEMA[field]
                    if dtype == pl.Boolean and not isinstance(value, bool):
                        assert value in ('T', 'F')
                        value = value == 'T'
                    elif dtype == pl.Int64:
                        value = int(value)
                    elif dtype == pl.Float64:
                        value = float(value)
                normalized[field] = value
            destination.append(normalized)
    assert_same_annotations(normalized_expected, normalized_actual, fields)
