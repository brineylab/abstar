# Copyright (c) 2025 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

"""Integrity checks for packaged germline databases."""

import json
from pathlib import Path

import pytest

from abstar.core.germline import get_germline_database_path

from .helpers import read_fasta_records


COUNT_SNAPSHOT_PATH = Path(__file__).parent / "data" / "germline_counts.json"
MMSEQS_SUFFIXES = (
    "",
    ".dbtype",
    ".index",
    ".lookup",
    ".source",
    "_h",
    "_h.dbtype",
    "_h.index",
)


def database_cases():
    if not COUNT_SNAPSHOT_PATH.is_file():
        return []
    counts = json.loads(COUNT_SNAPSHOT_PATH.read_text(encoding="utf-8"))
    return [
        pytest.param(receptor, database, segment, expected_count, id=case_id)
        for receptor, databases in counts.items()
        for database, segments in databases.items()
        for segment, expected_count in segments.items()
        for case_id in [f"{receptor}-{database}-{segment}"]
    ]


def test_reviewed_germline_count_snapshot_exists():
    assert COUNT_SNAPSHOT_PATH.is_file()


def test_read_fasta_records_rejects_duplicate_identifiers(tmp_path):
    fasta = tmp_path / "duplicates.fasta"
    fasta.write_text(
        ">IGHV1-1*01\nACGT\n>IGHV1-1*01\nTGCA\n",
        encoding="utf-8",
    )

    with pytest.raises(
        ValueError, match=r"duplicate FASTA identifier 'IGHV1-1\*01'"
    ):
        read_fasta_records(fasta)


@pytest.mark.parametrize(
    "receptor,database,segment,expected_count", database_cases()
)
def test_packaged_germline_database_integrity(
    receptor, database, segment, expected_count, monkeypatch, tmp_path
):
    monkeypatch.setenv("HOME", str(tmp_path))
    database_path = Path(
        get_germline_database_path(germdb_name=database, receptor=receptor)
    )
    ungapped = read_fasta_records(database_path / "ungapped" / f"{segment}.fasta")
    gapped = read_fasta_records(database_path / "imgt_gapped" / f"{segment}.fasta")

    assert len(ungapped) == expected_count
    assert set(gapped) == set(ungapped)
    assert len(gapped) == len(set(gapped))
    assert {key: value.replace(".", "") for key, value in gapped.items()} == ungapped
    assert (database_path / "manifest.txt").is_file()
    for suffix in MMSEQS_SUFFIXES:
        assert (database_path / "mmseqs" / f"{segment}{suffix}").is_file()

    expected_prefix = "IG" if receptor == "bcr" else "TR"
    assert all(identifier.startswith(expected_prefix) for identifier in ungapped)
    if segment == "d":
        assert {identifier[:3] for identifier in ungapped} <= {"IGH", "TRB", "TRD"}
