# Copyright (c) 2024 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

import os
from types import SimpleNamespace

import pytest
from abutils import Sequence

from ..annotation.antibody import Antibody
from ..annotation.germline import (
    get_germline,
    get_germline_database_path,
    process_cgene_alignment,
)

# ----------------------------
#      DATABASE PATHS
# ----------------------------


def test_get_germline_database_path():
    path = get_germline_database_path(germdb_name="human", receptor="bcr")
    assert path


@pytest.mark.parametrize("receptor", ["bcr", "tcr"])
def test_get_germline_database_path_prefers_custom_database(tmp_path, monkeypatch, receptor):
    monkeypatch.setenv("HOME", str(tmp_path))
    custom_db_path = (
        tmp_path / ".abstar" / "germline_dbs" / receptor / "custom_germdb"
    )
    custom_db_path.mkdir(parents=True)

    path = get_germline_database_path(germdb_name="custom_germdb", receptor=receptor)

    assert path == os.fspath(custom_db_path)


@pytest.mark.xfail(
    reason="a database named NotAGermlineDatabase does not exist",
    raises=FileNotFoundError,
)
def test_get_germline_database_name_invalid():
    path = get_germline_database_path(
        germdb_name="NotAGermlineDatabase", receptor="bcr"
    )
    assert path


@pytest.mark.xfail(
    reason="receptor type 'abc' is invalid (must be 'bcr' or 'tcr')",
    raises=ValueError,
)
def test_get_germline_database_receptor_invalid():
    path = get_germline_database_path(germdb_name="human", receptor="abc")
    assert path


# ----------------------------
#      ALIGNMENT HELPERS
# ----------------------------


@pytest.mark.parametrize(
    "semiglobal_query_begin,semiglobal_target_begin,local_query_begin,local_query_end",
    [
        (6, 2, 3, 11),
        (2, 5, 4, 10),
    ],
)
def test_process_cgene_alignment_uses_local_query_span_for_end(
    semiglobal_query_begin,
    semiglobal_target_begin,
    local_query_begin,
    local_query_end,
):
    oriented_input = "A" * 200
    semiglobal_aln = SimpleNamespace(
        query_begin=semiglobal_query_begin,
        target_begin=semiglobal_target_begin,
        target="A" * 200,
    )
    local_aln = SimpleNamespace(
        score=10.0,
        query_begin=local_query_begin,
        query_end=local_query_end,
        target_begin=1,
        target_end=20,
    )
    ab = Antibody()

    process_cgene_alignment(
        oriented_input=oriented_input,
        j_sequence_end=10,
        semiglobal_aln=semiglobal_aln,
        local_aln=local_aln,
        ab=ab,
    )

    expected_span = (local_query_end - local_query_begin) + 1
    assert ab.c_sequence_end - ab.c_sequence_start == expected_span
    assert len(ab.c_sequence) == expected_span


# ----------------------------
#        GET GERMLINE
# ----------------------------


def test_get_single_germline():
    germ = get_germline(
        germline_gene="IGHV1-2*02",
        germdb_name="human",
        receptor="bcr",
        exact_match=True,
    )
    assert isinstance(germ, Sequence)
    assert germ.id == "IGHV1-2*02"


def test_get_multiple_germlines():
    germs = get_germline(
        germline_gene="IGHV1-2",
        germdb_name="human",
        receptor="bcr",
        exact_match=False,
    )
    assert len(germs) >= 2
    assert all([isinstance(germ, Sequence) for germ in germs])
    assert all([germ.id.startswith("IGHV1-2") for germ in germs])


def test_get_single_germline_tcr():
    germ = get_germline(
        germline_gene="TRAV1-1*01",
        germdb_name="human",
        receptor="tcr",
        exact_match=True,
    )
    assert isinstance(germ, Sequence)
    assert germ.id == "TRAV1-1*01"


def test_get_multiple_germlines_tcr():
    germs = get_germline(
        germline_gene="TRAV1-1",
        germdb_name="human",
        receptor="tcr",
        exact_match=False,
    )
    assert len(germs) >= 2
    assert all([isinstance(germ, Sequence) for germ in germs])
    assert all([germ.id.startswith("TRAV1-1") for germ in germs])


@pytest.mark.xfail(
    reason="gene 'IGHV1-1*01' does not exist in the human bcr database",
    raises=ValueError,
)
def test_get_single_germline_nonexistent():
    germ = get_germline(
        germline_gene="IGHV1-1*01",
        germdb_name="human",
        receptor="bcr",
        exact_match=True,
    )
    assert isinstance(germ, Sequence)
    assert germ.id == "IGHV1-1*01"


@pytest.mark.xfail(
    reason="gene 'IGHV1-2' does not identically match anything in the human bcr database",
    raises=ValueError,
)
def test_get_single_germline_nonunique():
    germ = get_germline(
        germline_gene="IGHV1-2",
        germdb_name="human",
        receptor="bcr",
        exact_match=True,
    )
    assert isinstance(germ, Sequence)
    assert germ.id == "IGHV1-2"
