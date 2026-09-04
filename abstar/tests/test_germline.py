# Copyright (c) 2024 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

from types import SimpleNamespace

import pytest
from abutils import Sequence

from ..annotation.antibody import Antibody
from ..annotation.germline import (
    get_germline,
    get_germline_database_path,
    process_cgene_alignment,
    process_dgene_alignment,
    process_jgene_alignment,
    reassign_dgene,
)


def _alignment(**kwargs):
    defaults = {
        "query": "ABCDEFGHIJKL",
        "target": "abcdefghijkl",
        "query_begin": 0,
        "query_end": 0,
        "target_begin": 0,
        "target_end": 0,
        "score": 10,
    }
    return SimpleNamespace(**(defaults | kwargs))


def test_j_germline_end_uses_semiglobal_target_coordinate_once():
    ab = Antibody(sequence_id="j-coordinate")
    semiglobal = _alignment(query_begin=2, query_end=8, target_begin=3, target_end=9)
    local = _alignment(query_begin=1, query_end=4, target_begin=2, target_end=5)

    process_jgene_alignment("X" * 30, 10, semiglobal, local, ab)

    assert ab.j_germline_start == 5
    assert ab.j_germline_end == 10
    assert ab.j_germline == "fghij"


def test_d_germline_end_uses_target_span_when_alignment_has_indel():
    ab = Antibody(sequence_id="d-coordinate", d_call="IGHD1*01")
    local = _alignment(query_begin=1, query_end=3, target_begin=2, target_end=6)

    process_dgene_alignment("X" * 30, 10, local, ab)

    assert (ab.d_sequence_start, ab.d_sequence_end) == (11, 14)
    assert (ab.d_germline_start, ab.d_germline_end) == (2, 7)
    assert ab.d_germline == "cdefg"


def test_c_sequence_end_does_not_double_count_query_offset():
    ab = Antibody(sequence_id="c-coordinate")
    semiglobal = _alignment(query_begin=2, query_end=8, target_begin=3, target_end=9)
    local = _alignment(query_begin=1, query_end=4, target_begin=2, target_end=5)

    process_cgene_alignment("X" * 30, 10, semiglobal, local, ab)

    assert (ab.c_sequence_start, ab.c_sequence_end) == (13, 17)
    assert len(ab.c_sequence) == 4

# ----------------------------
#      DATABASE PATHS
# ----------------------------


def test_get_germline_database_path():
    path = get_germline_database_path(germdb_name="human", receptor="bcr")
    assert path


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


@pytest.mark.parametrize(
    "locus,sequence,expected_prefix",
    [
        ("TRB", "GGGACAGGGGGC", "TRBD"),
        ("TRD", "ACTGGGGGATACG", "TRDD"),
    ],
)
def test_reassign_dgene_uses_tcr_locus_database(locus, sequence, expected_prefix):
    alignment = reassign_dgene(
        sequence=sequence,
        germdb_name="human",
        locus=locus,
        receptor="tcr",
    )

    assert alignment is not None
    assert alignment.target.id.startswith(expected_prefix)


@pytest.mark.parametrize("locus", ["TRA", "TRG", "IGK", "IGL"])
def test_reassign_dgene_skips_loci_without_d_genes(locus):
    assert (
        reassign_dgene(
            sequence="GGGACAGGGGGC",
            germdb_name="human",
            locus=locus,
            receptor="tcr" if locus.startswith("TR") else "bcr",
        )
        is None
    )


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
