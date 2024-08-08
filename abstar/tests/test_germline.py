# Copyright (c) 2024 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

import pytest
from abutils import Sequence

from ..annotation.germline import get_germline, get_germline_database_path

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
