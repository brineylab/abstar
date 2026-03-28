# Copyright (c) 2024 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

import pytest
from abutils import Sequence

from abstar.annotation.germline_alignment import get_germline, get_germline_database_path

# ----------------------------
#      DATABASE PATHS
# ----------------------------


def test_get_germline_database_path():
    path = get_germline_database_path(germdb_name="human", receptor="bcr")
    assert path


def test_get_germline_database_name_invalid():
    with pytest.raises(FileNotFoundError):
        get_germline_database_path(
            germdb_name="NotAGermlineDatabase", receptor="bcr"
        )


def test_get_germline_database_receptor_invalid():
    with pytest.raises(ValueError):
        get_germline_database_path(germdb_name="human", receptor="abc")


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


def test_get_single_germline_nonexistent():
    with pytest.raises(ValueError):
        get_germline(
            germline_gene="IGHV1-1*01",
            germdb_name="human",
            receptor="bcr",
            exact_match=True,
        )


def test_get_single_germline_nonunique():
    with pytest.raises(ValueError):
        get_germline(
            germline_gene="IGHV1-2",
            germdb_name="human",
            receptor="bcr",
            exact_match=True,
        )
