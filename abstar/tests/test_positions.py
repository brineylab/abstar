# Copyright (c) 2024 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

import pytest

from ..annotation.positions import (
    get_gapped_position_from_raw,
    get_raw_position_from_aligned,
    get_raw_position_from_gapped,
)


@pytest.fixture
def gapped_sequence():
    return "ATGCATGC....ATGCATGC"


@pytest.fixture
def gapped_sequence_leading_gap():
    return "....ATGCATGCATGC"


@pytest.fixture
def gapped_sequence_trailing_gap():
    return "ATGCATGCATGC...."


@pytest.fixture
def gapped_sequence_multiple_gaps():
    return "ATGC....ATGC....ATGC"


@pytest.fixture
def aligned_sequence_no_gaps():
    return "ATGCATGCATGC"


@pytest.fixture
def aligned_sequence_with_gap():
    return "ATGC----ATGC"


@pytest.fixture
def aligned_sequence_leading_gap():
    return "----ATGCATGC"


@pytest.fixture
def aligned_sequence_trailing_gap():
    return "ATGCATGC----"


@pytest.fixture
def aligned_sequence_multiple_gaps():
    return "ATGC--GC--GC"


# ----------------------------
#      GAPPED FROM RAW
# ----------------------------


def test_get_gapped_position_from_raw(gapped_sequence):
    assert get_gapped_position_from_raw(1, gapped_sequence) == 1
    assert get_gapped_position_from_raw(2, gapped_sequence) == 2
    assert get_gapped_position_from_raw(3, gapped_sequence) == 3
    assert get_gapped_position_from_raw(4, gapped_sequence) == 4
    assert get_gapped_position_from_raw(5, gapped_sequence) == 5
    assert get_gapped_position_from_raw(6, gapped_sequence) == 6
    assert get_gapped_position_from_raw(7, gapped_sequence) == 7
    assert get_gapped_position_from_raw(8, gapped_sequence) == 8
    assert get_gapped_position_from_raw(9, gapped_sequence) == 13
    assert get_gapped_position_from_raw(10, gapped_sequence) == 14
    assert get_gapped_position_from_raw(11, gapped_sequence) == 15
    assert get_gapped_position_from_raw(12, gapped_sequence) == 16
    assert get_gapped_position_from_raw(13, gapped_sequence) == 17
    assert get_gapped_position_from_raw(14, gapped_sequence) == 18
    assert get_gapped_position_from_raw(15, gapped_sequence) == 19
    assert get_gapped_position_from_raw(16, gapped_sequence) == 20


def test_get_gapped_position_from_raw_leading_gap(gapped_sequence_leading_gap):
    assert get_gapped_position_from_raw(1, gapped_sequence_leading_gap) == 5
    assert get_gapped_position_from_raw(2, gapped_sequence_leading_gap) == 6
    assert get_gapped_position_from_raw(3, gapped_sequence_leading_gap) == 7
    assert get_gapped_position_from_raw(4, gapped_sequence_leading_gap) == 8
    assert get_gapped_position_from_raw(5, gapped_sequence_leading_gap) == 9
    assert get_gapped_position_from_raw(6, gapped_sequence_leading_gap) == 10
    assert get_gapped_position_from_raw(7, gapped_sequence_leading_gap) == 11
    assert get_gapped_position_from_raw(8, gapped_sequence_leading_gap) == 12
    assert get_gapped_position_from_raw(9, gapped_sequence_leading_gap) == 13
    assert get_gapped_position_from_raw(10, gapped_sequence_leading_gap) == 14
    assert get_gapped_position_from_raw(11, gapped_sequence_leading_gap) == 15
    assert get_gapped_position_from_raw(12, gapped_sequence_leading_gap) == 16


def test_get_gapped_position_from_raw_trailing_gap(gapped_sequence_trailing_gap):
    assert get_gapped_position_from_raw(1, gapped_sequence_trailing_gap) == 1
    assert get_gapped_position_from_raw(2, gapped_sequence_trailing_gap) == 2
    assert get_gapped_position_from_raw(3, gapped_sequence_trailing_gap) == 3
    assert get_gapped_position_from_raw(4, gapped_sequence_trailing_gap) == 4
    assert get_gapped_position_from_raw(5, gapped_sequence_trailing_gap) == 5
    assert get_gapped_position_from_raw(6, gapped_sequence_trailing_gap) == 6
    assert get_gapped_position_from_raw(7, gapped_sequence_trailing_gap) == 7
    assert get_gapped_position_from_raw(8, gapped_sequence_trailing_gap) == 8
    assert get_gapped_position_from_raw(9, gapped_sequence_trailing_gap) == 9
    assert get_gapped_position_from_raw(10, gapped_sequence_trailing_gap) == 10
    assert get_gapped_position_from_raw(11, gapped_sequence_trailing_gap) == 11
    assert get_gapped_position_from_raw(12, gapped_sequence_trailing_gap) == 12


def test_get_gapped_position_from_raw_multiple_gaps(gapped_sequence_multiple_gaps):
    assert get_gapped_position_from_raw(1, gapped_sequence_multiple_gaps) == 1
    assert get_gapped_position_from_raw(2, gapped_sequence_multiple_gaps) == 2
    assert get_gapped_position_from_raw(3, gapped_sequence_multiple_gaps) == 3
    assert get_gapped_position_from_raw(4, gapped_sequence_multiple_gaps) == 4
    assert get_gapped_position_from_raw(5, gapped_sequence_multiple_gaps) == 9
    assert get_gapped_position_from_raw(6, gapped_sequence_multiple_gaps) == 10
    assert get_gapped_position_from_raw(7, gapped_sequence_multiple_gaps) == 11
    assert get_gapped_position_from_raw(8, gapped_sequence_multiple_gaps) == 12
    assert get_gapped_position_from_raw(9, gapped_sequence_multiple_gaps) == 17
    assert get_gapped_position_from_raw(10, gapped_sequence_multiple_gaps) == 18
    assert get_gapped_position_from_raw(11, gapped_sequence_multiple_gaps) == 19
    assert get_gapped_position_from_raw(12, gapped_sequence_multiple_gaps) == 20


# ----------------------------
#      RAW FROM GAPPED
# ----------------------------


def test_get_raw_position_from_gapped(gapped_sequence):
    assert get_raw_position_from_gapped(1, gapped_sequence) == 1
    assert get_raw_position_from_gapped(2, gapped_sequence) == 2
    assert get_raw_position_from_gapped(3, gapped_sequence) == 3
    assert get_raw_position_from_gapped(4, gapped_sequence) == 4
    assert get_raw_position_from_gapped(5, gapped_sequence) == 5
    assert get_raw_position_from_gapped(6, gapped_sequence) == 6
    assert get_raw_position_from_gapped(7, gapped_sequence) == 7
    assert get_raw_position_from_gapped(8, gapped_sequence) == 8
    assert get_raw_position_from_gapped(9, gapped_sequence) == 8
    assert get_raw_position_from_gapped(10, gapped_sequence) == 8
    assert get_raw_position_from_gapped(11, gapped_sequence) == 8
    assert get_raw_position_from_gapped(12, gapped_sequence) == 8
    assert get_raw_position_from_gapped(13, gapped_sequence) == 9
    assert get_raw_position_from_gapped(14, gapped_sequence) == 10
    assert get_raw_position_from_gapped(15, gapped_sequence) == 11
    assert get_raw_position_from_gapped(16, gapped_sequence) == 12


def test_get_raw_position_from_gapped_with_sequence_start(gapped_sequence):
    # test sequence_start = 1 at multiple positions
    assert get_raw_position_from_gapped(1, gapped_sequence, sequence_start=1) == 0
    assert get_raw_position_from_gapped(2, gapped_sequence, sequence_start=1) == 1
    assert get_raw_position_from_gapped(3, gapped_sequence, sequence_start=1) == 2
    assert get_raw_position_from_gapped(4, gapped_sequence, sequence_start=1) == 3
    assert get_raw_position_from_gapped(5, gapped_sequence, sequence_start=1) == 4
    assert get_raw_position_from_gapped(6, gapped_sequence, sequence_start=1) == 5
    assert get_raw_position_from_gapped(7, gapped_sequence, sequence_start=1) == 6
    assert get_raw_position_from_gapped(8, gapped_sequence, sequence_start=1) == 7
    assert get_raw_position_from_gapped(9, gapped_sequence, sequence_start=1) == 7
    assert get_raw_position_from_gapped(10, gapped_sequence, sequence_start=1) == 7
    assert get_raw_position_from_gapped(11, gapped_sequence, sequence_start=1) == 7
    assert get_raw_position_from_gapped(12, gapped_sequence, sequence_start=1) == 7
    assert get_raw_position_from_gapped(13, gapped_sequence, sequence_start=1) == 8
    assert get_raw_position_from_gapped(14, gapped_sequence, sequence_start=1) == 9
    assert get_raw_position_from_gapped(15, gapped_sequence, sequence_start=1) == 10
    assert get_raw_position_from_gapped(16, gapped_sequence, sequence_start=1) == 11
    # test increasign sequence_start values
    assert get_raw_position_from_gapped(8, gapped_sequence, sequence_start=2) == 6
    assert get_raw_position_from_gapped(8, gapped_sequence, sequence_start=3) == 5
    assert get_raw_position_from_gapped(8, gapped_sequence, sequence_start=4) == 4
    # test increasing values at a gapped position
    assert get_raw_position_from_gapped(12, gapped_sequence, sequence_start=2) == 6
    assert get_raw_position_from_gapped(12, gapped_sequence, sequence_start=3) == 5
    assert get_raw_position_from_gapped(12, gapped_sequence, sequence_start=4) == 4


def test_get_raw_position_from_gapped_with_germline_start(gapped_sequence):
    # test germline_start = 1 at multiple positions
    assert get_raw_position_from_gapped(1, gapped_sequence, germline_start=1) == 0
    assert get_raw_position_from_gapped(2, gapped_sequence, germline_start=1) == 1
    assert get_raw_position_from_gapped(3, gapped_sequence, germline_start=1) == 2
    assert get_raw_position_from_gapped(4, gapped_sequence, germline_start=1) == 3
    assert get_raw_position_from_gapped(5, gapped_sequence, germline_start=1) == 4
    assert get_raw_position_from_gapped(6, gapped_sequence, germline_start=1) == 5
    assert get_raw_position_from_gapped(7, gapped_sequence, germline_start=1) == 6
    assert get_raw_position_from_gapped(8, gapped_sequence, germline_start=1) == 7
    assert get_raw_position_from_gapped(9, gapped_sequence, germline_start=1) == 7
    assert get_raw_position_from_gapped(10, gapped_sequence, germline_start=1) == 7
    assert get_raw_position_from_gapped(11, gapped_sequence, germline_start=1) == 7
    assert get_raw_position_from_gapped(12, gapped_sequence, germline_start=1) == 7
    assert get_raw_position_from_gapped(13, gapped_sequence, germline_start=1) == 8
    assert get_raw_position_from_gapped(14, gapped_sequence, germline_start=1) == 9
    assert get_raw_position_from_gapped(15, gapped_sequence, germline_start=1) == 10
    assert get_raw_position_from_gapped(16, gapped_sequence, germline_start=1) == 11
    # test increasign germline_start values
    assert get_raw_position_from_gapped(8, gapped_sequence, germline_start=2) == 6
    assert get_raw_position_from_gapped(8, gapped_sequence, germline_start=3) == 5
    assert get_raw_position_from_gapped(8, gapped_sequence, germline_start=4) == 4
    # test increasing values at a gapped position
    assert get_raw_position_from_gapped(12, gapped_sequence, germline_start=2) == 6
    assert get_raw_position_from_gapped(12, gapped_sequence, germline_start=3) == 5
    assert get_raw_position_from_gapped(12, gapped_sequence, germline_start=4) == 4


# ----------------------------
#      RAW FROM ALIGNED
# ----------------------------


def test_get_raw_position_from_aligned_no_gaps(aligned_sequence_no_gaps):
    assert (
        get_raw_position_from_aligned(
            position=1,
            aligned_sequence=aligned_sequence_no_gaps,
            aligned_reference=aligned_sequence_no_gaps,
        )
        == 1
    )
    assert (
        get_raw_position_from_aligned(
            position=2,
            aligned_sequence=aligned_sequence_no_gaps,
            aligned_reference=aligned_sequence_no_gaps,
        )
        == 2
    )
    assert (
        get_raw_position_from_aligned(
            position=3,
            aligned_sequence=aligned_sequence_no_gaps,
            aligned_reference=aligned_sequence_no_gaps,
        )
        == 3
    )
    assert (
        get_raw_position_from_aligned(
            position=4,
            aligned_sequence=aligned_sequence_no_gaps,
            aligned_reference=aligned_sequence_no_gaps,
        )
        == 4
    )
    assert (
        get_raw_position_from_aligned(
            position=5,
            aligned_sequence=aligned_sequence_no_gaps,
            aligned_reference=aligned_sequence_no_gaps,
        )
        == 5
    )
    assert (
        get_raw_position_from_aligned(
            position=6,
            aligned_sequence=aligned_sequence_no_gaps,
            aligned_reference=aligned_sequence_no_gaps,
        )
        == 6
    )
    assert (
        get_raw_position_from_aligned(
            position=7,
            aligned_sequence=aligned_sequence_no_gaps,
            aligned_reference=aligned_sequence_no_gaps,
        )
        == 7
    )
    assert (
        get_raw_position_from_aligned(
            position=8,
            aligned_sequence=aligned_sequence_no_gaps,
            aligned_reference=aligned_sequence_no_gaps,
        )
        == 8
    )
    assert (
        get_raw_position_from_aligned(
            position=9,
            aligned_sequence=aligned_sequence_no_gaps,
            aligned_reference=aligned_sequence_no_gaps,
        )
        == 9
    )
    assert (
        get_raw_position_from_aligned(
            position=10,
            aligned_sequence=aligned_sequence_no_gaps,
            aligned_reference=aligned_sequence_no_gaps,
        )
        == 10
    )
    assert (
        get_raw_position_from_aligned(
            position=11,
            aligned_sequence=aligned_sequence_no_gaps,
            aligned_reference=aligned_sequence_no_gaps,
        )
        == 11
    )
    assert (
        get_raw_position_from_aligned(
            position=12,
            aligned_sequence=aligned_sequence_no_gaps,
            aligned_reference=aligned_sequence_no_gaps,
        )
        == 12
    )


def test_get_raw_position_from_aligned_with_reference_gap(
    aligned_sequence_no_gaps, aligned_sequence_with_gap
):
    assert (
        get_raw_position_from_aligned(
            position=1,
            aligned_sequence=aligned_sequence_no_gaps,
            aligned_reference=aligned_sequence_with_gap,
        )
        == 1
    )
    assert (
        get_raw_position_from_aligned(
            position=2,
            aligned_sequence=aligned_sequence_no_gaps,
            aligned_reference=aligned_sequence_with_gap,
        )
        == 2
    )
    assert (
        get_raw_position_from_aligned(
            position=3,
            aligned_sequence=aligned_sequence_no_gaps,
            aligned_reference=aligned_sequence_with_gap,
        )
        == 3
    )
    assert (
        get_raw_position_from_aligned(
            position=4,
            aligned_sequence=aligned_sequence_no_gaps,
            aligned_reference=aligned_sequence_with_gap,
        )
        == 4
    )
    assert (
        get_raw_position_from_aligned(
            position=5,
            aligned_sequence=aligned_sequence_no_gaps,
            aligned_reference=aligned_sequence_with_gap,
        )
        == 9
    )
    assert (
        get_raw_position_from_aligned(
            position=6,
            aligned_sequence=aligned_sequence_no_gaps,
            aligned_reference=aligned_sequence_with_gap,
        )
        == 10
    )
    assert (
        get_raw_position_from_aligned(
            position=7,
            aligned_sequence=aligned_sequence_no_gaps,
            aligned_reference=aligned_sequence_with_gap,
        )
        == 11
    )
    assert (
        get_raw_position_from_aligned(
            position=8,
            aligned_sequence=aligned_sequence_no_gaps,
            aligned_reference=aligned_sequence_with_gap,
        )
        == 12
    )


def test_get_raw_position_from_aligned_with_trailing_reference_gap(
    aligned_sequence_no_gaps, aligned_sequence_trailing_gap
):
    assert (
        get_raw_position_from_aligned(
            position=1,
            aligned_sequence=aligned_sequence_no_gaps,
            aligned_reference=aligned_sequence_trailing_gap,
        )
        == 1
    )
    assert (
        get_raw_position_from_aligned(
            position=2,
            aligned_sequence=aligned_sequence_no_gaps,
            aligned_reference=aligned_sequence_trailing_gap,
        )
        == 2
    )
    assert (
        get_raw_position_from_aligned(
            position=3,
            aligned_sequence=aligned_sequence_no_gaps,
            aligned_reference=aligned_sequence_trailing_gap,
        )
        == 3
    )
    assert (
        get_raw_position_from_aligned(
            position=4,
            aligned_sequence=aligned_sequence_no_gaps,
            aligned_reference=aligned_sequence_trailing_gap,
        )
        == 4
    )
    assert (
        get_raw_position_from_aligned(
            position=5,
            aligned_sequence=aligned_sequence_no_gaps,
            aligned_reference=aligned_sequence_trailing_gap,
        )
        == 5
    )
    assert (
        get_raw_position_from_aligned(
            position=6,
            aligned_sequence=aligned_sequence_no_gaps,
            aligned_reference=aligned_sequence_trailing_gap,
        )
        == 6
    )
    assert (
        get_raw_position_from_aligned(
            position=7,
            aligned_sequence=aligned_sequence_no_gaps,
            aligned_reference=aligned_sequence_trailing_gap,
        )
        == 7
    )
    assert (
        get_raw_position_from_aligned(
            position=8,
            aligned_sequence=aligned_sequence_no_gaps,
            aligned_reference=aligned_sequence_trailing_gap,
        )
        == 8
    )


def test_get_raw_position_from_aligned_with_leading_reference_gap(
    aligned_sequence_no_gaps, aligned_sequence_leading_gap
):
    assert (
        get_raw_position_from_aligned(
            position=1,
            aligned_sequence=aligned_sequence_no_gaps,
            aligned_reference=aligned_sequence_leading_gap,
        )
        == 5
    )
    assert (
        get_raw_position_from_aligned(
            position=2,
            aligned_sequence=aligned_sequence_no_gaps,
            aligned_reference=aligned_sequence_leading_gap,
        )
        == 6
    )
    assert (
        get_raw_position_from_aligned(
            position=3,
            aligned_sequence=aligned_sequence_no_gaps,
            aligned_reference=aligned_sequence_leading_gap,
        )
        == 7
    )
    assert (
        get_raw_position_from_aligned(
            position=4,
            aligned_sequence=aligned_sequence_no_gaps,
            aligned_reference=aligned_sequence_leading_gap,
        )
        == 8
    )
    assert (
        get_raw_position_from_aligned(
            position=5,
            aligned_sequence=aligned_sequence_no_gaps,
            aligned_reference=aligned_sequence_leading_gap,
        )
        == 9
    )
    assert (
        get_raw_position_from_aligned(
            position=6,
            aligned_sequence=aligned_sequence_no_gaps,
            aligned_reference=aligned_sequence_leading_gap,
        )
        == 10
    )
    assert (
        get_raw_position_from_aligned(
            position=7,
            aligned_sequence=aligned_sequence_no_gaps,
            aligned_reference=aligned_sequence_leading_gap,
        )
        == 11
    )
    assert (
        get_raw_position_from_aligned(
            position=8,
            aligned_sequence=aligned_sequence_no_gaps,
            aligned_reference=aligned_sequence_leading_gap,
        )
        == 12
    )
