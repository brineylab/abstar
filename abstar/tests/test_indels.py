# Copyright (c) 2024 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

import pytest

from ..annotation.indels import annotate_deletions, annotate_insertions

# =============================================
#                  FIXTURES
# =============================================


@pytest.fixture
def aligned_sequence_no_indels():
    return "ATGCATGC"


@pytest.fixture
def aligned_germline_no_indels():
    return "ATGCATGC"


@pytest.fixture
def gapped_germline_no_indels():
    return "ATGCATGC"


@pytest.fixture
def aligned_sequence_single_insertion():
    """Sequence with a single insertion (GGG) at position 4"""
    return "ATGCGGGATGC"


@pytest.fixture
def aligned_germline_single_insertion():
    """Germline with a gap corresponding to the insertion in aligned_sequence_single_insertion"""
    return "ATGC---ATGC"


@pytest.fixture
def gapped_germline_single_insertion():
    """IMGT-gapped germline with periods at positions 3, 4, and 5"""
    return "AT.GCA...TGC"


@pytest.fixture
def aligned_sequence_multi_insertion():
    """Sequence with multiple insertions: (GGG) at position 4 and (AA) at position 8"""
    return "ATGCGGGATAAATGC"


@pytest.fixture
def aligned_germline_multi_insertion():
    """Germline with gaps corresponding to the insertions in aligned_sequence_multi_insertion"""
    return "ATGC---AT--ATGC"


@pytest.fixture
def gapped_germline_multi_insertion():
    """IMGT-gapped germline with appropriate gaps"""
    return "AT.GCA...T..GATGC"


@pytest.fixture
def aligned_sequence_frameshift_insertion():
    """Sequence with a frameshift insertion (A) at position 4"""
    return "ATGCAATGC"


@pytest.fixture
def aligned_germline_frameshift_insertion():
    """Germline with a gap corresponding to the frameshift insertion"""
    return "ATGC-ATGC"


@pytest.fixture
def gapped_germline_frameshift_insertion():
    """IMGT-gapped germline with a period at position 5"""
    return "AT.GCA.TGC"


@pytest.fixture
def aligned_sequence_single_deletion():
    """Sequence with a single deletion at position 5 (deletion of T)"""
    return "ATGCAGC"


@pytest.fixture
def aligned_germline_single_deletion():
    """Germline with a character at the deletion position"""
    return "ATGCATGC"


@pytest.fixture
def gapped_germline_single_deletion():
    """IMGT-gapped germline"""
    return "AT.GCATGC"


@pytest.fixture
def aligned_sequence_multi_deletion():
    """Sequence with multiple deletions: positions 5-6 (deletion of TG)"""
    return "ATGCA--C"


@pytest.fixture
def aligned_germline_multi_deletion():
    """Germline with characters at the deletion positions"""
    return "ATGCATGC"


@pytest.fixture
def gapped_germline_multi_deletion():
    """IMGT-gapped germline"""
    return "AT.GCATGC"


@pytest.fixture
def aligned_sequence_frameshift_deletion():
    """Sequence with a frameshift deletion (single nucleotide)"""
    return "ATGCA-GC"


@pytest.fixture
def aligned_germline_frameshift_deletion():
    """Germline with character at the frameshift deletion position"""
    return "ATGCATGC"


@pytest.fixture
def gapped_germline_frameshift_deletion():
    """IMGT-gapped germline"""
    return "AT.GCATGC"


# =============================================
#             INSERTION TESTS
# =============================================


def test_annotate_insertions_no_insertions(
    aligned_sequence_no_indels,
    aligned_germline_no_indels,
    gapped_germline_no_indels,
):
    result = annotate_insertions(
        aligned_sequence=aligned_sequence_no_indels,
        aligned_germline=aligned_germline_no_indels,
        gapped_germline=gapped_germline_no_indels,
        germline_start=0,
    )
    assert result == ""


def test_annotate_insertions_single_insertion(
    aligned_sequence_single_insertion,
    aligned_germline_single_insertion,
    gapped_germline_single_insertion,
):
    result = annotate_insertions(
        aligned_sequence=aligned_sequence_single_insertion,
        aligned_germline=aligned_germline_single_insertion,
        gapped_germline=gapped_germline_single_insertion,
        germline_start=0,
    )
    assert result == "5:3>GGG"


def test_annotate_insertions_multi_insertion(
    aligned_sequence_multi_insertion,
    aligned_germline_multi_insertion,
    gapped_germline_multi_insertion,
):
    result = annotate_insertions(
        aligned_sequence=aligned_sequence_multi_insertion,
        aligned_germline=aligned_germline_multi_insertion,
        gapped_germline=gapped_germline_multi_insertion,
        germline_start=0,
    )
    assert "5:3>GGG" in result
    assert "10:2>AA" in result
    assert "|" in result  # Separator for multiple insertions
    assert result.count("|") == 1  # Only one separator


def test_annotate_insertions_frameshift(
    aligned_sequence_frameshift_insertion,
    aligned_germline_frameshift_insertion,
    gapped_germline_frameshift_insertion,
):
    result = annotate_insertions(
        aligned_sequence=aligned_sequence_frameshift_insertion,
        aligned_germline=aligned_germline_frameshift_insertion,
        gapped_germline=gapped_germline_frameshift_insertion,
        germline_start=0,
    )
    assert result == "5:1>A!"  # Note the ! marking frameshift


def test_annotate_insertions_germline_offset(
    aligned_sequence_single_insertion,
    aligned_germline_single_insertion,
    gapped_germline_single_insertion,
):
    """Test with non-zero germline_start parameter"""
    result = annotate_insertions(
        aligned_sequence=aligned_sequence_single_insertion,
        aligned_germline=aligned_germline_single_insertion,
        gapped_germline="AAAAAAAAAA" + gapped_germline_single_insertion,
        germline_start=10,  # Offset the germline start
    )
    assert result == "15:3>GGG"  # Position should be offset by 10


# =============================================
#              DELETION TESTS
# =============================================


def test_annotate_deletions_no_deletions(
    aligned_sequence_no_indels,
    aligned_germline_no_indels,
    gapped_germline_no_indels,
):
    result = annotate_deletions(
        aligned_sequence=aligned_sequence_no_indels,
        aligned_germline=aligned_germline_no_indels,
        gapped_germline=gapped_germline_no_indels,
        germline_start=0,
    )
    assert result == ""


def test_annotate_deletions_single_deletion(
    aligned_sequence_single_deletion,
    aligned_germline_single_deletion,
    gapped_germline_single_deletion,
):
    # Assuming the sequence has a deletion at position 5 (T)
    result = annotate_deletions(
        aligned_sequence="ATGCA-GC",  # Explicit dash to mark deletion
        aligned_germline=aligned_germline_single_deletion,
        gapped_germline=gapped_germline_single_deletion,
        germline_start=0,
    )
    assert result == "7:1>T!"


def test_annotate_deletions_multi_deletion(
    aligned_germline_multi_deletion,
    gapped_germline_multi_deletion,
):
    # Using explicit sequence with dashes to ensure proper alignment
    result = annotate_deletions(
        aligned_sequence="ATGCA--C",  # Two consecutive deletions
        aligned_germline=aligned_germline_multi_deletion,
        gapped_germline=gapped_germline_multi_deletion,
        germline_start=0,
    )
    assert result == "7-8:2>TG!"


def test_annotate_deletions_frameshift(
    aligned_germline_frameshift_deletion,
    gapped_germline_frameshift_deletion,
):
    # Using explicit sequence with dashes to ensure proper alignment
    result = annotate_deletions(
        aligned_sequence="ATGCA-GC",  # Single deletion causing frameshift
        aligned_germline=aligned_germline_frameshift_deletion,
        gapped_germline=gapped_germline_frameshift_deletion,
        germline_start=0,
    )
    assert result == "7:1>T!"  # Note the ! marking frameshift


def test_annotate_deletions_germline_offset(
    aligned_germline_single_deletion,
    gapped_germline_single_deletion,
):
    """Test with non-zero germline_start parameter"""
    result = annotate_deletions(
        aligned_sequence="ATGCA-GC",  # Explicit dash to mark deletion
        aligned_germline=aligned_germline_single_deletion,
        gapped_germline="AAAAAAAAAA" + gapped_germline_single_deletion,
        germline_start=10,  # Offset the germline start
    )
    assert result == "17:1>T!"  # Position should be offset by 10


def test_annotate_deletions_in_frame_deletion(
    aligned_germline_multi_deletion,
    gapped_germline_multi_deletion,
):
    # Using a sequence with a 3-nucleotide deletion (in-frame)
    result = annotate_deletions(
        aligned_sequence="ATG---GC",  # 3-nucleotide deletion
        aligned_germline="ATGCATGC",  # Original germline
        gapped_germline=gapped_germline_multi_deletion,
        germline_start=0,
    )
    assert result == "5-7:3>CAT"  # Should be in-frame (no ! mark)


def test_annotate_multiple_scattered_deletions():
    """Test multiple deletions that are not adjacent"""
    result = annotate_deletions(
        aligned_sequence="A-GC-TGC",  # Deletions at positions 1 and 4
        aligned_germline="ATGCATGC",
        gapped_germline="ATGCATGC",
        germline_start=0,
    )
    assert "2:1>T!" in result  # First deletion
    assert "5:1>A!" in result  # Second deletion
    assert "|" in result  # Separator
    assert result.count("|") == 1  # One separator for two deletions


# =============================================
#           COMBINED/EDGE CASE TESTS
# =============================================


def test_annotate_insertions_with_empty_strings():
    """Test with empty strings to ensure graceful handling"""
    result = annotate_insertions(
        aligned_sequence="",
        aligned_germline="",
        gapped_germline="",
        germline_start=0,
    )
    assert result == ""


def test_annotate_deletions_with_empty_strings():
    """Test with empty strings to ensure graceful handling"""
    result = annotate_deletions(
        aligned_sequence="",
        aligned_germline="",
        gapped_germline="",
        germline_start=0,
    )
    assert result == ""


def test_annotate_insertions_complex_case():
    """Test a more complex case with multiple insertions of varying lengths"""
    result = annotate_insertions(
        aligned_sequence="ATGCAGGGCATAACC",
        aligned_germline="ATGCA---CAT--CC",
        gapped_germline="ATGCA...C..AT..CC",
        germline_start=0,
    )
    assert "5:3>GGG" in result
    assert "13:2>AA" in result
    assert result.count("|") == 1  # One separator for two insertions


def test_annotate_deletions_complex_case():
    """Test a more complex case with multiple deletions of varying lengths"""
    result = annotate_deletions(
        aligned_sequence="AT--A-G--C",
        aligned_germline="ATGCATGAGC",
        gapped_germline="ATGCATGAGC",
        germline_start=0,
    )
    assert "3-4:2>GC" in result  # First deletion (2 nucleotides)
    assert "6:1>T" in result  # Second deletion (1 nucleotide)
    assert "8-9:2>AG" in result  # Third deletion (2 nucleotides)
    assert result.count("|") == 2  # Two separators for three deletions
