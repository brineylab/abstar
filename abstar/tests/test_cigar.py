#!usr/env/python
# filename: test_cigar.py

import pytest

from ..utils.cigar import make_alignment_cigar, make_cigar


@pytest.fixture
def germline_segment():
    class Realignment:
        def __init__(self, aligned_query, aligned_target):
            self.aligned_query = aligned_query
            self.aligned_target = aligned_target

    class GermlineSegment:
        def __init__(self, realignment, germline_start, query_start):
            self.realignment = realignment
            self.germline_start = germline_start
            self.query_start = query_start

    realignment = Realignment(
        aligned_query="ATCG",
        aligned_target="ATCG",
    )
    germline_segment = GermlineSegment(
        realignment=realignment,
        germline_start=0,
        query_start=0,
    )
    return germline_segment


# ----------------------------
#          cigar
# ----------------------------


def test_make_cigar(germline_segment):
    # test make_cigar function
    cigar = make_cigar(
        query_start=germline_segment.query_start,
        germline_start=germline_segment.germline_start,
        aligned_query=germline_segment.realignment.aligned_query,
        aligned_germline=germline_segment.realignment.aligned_target,
    )
    assert cigar == "4="


def test_allmatch_make_alignment_cigar():
    # test make_alignment_cigar function with identical sequences
    query = "ATCG"
    target = "ATCG"
    cigar = make_alignment_cigar(query, target)
    assert cigar == "4="


def test_mismatch_make_alignment_cigar():
    # test make_alignment_cigar function with different sequences
    query = "ATCG"
    target = "AGCG"
    cigar = make_alignment_cigar(query, target)
    assert cigar == "1=1X2="


def test_insertion_make_alignment_cigar():
    # test make_alignment_cigar function with insertions and deletions
    query = "ATCG"
    target = "A-CG"
    cigar = make_alignment_cigar(query, target)
    assert cigar == "1=1I2="


def test_insertion_mismatch_make_alignment_cigar():
    # test make_alignment_cigar function with insertions and deletions
    query = "ATGG"
    target = "A-CG"
    cigar = make_alignment_cigar(query, target)
    assert cigar == "1=1I1X1="


def test_deletion_make_alignment_cigar():
    # test make_alignment_cigar function with insertions and deletions
    query = "A-CG"
    target = "ATCG"
    cigar = make_alignment_cigar(query, target)
    assert cigar == "1=1D2="


def test_deletion_mismatch_make_alignment_cigar():
    # test make_alignment_cigar function with insertions and deletions
    query = "A-GG"
    target = "ATCG"
    cigar = make_alignment_cigar(query, target)
    assert cigar == "1=1D1X1="
