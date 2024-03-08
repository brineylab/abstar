import pytest

from ..utils.cigar import make_cigar, make_alignment_cigar


@pytest.fixture
def germline_segment():
    class Realignment:
        def __init__(self, aligned_query, aligned_target):
            self.aligned_query = aligned_query
            self.aligned_target = aligned_target

    class GermlineSegment:
        def __init__(
            self, realignment, germline_start, germline_end, query_start, query_end
        ):
            self.realignment = realignment
            self.germline_start = germline_start
            self.germline_end = germline_end
            self.query_start = query_start
            self.query_end = query_end

    realignment = Realignment(
        aligned_query="ATCG",
        aligned_target="ATCG",
    )
    germline_segment = GermlineSegment(
        realignment=realignment,
        germline_start=0,
        germline_end=4,
        query_start=0,
        query_end=4,
    )
    return germline_segment


# ----------------------------
#          cigar
# ----------------------------


def test_make_cigar(germline_segment):
    # test make_cigar function
    cigar = make_cigar(germline_segment)
    assert cigar == "4M"


def test_make_alignment_cigar():
    # test make_alignment_cigar function with identical sequences
    query = "ATCG"
    target = "ATCG"
    cigar = make_alignment_cigar(query, target)
    assert cigar == "4M"

    # test make_alignment_cigar function with different sequences
    query = "ATCG"
    target = "AGCG"
    cigar = make_alignment_cigar(query, target)
    assert cigar == "1M1D2M"

    # test make_alignment_cigar function with insertions and deletions
    query = "ATCG"
    target = "A-TG"
    cigar = make_alignment_cigar(query, target)
    assert cigar == "1M1D1M1I1M"
