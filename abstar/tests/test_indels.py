from ..utils.indels import (
    Deletion,
    Indel,
    Insertion,
    _annotate_deletion,
    _annotate_insertion,
    _fix_frameshift_deletion,
    _fix_frameshift_insertion,
    find_deletions,
    find_insertions,
)


class GermlineSegment:
    def __init__(
        self,
        realignment,
        germline_start,
        germline_end,
        query_start,
        query_end,
        query_sequence,
        germline_sequence,
    ):
        self.realignment = realignment
        self.germline_start = germline_start
        self.germline_end = germline_end
        self.query_start = query_start
        self.query_end = query_end
        self.query_sequence = query_sequence
        self.germline_sequence = germline_sequence


class Sequence:
    def __init__(self, sequence):
        self.sequence = sequence


class Antibody:
    def __init__(self, oriented_input_sequence):
        self.oriented_input = Sequence(oriented_input_sequence)


# ----------------------------
#        Indel class
# ----------------------------


def test_indel_class_init():
    # test initialization with all fields
    indel = Indel(
        {
            "len": 3,
            "pos": 10,
            "seq": "ATC",
            "fixed": True,
            "in frame": "yes",
        }
    )
    assert indel.length == 3
    assert indel.raw_position == 10
    assert indel.sequence == "ATC"
    assert indel.fixed is True
    assert indel.in_frame is True
    assert indel.imgt_position is None
    assert indel.imgt_codon is None


def test_indel_class_init_missing_fields():
    # test initialization with some fields missing
    indel = Indel(
        {
            "len": 3,
            "seq": "ATC",
        }
    )
    assert indel.length == 3
    assert indel.raw_position is None
    assert indel.sequence == "ATC"
    assert indel.fixed is None
    assert indel.in_frame is None
    assert indel.imgt_position is None
    assert indel.imgt_codon is None


def test_indel_class_contains():
    indel = Indel(
        {
            "len": 3,
            "pos": 10,
            "seq": "ATC",
            "fixed": True,
            "in frame": "yes",
        }
    )
    assert "len" in indel
    assert "pos" in indel
    assert "seq" in indel
    assert "fixed" in indel
    assert "in frame" in indel
    assert "imgt_position" not in indel
    assert "imgt_codon" not in indel


def test_indel_class_getitem():
    indel = Indel(
        {
            "len": 3,
            "pos": 10,
            "seq": "ATC",
            "fixed": True,
            "in frame": "yes",
        }
    )
    assert indel["len"] == 3
    assert indel["pos"] == 10
    assert indel["seq"] == "ATC"
    assert indel["fixed"] is True
    assert indel["in frame"] == "yes"
    assert indel["imgt_position"] is None
    assert indel["imgt_codon"] is None

    # test getting a non-existent key
    assert indel["foo"] is None


def test_indel_class_get_frame():
    # test getting in-frame status with "in frame" field
    indel = Indel(
        {
            "len": 3,
            "pos": 10,
            "seq": "ATC",
            "fixed": True,
            "in frame": "yes",
        }
    )
    assert indel._get_frame() is True


def test_indel_class_get_frame_boolean_true():
    # test getting in-frame status with boolean "in frame" field
    indel = Indel(
        {
            "len": 3,
            "pos": 10,
            "seq": "ATC",
            "fixed": True,
            "in frame": True,
        }
    )
    assert indel._get_frame() is True


def test_indel_class_get_frame_boolean_false():
    # test getting in-frame status with boolean "in frame" field
    indel = Indel(
        {
            "len": 3,
            "pos": 10,
            "seq": "ATC",
            "fixed": True,
            "in frame": False,
        }
    )
    assert indel._get_frame() is False


def test_indel_class_get_frame_missing():
    # test getting in-frame status with missing "in frame" field
    indel = Indel(
        {
            "len": 3,
            "pos": 10,
            "seq": "ATC",
            "fixed": True,
        }
    )
    assert indel._get_frame() is None


# ----------------------------
#      Insertion class
# ----------------------------


def test_insertion_class_init():
    # test initialization with all fields
    indel = Insertion(
        {
            "len": 3,
            "pos": 10,
            "seq": "ATC",
            "fixed": True,
            "in frame": "yes",
        }
    )
    assert indel.length == 3
    assert indel.raw_position == 10
    assert indel.sequence == "ATC"
    assert indel.fixed is True
    assert indel.in_frame is True
    assert indel.imgt_position is None
    assert indel.imgt_codon is None
    assert indel.type == "insertion"


def test_insertion_class_init_missing_fields():
    # test initialization with some fields missing
    indel = Insertion(
        {
            "len": 3,
            "seq": "ATC",
        }
    )
    assert indel.length == 3
    assert indel.raw_position is None
    assert indel.sequence == "ATC"
    assert indel.fixed is None
    assert indel.in_frame is None
    assert indel.imgt_position is None
    assert indel.imgt_codon is None
    assert indel.type == "insertion"


def test_insertion_class_imgt_formatting_inframe():
    # test imgt_formatted property with in-frame insertion
    indel = Insertion(
        {
            "len": 3,
            "pos": 10,
            "seq": "ATC",
            # "fixed": True,
            "in frame": True,
        }
    )
    indel.imgt_position = 11
    assert indel.imgt_formatted == "11^12>ins^atc"


def test_insertion_class_imgt_formatting_frameshift():
    # test imgt_formatted property with out-of-frame insertion
    indel = Insertion(
        {
            "len": 1,
            "pos": 10,
            "seq": "A",
            "fixed": True,
            "in frame": False,
        }
    )
    indel.imgt_position = 11
    assert indel.imgt_formatted == "11^12>ins^a#"


def test_insertion_class_abstar_formatting_inframe():
    # test abstar_formatted property with in-frame insertion
    indel = Insertion(
        {
            "len": 3,
            "pos": 10,
            "seq": "ATC",
            "fixed": True,
            "in frame": True,
        }
    )
    indel.imgt_position = 11
    assert indel.abstar_formatted == "11:3>ATC"


def test_insertion_class_abstar_formatting_frameshift():
    # test abstar_formatted property with out-of-frame insertion
    indel = Insertion(
        {
            "len": 1,
            "pos": 10,
            "seq": "A",
            "fixed": True,
            "in frame": False,
        }
    )
    indel.imgt_position = 11
    assert indel.abstar_formatted == "11:1!>A"


def test_insertion_class_json_formatting_inframe():
    # test json_formatted property with in-frame insertion
    indel = Insertion(
        {
            "len": 3,
            "pos": 10,
            "seq": "ATC",
            "fixed": True,
            "in frame": True,
        }
    )
    indel.imgt_position = 11
    assert indel.json_formatted == {
        "in_frame": "yes",
        "length": 3,
        "sequence": "ATC",
        "position": "11",
        "codon": None,
    }


def test_insertion_class_json_formatting_frameshift():
    # test json_formatted property with out-of-frame insertion
    indel = Insertion(
        {
            "len": 1,
            "pos": 10,
            "seq": "A",
            "fixed": True,
            "in frame": False,
        }
    )
    indel.imgt_position = 11
    assert indel.json_formatted == {
        "in_frame": "no",
        "length": 1,
        "sequence": "A",
        "position": "11",
        "codon": None,
    }


# ----------------------------
#      Deletion class
# ----------------------------


def test_deletion_class_init():
    # test initialization with all fields
    indel = Deletion(
        {
            "len": 3,
            "pos": 10,
            "seq": "ATC",
            "fixed": True,
            "in frame": "yes",
        }
    )
    assert indel.length == 3
    assert indel.raw_position == 10
    assert indel.sequence == "ATC"
    assert indel.fixed is True
    assert indel.in_frame == "yes"
    assert indel.imgt_position is None
    assert indel.imgt_codon is None
    assert indel.type == "deletion"

    # test initialization with some fields missing
    indel = Deletion(
        {
            "len": 3,
            "seq": "ATC",
        }
    )
    assert indel.length == 3
    assert indel.raw_position is None
    assert indel.sequence == "ATC"
    assert indel.fixed is None
    assert indel.in_frame is None
    assert indel.imgt_position is None
    assert indel.imgt_codon is None
    assert indel.type == "deletion"


def test_imgt_formatted():
    # test imgt_formatted property with single nucleotide deletion
    indel = Deletion(
        {
            "len": 1,
            "pos": 10,
            "seq": "A",
            "fixed": True,
            "in frame": True,
        }
    )
    assert indel.imgt_formatted == "a10del#"

    # test imgt_formatted property with multi-nucleotide deletion
    indel = Deletion(
        {
            "len": 3,
            "pos": 10,
            "seq": "ATC",
            "fixed": True,
            "in frame": False,
        }
    )
    assert indel.imgt_formatted == "a10-t12del(3nt)#"


def test_abstar_formatted():
    # test abstar_formatted property with single nucleotide deletion
    indel = Deletion(
        {
            "len": 1,
            "pos": 10,
            "seq": "A",
            "fixed": True,
            "in frame": False,
        }
    )
    assert indel.abstar_formatted == "10:1>!A"

    # test abstar_formatted property with multi-nucleotide deletion
    indel = Deletion(
        {
            "len": 3,
            "pos": 10,
            "seq": "ATC",
            "fixed": True,
            "in frame": True,
        }
    )
    assert indel.abstar_formatted == "10-12:3>ATC"


def test_json_formatted():
    # test json_formatted property with single nucleotide deletion
    indel = Deletion(
        {
            "len": 1,
            "pos": 10,
            "seq": "A",
            "fixed": True,
            "in frame": True,
        }
    )
    assert indel.json_formatted == {
        "in_frame": "yes",
        "length": 1,
        "sequence": "A",
        "position": "10",
        "codon": None,
    }

    # test json_formatted property with multi-nucleotide deletion
    indel = Deletion(
        {
            "len": 3,
            "pos": 10,
            "seq": "ATC",
            "fixed": True,
            "in frame": False,
        }
    )
    assert indel.json_formatted == {
        "in_frame": "no",
        "length": 3,
        "sequence": "ATC",
        "position": "10-12",
        "codon": None,
    }


# ----------------------------
#     find_insertions
# ----------------------------


def test_find_insertions():
    # test with no insertions
    segment = GermlineSegment(
        realignment=None,
        germline_start=0,
        germline_end=4,
        query_start=0,
        query_end=4,
        query_sequence="ATCG",
        germline_sequence="ATCG",
    )
    antibody = Antibody(segment.query_sequence)
    insertions = find_insertions(antibody, segment)
    assert insertions is None

    # test with a single codon-length insertion
    segment = GermlineSegment(
        realignment=None,
        germline_start=0,
        germline_end=4,
        query_start=0,
        query_end=7,
        query_sequence="ATCGTCTA",
        germline_sequence="ATCGA",
    )
    antibody = Antibody(segment.query_sequence)
    insertions = find_insertions(antibody, segment)
    assert len(insertions) == 1
    assert insertions[0].raw_position == 4
    assert insertions[0].length == 3
    assert insertions[0].sequence == "TCT"
    assert insertions[0].fixed is False
    assert insertions[0].in_frame == "yes"

    # test with a single frameshift insertion
    segment = GermlineSegment(
        realignment=None,
        germline_start=0,
        germline_end=4,
        query_start=0,
        query_end=8,
        query_sequence="ATCGTCTGA",
        germline_sequence="ATCGTTGA",
    )
    antibody = Antibody(segment.query_sequence)
    insertions = find_insertions(antibody, segment)
    assert len(insertions) == 1
    assert insertions[0].raw_position == 4
    assert insertions[0].length == 3
    assert insertions[0].sequence == "TCT"
    assert insertions[0].fixed is True
    assert insertions[0].in_frame == "no"


def test_annotate_insertion():
    # test with a codon-length insertion
    insertion = _annotate_insertion(None, 4, 3, "TCT")
    assert insertion.raw_position == 4
    assert insertion.length == 3
    assert insertion.sequence == "TCT"
    assert insertion.fixed is False
    assert insertion.in_frame == "yes"

    # test with a non-codon-length insertion
    insertion = _annotate_insertion(None, 4, 2, "TC")
    assert insertion.raw_position == 4
    assert insertion.length == 2
    assert insertion.sequence == "TC"
    assert insertion.fixed is False
    assert insertion.in_frame == "no"


def test_fix_frameshift_insertion():
    # test with a frameshift deletion
    segment = GermlineSegment(
        realignment=None,
        germline_start=0,
        germline_end=4,
        query_start=0,
        query_end=2,
        query_sequence="ATACG",
        germline_sequence="ATCG",
    )
    antibody = Antibody(segment.query_sequence)
    _fix_frameshift_insertion(antibody, segment, 2, 3)
    assert segment.query_alignment == "ATCG"
    assert segment.alignment_midline == "||||"
    assert antibody.oriented_input.sequence == "ATCG"


# ----------------------------
#    find_deletions
# ----------------------------


def test_find_deletions():
    # test with no deletions
    segment = GermlineSegment(
        realignment=None,
        germline_start=0,
        germline_end=4,
        query_start=0,
        query_end=4,
        query_sequence="ATCG",
        germline_sequence="ATCG",
    )
    antibody = Antibody(segment.query_sequence)
    deletions = find_deletions(antibody, segment)
    assert deletions is None

    # test with a single codon-length deletion
    segment = GermlineSegment(
        realignment=None,
        germline_start=0,
        germline_end=4,
        query_start=0,
        query_end=3,
        query_sequence="ATC",
        germline_sequence="ATCG",
    )
    antibody = Antibody(segment.query_sequence)
    deletions = find_deletions(antibody, segment)
    assert len(deletions) == 1
    assert deletions[0].raw_position == 3
    assert deletions[0].length == 1
    assert deletions[0].sequence == "G"
    assert deletions[0].fixed is False
    assert deletions[0].in_frame == "yes"

    # test with a single frameshift deletion
    segment = GermlineSegment(
        realignment=None,
        germline_start=0,
        germline_end=4,
        query_start=0,
        query_end=2,
        query_sequence="AT",
        germline_sequence="ATCG",
    )
    antibody = Antibody(segment.query_sequence)
    deletions = find_deletions(antibody, segment)
    assert len(deletions) == 1
    assert deletions[0].raw_position == 2
    assert deletions[0].length == 2
    assert deletions[0].sequence == "CG"
    assert deletions[0].fixed is True
    assert deletions[0].in_frame == "no"


def test_annotate_deletion():
    # test with a codon-length deletion
    deletion = _annotate_deletion(None, 3, 1, "G")
    assert deletion.raw_position == 3
    assert deletion.length == 1
    assert deletion.sequence == "G"
    assert deletion.fixed is False
    assert deletion.in_frame == "yes"

    # test with a frameshift deletion
    deletion = _annotate_deletion(None, 2, 2, "CG", fixed=True)
    assert deletion.raw_position == 2
    assert deletion.length == 2
    assert deletion.sequence == "CG"
    assert deletion.fixed is True
    assert deletion.in_frame == "no"


def test_fix_frameshift_deletion():
    # test with a frameshift deletion
    segment = GermlineSegment(
        realignment=None,
        germline_start=0,
        germline_end=4,
        query_start=0,
        query_end=2,
        query_sequence="ATG",
        germline_sequence="ATCG",
    )
    antibody = Antibody(segment.query_sequence)
    _fix_frameshift_deletion(antibody, segment, 2, 3)
    assert segment.query_alignment == "ATCG"
    assert segment.alignment_midline == "||||"
    assert antibody.oriented_input.sequence == "ATCG"
