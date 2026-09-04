# Copyright (c) 2025 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

import os
import pytest
from abutils import Sequence

from ..annotation.umi import UMI, parse_umis
from ..core.abstar import run


def test_umi_slice_positive_length():
    u = UMI(sequence=Sequence("ACGTACGTACGT"), pattern=None, length=4)
    assert u.umi == "ACGT"
    assert u.num_mismatches == 0


def test_umi_slice_negative_length():
    u = UMI(sequence=Sequence("ACGTACGTACGT"), pattern=None, length=-4)
    assert u.umi == "ACGT"


def test_umi_leading_pattern_extract():
    # pattern: ATGC[UMI], expect 6 bases after leading pattern
    seq = Sequence("ATGC" + "NNNNNN" + "TTTT")
    u = UMI(sequence=seq, pattern="ATGC[UMI]", length=6)
    assert u.leading == "ATGC"
    assert u.trailing is None
    assert u.umi == "NNNNNN"


def test_umi_trailing_pattern_extract():
    # pattern: [UMI]ATGC, extract 6 bases before trailing pattern
    seq = Sequence("NNNNNN" + "ATGC")
    u = UMI(sequence=seq, pattern="[UMI]ATGC", length=6)
    assert u.trailing == "ATGC"
    assert u.leading is None
    assert u.umi == "NNNNNN"


def test_parse_umis_single_sequence_no_pattern():
    umi = parse_umis(sequences=Sequence("ACGTACGT"), pattern=None, length=3)
    assert umi == "ACG"


def test_parse_umis_list_sequences_with_pattern_and_mismatch():
    seqs = [Sequence("ATGC" + "AAAAAA" + "TTTT"), Sequence("ATGC" + "CCCCCC" + "TTTT")]
    out = parse_umis(
        sequences=seqs,
        pattern="ATGC[UMI]",
        length=6,
        allowed_mismatches=0,
        output_file=None,
    )
    assert isinstance(out, list)
    assert all(isinstance(s, Sequence) for s in out)
    assert all("umi" in s for s in out)
    assert out[0]["umi"] == "AAAAAA"
    assert out[1]["umi"] == "CCCCCC"


def test_parse_umis_file_roundtrip(tmp_path):
    # create temp fasta file
    input_path = tmp_path / "umis.fasta"
    with open(input_path, "w") as f:
        f.write(">s1\nATGC" + "AAAA" + "TT\n")
        f.write(">s2\nATGC" + "GGGG" + "TT\n")
    out = parse_umis(
        sequences=str(input_path),
        pattern="ATGC[UMI]",
        length=4,
        allowed_mismatches=0,
        fmt="fasta",
    )
    original = input_path.read_text()
    # A sibling output is created; the source file is preserved.
    assert os.path.exists(out)
    assert out != str(input_path)
    assert input_path.read_text() == original
    with open(out) as f:
        data = f.read()
    assert "s1_AAAA" in data
    assert "s2_GGGG" in data


def test_parse_umis_builtin_pattern_defaults():
    # The built-in permits two mismatches when no explicit override is supplied.
    seq = Sequence("ACGTACGTACGT" + "TCAGCGGGAAGACACC")
    umi = parse_umis(sequences=seq, pattern="smartseq-human-bcr", length=None)
    assert umi == "ACGTACGTACGT"


def test_pattern_without_length_infers_umi_before_trailing_anchor():
    umi = parse_umis(
        Sequence("ACGTAC" + "TTGGCC"),
        pattern="[UMI]TTGGCC",
        length=None,
        allowed_mismatches=0,
    )

    assert umi == "ACGTAC"


def test_pattern_ending_in_umi_requires_length():
    with pytest.raises(ValueError, match="length is required"):
        parse_umis(Sequence("ATGCACGT"), pattern="ATGC[UMI]", length=None)


def test_pattern_and_length_cannot_both_be_missing():
    with pytest.raises(ValueError, match="Either pattern or length"):
        parse_umis(Sequence("ATGC"), pattern=None, length=None)


def test_negative_pattern_length_uses_absolute_slice_after_reverse_complement():
    # Reverse-complementing the input produces ATGC + AAAA at its 5' end.
    umi = parse_umis(
        Sequence("GGGGTTTTGCAT"),
        pattern="ATGC[UMI]",
        length=-4,
        allowed_mismatches=0,
    )

    assert umi == "AAAA"


def test_pattern_search_is_limited_to_sequence_end():
    sequence = Sequence("A" * 60 + "ATGC" + "CCCC")

    umi = parse_umis(
        sequence,
        pattern="ATGC[UMI]",
        length=4,
        allowed_mismatches=0,
    )

    assert umi is None


def test_end_search_window_can_be_extended():
    sequence = Sequence("A" * 60 + "ATGC" + "CCCC")

    umi = parse_umis(
        sequence,
        pattern="ATGC[UMI]",
        length=4,
        allowed_mismatches=0,
        extra_length_for_alignment=60,
    )

    assert umi == "CCCC"


def test_iterable_retains_record_without_detected_umi():
    matched = Sequence("ATGCAAAA", id="matched")
    unmatched = Sequence("TTTTCCCC", id="unmatched")

    output = parse_umis(
        [matched, unmatched],
        pattern="ATGC[UMI]",
        length=4,
        allowed_mismatches=0,
    )

    assert [sequence.id for sequence in output] == ["matched", "unmatched"]
    assert output[0]["umi"] == "AAAA"
    assert output[1]["umi"] is None


def test_file_retains_record_without_detected_umi(tmp_path):
    input_path = tmp_path / "input.fasta"
    input_path.write_text(">matched\nATGCAAAA\n>unmatched\nTTTTCCCC\n")

    output_path = parse_umis(
        str(input_path),
        pattern="ATGC[UMI]",
        length=4,
        allowed_mismatches=0,
    )

    output = open(output_path).read()
    assert ">matched_AAAA" in output
    assert ">unmatched\n" in output


def test_file_rejects_explicit_in_place_output(tmp_path):
    input_path = tmp_path / "input.fasta"
    input_path.write_text(">read\nATGCAAAA\n")

    with pytest.raises(ValueError, match="must not replace"):
        parse_umis(
            str(input_path),
            output_file=str(input_path),
            pattern="ATGC[UMI]",
            length=4,
        )


@pytest.mark.e2e
def test_annotation_pipeline_retains_sequence_without_umi(single_hc_sequence):
    result = run(
        single_hc_sequence,
        umi_pattern="TTTTTTTT[UMI]",
        umi_length=4,
        n_processes=1,
    )

    assert isinstance(result, Sequence)
    assert result.id == single_hc_sequence.id
    assert result["umi"] is None
