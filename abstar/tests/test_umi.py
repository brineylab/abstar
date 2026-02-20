# Copyright (c) 2025 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

import os
import tempfile

import pytest
from abutils import Sequence

from ..annotation.umi import UMI, parse_umis


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
    # output should overwrite input when output_file is None
    assert os.path.exists(out)
    with open(out) as f:
        data = f.read()
    assert "s1_AAAA" in data
    assert "s2_GGGG" in data


def test_parse_umis_builtin_pattern_defaults():
    # builtin uses multiple trailing patterns and fixed length; just sanity-check it runs
    seq = Sequence("N" * 20 + "TCAGCGGGAAGACATT")
    umi = parse_umis(sequences=seq, pattern="smartseq-human-bcr", length=None)
    # UMI may be empty if align fails with default mismatch, but should be str or None
    assert umi is None or isinstance(umi, str)


def test_parse_umis_builtin_uses_builtin_mismatch_default():
    seq = Sequence("ACGTACGTACGT" + "TCAGCGGGAAGTCGTT")
    umi_default = parse_umis(sequences=seq, pattern="smartseq-human-bcr", length=None)
    umi_strict = parse_umis(
        sequences=seq,
        pattern="smartseq-human-bcr",
        length=None,
        allowed_mismatches=1,
    )

    assert umi_default == "ACGTACGTACGT"
    assert umi_strict is None


def test_umi_flanked_pattern_without_length_does_not_error():
    seq = Sequence("ATGC" + "NNNN" + "TTTT")
    u = UMI(sequence=seq, pattern="ATGC[UMI]TTTT", length=None)
    assert u.umi == "NNNN"
