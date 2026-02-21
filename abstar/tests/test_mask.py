# Copyright (c) 2025 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

import pytest

from ..annotation.antibody import Antibody
from ..annotation.mask import (
    generate_cdr_mask,
    generate_gene_segment_mask,
    generate_nongermline_mask,
)


@pytest.fixture
def minimal_ab():
    # Construct a minimal Antibody with regions and masks
    ab = Antibody(sequence_id="mask_test")
    # Set nucleotide regions
    ab.fwr1 = "AAAA"
    ab.cdr1 = "BBBBB"
    ab.fwr2 = "CCCC"
    ab.cdr2 = "DDD"
    ab.fwr3 = "EEEEEE"
    ab.cdr3 = "FFFFF"
    ab.fwr4 = "GGG"
    ab.sequence = ab.fwr1 + ab.cdr1 + ab.fwr2 + ab.cdr2 + ab.fwr3 + ab.cdr3 + ab.fwr4
    # Amino-acid regions mirror for simplicity
    ab.fwr1_aa = "A" * len(ab.fwr1)
    ab.cdr1_aa = "B" * len(ab.cdr1)
    ab.fwr2_aa = "C" * len(ab.fwr2)
    ab.cdr2_aa = "D" * len(ab.cdr2)
    ab.fwr3_aa = "E" * len(ab.fwr3)
    ab.cdr3_aa = "F" * len(ab.cdr3)
    ab.fwr4_aa = "G" * len(ab.fwr4)
    # Full AA sequence
    ab.sequence_aa = (
        ab.fwr1_aa
        + ab.cdr1_aa
        + ab.fwr2_aa
        + ab.cdr2_aa
        + ab.fwr3_aa
        + ab.cdr3_aa
        + ab.fwr4_aa
    )
    # CDR3 subregions
    ab.cdr3_v = "E" * 2
    ab.cdr3_n1 = "N" * 1
    ab.cdr3_d = None
    ab.cdr3_n2 = None
    ab.cdr3_j = "E" * 2
    ab.cdr3_v_aa = ab.cdr3_v
    ab.cdr3_n1_aa = ab.cdr3_n1
    ab.cdr3_d_aa = None
    ab.cdr3_n2_aa = None
    ab.cdr3_j_aa = ab.cdr3_j
    # Gene segment mask requires these
    ab.v_sequence = ab.fwr1 + ab.cdr1 + ab.fwr2 + ab.cdr2 + ab.fwr3
    ab.v_sequence_aa = ab.fwr1_aa + ab.cdr1_aa + ab.fwr2_aa + ab.cdr2_aa + ab.fwr3_aa
    ab.j_sequence = ab.cdr3_j + ab.fwr4
    ab.j_sequence_aa = ab.cdr3_j_aa + ab.fwr4_aa
    return ab


def test_generate_cdr_mask_nt(minimal_ab):
    mask = generate_cdr_mask(minimal_ab, aa=False, as_string=True)
    # 0s for FWRs, 1/2/3 for CDRs
    expected = (
        "0" * len(minimal_ab.fwr1)
        + "1" * len(minimal_ab.cdr1)
        + "0" * len(minimal_ab.fwr2)
        + "2" * len(minimal_ab.cdr2)
        + "0" * len(minimal_ab.fwr3)
        + "3" * len(minimal_ab.cdr3)
        + "0" * len(minimal_ab.fwr4)
    )
    assert mask == expected


def test_generate_cdr_mask_aa(minimal_ab):
    mask = generate_cdr_mask(minimal_ab, aa=True, as_string=True)
    expected = (
        "0" * len(minimal_ab.fwr1_aa)
        + "1" * len(minimal_ab.cdr1_aa)
        + "0" * len(minimal_ab.fwr2_aa)
        + "2" * len(minimal_ab.cdr2_aa)
        + "0" * len(minimal_ab.fwr3_aa)
        + "3" * len(minimal_ab.cdr3_aa)
        + "0" * len(minimal_ab.fwr4_aa)
    )
    assert mask == expected


def test_generate_gene_segment_mask_nt(minimal_ab):
    # With no D-call, mask should be V, then N, then J through FWR4
    seg = generate_gene_segment_mask(minimal_ab, aa=False, as_string=True)
    # pre-CDR3 length is sequence up to start of cdr3
    cdr3_start = minimal_ab.sequence.find(minimal_ab.cdr3)
    expected = (
        "V" * cdr3_start
        + "V" * len(minimal_ab.cdr3_v)
        + "N" * len(minimal_ab.cdr3_n1)
        + ("D" * len(minimal_ab.cdr3_d) if minimal_ab.cdr3_d else "")
        + ("N" * len(minimal_ab.cdr3_n2) if minimal_ab.cdr3_n2 else "")
        + "J" * len(minimal_ab.cdr3_j)
        + "J" * len(minimal_ab.fwr4)
    )
    assert seg == expected


def test_generate_gene_segment_mask_aa(minimal_ab):
    seg = generate_gene_segment_mask(minimal_ab, aa=True, as_string=True)
    cdr3_start = (
        minimal_ab.sequence_aa.find(minimal_ab.cdr3_aa)
        if minimal_ab.sequence_aa
        else minimal_ab.sequence.find(minimal_ab.cdr3)
    )
    expected = (
        "V" * cdr3_start
        + "V" * len(minimal_ab.cdr3_v_aa)
        + "N" * len(minimal_ab.cdr3_n1_aa)
        + ("D" * len(minimal_ab.cdr3_d_aa) if minimal_ab.cdr3_d_aa else "")
        + ("N" * len(minimal_ab.cdr3_n2_aa) if minimal_ab.cdr3_n2_aa else "")
        + "J" * len(minimal_ab.cdr3_j_aa)
        + "J" * len(minimal_ab.fwr4_aa)
    )
    assert seg == expected


def test_generate_gene_segment_mask_nt_raises_when_cdr3_not_found(minimal_ab):
    minimal_ab.cdr3 = "ZZZZ"
    with pytest.raises(ValueError, match="CDR3"):
        generate_gene_segment_mask(minimal_ab, aa=False, as_string=True)


def test_generate_gene_segment_mask_aa_raises_when_cdr3_not_found(minimal_ab):
    minimal_ab.cdr3_aa = "ZZZZ"
    with pytest.raises(ValueError, match="CDR3"):
        generate_gene_segment_mask(minimal_ab, aa=True, as_string=True)


def test_generate_gene_segment_mask_nt_raises_on_length_mismatch(minimal_ab):
    minimal_ab.j_sequence = ""
    with pytest.raises(ValueError, match="length"):
        generate_gene_segment_mask(minimal_ab, aa=False, as_string=True)


def test_generate_nongermline_mask_nt(minimal_ab):
    # Provide simplistic alignments and a basic gene segment mask
    minimal_ab.gene_segment_mask = generate_gene_segment_mask(
        minimal_ab, aa=False, as_string=True
    )
    # aligned sequence = same as sequence, aligned germline = same, expect all germline (0) except N segments (1)
    minimal_ab.sequence_alignment = minimal_ab.sequence
    minimal_ab.germline_alignment = minimal_ab.sequence
    mask = generate_nongermline_mask(minimal_ab, aa=False, as_string=True)
    expected = "".join("1" if c == "N" else "0" for c in minimal_ab.gene_segment_mask)
    assert mask == expected


def test_generate_nongermline_mask_nt_raises_on_length_mismatch(minimal_ab):
    minimal_ab.gene_segment_mask = "V" * (len(minimal_ab.sequence) + 1)
    minimal_ab.sequence_alignment = minimal_ab.sequence
    minimal_ab.germline_alignment = minimal_ab.sequence
    with pytest.raises(ValueError, match="length"):
        generate_nongermline_mask(minimal_ab, aa=False, as_string=True)


def test_generate_nongermline_mask_aa(minimal_ab):
    minimal_ab.gene_segment_mask_aa = generate_gene_segment_mask(
        minimal_ab, aa=True, as_string=True
    )
    minimal_ab.sequence_alignment_aa = (
        minimal_ab.fwr1_aa
        + minimal_ab.cdr1_aa
        + minimal_ab.fwr2_aa
        + minimal_ab.cdr2_aa
        + minimal_ab.fwr3_aa
        + minimal_ab.cdr3_aa
        + minimal_ab.fwr4_aa
    )
    minimal_ab.germline_alignment_aa = minimal_ab.sequence_alignment_aa
    mask = generate_nongermline_mask(minimal_ab, aa=True, as_string=True)
    expected = "".join(
        "1" if c == "N" else "0" for c in minimal_ab.gene_segment_mask_aa
    )
    assert mask == expected
