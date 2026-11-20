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


@pytest.fixture
def assembled_ab():
    # The CDR/FWR annotations deliberately remain unset: segment masks follow
    # the complete assembled sequence, including bases outside those regions.
    return Antibody(
        sequence_id="assembled_mask", frame=1,
        v_sequence="AAAA", np1="CC", d_sequence="GGGG", np2="T",
        j_sequence="CCCC", d_call="IGHD1-14*01",
        sequence="AAAACCGGGGTCCCC", sequence_aa="KTGVP",
    )


def test_generate_gene_segment_mask_nt(assembled_ab):
    assert generate_gene_segment_mask(assembled_ab) == "VVVVNNDDDDNJJJJ"
    assert generate_gene_segment_mask(assembled_ab, as_string=False) == list("VVVVNNDDDDNJJJJ")


@pytest.mark.parametrize('frame,translation,expected', [
    (1, 'KTGVP', 'VNDNJ'),
    (2, 'KPGS', 'VNDN'),
    (3, 'NRGP', 'NNNJ'),
])
def test_generate_gene_segment_mask_aa(assembled_ab, frame, translation, expected):
    assembled_ab.frame = frame
    assembled_ab.sequence_aa = translation
    assert generate_gene_segment_mask(assembled_ab, aa=True) == expected
    assert generate_gene_segment_mask(assembled_ab, aa=True, as_string=False) == list(expected)


@pytest.mark.parametrize('v,np1,d,np2,j,sequence,translation,nt_mask,aa_mask', [
    ('AAACC', 'G', None, None, 'TTTGGA', 'AAACCGTTTGGA', 'KPFG', 'VVVVVNJJJJJJ', 'VNJJ'),
    ('AAAA', '', 'GGG', '', 'TTTTT', 'AAAAGGGTTTTT', 'KRVF', 'VVVVDDDJJJJJ', 'VNNJ'),
    ('AAA', '', None, None, 'TTT', 'AAATTT', 'KF', 'VVVJJJ', 'VJ'),
    ('TGTTGT', '', None, None, 'TTT', 'TGTTGTTTT', 'CCF', 'VVVVVVJJJ', 'VVJ'),
])
def test_gene_segment_mask_uses_complete_assembled_spans(
    v, np1, d, np2, j, sequence, translation, nt_mask, aa_mask,
):
    ab = Antibody(v_sequence=v, np1=np1, d_sequence=d, np2=np2,
                  j_sequence=j, sequence=sequence, sequence_aa=translation, frame=1)
    assert generate_gene_segment_mask(ab) == nt_mask
    assert generate_gene_segment_mask(ab, aa=True) == aa_mask


@pytest.mark.parametrize('aa,field,value', [
    (False, 'sequence', 'AAAACCGGGGTCCC'),
    (True, 'sequence_aa', 'KTGV'),
    (True, 'sequence_aa', 'KTGVPX'),
])
def test_gene_segment_mask_rejects_inconsistent_assembly(assembled_ab, aa, field, value):
    setattr(assembled_ab, field, value)
    with pytest.raises(ValueError, match="mask length"):
        generate_gene_segment_mask(assembled_ab, aa=aa)


def test_generate_nongermline_mask_nt(assembled_ab):
    assembled_ab.gene_segment_mask = generate_gene_segment_mask(assembled_ab)
    assembled_ab.sequence_alignment = assembled_ab.sequence
    assembled_ab.germline_alignment = assembled_ab.sequence
    assert generate_nongermline_mask(assembled_ab) == "000011000010000"


def test_generate_nongermline_mask_aa(assembled_ab):
    assembled_ab.gene_segment_mask_aa = generate_gene_segment_mask(assembled_ab, aa=True)
    assembled_ab.sequence_alignment_aa = assembled_ab.sequence_aa
    assembled_ab.germline_alignment_aa = assembled_ab.sequence_aa
    assert generate_nongermline_mask(assembled_ab, aa=True) == "01010"


def test_gene_segment_mask_does_not_search_for_repeated_cdr3(minimal_ab):
    minimal_ab.fwr1 = minimal_ab.cdr3 + "A"
    minimal_ab.sequence = (
        minimal_ab.fwr1
        + minimal_ab.cdr1
        + minimal_ab.fwr2
        + minimal_ab.cdr2
        + minimal_ab.fwr3
        + minimal_ab.cdr3
        + minimal_ab.fwr4
    )
    minimal_ab.v_sequence = (
        minimal_ab.fwr1
        + minimal_ab.cdr1
        + minimal_ab.fwr2
        + minimal_ab.cdr2
        + minimal_ab.fwr3
        + minimal_ab.cdr3_v
    )
    minimal_ab.np1 = "F"

    mask = generate_gene_segment_mask(minimal_ab, aa=False, as_string=True)

    expected_cdr3_start = sum(
        len(region)
        for region in (
            minimal_ab.fwr1,
            minimal_ab.cdr1,
            minimal_ab.fwr2,
            minimal_ab.cdr2,
            minimal_ab.fwr3,
        )
    )
    assert mask[:expected_cdr3_start] == "V" * expected_cdr3_start


def test_nongermline_mask_handles_terminal_deletion_after_mask_is_consumed(minimal_ab):
    minimal_ab.sequence_alignment = "A-"
    minimal_ab.germline_alignment = "AT"
    minimal_ab.gene_segment_mask = "V"

    assert generate_nongermline_mask(minimal_ab) == "0"


def test_nongermline_mask_rejects_short_segment_mask(minimal_ab):
    minimal_ab.sequence_alignment = "AA"
    minimal_ab.germline_alignment = "AA"
    minimal_ab.gene_segment_mask = "V"

    with pytest.raises(ValueError, match="mask length"):
        generate_nongermline_mask(minimal_ab)
