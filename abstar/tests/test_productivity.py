# Copyright (c) 2024 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT


import pytest

from ..annotation.antibody import Antibody
from ..annotation.productivity import assess_productivity


@pytest.fixture
def productive_antibody():
    return Antibody(
        sequence="CAGGTTCAGCTGGTGCAGTCTGGAGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGTTACACCTTTAGCAGTTTTGGTATCAGCTGGGTGCGACAGGCCCCTGGGCAAGGGCTTGAGTGGCTGGGATGGAGCAGCACTGACAATGGTAACACAAACTATGCACAGAAGTTCCAGGGCAGAGTCACCATGACCACAGACACATCCACGACCACAGCCTACATGGAGCTGAGGAGCCTAAGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGATTTAGGGCGGTGTACCAATACCGGGTGCTATCGCAACTGGTTCGACCCCTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG",
        sequence_aa="QVQLVQSGAEVKKPGASVKVSCKASGYTFSSFGISWVRQAPGQGLEWLGWSSTDNGNTNYAQKFQGRVTMTTDTSTTTAYMELRSLRSDDTAVYYCARDLGRCTNTGCYRNWFDPWGQGTLVTVSS",
        v_call="IGHV1-18*01",
        j_call="IGHJ5*01",
        junction_aa="CARDLGRCTNTGCYRNWFDPW",
    )


@pytest.fixture
def stop_codon_antibody():
    return Antibody(
        sequence="CAGGTTCAGCTGGTGCAGTCTGGAGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGTTACACCTTTAGCAGTTTTGGTATCAGCTGGGTGCGACAGGCCCCTGGGCAAGGGCTTGAGTGGCTGGGATGGAGCAGCACTGACAATGGTAACACAAACTATGCACAGAAGTTCCAGGGCAGAGTCACCATGACCACAGACACATCCACGACCACAGCCTACATGGAGCTGAGGAGCCTAAGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGATTTAGGGCGGTGTACCAATACCGGGTGCTATCGCAACTGGTTCGACCCCTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG",
        sequence_aa="QVQLVQSGAEVKKPGASVKVSCKASGYTFSSFGISWVRQAPGQGLEWLGWSSTDNGNTNYAQKFQGRVT*TTDTSTTTAYMELRSLRSDDTAVYYCARDLGRCTNTGCYRNWFDPWGQGTLVTVSS",
        v_call="IGHV1-18*01",
        j_call="IGHJ5*01",
        junction_aa="CARDLGRCTNTGCYRNWFDPW",
    )


@pytest.fixture
def missing_conserved_cysteine_antibody():
    return Antibody(
        sequence="CAGGTTCAGCTGGTGCAGTCTGGAGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGTTACACCTTTAGCAGTTTTGGTATCAGCTGGGTGCGACAGGCCCCTGGGCAAGGGCTTGAGTGGCTGGGATGGAGCAGCACTGACAATGGTAACACAAACTATGCACAGAAGTTCCAGGGCAGAGTCACCATGACCACAGACACATCCACGACCACAGCCTACATGGAGCTGAGGAGCCTAAGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGATTTAGGGCGGTGTACCAATACCGGGTGCTATCGCAACTGGTTCGACCCCTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG",
        sequence_aa="QVQLVQSGAEVKKPGASVKVSCKASGYTFSSFGISWVRQAPGQGLEWLGWSSTDNGNTNYAQKFQGRVTMTTDTSTTTAYMELRSLRSDDTAVYYAARDLGRCTNTGCYRNWFDPWGQGTLVTVSS",
        v_call="IGHV1-18*01",
        j_call="IGHJ5*01",
        junction_aa="AARDLGRCTNTGCYRNWFDPW",
    )


@pytest.fixture
def missing_conserved_tryptophan_antibody():
    return Antibody(
        sequence="CAGGTTCAGCTGGTGCAGTCTGGAGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGTTACACCTTTAGCAGTTTTGGTATCAGCTGGGTGCGACAGGCCCCTGGGCAAGGGCTTGAGTGGCTGGGATGGAGCAGCACTGACAATGGTAACACAAACTATGCACAGAAGTTCCAGGGCAGAGTCACCATGACCACAGACACATCCACGACCACAGCCTACATGGAGCTGAGGAGCCTAAGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGATTTAGGGCGGTGTACCAATACCGGGTGCTATCGCAACTGGTTCGACCCCTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG",
        sequence_aa="QVQLVQSGAEVKKPGASVKVSCKASGYTFSSFGISWVRQAPGQGLEWLGWSSTDNGNTNYAQKFQGRVTMTTDTSTTTAYMELRSLRSDDTAVYYCARDLGRCTNTGCYRNWFDPAGQGTLVTVSS",
        v_call="IGHV1-18*01",
        j_call="IGHJ5*01",
        junction_aa="CARDLGRCTNTGCYRNWFDPA",
    )


@pytest.fixture
def locus_mismatch_antibody():
    return Antibody(
        sequence="CAGGTTCAGCTGGTGCAGTCTGGAGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGTTACACCTTTAGCAGTTTTGGTATCAGCTGGGTGCGACAGGCCCCTGGGCAAGGGCTTGAGTGGCTGGGATGGAGCAGCACTGACAATGGTAACACAAACTATGCACAGAAGTTCCAGGGCAGAGTCACCATGACCACAGACACATCCACGACCACAGCCTACATGGAGCTGAGGAGCCTAAGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGATTTAGGGCGGTGTACCAATACCGGGTGCTATCGCAACTGGTTCGACCCCTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG",
        sequence_aa="QVQLVQSGAEVKKPGASVKVSCKASGYTFSSFGISWVRQAPGQGLEWLGWSSTDNGNTNYAQKFQGRVTMTTDTSTTTAYMELRSLRSDDTAVYYCARDLGRCTNTGCYRNWFDPWGQGTLVTVSS",
        v_call="IGHV1-18*01",
        j_call="IGKJ5*01",
        junction_aa="CARDLGRCTNTGCYRNWFDPW",
    )


@pytest.fixture
def ambiguous_nucleotide_antibody():
    return Antibody(
        sequence="CAGGTTCAGCTGGTGCAGTCTGGAGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGTTACACCTTTAGCAGTTTTGGTATCAGCTGGGTGCGACAGGCCCCTGGGCAAGGGCNTGAGTGGCTGGGATGGAGCAGCACTGACAATGGTAACACAAACTATGCACAGAAGTTCCAGGGCAGAGTCACCATGACCACAGACACATCCACGACCACAGCCTACATGGAGCTGAGGAGCCTAAGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGATTTAGGGCGGTGTACCAATACCGGGTGCTATCGCAACTGGTTCGACCCCTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG",
        sequence_aa="QVQLVQSGAEVKKPGASVKVSCKASGYTFSSFGISWVRQAPGQGLEWLGWSSTDNGNTNYAQKFQGRVTMTTDTSTTTAYMELRSLRSDDTAVYYCARDLGRCTNTGCYRNWFDPWGQGTLVTVSS",
        v_call="IGHV1-18*01",
        j_call="IGHJ5*01",
        junction_aa="CARDLGRCTNTGCYRNWFDPW",
    )


def test_productive_antibody(productive_antibody):
    ab = assess_productivity(productive_antibody)
    assert ab.productive
    assert not ab.productivity_issues


def test_stop_codon_antibody(stop_codon_antibody):
    ab = assess_productivity(stop_codon_antibody)
    assert not ab.productive
    assert "stop codon" in ab.productivity_issues


def test_missing_conserved_cysteine_antibody(missing_conserved_cysteine_antibody):
    ab = assess_productivity(missing_conserved_cysteine_antibody)
    assert not ab.productive
    assert "conserved C" in ab.productivity_issues


def test_missing_conserved_tryptophan_antibody(missing_conserved_tryptophan_antibody):
    ab = assess_productivity(missing_conserved_tryptophan_antibody)
    assert not ab.productive
    assert "conserved W" in ab.productivity_issues


def test_locus_mismatch_antibody(locus_mismatch_antibody):
    ab = assess_productivity(locus_mismatch_antibody)
    assert not ab.productive
    assert "locus mismatch" in ab.productivity_issues


def test_ambiguous_nucleotide_antibody(ambiguous_nucleotide_antibody):
    ab = assess_productivity(ambiguous_nucleotide_antibody)
    assert not ab.productive
    assert "ambiguous nucleotide" in ab.productivity_issues


@pytest.mark.parametrize("ambiguous_base", ["r", "Y", "u", "-"])
def test_all_non_acgt_bases_are_ambiguous(productive_antibody, ambiguous_base):
    productive_antibody.sequence = productive_antibody.sequence[:10] + ambiguous_base

    ab = assess_productivity(productive_antibody)

    assert not ab.productive
    assert "ambiguous nucleotide(s)" in ab.productivity_issues


@pytest.mark.parametrize("junction,junction_aa", [("", ""), ("TGT", "C")])
def test_empty_or_truncated_junction_is_nonproductive(
    productive_antibody, junction, junction_aa
):
    productive_antibody.junction = junction
    productive_antibody.junction_aa = junction_aa

    ab = assess_productivity(productive_antibody)

    assert not ab.productive
    assert "missing or truncated junction" in ab.productivity_issues
    assert ab.vj_in_frame is False


def test_junction_length_must_be_a_multiple_of_three(productive_antibody):
    productive_antibody.junction = "TGTAAAT"

    ab = assess_productivity(productive_antibody)

    assert not ab.productive
    assert "junction length is not a multiple of 3" in ab.productivity_issues
    assert ab.vj_in_frame is False


def test_junction_must_start_in_v_reading_frame(productive_antibody):
    productive_antibody.junction = "TGTAAATGG"
    productive_antibody.frame = 1
    productive_antibody.junction_start = 1

    ab = assess_productivity(productive_antibody)

    assert not ab.productive
    assert "V/J junction is out of frame" in ab.productivity_issues
    assert ab.vj_in_frame is False


@pytest.mark.parametrize("locus", ["TRA", "TRB", "TRD", "TRG"])
def test_tcr_junction_requires_terminal_phenylalanine(locus):
    ab = Antibody(
        sequence="TGTGCTTTT",
        sequence_aa="CAF",
        v_call=f"{locus}V1*01",
        j_call=f"{locus}J1*01",
        locus=locus,
        junction="TGTGCTTTT",
        junction_aa="CAF",
        frame=1,
    )
    ab.junction_start = 0

    assessed = assess_productivity(ab)

    assert assessed.productive
    assert assessed.vj_in_frame is True


def test_tcr_junction_rejects_terminal_tryptophan():
    ab = Antibody(
        sequence="TGTGCTTGG",
        sequence_aa="CAW",
        v_call="TRAV1*01",
        j_call="TRAJ1*01",
        locus="TRA",
        junction="TGTGCTTGG",
        junction_aa="CAW",
        frame=1,
    )
    ab.junction_start = 0

    assessed = assess_productivity(ab)

    assert not assessed.productive
    assert "junction does not end with conserved F" in assessed.productivity_issues
