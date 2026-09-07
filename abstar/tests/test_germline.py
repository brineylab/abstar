# Copyright (c) 2024 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

from types import SimpleNamespace

import pytest
from abutils import Sequence

from ..annotation.antibody import Antibody
from ..annotation.germline import (
    get_germline,
    get_germline_database_path,
    process_cgene_alignment,
    process_dgene_alignment,
    process_jgene_alignment,
    reassign_dgene,
)


@pytest.mark.parametrize('index', range(36))
@pytest.mark.parametrize('segment', ['v', 'j'])
def test_real_bcr_full_query_boundaries_match_authenticated_traces(index, segment):
    from abstar.tests.corpus import load_real_bcr_cases
    from abstar.annotation.germline import realign_germline, VJ_BOUNDARY_PARAMS

    case = load_real_bcr_cases()[index]
    trace = case.source['alignment'][segment]
    origin = 0 if segment == 'v' else case.expected['junction_start'] + 3
    _, local = realign_germline(
        case.sequence[origin:], trace['reference'] + '__homo_sapiens', 'human',
        receptor='bcr', local_full_query=True, local_aln_params=VJ_BOUNDARY_PARAMS,
    )
    assert (origin + local.query_begin, origin + local.query_end + 1,
            local.target_begin, local.target_end + 1) == (
        trace['query_start'], trace['query_end'],
        trace['germline_start'], trace['germline_end'],
    )
    assert local.score == trace['score']
    assert local.aligned_query == trace['query_aligned']
    assert local.aligned_target == trace['germline_aligned']


def test_real_bcr_equal_score_j_repeat_selects_earlier_query_endpoint():
    import abutils
    from abstar.tests.corpus import load_real_bcr_cases
    from abstar.annotation.germline import VJ_BOUNDARY_PARAMS

    case = next(c for c in load_real_bcr_cases()
                if c.sequence_id == 'CCATGTCCAGTCTTCC-1_contig_1')
    germline = get_germline('IGHJ4*02', 'human', receptor='bcr', exact_match=True)
    primary = abutils.tl.local_alignment(case.sequence[412:487], germline,
                                         **VJ_BOUNDARY_PARAMS)
    secondary = abutils.tl.local_alignment(case.sequence[487:], germline,
                                           **VJ_BOUNDARY_PARAMS)
    combined = abutils.tl.local_alignment(case.sequence[412:], germline,
                                          **VJ_BOUNDARY_PARAMS)
    assert primary.score == secondary.score == combined.score == 65
    assert (412 + primary.query_begin, 413 + primary.query_end) == (452, 487)
    assert (487 + secondary.query_begin, 488 + secondary.query_end) == (493, 528)
    assert (412 + combined.query_begin, 413 + combined.query_end) == (452, 487)


def _alignment(**kwargs):
    defaults = {
        "query": "ABCDEFGHIJKL",
        "target": "abcdefghijkl",
        "query_begin": 0,
        "query_end": 0,
        "target_begin": 0,
        "target_end": 0,
        "score": 10,
    }
    return SimpleNamespace(**(defaults | kwargs))


@pytest.mark.parametrize('nt_start,frame,aa_start', [
    (0, 1, 0), (1, 3, 1), (2, 2, 1),
    (3, 1, 1), (4, 3, 2), (137, 2, 46),
])
def test_translated_reference_origin_accounts_for_partial_first_codon(nt_start, frame, aa_start):
    from abstar.annotation.germline import translated_reference_start

    assert translated_reference_start(nt_start, frame) == aa_start


def test_j_germline_end_uses_semiglobal_target_coordinate_once():
    ab = Antibody(sequence_id="j-coordinate")
    semiglobal = _alignment(query_begin=2, query_end=8, target_begin=3, target_end=9)
    local = _alignment(query_begin=1, query_end=4, target_begin=2, target_end=5)

    process_jgene_alignment("X" * 30, 10, semiglobal, local, ab)

    assert ab.j_germline_start == 5
    assert ab.j_germline_end == 10
    assert ab.j_germline == "fghij"


def test_d_germline_end_uses_target_span_when_alignment_has_indel():
    ab = Antibody(sequence_id="d-coordinate", d_call="IGHD1*01")
    local = _alignment(query_begin=1, query_end=3, target_begin=2, target_end=6)

    process_dgene_alignment("X" * 30, 10, local, ab)

    assert (ab.d_sequence_start, ab.d_sequence_end) == (11, 14)
    assert (ab.d_germline_start, ab.d_germline_end) == (2, 7)
    assert ab.d_germline == "cdefg"


def test_c_sequence_end_does_not_double_count_query_offset():
    ab = Antibody(sequence_id="c-coordinate")
    semiglobal = _alignment(query_begin=2, query_end=8, target_begin=3, target_end=9)
    local = _alignment(query_begin=1, query_end=4, target_begin=2, target_end=5)

    process_cgene_alignment("X" * 30, 10, semiglobal, local, ab)

    assert (ab.c_sequence_start, ab.c_sequence_end) == (13, 17)
    assert len(ab.c_sequence) == 4

# ----------------------------
#      DATABASE PATHS
# ----------------------------


def test_get_germline_database_path():
    path = get_germline_database_path(germdb_name="human", receptor="bcr")
    assert path


def test_get_germline_database_name_invalid():
    with pytest.raises(FileNotFoundError, match="NotAGermlineDatabase"):
        get_germline_database_path("NotAGermlineDatabase", receptor="bcr")


def test_get_germline_database_receptor_invalid():
    with pytest.raises(ValueError, match="receptor"):
        get_germline_database_path("human", receptor="abc")


# ----------------------------
#        GET GERMLINE
# ----------------------------


def test_get_single_germline():
    germ = get_germline(
        germline_gene="IGHV1-2*02",
        germdb_name="human",
        receptor="bcr",
        exact_match=True,
    )
    assert isinstance(germ, Sequence)
    assert germ.id == "IGHV1-2*02"


def test_get_multiple_germlines():
    germs = get_germline(
        germline_gene="IGHV1-2",
        germdb_name="human",
        receptor="bcr",
        exact_match=False,
    )
    assert len(germs) >= 2
    assert all([isinstance(germ, Sequence) for germ in germs])
    assert all([germ.id.startswith("IGHV1-2") for germ in germs])


def test_get_single_germline_tcr():
    germ = get_germline(
        germline_gene="TRAV1-1*01",
        germdb_name="human",
        receptor="tcr",
        exact_match=True,
    )
    assert isinstance(germ, Sequence)
    assert germ.id == "TRAV1-1*01"


def test_get_multiple_germlines_tcr():
    germs = get_germline(
        germline_gene="TRAV1-1",
        germdb_name="human",
        receptor="tcr",
        exact_match=False,
    )
    assert len(germs) >= 2
    assert all([isinstance(germ, Sequence) for germ in germs])
    assert all([germ.id.startswith("TRAV1-1") for germ in germs])


@pytest.mark.parametrize(
    "locus,sequence,expected_prefix",
    [
        ("TRB", "GGGACAGGGGGC", "TRBD"),
        ("TRD", "ACTGGGGGATACG", "TRDD"),
    ],
)
def test_reassign_dgene_uses_tcr_locus_database(locus, sequence, expected_prefix):
    alignment = reassign_dgene(
        sequence=sequence,
        germdb_name="human",
        locus=locus,
        receptor="tcr",
    )

    assert alignment is not None
    assert alignment.target.id.startswith(expected_prefix)


@pytest.mark.parametrize("locus", ["TRA", "TRG", "IGK", "IGL"])
def test_reassign_dgene_skips_loci_without_d_genes(locus):
    assert (
        reassign_dgene(
            sequence="GGGACAGGGGGC",
            germdb_name="human",
            locus=locus,
            receptor="tcr" if locus.startswith("TR") else "bcr",
        )
        is None
    )


def test_get_single_germline_nonexistent():
    with pytest.raises(ValueError, match=r"IGHV1-1\*01"):
        get_germline(
            germline_gene="IGHV1-1*01",
            germdb_name="human",
            receptor="bcr",
            exact_match=True,
        )


def test_get_single_germline_nonunique():
    with pytest.raises(ValueError, match="IGHV1-2"):
        get_germline(
            germline_gene="IGHV1-2",
            germdb_name="human",
            receptor="bcr",
            exact_match=True,
        )
