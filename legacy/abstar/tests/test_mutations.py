# Copyright (c) 2024 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

import pytest

from ..annotation.antibody import Antibody
from ..annotation.mutations import annotate_mutations


@pytest.fixture
def germline():
    return "ATGCATGC"


@pytest.fixture
def aligned_sequence_with_insertion():
    return "ATGCAAAATGC"


@pytest.fixture
def aligned_germline_with_insertion():
    return "ATGC---ATGC"


@pytest.fixture
def gapped_germline():
    return "ATG.CAT...GC"


@pytest.fixture
def sequence_single_mutation_5at():
    # mutation is A>T at raw position 5
    return "ATGCTTGC"


@pytest.fixture
def sequence_single_mutation_7gc():
    # mutation is G>C at raw position 7
    return "ATGCATCC"


@pytest.fixture
def sequence_multi_mutation_5at_7gc():
    # mutations (raw positions) are 5:A>T and 7:G>C
    return "ATGCTTCC"


def test_annotate_mutations_no_mutations(germline):
    mutations = annotate_mutations(
        aligned_sequence=germline,
        aligned_germline=germline,
        gapped_germline=germline,
        germline_start=0,
    )
    assert len(mutations) == 0


def test_annotate_mutations_single_mutation_5at(germline, sequence_single_mutation_5at):
    mutations = annotate_mutations(
        aligned_sequence=sequence_single_mutation_5at,
        aligned_germline=germline,
        gapped_germline=germline,
        germline_start=0,
    )
    assert "5:A>T" in mutations
    assert len(mutations) == 1


def test_annotate_mutations_single_mutation_7gc(germline, sequence_single_mutation_7gc):
    mutations = annotate_mutations(
        aligned_sequence=sequence_single_mutation_7gc,
        aligned_germline=germline,
        gapped_germline=germline,
        germline_start=0,
    )
    assert "7:G>C" in mutations
    assert len(mutations) == 1


def test_annotate_mutations_single_mutation_single_germline_gap(
    germline, sequence_single_mutation_5at, gapped_germline
):
    mutations = annotate_mutations(
        aligned_sequence=sequence_single_mutation_5at,
        aligned_germline=germline,
        gapped_germline=gapped_germline,
        germline_start=0,
    )
    assert "6:A>T" in mutations
    assert len(mutations) == 1


def test_annotate_mutations_single_mutation_multi_germline_gap(
    germline, sequence_single_mutation_7gc, gapped_germline
):
    mutations = annotate_mutations(
        aligned_sequence=sequence_single_mutation_7gc,
        aligned_germline=germline,
        gapped_germline=gapped_germline,
        germline_start=0,
    )
    assert "11:G>C" in mutations
    assert len(mutations) == 1


def test_annotate_mutations_multi_mutation(germline, sequence_multi_mutation_5at_7gc):
    mutations = annotate_mutations(
        aligned_sequence=sequence_multi_mutation_5at_7gc,
        aligned_germline=germline,
        gapped_germline=germline,
        germline_start=0,
    )
    assert "5:A>T" in mutations
    assert "7:G>C" in mutations
    assert len(mutations) == 2


def test_annotate_mutations_multi_mutation_germline_gaps(
    germline, sequence_multi_mutation_5at_7gc, gapped_germline
):
    mutations = annotate_mutations(
        aligned_sequence=sequence_multi_mutation_5at_7gc,
        aligned_germline=germline,
        gapped_germline=gapped_germline,
        germline_start=0,
    )
    assert "6:A>T" in mutations
    assert "11:G>C" in mutations
    assert len(mutations) == 2


def test_annotate_mutations_single_mutation_5at_aa(
    germline, sequence_single_mutation_5at
):
    mutations = annotate_mutations(
        aligned_sequence=sequence_single_mutation_5at,
        aligned_germline=germline,
        gapped_germline=germline,
        germline_start=0,
    )
    assert "5:A>T" in mutations
    assert len(mutations) == 1
