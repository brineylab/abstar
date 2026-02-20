# Copyright (c) 2025 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

"""
Tests for the `as_dataframe` option in the `run()` function.
"""

import os

import polars as pl
import pytest
from abutils import Sequence

from ..core.abstar import run


# =============================================
#                  FIXTURES
# =============================================


@pytest.fixture
def test_data_path():
    """Path to the test FASTA file."""
    return os.path.join(
        os.path.dirname(__file__),
        "..",
        "test_data",
        "test_hiv_bnab_hcs.fasta",
    )


@pytest.fixture
def single_sequence():
    """A single antibody sequence for testing."""
    # 10E8 heavy chain sequence from test data
    return Sequence(
        "GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTGAAGCCTGGAGGATCCCTTAGACTCTCATGTTCAGCCTCTGGTTTCGACTTCGATAACGCCTGGATGACTTGGGTCCGCCAGCCTCCAGGGAAGGGCCTCGAATGGGTTGGTCGTATTACGGGTCCAGGTGAAGGTTGGTCAGTGGACTATGCTGCACCCGTGGAAGGCAGATTTACCATCTCGAGACTCAATTCAATAAATTTCTTATATTTGGAGATGAACAATTTAAGAATGGAAGACTCAGGCCTTTACTTCTGTGCCCGCACGGGAAAATATTATGATTTTTGGAGTGGCTATCCGCCGGGAGAAGAATACTTCCAAGACTGGGGCCGGGGCACCCTGGTCACCGTCTCCTCA",
        id="10E8",
    )


@pytest.fixture
def multiple_sequences():
    """A list of multiple antibody sequences for testing."""
    sequences = [
        Sequence(
            "GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTGAAGCCTGGAGGATCCCTTAGACTCTCATGTTCAGCCTCTGGTTTCGACTTCGATAACGCCTGGATGACTTGGGTCCGCCAGCCTCCAGGGAAGGGCCTCGAATGGGTTGGTCGTATTACGGGTCCAGGTGAAGGTTGGTCAGTGGACTATGCTGCACCCGTGGAAGGCAGATTTACCATCTCGAGACTCAATTCAATAAATTTCTTATATTTGGAGATGAACAATTTAAGAATGGAAGACTCAGGCCTTTACTTCTGTGCCCGCACGGGAAAATATTATGATTTTTGGAGTGGCTATCCGCCGGGAGAAGAATACTTCCAAGACTGGGGCCGGGGCACCCTGGTCACCGTCTCCTCA",
            id="10E8",
        ),
        Sequence(
            "CAGGGTCAACTAGTCCAGTCTGGAGGTGAATTGAAGAAGCCTGGGGCCTCGGTGAAGATTTCCTGTAAGACCTCGGGTTATAGATTTAGTTTCTATCATATTAATTGGATTCGACAACTAGTAGGGCGCGGACCTGAGTGGATGGGCTGGATCAGCCCTTACAACGGAGGCACAAACCTCGCACCTGAGTTGCGAGGCAGACTCGTGCTAACCACAGAGAGAGAGGTCGTGGACACCATGACCCTGTCCACGGGCACAGCCCACATGGAACTAAGGAACCTAAGATCTGACGACACGGGCATCTACTTCTGTGCAAAGGGCCTCTTGCGCGACGGTTCGTCGACGTGGCTGCCTCATTTGTGGGGCCAGGGAACCCTGCTCACCGTCTCGTCA",
            id="10J4",
        ),
        Sequence(
            "CAGGGTCAACTAGTCCAGTCTGGAGGTGAATTGAAGAAGCCTGGGGCCTCGGTGAAGATTTCCTGTAAGACCTCGGGTTATAGATTTAGTTTCTATCATATTAATTGGATTCGACAAGTAATAGGGCGCGGACCTGAGTGGATGGGCTGGATCAGCCCTTACAGCGGAGGCACAAACCTCGCACCTGAGTTCCGAGGCAGACTCGTGCTGACCACAGAGAGAGAGGTCGTGGACACCATGACCCTGTCCACGGGCACAGCCCACATGGAACTGAGGAACCTAAAATCTGACGACACGGGCATCTACTTCTGTGCAAAGGGCCTCTTGCGCGACGGTTCGTCGACGTGGCTGCCTCATTTGTGGGGCCAGGGAACCCTGCTCACCGTCTCGTCA",
            id="10M6",
        ),
    ]
    return sequences


# =============================================
#              DATAFRAME RETURN TESTS
# =============================================


def test_as_dataframe_returns_polars_dataframe(single_sequence):
    """Test that run() returns a polars DataFrame when as_dataframe=True."""
    result = run(single_sequence, as_dataframe=True)
    assert isinstance(result, pl.DataFrame), "Expected result to be a polars DataFrame"


def test_as_dataframe_false_returns_sequence(single_sequence):
    """Test that run() returns a Sequence object when as_dataframe=False (default)."""
    result = run(single_sequence, as_dataframe=False)
    assert isinstance(result, Sequence), "Expected result to be a Sequence object"


def test_as_dataframe_default_returns_sequence(single_sequence):
    """Test that the default behavior returns a Sequence object."""
    result = run(single_sequence)
    assert isinstance(result, Sequence), "Expected default result to be a Sequence object"


# =============================================
#           DATAFRAME STRUCTURE TESTS
# =============================================


def test_dataframe_has_expected_columns(single_sequence):
    """Test that the DataFrame contains expected core columns from the output schema."""
    result = run(single_sequence, as_dataframe=True)

    # Check for essential columns that should always be present
    expected_columns = [
        "sequence_id",
        "sequence",
        "v_gene",
        "j_gene",
        "productive",
        "locus",
    ]
    for col in expected_columns:
        assert col in result.columns, f"Expected column '{col}' not found in DataFrame"


def test_dataframe_has_sequence_id(single_sequence):
    """Test that sequence_id is populated in the DataFrame (not null/empty)."""
    result = run(single_sequence, as_dataframe=True)
    seq_id = result["sequence_id"][0]
    assert seq_id is not None, "sequence_id should not be None"
    assert len(seq_id) > 0, "sequence_id should not be empty"


# =============================================
#            ROW COUNT TESTS
# =============================================


def test_dataframe_row_count_single_sequence(single_sequence):
    """Test that DataFrame has correct row count for a single sequence."""
    result = run(single_sequence, as_dataframe=True)
    assert result.height == 1, "Expected DataFrame to have 1 row for single sequence input"


def test_dataframe_row_count_multiple_sequences(multiple_sequences):
    """Test that DataFrame has correct row count for multiple sequences."""
    result = run(multiple_sequences, as_dataframe=True)
    assert result.height == 3, "Expected DataFrame to have 3 rows for 3 sequence inputs"


# =============================================
#           MULTIPLE SEQUENCE TESTS
# =============================================


def test_multiple_sequences_as_dataframe(multiple_sequences):
    """Test that multiple sequences return a single DataFrame with multiple rows."""
    result = run(multiple_sequences, as_dataframe=True)

    assert isinstance(result, pl.DataFrame), "Expected result to be a polars DataFrame"
    assert result.height == 3, "Expected 3 rows in DataFrame"

    # Check that sequence_id column has values
    sequence_ids = result["sequence_id"].to_list()
    assert len(sequence_ids) == 3, "Expected 3 sequence IDs"


def test_multiple_sequences_default_returns_list(multiple_sequences):
    """Test that multiple sequences return a list of Sequence objects by default."""
    result = run(multiple_sequences, as_dataframe=False)

    assert isinstance(result, list), "Expected result to be a list"
    assert len(result) == 3, "Expected 3 Sequence objects"
    assert all(
        isinstance(s, Sequence) for s in result
    ), "All items should be Sequence objects"


# =============================================
#             DATA INTEGRITY TESTS
# =============================================


def test_dataframe_contains_annotation_data(single_sequence):
    """Test that the DataFrame contains actual annotation data, not just empty columns."""
    result = run(single_sequence, as_dataframe=True)

    # sequence should not be null
    assert result["sequence"][0] is not None, "sequence should not be null"

    # v_gene should be assigned (10E8 is a well-known antibody)
    v_gene = result["v_gene"][0]
    assert v_gene is not None, "v_gene should be assigned for a valid antibody sequence"

    # locus should be IGH for heavy chain
    locus = result["locus"][0]
    assert locus == "IGH", "locus should be IGH for heavy chain sequence"
