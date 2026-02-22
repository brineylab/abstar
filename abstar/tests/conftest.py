# Copyright (c) 2025 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

"""
Shared pytest fixtures for abstar tests.
"""

import csv
import os
import pytest
from abutils import Sequence


# Path constants
TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "test_data")

# Required columns in ground truth CSV
_GROUND_TRUTH_COLUMNS = [
    "sequence_id",
    "sequence_input",
    "fwr1",
    "cdr1",
    "fwr2",
    "cdr2",
    "fwr3",
    "cdr3",
    "fwr4",
    "fwr1_aa",
    "cdr1_aa",
    "fwr2_aa",
    "cdr2_aa",
    "fwr3_aa",
    "cdr3_aa",
    "fwr4_aa",
]


@pytest.fixture
def test_data_path():
    """Return path to test data directory."""
    return TEST_DATA_DIR


@pytest.fixture
def hiv_bnab_hc_path():
    """Path to HIV bnAb heavy chain test file."""
    return os.path.join(TEST_DATA_DIR, "test_hiv_bnab_hcs.fasta")


@pytest.fixture
def hiv_bnab_lc_path():
    """Path to HIV bnAb light chain test file."""
    return os.path.join(TEST_DATA_DIR, "test_hiv_bnab_lcs.fasta")


@pytest.fixture
def fastq_test_path():
    """Path to FASTQ test file."""
    return os.path.join(TEST_DATA_DIR, "test.fastq")


@pytest.fixture
def single_hc_sequence():
    """10E8 heavy chain - well-characterized HIV bnAb."""
    return Sequence(
        "GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTGAAGCCTGGAGGATCCCTTAGACTCTCATGTTCAGCCTCTGGTTTCGACTTCGATAACGCCTGGATGACTTGGGTCCGCCAGCCTCCAGGGAAGGGCCTCGAATGGGTTGGTCGTATTACGGGTCCAGGTGAAGGTTGGTCAGTGGACTATGCTGCACCCGTGGAAGGCAGATTTACCATCTCGAGACTCAATTCAATAAATTTCTTATATTTGGAGATGAACAATTTAAGAATGGAAGACTCAGGCCTTTACTTCTGTGCCCGCACGGGAAAATATTATGATTTTTGGAGTGGCTATCCGCCGGGAGAAGAATACTTCCAAGACTGGGGCCGGGGCACCCTGGTCACCGTCTCCTCA",
        id="10E8",
    )


@pytest.fixture
def single_lc_sequence():
    """Light chain sequence for testing."""
    return Sequence(
        "GACATCCAGATGACCCAGTCTCCATCCTCACTGTCTGCATCTGTAGGAGACAGAGTCACCATCACTTGTCGGGCGAGTCAGGGTATTAGCAGCTGGTTAGCCTGGTATCAGCAGAAACCAGGGAAAGCCCCTAAGCTCCTGATCTATGCTGCATCCAGTTTGCAAAGTGGGGTCCCATCAAGGTTCAGCGGCAGTGGATCTGGGACAGATTTCACTCTCACCATCAGCAGCCTGCAGCCTGAAGATTTTGCAACTTACTATTGTCAACAGGCTAACAGTTTCCCGCTCACTTTCGGCGGAGGGACCAAGGTGGAGATCAAACGA",
        id="test_lc",
    )


@pytest.fixture
def multiple_hc_sequences():
    """Multiple heavy chain sequences for testing."""
    return [
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


@pytest.fixture
def temp_directories(tmp_path):
    """Create temporary output and log directories."""
    output_dir = tmp_path / "output"
    log_dir = tmp_path / "logs"
    output_dir.mkdir()
    log_dir.mkdir()
    return str(output_dir), str(log_dir)


@pytest.fixture
def small_fasta_file(tmp_path, single_hc_sequence):
    """Create a small FASTA file with a single sequence."""
    fasta_path = tmp_path / "test_sequences.fasta"
    with open(fasta_path, "w") as f:
        f.write(f">{single_hc_sequence.id}\n{single_hc_sequence.sequence}\n")
    return str(fasta_path)


@pytest.fixture
def multi_sequence_fasta_file(tmp_path, multiple_hc_sequences):
    """Create a FASTA file with multiple sequences."""
    fasta_path = tmp_path / "multi_sequences.fasta"
    with open(fasta_path, "w") as f:
        for seq in multiple_hc_sequences:
            f.write(f">{seq.id}\n{seq.sequence}\n")
    return str(fasta_path)


# ---------------------------------------------------------------------------
# Ground truth fixtures for edge case tests
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def ground_truth_sequences():
    """Load all 50 ground truth sequences with their known regions.

    Returns a list of dicts parsed from ``test_50.csv``.  Each dict maps column
    names to string values (e.g. ``row["fwr1"]``, ``row["cdr3_aa"]``).

    The CSV is opened with ``encoding="utf-8-sig"`` to handle a possible
    UTF-8 BOM in the header row.
    """
    csv_path = os.path.join(TEST_DATA_DIR, "test_50.csv")
    with open(csv_path, encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh)
        rows = list(reader)
    # Validate required columns
    missing = set(_GROUND_TRUTH_COLUMNS) - set(rows[0].keys())
    assert not missing, f"Ground truth CSV is missing columns: {missing}"
    assert len(rows) == 50, f"Expected 50 ground truth rows, got {len(rows)}"
    return rows


@pytest.fixture(scope="module")
def heavy_chain_ground_truth(ground_truth_sequences):
    """Filter ground truth to heavy chain sequences only."""
    return [r for r in ground_truth_sequences if r["sequence_id"].endswith("_heavy")]


@pytest.fixture(scope="module")
def light_chain_ground_truth(ground_truth_sequences):
    """Filter ground truth to light chain sequences only."""
    return [r for r in ground_truth_sequences if r["sequence_id"].endswith("_light")]
