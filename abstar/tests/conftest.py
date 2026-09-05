# Copyright (c) 2025 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

"""
Shared pytest fixtures for abstar tests.
"""

import os
import pytest
from abutils import Sequence


# Path constants
TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "test_data")


@pytest.fixture
def test_data_path():
    """Return path to test data directory."""
    return TEST_DATA_DIR


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


@pytest.fixture
def real_bcr_cases():
    """Fresh immutable published BCR cases; loading never invokes annotation."""
    from abstar.tests.corpus import load_real_bcr_cases
    return load_real_bcr_cases()


@pytest.fixture
def pilot_loss_cases():
    """Fresh cases for the eight original dataset 1279068 annotation losses."""
    from abstar.tests.corpus import PILOT_LOSS_IDS, load_real_bcr_cases
    return tuple(case for case in load_real_bcr_cases()
                 if case.dataset == '1279068' and case.sequence_id in PILOT_LOSS_IDS)


@pytest.fixture(scope="module")
def public_bcr_cases():
    """Authenticated concordant IGH, IGK and IGL controls, in input order."""
    from abstar.tests.corpus import load_real_bcr_cases

    cases = {(case.dataset, case.sequence_id): case for case in load_real_bcr_cases()}
    return tuple(cases[key] for key in (
        ("1287199", "GCTGCGAGTCCTGCTT-1_contig_1"),
        ("1287191", "AAAGCAACAATCGAAA-1_contig_3"),
        ("1287157", "CTCATTAAGGATCGCA-1_contig_1"),
    ))


@pytest.fixture
def public_bcr_inputs(public_bcr_cases, tmp_path):
    """Fresh iterables and equivalent files; filenames define the same order."""
    def sequence_list():
        return [case.as_sequence() for case in public_bcr_cases]

    def sequence_iterator():
        return iter(sequence_list())

    def sequence_generator():
        return (case.as_sequence() for case in public_bcr_cases)

    fasta = tmp_path / "controls.fasta"
    fastq = tmp_path / "controls.fastq"
    fasta.write_text("".join(f">{case.sequence_id}\n{case.sequence}\n"
                             for case in public_bcr_cases))
    fastq.write_text("".join(f"@{case.sequence_id}\n{case.sequence}\n+\n"
                             f"{'I' * len(case.sequence)}\n"
                             for case in public_bcr_cases))
    flat = tmp_path / "flat"
    nested = tmp_path / "nested"
    # Create files out of input order and use 2/10 to distinguish natural
    # from lexical sorting. The shared parent must survive copying too.
    for index in (2, 0, 1):
        ordinal = (1, 2, 10)[index]
        case = public_bcr_cases[index]
        for path in (flat / f"sample{ordinal}.fasta",
                     nested / "batch" / f"sample{ordinal}" / "reads.fasta"):
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(f">{case.sequence_id}\n{case.sequence}\n")
    return {
        "list": sequence_list,
        "iterator": sequence_iterator,
        "generator": sequence_generator,
        "fasta": lambda: str(fasta),
        "fastq": lambda: str(fastq),
        "flat_directory": lambda: str(flat),
        "nested_directory": lambda: str(nested),
    }
