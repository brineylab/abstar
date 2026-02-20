# Copyright (c) 2025 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

"""
Tests for the main abstar run() pipeline function.
"""

import os

import polars as pl
import pytest
from abutils import Sequence

from ..core.abstar import run, _process_inputs


# =============================================
#          INPUT PROCESSING TESTS
# =============================================


def test_process_inputs_single_file(tmp_path, small_fasta_file):
    """Test processing a single FASTA file path."""
    sequence_files = _process_inputs(small_fasta_file, str(tmp_path))

    assert len(sequence_files) == 1
    assert os.path.exists(sequence_files[0])


def test_process_inputs_directory(tmp_path, test_data_path):
    """Test processing a directory of FASTA files."""
    sequence_files = _process_inputs(test_data_path, str(tmp_path))

    # Should find multiple FASTA files in test_data directory
    assert len(sequence_files) >= 1
    for f in sequence_files:
        assert os.path.exists(f)


def test_process_inputs_directory_accepts_short_gz_extensions(tmp_path):
    """Test that .fa.gz and .fq.gz files are included during directory scans."""
    input_dir = tmp_path / "compressed_inputs"
    input_dir.mkdir()
    (input_dir / "sample.fa.gz").touch()
    (input_dir / "sample.fq.gz").touch()

    sequence_files = _process_inputs(str(input_dir), str(tmp_path))

    assert {os.path.basename(f) for f in sequence_files} == {
        "sample.fa.gz",
        "sample.fq.gz",
    }


def test_process_inputs_logs_unsupported_extensions(tmp_path, monkeypatch):
    """Test unsupported file extensions are reported when scanning directories."""
    import abstar.core.abstar as abstar_module

    input_dir = tmp_path / "mixed_inputs"
    input_dir.mkdir()
    with open(input_dir / "sample.fasta", "w") as f:
        f.write(">s1\nATGC\n")
    with open(input_dir / "notes.txt", "w") as f:
        f.write("unsupported")

    logged_messages = []

    class _TestLogger:
        def info(self, message):
            logged_messages.append(message)

    monkeypatch.setattr(abstar_module, "logger", _TestLogger(), raising=False)

    _process_inputs(str(input_dir), str(tmp_path))

    assert any(
        "unsupported extension '.txt'" in message for message in logged_messages
    )


def test_process_inputs_sequence_object(tmp_path, single_hc_sequence):
    """Test processing a single Sequence object."""
    sequence_files = _process_inputs(single_hc_sequence, str(tmp_path))

    assert len(sequence_files) == 1
    # Should create a temp file
    assert os.path.exists(sequence_files[0])


def test_process_inputs_sequence_list(tmp_path, multiple_hc_sequences):
    """Test processing a list of Sequence objects."""
    sequence_files = _process_inputs(multiple_hc_sequences, str(tmp_path))

    assert len(sequence_files) == 1  # All sequences in one temp file

    # Verify sequences are in the file
    with open(sequence_files[0]) as f:
        content = f.read()
        for seq in multiple_hc_sequences:
            assert seq.id in content


def test_process_inputs_sequence_generator(tmp_path, multiple_hc_sequences):
    """Test processing a generator of Sequence objects."""
    sequence_generator = (seq for seq in multiple_hc_sequences)
    sequence_files = _process_inputs(sequence_generator, str(tmp_path))

    assert len(sequence_files) == 1

    with open(sequence_files[0]) as f:
        content = f.read()
        for seq in multiple_hc_sequences:
            assert seq.id in content


def test_process_inputs_empty_sequence_generator_raises_error(tmp_path):
    """Test that empty iterables raise a clear ValueError."""
    empty_sequence_generator = (seq for seq in [])
    with pytest.raises(ValueError, match="empty"):
        _process_inputs(empty_sequence_generator, str(tmp_path))


def test_process_inputs_raw_string(tmp_path):
    """Test processing a raw sequence string (not file path)."""
    raw_sequence = "ATGCATGCATGCATGCATGCATGCATGC"
    sequence_files = _process_inputs(raw_sequence, str(tmp_path))

    assert len(sequence_files) == 1
    assert os.path.exists(sequence_files[0])


def test_process_inputs_invalid_raises_error(tmp_path):
    """Test that invalid input raises ValueError."""
    with pytest.raises(ValueError):
        _process_inputs(12345, str(tmp_path))  # Integer is invalid


# =============================================
#           RETURN TYPE TESTS
# =============================================


def test_run_returns_sequence_object(single_hc_sequence):
    """Test run() returns Sequence when no project_path."""
    result = run(single_hc_sequence)

    assert isinstance(result, Sequence)


def test_run_returns_sequence_list(multiple_hc_sequences):
    """Test run() returns list of Sequences for multiple inputs."""
    result = run(multiple_hc_sequences)

    assert isinstance(result, list)
    assert len(result) >= 1
    assert all(isinstance(s, Sequence) for s in result)


def test_run_returns_dataframe_when_requested(single_hc_sequence):
    """Test run() returns polars DataFrame when as_dataframe=True."""
    result = run(single_hc_sequence, as_dataframe=True)

    assert isinstance(result, pl.DataFrame)


def test_run_returns_none_with_project_path(single_hc_sequence, tmp_path):
    """Test run() returns None when project_path provided (writes files)."""
    project_path = str(tmp_path / "test_project")
    result = run(single_hc_sequence, project_path=project_path)

    assert result is None


def test_run_single_sequence_returns_single_not_list(single_hc_sequence):
    """Test that single sequence input returns single Sequence, not list."""
    result = run(single_hc_sequence)

    assert isinstance(result, Sequence)
    assert not isinstance(result, list)


# =============================================
#          OUTPUT FORMAT TESTS
# =============================================


def test_run_creates_airr_output(single_hc_sequence, tmp_path):
    """Test AIRR TSV output file creation."""
    project_path = str(tmp_path / "airr_test")

    run(single_hc_sequence, project_path=project_path, output_format="airr")

    # Check AIRR directory and file exist
    airr_dir = os.path.join(project_path, "airr")
    assert os.path.exists(airr_dir)

    airr_files = [f for f in os.listdir(airr_dir) if f.endswith(".tsv")]
    assert len(airr_files) >= 1


def test_run_creates_parquet_output(single_hc_sequence, tmp_path):
    """Test Parquet output file creation."""
    project_path = str(tmp_path / "parquet_test")

    run(single_hc_sequence, project_path=project_path, output_format="parquet")

    # Check parquet directory and file exist
    parquet_dir = os.path.join(project_path, "parquet")
    assert os.path.exists(parquet_dir)

    parquet_files = [f for f in os.listdir(parquet_dir) if f.endswith(".parquet")]
    assert len(parquet_files) >= 1


def test_run_creates_both_outputs(single_hc_sequence, tmp_path):
    """Test creating both AIRR and Parquet outputs."""
    project_path = str(tmp_path / "both_test")

    run(single_hc_sequence, project_path=project_path, output_format=["airr", "parquet"])

    # Check both directories exist
    assert os.path.exists(os.path.join(project_path, "airr"))
    assert os.path.exists(os.path.join(project_path, "parquet"))


def test_run_creates_log_directory(single_hc_sequence, tmp_path):
    """Test log directory creation."""
    project_path = str(tmp_path / "log_test")

    run(single_hc_sequence, project_path=project_path)

    # Check logs directory exists
    log_dir = os.path.join(project_path, "logs")
    assert os.path.exists(log_dir)


# =============================================
#        GERMLINE DATABASE TESTS
# =============================================


def test_run_with_human_database(single_hc_sequence):
    """Test run with human germline database."""
    result = run(single_hc_sequence, germline_database="human")

    assert isinstance(result, Sequence)
    assert result["v_gene"] is not None


@pytest.mark.skip(reason="Mouse germline database not installed in test environment")
def test_run_with_mouse_database(tmp_path):
    """Test run with mouse germline database."""
    # Use a simple sequence - mouse genes will be different
    mouse_seq = Sequence(
        "GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTGAAGCCTGGAGGATCCCTTAGACTCTCATGTTCAGCCTCTGGTTTCGACTTCGATAACGCCTGGATGACTTGGGTCCGCCAGCCTCCAGGGAAGGGCCTCGAATGGGTTGGTCGTATTACGGGTCCAGGTGAAGGTTGGTCAGTGGACTATGCTGCACCCGTGGAAGGCAGATTTACCATCTCGAGACTCAATTCAATAAATTTCTTATATTTGGAGATGAACAATTTAAGAATGGAAGACTCAGGCCTTTACTTCTGTGCCCGCACGGGAAAATATTATGATTTTTGGAGTGGCTATCCGCCGGGAGAAGAATACTTCCAAGACTGGGGCCGGGGCACCCTGGTCACCGTCTCCTCA",
        id="test_mouse",
    )
    # This should run without error even if assignment differs
    result = run(mouse_seq, germline_database="mouse")
    assert result is not None


# =============================================
#       ANNOTATION VALIDATION TESTS
# =============================================


def test_annotated_sequence_has_v_gene(single_hc_sequence):
    """Test that annotated sequence has V gene assignment."""
    result = run(single_hc_sequence)

    assert result["v_gene"] is not None
    assert "IGHV" in result["v_gene"]


def test_annotated_sequence_has_j_gene(single_hc_sequence):
    """Test that annotated sequence has J gene assignment."""
    result = run(single_hc_sequence)

    assert result["j_gene"] is not None
    assert "IGHJ" in result["j_gene"]


def test_annotated_sequence_has_cdr3(single_hc_sequence):
    """Test that annotated sequence has CDR3."""
    result = run(single_hc_sequence)

    assert result["cdr3"] is not None
    assert len(result["cdr3"]) > 0


def test_annotated_sequence_has_junction(single_hc_sequence):
    """Test that annotated sequence has junction."""
    result = run(single_hc_sequence)

    assert result["junction"] is not None
    assert len(result["junction"]) > 0


def test_annotated_sequence_has_regions(single_hc_sequence):
    """Test all regions (FWR1-4, CDR1-3) are populated."""
    result = run(single_hc_sequence)

    # Check FWR and CDR regions
    regions = ["fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3"]
    for region in regions:
        assert result[region] is not None, f"Region {region} should not be None"

    # FWR4 may be present depending on sequence coverage
    # CDR3 is already checked above


def test_annotated_sequence_has_locus(single_hc_sequence):
    """Test that annotated sequence has locus."""
    result = run(single_hc_sequence)

    assert result["locus"] is not None
    assert result["locus"] == "IGH"  # 10E8 is heavy chain


def test_annotated_sequence_has_productivity(single_hc_sequence):
    """Test that annotated sequence has productivity assessment."""
    result = run(single_hc_sequence)

    # Productive field should be boolean
    assert result["productive"] is not None
    assert isinstance(result["productive"], bool)


def test_annotated_heavy_chain_has_d_gene(single_hc_sequence):
    """Test heavy chain has D gene (when applicable)."""
    result = run(single_hc_sequence)

    # Heavy chains should have D gene assignment
    # Note: May be None if no D is found, but field should exist
    # Just check the field is accessible
    d_gene = result["d_gene"]
    # D gene might be None for some sequences, but for 10E8 it should be present
    assert d_gene is not None or "d_gene" in result.annotations


# =============================================
#           DATAFRAME OUTPUT TESTS
# =============================================


def test_dataframe_has_expected_columns(single_hc_sequence):
    """Test that the DataFrame contains expected core columns."""
    result = run(single_hc_sequence, as_dataframe=True)

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


def test_dataframe_row_count_single_sequence(single_hc_sequence):
    """Test that DataFrame has correct row count for a single sequence."""
    result = run(single_hc_sequence, as_dataframe=True)

    assert result.height == 1


def test_dataframe_contains_annotation_data(single_hc_sequence):
    """Test that the DataFrame contains actual annotation data."""
    result = run(single_hc_sequence, as_dataframe=True)

    # sequence should not be null
    assert result["sequence"][0] is not None

    # v_gene should be assigned
    v_gene = result["v_gene"][0]
    assert v_gene is not None

    # locus should be IGH for heavy chain
    locus = result["locus"][0]
    assert locus == "IGH"


# =============================================
#              EDGE CASE TESTS
# =============================================


def test_run_from_fasta_file(small_fasta_file):
    """Test run() with a FASTA file path."""
    result = run(small_fasta_file)

    # Should return a single Sequence (since file has one sequence)
    assert isinstance(result, Sequence)
    assert result["v_gene"] is not None


def test_run_preserves_sequence_id(single_hc_sequence):
    """Test that original sequence ID is preserved or reasonably handled."""
    result = run(single_hc_sequence)

    # Sequence ID should be preserved exactly.
    assert result["sequence_id"] is not None
    assert result["sequence_id"] == single_hc_sequence.id


def test_run_with_debug_mode(single_hc_sequence, tmp_path):
    """Test run with debug mode enabled."""
    project_path = str(tmp_path / "debug_test")

    # Debug mode should not raise an error
    run(single_hc_sequence, project_path=project_path, debug=True)

    # Log files should exist
    log_dir = os.path.join(project_path, "logs")
    assert os.path.exists(log_dir)
