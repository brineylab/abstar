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

from ..core.abstar import (
    _copy_inputs_to_project,
    _get_sample_names,
    _process_inputs,
    run,
)


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


def test_process_inputs_generator_is_consumed_once(tmp_path, multiple_hc_sequences):
    """A one-shot iterable must retain every sequence when written to FASTA."""
    sequence_generator = (sequence for sequence in multiple_hc_sequences)

    sequence_files = _process_inputs(sequence_generator, str(tmp_path))

    with open(sequence_files[0]) as fasta_file:
        content = fasta_file.read()
    assert [sequence.id for sequence in multiple_hc_sequences] == [
        line[1:] for line in content.splitlines() if line.startswith(">")
    ]


def test_process_inputs_empty_directory_raises_error(tmp_path):
    with pytest.raises(ValueError, match="No supported FASTA or FASTQ"):
        _process_inputs(str(tmp_path), str(tmp_path))


def test_sample_names_are_unique_for_nested_and_extension_collisions(tmp_path):
    paths = [
        tmp_path / "a" / "sample.fasta",
        tmp_path / "b" / "sample.fasta",
        tmp_path / "b" / "sample.fa",
    ]
    for path in paths:
        path.parent.mkdir(exist_ok=True)
        path.touch()

    names = _get_sample_names([str(path) for path in paths])

    assert len(set(names.values())) == 3
    assert names[str(paths[0])] == "a__sample"


def test_copy_inputs_preserves_relative_directories(tmp_path):
    source = tmp_path / "source"
    project = tmp_path / "project"
    paths = [source / "a" / "sample.fasta", source / "b" / "sample.fasta"]
    for path in paths:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(">id\nACGT\n")

    _copy_inputs_to_project([str(path) for path in paths], str(project))

    assert (project / "input" / "a" / "sample.fasta").is_file()
    assert (project / "input" / "b" / "sample.fasta").is_file()


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


@pytest.mark.e2e
def test_run_returns_sequence_object(single_hc_sequence):
    """Test run() returns Sequence when no project_path."""
    result = run(single_hc_sequence)

    assert isinstance(result, Sequence)


@pytest.mark.e2e
def test_run_returns_sequence_list(multiple_hc_sequences):
    """Test run() returns list of Sequences for multiple inputs."""
    result = run(multiple_hc_sequences)

    assert isinstance(result, list)
    assert len(result) == 3
    assert all(isinstance(s, Sequence) for s in result)
    assert [sequence.id for sequence in result] == ["10E8", "10J4", "10M6"]


@pytest.mark.e2e
def test_run_returns_dataframe_when_requested(single_hc_sequence):
    """Test run() returns polars DataFrame when as_dataframe=True."""
    result = run(single_hc_sequence, as_dataframe=True)

    assert isinstance(result, pl.DataFrame)


@pytest.mark.e2e
def test_run_accumulates_dataframes_from_multiple_files(single_hc_sequence, tmp_path):
    input_dir = tmp_path / "inputs"
    for directory, suffix, sequence_id in (
        ("a", "fasta", "first"),
        ("b", "fa", "second"),
    ):
        path = input_dir / directory / f"sample.{suffix}"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(f">{sequence_id}\n{single_hc_sequence.sequence}\n")

    result = run(str(input_dir), as_dataframe=True, n_processes=1)

    assert result.height == 2
    assert result.get_column("sequence_id").to_list() == ["first", "second"]


@pytest.mark.e2e
def test_multi_record_input_with_one_annotation_returns_list(single_hc_sequence):
    unassignable = Sequence("NNNNNNNNNNNNNNNN", id="unassignable")

    result = run([single_hc_sequence, unassignable], n_processes=1)

    assert isinstance(result, list)
    assert [sequence.id for sequence in result] == [single_hc_sequence.id]


@pytest.mark.e2e
def test_five_prime_truncated_read_is_annotated_as_partial(single_hc_sequence):
    truncated = Sequence(
        single_hc_sequence.sequence[90:], id="five_prime_truncated"
    )

    result = run(truncated, n_processes=1)

    assert isinstance(result, Sequence)
    assert result.id == "five_prime_truncated"
    assert result["fwr1"] == ""
    assert result["cdr1"] == "AACGCCTGG"
    assert result["fwr2"] == "ATGACTTGGGTCCGCCAGCCTCCAGGGAAGGGCCTCGAATGGGTTGGTCGT"


@pytest.mark.e2e
def test_run_returns_none_with_project_path(single_hc_sequence, tmp_path):
    """Test run() returns None when project_path provided (writes files)."""
    project_path = str(tmp_path / "test_project")
    result = run(single_hc_sequence, project_path=project_path)

    assert result is None


@pytest.mark.e2e
def test_run_single_sequence_returns_single_not_list(single_hc_sequence):
    """Test that single sequence input returns single Sequence, not list."""
    result = run(single_hc_sequence)

    assert isinstance(result, Sequence)
    assert not isinstance(result, list)


# =============================================
#          OUTPUT FORMAT TESTS
# =============================================


@pytest.mark.e2e
def test_run_rejects_unsupported_output_before_creating_project(
    single_hc_sequence, tmp_path
):
    project_path = tmp_path / "invalid_output"

    with pytest.raises(ValueError, match="Unsupported output format"):
        run(
            single_hc_sequence,
            project_path=str(project_path),
            output_format="csv",
        )

    assert not project_path.exists()


@pytest.mark.e2e
def test_run_creates_airr_output(single_hc_sequence, tmp_path):
    """Test AIRR TSV output file creation."""
    project_path = str(tmp_path / "airr_test")

    run(single_hc_sequence, project_path=project_path, output_format="airr")

    # Check AIRR directory and file exist
    airr_dir = os.path.join(project_path, "airr")
    assert os.path.exists(airr_dir)

    airr_files = [f for f in os.listdir(airr_dir) if f.endswith(".tsv")]
    assert len(airr_files) >= 1
    output_df = pl.read_csv(
        os.path.join(airr_dir, airr_files[0]),
        separator="\t",
        columns=["sequence_id"],
        schema_overrides={"sequence_id": pl.String},
    )
    assert output_df.height == 1


@pytest.mark.e2e
def test_run_creates_parquet_output(single_hc_sequence, tmp_path):
    """Test Parquet output file creation."""
    project_path = str(tmp_path / "parquet_test")

    run(single_hc_sequence, project_path=project_path, output_format="parquet")

    # Check parquet directory and file exist
    parquet_dir = os.path.join(project_path, "parquet")
    assert os.path.exists(parquet_dir)

    parquet_files = [f for f in os.listdir(parquet_dir) if f.endswith(".parquet")]
    assert len(parquet_files) >= 1


@pytest.mark.e2e
def test_run_creates_both_outputs(single_hc_sequence, tmp_path):
    """Test creating both AIRR and Parquet outputs."""
    project_path = str(tmp_path / "both_test")

    run(single_hc_sequence, project_path=project_path, output_format=["airr", "parquet"])

    # Check both directories exist
    assert os.path.exists(os.path.join(project_path, "airr"))
    assert os.path.exists(os.path.join(project_path, "parquet"))


@pytest.mark.e2e
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


@pytest.mark.e2e
def test_run_with_human_database(single_hc_sequence):
    """Test run with human germline database."""
    result = run(single_hc_sequence, germline_database="human")

    assert isinstance(result, Sequence)
    assert result["v_gene"] is not None


# =============================================
#       ANNOTATION VALIDATION TESTS
# =============================================


@pytest.mark.e2e
def test_annotated_sequence_has_v_gene(single_hc_sequence):
    """Test that annotated sequence has V gene assignment."""
    result = run(single_hc_sequence)

    assert result["v_gene"] is not None
    assert "IGHV" in result["v_gene"]


@pytest.mark.e2e
def test_annotated_sequence_has_j_gene(single_hc_sequence):
    """Test that annotated sequence has J gene assignment."""
    result = run(single_hc_sequence)

    assert result["j_gene"] is not None
    assert "IGHJ" in result["j_gene"]


@pytest.mark.e2e
def test_annotated_sequence_has_cdr3(single_hc_sequence):
    """Test that annotated sequence has CDR3."""
    result = run(single_hc_sequence)

    assert result["cdr3"] is not None
    assert len(result["cdr3"]) > 0


@pytest.mark.e2e
def test_annotated_sequence_has_junction(single_hc_sequence):
    """Test that annotated sequence has junction."""
    result = run(single_hc_sequence)

    assert result["junction"] is not None
    assert len(result["junction"]) > 0


@pytest.mark.e2e
def test_annotated_sequence_has_regions(single_hc_sequence):
    """Test all regions (FWR1-4, CDR1-3) are populated."""
    result = run(single_hc_sequence)

    # Check FWR and CDR regions
    regions = ["fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3"]
    for region in regions:
        assert result[region] is not None, f"Region {region} should not be None"

    # FWR4 may be present depending on sequence coverage
    # CDR3 is already checked above


@pytest.mark.e2e
def test_annotated_sequence_has_locus(single_hc_sequence):
    """Test that annotated sequence has locus."""
    result = run(single_hc_sequence)

    assert result["locus"] is not None
    assert result["locus"] == "IGH"  # 10E8 is heavy chain


@pytest.mark.e2e
def test_annotated_sequence_has_productivity(single_hc_sequence):
    """Test that annotated sequence has productivity assessment."""
    result = run(single_hc_sequence)

    # Productive field should be boolean
    assert result["productive"] is not None
    assert isinstance(result["productive"], bool)


@pytest.mark.e2e
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


@pytest.mark.e2e
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


@pytest.mark.e2e
def test_dataframe_row_count_single_sequence(single_hc_sequence):
    """Test that DataFrame has correct row count for a single sequence."""
    result = run(single_hc_sequence, as_dataframe=True)

    assert result.height == 1


@pytest.mark.e2e
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


@pytest.mark.e2e
def test_run_from_fasta_file(small_fasta_file):
    """Test run() with a FASTA file path."""
    result = run(small_fasta_file)

    # Should return a single Sequence (since file has one sequence)
    assert isinstance(result, Sequence)
    assert result["v_gene"] is not None


@pytest.mark.e2e
def test_run_preserves_sequence_id(single_hc_sequence):
    """Test that original sequence ID is preserved or reasonably handled."""
    result = run(single_hc_sequence)

    # Sequence ID should be present (may be transformed by internal processing)
    assert result["sequence_id"] is not None
    assert len(result["sequence_id"]) > 0


@pytest.mark.e2e
def test_run_with_debug_mode(single_hc_sequence, tmp_path):
    """Test run with debug mode enabled."""
    project_path = str(tmp_path / "debug_test")

    # Debug mode should not raise an error
    run(single_hc_sequence, project_path=project_path, debug=True)

    # Log files should exist
    log_dir = os.path.join(project_path, "logs")
    assert os.path.exists(log_dir)
