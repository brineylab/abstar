# Copyright (c) 2025 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

"""
Tests for the CLI entry point (scripts/abstar.py).
"""

import os

import pytest
from click.testing import CliRunner

from abstar.scripts.abstar import cli

# =============================================
#              FIXTURES
# =============================================


@pytest.fixture
def runner():
    """Click CLI test runner."""
    return CliRunner()


@pytest.fixture
def small_fasta(tmp_path):
    """Create a small FASTA file for CLI testing."""
    fasta = tmp_path / "test.fasta"
    fasta.write_text(
        ">10E8\n"
        "GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTGAAGCCTGGAGGATCCCTTAGACTCTCATGTTCAGCCTCTGGTTTCGA"
        "CTTCGATAACGCCTGGATGACTTGGGTCCGCCAGCCTCCAGGGAAGGGCCTCGAATGGGTTGGTCGTATTACGGGTCCAGGTG"
        "AAGGTTGGTCAGTGGACTATGCTGCACCCGTGGAAGGCAGATTTACCATCTCGAGACTCAATTCAATAAATTTCTTATATTTGG"
        "AGATGAACAATTTAAGAATGGAAGACTCAGGCCTTTACTTCTGTGCCCGCACGGGAAAATATTATGATTTTTGGAGTGGCTATC"
        "CGCCGGGAGAAGAATACTTCCAAGACTGGGGCCGGGGCACCCTGGTCACCGTCTCCTCA\n"
    )
    return str(fasta)


# =============================================
#           CLI HELP TESTS
# =============================================


def test_cli_help(runner):
    """Test CLI top-level help."""
    result = runner.invoke(cli, ["--help"])
    assert result.exit_code == 0
    assert "Usage" in result.output


def test_run_help(runner):
    """Test 'run' command help text."""
    result = runner.invoke(cli, ["run", "--help"])
    assert result.exit_code == 0
    assert "INPUT_PATH" in result.output
    assert "PROJECT_PATH" in result.output
    assert "--germline_database" in result.output
    assert "--receptor" in result.output
    assert "--output_format" in result.output


def test_build_germline_database_help(runner):
    """Test 'build_germline_database' command help text."""
    result = runner.invoke(cli, ["build_germline_database", "--help"])
    assert result.exit_code == 0
    assert "NAME" in result.output
    assert "--fasta" in result.output
    assert "--json" in result.output
    assert "--receptor" in result.output


# =============================================
#         CLI RUN COMMAND TESTS
# =============================================


def test_run_command_with_output(runner, small_fasta, tmp_path):
    """Test the run command produces output files."""
    project_path = str(tmp_path / "project")
    result = runner.invoke(
        cli,
        [
            "run",
            small_fasta,
            project_path,
            "--output_format",
            "airr",
            "--quiet",
        ],
    )
    assert result.exit_code == 0, f"CLI failed: {result.output}"
    assert os.path.exists(os.path.join(project_path, "airr"))


def test_run_command_parquet_output(runner, small_fasta, tmp_path):
    """Test run command with parquet output format."""
    project_path = str(tmp_path / "project")
    result = runner.invoke(
        cli,
        [
            "run",
            small_fasta,
            project_path,
            "--output_format",
            "parquet",
            "--quiet",
        ],
    )
    assert result.exit_code == 0, f"CLI failed: {result.output}"
    assert os.path.exists(os.path.join(project_path, "parquet"))


def test_run_command_both_outputs(runner, small_fasta, tmp_path):
    """Test run command with both airr and parquet output."""
    project_path = str(tmp_path / "project")
    result = runner.invoke(
        cli,
        [
            "run",
            small_fasta,
            project_path,
            "-o",
            "airr",
            "-o",
            "parquet",
            "--quiet",
        ],
    )
    assert result.exit_code == 0, f"CLI failed: {result.output}"
    assert os.path.exists(os.path.join(project_path, "airr"))
    assert os.path.exists(os.path.join(project_path, "parquet"))


def test_run_command_invalid_receptor(runner, small_fasta, tmp_path):
    """Test run command with invalid receptor type."""
    project_path = str(tmp_path / "project")
    result = runner.invoke(
        cli,
        [
            "run",
            small_fasta,
            project_path,
            "--receptor",
            "invalid",
            "--quiet",
        ],
    )
    assert result.exit_code != 0


def test_run_command_missing_args(runner):
    """Test run command with missing required arguments."""
    result = runner.invoke(cli, ["run"])
    assert result.exit_code != 0


def test_run_command_tcr_receptor(runner, small_fasta, tmp_path):
    """Test run command with TCR receptor type."""
    project_path = str(tmp_path / "project")
    result = runner.invoke(
        cli,
        [
            "run",
            small_fasta,
            project_path,
            "--receptor",
            "tcr",
            "--quiet",
        ],
    )
    # Should run without error even if sequence is BCR
    assert result.exit_code == 0, f"CLI failed: {result.output}"
