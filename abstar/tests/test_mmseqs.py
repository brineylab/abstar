# Copyright (c) 2025 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

"""
Tests for the MMseqs2 germline assigner.
"""

import os

import polars as pl
import pytest

from ..assigners.mmseqs import MMseqs


# =============================================
#              FIXTURES
# =============================================


@pytest.fixture
def mmseqs_instance(temp_directories):
    """Create an MMseqs assigner instance."""
    output_dir, log_dir = temp_directories
    return MMseqs(
        output_directory=output_dir,
        log_directory=log_dir,
        germdb_name="human",
        receptor="bcr",
        chunksize=1000,  # Small chunksize for testing
    )


@pytest.fixture
def mmseqs_tcr_instance(temp_directories):
    """Create an MMseqs assigner instance for TCR."""
    output_dir, log_dir = temp_directories
    return MMseqs(
        output_directory=output_dir,
        log_directory=log_dir,
        germdb_name="human",
        receptor="tcr",
        chunksize=1000,
    )


# =============================================
#          INITIALIZATION TESTS
# =============================================


def test_mmseqs_initialization(temp_directories):
    """Test that MMseqs initializes with correct attributes."""
    output_dir, log_dir = temp_directories
    mmseqs = MMseqs(
        output_directory=output_dir,
        log_directory=log_dir,
        germdb_name="human",
        receptor="bcr",
    )

    assert mmseqs.output_directory == output_dir
    assert mmseqs.log_directory == log_dir
    assert mmseqs.germdb_name == "human"
    assert mmseqs.receptor == "bcr"
    assert mmseqs.germdb_path is not None
    assert os.path.exists(mmseqs.germdb_path)
    assert mmseqs.to_delete == []


def test_mmseqs_initialization_with_threads(temp_directories):
    """Test MMseqs initialization with custom thread count."""
    output_dir, log_dir = temp_directories
    mmseqs = MMseqs(
        output_directory=output_dir,
        log_directory=log_dir,
        germdb_name="human",
        receptor="bcr",
        threads=4,
    )

    assert mmseqs.threads == 4


def test_mmseqs_initialization_with_chunksize(temp_directories):
    """Test MMseqs initialization with custom chunksize."""
    output_dir, log_dir = temp_directories
    mmseqs = MMseqs(
        output_directory=output_dir,
        log_directory=log_dir,
        germdb_name="human",
        receptor="bcr",
        chunksize=500,
    )

    assert mmseqs.chunksize == 500


def test_mmseqs_initialization_tcr(temp_directories):
    """Test MMseqs initialization for TCR receptor."""
    output_dir, log_dir = temp_directories
    mmseqs = MMseqs(
        output_directory=output_dir,
        log_directory=log_dir,
        germdb_name="human",
        receptor="tcr",
    )

    assert mmseqs.receptor == "tcr"
    assert mmseqs.germdb_path is not None


@pytest.mark.skip(reason="Mouse germline database not installed in test environment")
def test_mmseqs_initialization_mouse_database(temp_directories):
    """Test MMseqs initialization with mouse germline database."""
    output_dir, log_dir = temp_directories
    mmseqs = MMseqs(
        output_directory=output_dir,
        log_directory=log_dir,
        germdb_name="mouse",
        receptor="bcr",
    )

    assert mmseqs.germdb_name == "mouse"
    assert mmseqs.germdb_path is not None


# =============================================
#        INPUT PREPARATION TESTS
# =============================================


def test_prepare_input_files_fasta(mmseqs_instance, small_fasta_file):
    """Test prepare_input_files with FASTA input."""
    # Set sample_name which is normally set in __call__
    mmseqs_instance.sample_name = "test_sequences"

    fasta_paths, tsv_paths, sequence_count = mmseqs_instance.prepare_input_files(
        small_fasta_file, chunksize=1000
    )

    # Check return values
    assert len(fasta_paths) == 1
    assert len(tsv_paths) == 1
    assert sequence_count == 1

    # Check FASTA file exists and has content
    assert os.path.exists(fasta_paths[0])
    with open(fasta_paths[0]) as f:
        content = f.read()
        assert ">10E8" in content

    # Check TSV file has correct columns
    assert os.path.exists(tsv_paths[0])
    df = pl.read_csv(tsv_paths[0], separator="\t")
    assert "sequence_id" in df.columns
    assert "sequence_input" in df.columns
    assert "quality" in df.columns


def test_prepare_input_files_multiple_sequences(mmseqs_instance, multi_sequence_fasta_file):
    """Test prepare_input_files with multiple sequences."""
    mmseqs_instance.sample_name = "multi_sequences"

    fasta_paths, tsv_paths, sequence_count = mmseqs_instance.prepare_input_files(
        multi_sequence_fasta_file, chunksize=1000
    )

    assert sequence_count == 3

    # Check all sequences are in the TSV
    df = pl.read_csv(tsv_paths[0], separator="\t")
    assert df.height == 3


def test_prepare_input_files_fastq(mmseqs_instance, fastq_test_path):
    """Test prepare_input_files with FASTQ input (quality scores preserved)."""
    mmseqs_instance.sample_name = "test"

    fasta_paths, tsv_paths, sequence_count = mmseqs_instance.prepare_input_files(
        fastq_test_path, chunksize=1000
    )

    # Check that quality scores are present in TSV
    df = pl.read_csv(tsv_paths[0], separator="\t")
    quality_col = df["quality"]
    # At least some quality scores should be non-empty
    non_empty_quality = quality_col.filter(quality_col.str.len_chars() > 0)
    assert len(non_empty_quality) > 0


# =============================================
#       QUERY FASTA BUILDER TESTS
# =============================================


def test_build_jquery_fasta(mmseqs_instance, tmp_path):
    """Test J-gene query FASTA generation from V results."""
    # Create mock V result DataFrame
    vresult_df = pl.DataFrame({
        "v_query": ["seq1", "seq2"],
        "v_qstart": [1, 1],
        "v_qend": [100, 80],
        "v_qseq": [
            "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC",
            "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"
        ],
    })

    jquery_path = str(tmp_path / "jquery.fasta")
    mmseqs_instance.build_jquery_fasta(vresult_df, jquery_path)

    # Check file was created
    assert os.path.exists(jquery_path)

    # Check content
    with open(jquery_path) as f:
        content = f.read()
        # Should have sequences that follow the V alignment
        assert ">seq1" in content or ">seq2" in content


def test_build_jquery_fasta_filters_short_sequences(mmseqs_instance, tmp_path):
    """Test that J-gene query FASTA excludes sequences < 14 nt."""
    # Create mock V result where remaining sequence is too short
    vresult_df = pl.DataFrame({
        "v_query": ["short_seq"],
        "v_qstart": [1],
        "v_qend": [10],  # Only 10 nt remain after V end
        "v_qseq": ["ATGCATGCATGC"],  # 12 nt total, only 2 remain after position 10
    })

    jquery_path = str(tmp_path / "jquery_short.fasta")
    mmseqs_instance.build_jquery_fasta(vresult_df, jquery_path)

    # File should be empty or have no sequences
    with open(jquery_path) as f:
        content = f.read().strip()
        # Short sequences should be filtered out
        assert content == "" or ">" not in content


def test_build_dquery_fasta_heavy_chain_only(mmseqs_instance, tmp_path):
    """Test D-gene query FASTA only includes IGH/TRA/TRD."""
    # Create mock VJ result with both heavy and light chain
    vjresult_df = pl.DataFrame({
        "v_query": ["heavy_seq", "light_seq"],
        "v_call": ["IGHV3-23*01", "IGKV1-39*01"],  # IGH is heavy, IGK is light
        "j_qstart": [50, 40],
        "j_qend": [100, 80],
        "j_qseq": [
            "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC",
            "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"
        ],
    })

    dquery_path = str(tmp_path / "dquery.fasta")
    mmseqs_instance.build_dquery_fasta(vjresult_df, dquery_path)

    # Check content - should only have heavy chain sequence
    with open(dquery_path) as f:
        content = f.read()
        # Heavy chain should be included
        if content.strip():  # If not empty
            assert "heavy_seq" in content or content.count(">") >= 1
            # Light chain should NOT be included
            assert "light_seq" not in content


def test_build_cquery_fasta(mmseqs_instance, tmp_path):
    """Test C-gene query FASTA generation from VJ results."""
    vjresult_df = pl.DataFrame({
        "v_query": ["seq1"],
        "j_qstart": [1],
        "j_qend": [50],
        "j_qseq": [
            "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"
        ],
    })

    cquery_path = str(tmp_path / "cquery.fasta")
    mmseqs_instance.build_cquery_fasta(vjresult_df, cquery_path)

    # Check file was created
    assert os.path.exists(cquery_path)


# =============================================
#      MMSEQS RESULT PARSING TESTS
# =============================================


def test_read_mmseqs_results_empty_file_has_stable_schema(mmseqs_instance, tmp_path):
    """Test empty MMseqs outputs produce an empty frame with expected schema."""
    dresult_path = tmp_path / "empty_dresult.tsv"
    dresult_path.write_text("")

    dresult_df = mmseqs_instance._read_mmseqs_results(str(dresult_path), prefix="d")

    assert dresult_df.is_empty()
    assert dresult_df.schema["d_query"] == pl.String
    assert dresult_df.schema["d_call"] == pl.String
    assert dresult_df.schema["d_support"] == pl.Float64


def test_read_mmseqs_results_keeps_best_hit_per_query(mmseqs_instance, tmp_path):
    """Test MMseqs parser keeps highest nident assignment for duplicate queries."""
    dresult_path = tmp_path / "dresult.tsv"
    dresult_path.write_text(
        "\n".join(
            [
                "query\ttarget\tevalue\tqstart\tqend\tqseq\tnident",
                "seq1\tIGHD1-1*01\t1e-5\t1\t10\tACGTACGTAA\t8",
                "seq1\tIGHD2-2*01\t1e-5\t1\t10\tACGTACGTAA\t10",
                "seq2\tIGHD3-3*01\t1e-5\t2\t11\tCGTACGTAAC\t9",
            ]
        )
    )

    dresult_df = mmseqs_instance._read_mmseqs_results(str(dresult_path), prefix="d")

    assert dresult_df.height == 2
    assert dresult_df.schema["d_query"] == pl.String
    assert (
        dresult_df.filter(pl.col("d_query") == "seq1").select("d_call").item()
        == "IGHD2-2*01"
    )


# =============================================
#        FULL ASSIGNMENT TESTS
# =============================================


def test_mmseqs_call_returns_parquet_path_and_count(mmseqs_instance, small_fasta_file):
    """Test __call__ method returns correct path and count."""
    assigned_path, sequence_count = mmseqs_instance(small_fasta_file)

    # Check return values
    assert assigned_path is not None
    assert os.path.exists(assigned_path)
    assert assigned_path.endswith(".parquet")
    assert sequence_count == 1


def test_mmseqs_call_parquet_has_required_columns(mmseqs_instance, small_fasta_file):
    """Test that output parquet has required columns."""
    assigned_path, _ = mmseqs_instance(small_fasta_file)

    df = pl.read_parquet(assigned_path)

    # Check required columns exist
    required_columns = [
        "sequence_id",
        "sequence_input",
        "quality",
        "rev_comp",
        "v_call",
        "v_support",
        "d_call",
        "d_support",
        "j_call",
        "j_support",
        "c_call",
        "c_support",
    ]

    for col in required_columns:
        assert col in df.columns, f"Missing required column: {col}"


def test_mmseqs_assigns_v_gene(mmseqs_instance, small_fasta_file):
    """Test that V gene is assigned for a valid antibody sequence."""
    assigned_path, _ = mmseqs_instance(small_fasta_file)

    df = pl.read_parquet(assigned_path)

    # Check V gene assignment
    v_calls = df["v_call"].to_list()
    assert len(v_calls) > 0
    assert v_calls[0] is not None
    assert "IGHV" in v_calls[0]  # 10E8 is a heavy chain


def test_mmseqs_assigns_j_gene(mmseqs_instance, small_fasta_file):
    """Test that J gene is assigned for a valid antibody sequence."""
    assigned_path, _ = mmseqs_instance(small_fasta_file)

    df = pl.read_parquet(assigned_path)

    # Check J gene assignment
    j_calls = df["j_call"].to_list()
    assert len(j_calls) > 0
    assert j_calls[0] is not None
    assert "IGHJ" in j_calls[0]


def test_mmseqs_preserves_scientific_notation_like_sequence_ids(
    mmseqs_instance, single_hc_sequence, tmp_path
):
    """Sequence IDs like 10E8 should be preserved exactly through assignment joins."""
    input_fasta = tmp_path / "scientific_like_id.fasta"
    with open(input_fasta, "w") as f:
        f.write(f">{single_hc_sequence.id}\n{single_hc_sequence.sequence}\n")

    assigned_path, _ = mmseqs_instance(str(input_fasta))
    df = pl.read_parquet(assigned_path)

    assert df.height >= 1
    sequence_ids = df["sequence_id"].to_list()
    assert single_hc_sequence.id in sequence_ids
    assert "1000000000.0" not in sequence_ids


def test_mmseqs_multiple_sequences(mmseqs_instance, multi_sequence_fasta_file):
    """Test assignment of multiple sequences."""
    assigned_path, sequence_count = mmseqs_instance(multi_sequence_fasta_file)

    assert sequence_count == 3

    df = pl.read_parquet(assigned_path)
    # Some or all sequences should be successfully assigned
    assert df.height >= 1


# =============================================
#            CLEANUP TESTS
# =============================================


def test_cleanup_removes_temp_files(mmseqs_instance, small_fasta_file):
    """Test that cleanup() removes files in to_delete list."""
    # Run assignment to populate to_delete
    mmseqs_instance(small_fasta_file)

    # Manually add a file to to_delete for testing
    test_file = os.path.join(mmseqs_instance.output_directory, "test_cleanup.txt")
    with open(test_file, "w") as f:
        f.write("test")
    mmseqs_instance.to_delete.append(test_file)

    # Verify file exists
    assert os.path.exists(test_file)

    # Run cleanup
    mmseqs_instance.cleanup()

    # Verify file was deleted
    assert not os.path.exists(test_file)
    assert mmseqs_instance.to_delete == []


# =============================================
#          OUTPUT ASSEMBLY TESTS
# =============================================


def test_assemble_output_files(mmseqs_instance, tmp_path):
    """Test merging multiple parquet/CSV files into single outputs."""
    # Create mock assigned parquet files
    df1 = pl.DataFrame({
        "sequence_id": ["seq1"],
        "sequence_input": ["ATGC"],
        "quality": [""],
        "rev_comp": [False],
        "v_call": ["IGHV3-23*01"],
        "v_support": [1e-10],
        "d_call": [None],
        "d_support": [None],
        "j_call": ["IGHJ4*02"],
        "j_support": [1e-5],
        "c_call": [None],
        "c_support": [None],
    })

    df2 = pl.DataFrame({
        "sequence_id": ["seq2"],
        "sequence_input": ["GCTA"],
        "quality": [""],
        "rev_comp": [False],
        "v_call": ["IGHV1-2*01"],
        "v_support": [1e-12],
        "d_call": [None],
        "d_support": [None],
        "j_call": ["IGHJ6*01"],
        "j_support": [1e-6],
        "c_call": [None],
        "c_support": [None],
    })

    assigned_path1 = str(tmp_path / "assigned1.parquet")
    assigned_path2 = str(tmp_path / "assigned2.parquet")
    df1.write_parquet(assigned_path1)
    df2.write_parquet(assigned_path2)

    # Create mock unassigned CSV files
    unassigned_df = pl.DataFrame({
        "sequence_id": ["unassigned1"],
        "sequence_input": ["NNNN"],
    })
    unassigned_path1 = str(tmp_path / "unassigned1.csv")
    unassigned_path2 = str(tmp_path / "unassigned2.csv")
    unassigned_df.write_csv(unassigned_path1)
    unassigned_df.write_csv(unassigned_path2)

    # Set sample name for output
    mmseqs_instance.sample_name = "test_assembly"

    # Assemble outputs
    assembled_path = mmseqs_instance.assemble_output_files(
        [assigned_path1, assigned_path2],
        [unassigned_path1, unassigned_path2]
    )

    # Check assembled file
    assert os.path.exists(assembled_path)
    assembled_df = pl.read_parquet(assembled_path)
    assert assembled_df.height == 2  # Both sequences combined
