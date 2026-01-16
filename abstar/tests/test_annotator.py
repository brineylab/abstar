# Copyright (c) 2025 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

"""
Tests for the annotation orchestrator module.
"""

import os

import polars as pl
import pytest

from ..annotation.annotator import annotate, annotate_single_sequence
from ..annotation.antibody import Antibody
from ..core.abstar import run


# =============================================
#              FIXTURES
# =============================================


@pytest.fixture
def mock_assignment_parquet(tmp_path):
    """Create a mock assignment parquet file with V/J assignments."""
    # Use a real sequence that will successfully annotate
    df = pl.DataFrame({
        "sequence_id": ["10E8_test"],
        "sequence_input": [
            "GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTGAAGCCTGGAGGATCCCTTAGACTCTCATGTTCAGCCTCTGGTTTCGACTTCGATAACGCCTGGATGACTTGGGTCCGCCAGCCTCCAGGGAAGGGCCTCGAATGGGTTGGTCGTATTACGGGTCCAGGTGAAGGTTGGTCAGTGGACTATGCTGCACCCGTGGAAGGCAGATTTACCATCTCGAGACTCAATTCAATAAATTTCTTATATTTGGAGATGAACAATTTAAGAATGGAAGACTCAGGCCTTTACTTCTGTGCCCGCACGGGAAAATATTATGATTTTTGGAGTGGCTATCCGCCGGGAGAAGAATACTTCCAAGACTGGGGCCGGGGCACCCTGGTCACCGTCTCCTCA"
        ],
        "quality": [""],
        "rev_comp": [False],
        "v_call": ["IGHV3-23*01"],
        "v_support": [1e-50],
        "d_call": ["IGHD3-10*01"],
        "d_support": [1e-5],
        "j_call": ["IGHJ4*02"],
        "j_support": [1e-20],
        "c_call": [None],
        "c_support": [None],
    })

    parquet_path = str(tmp_path / "assignment.parquet")
    df.write_parquet(parquet_path)
    return parquet_path


@pytest.fixture
def mock_light_chain_parquet(tmp_path):
    """Create a mock assignment parquet file for light chain."""
    df = pl.DataFrame({
        "sequence_id": ["lc_test"],
        "sequence_input": [
            "GACATCCAGATGACCCAGTCTCCATCCTCACTGTCTGCATCTGTAGGAGACAGAGTCACCATCACTTGTCGGGCGAGTCAGGGTATTAGCAGCTGGTTAGCCTGGTATCAGCAGAAACCAGGGAAAGCCCCTAAGCTCCTGATCTATGCTGCATCCAGTTTGCAAAGTGGGGTCCCATCAAGGTTCAGCGGCAGTGGATCTGGGACAGATTTCACTCTCACCATCAGCAGCCTGCAGCCTGAAGATTTTGCAACTTACTATTGTCAACAGGCTAACAGTTTCCCGCTCACTTTCGGCGGAGGGACCAAGGTGGAGATCAAACGA"
        ],
        "quality": [""],
        "rev_comp": [False],
        "v_call": ["IGKV1-39*01"],
        "v_support": [1e-40],
        "d_call": [None],
        "d_support": [None],
        "j_call": ["IGKJ1*01"],
        "j_support": [1e-15],
        "c_call": [None],
        "c_support": [None],
    })

    parquet_path = str(tmp_path / "lc_assignment.parquet")
    df.write_parquet(parquet_path)
    return parquet_path


@pytest.fixture
def annotated_heavy_chain_result(single_hc_sequence):
    """Get a fully annotated heavy chain sequence using run()."""
    return run(single_hc_sequence)


# =============================================
#        annotate() FUNCTION TESTS
# =============================================


def test_annotate_returns_correct_files(mock_assignment_parquet, tmp_path):
    """Test annotate() returns (output_file, failed_logfile, succeeded_logfile)."""
    output_dir = str(tmp_path / "output")
    log_dir = str(tmp_path / "logs")
    os.makedirs(output_dir)
    os.makedirs(log_dir)

    output_file, failed, succeeded = annotate(
        input_file=mock_assignment_parquet,
        output_directory=output_dir,
        germline_database="human",
        log_directory=log_dir,
    )

    assert output_file is not None
    assert os.path.exists(output_file)
    assert output_file.endswith(".parquet")

    # Failed log file should exist (may be empty)
    assert failed is not None
    assert os.path.exists(failed)


def test_annotate_creates_parquet_output(mock_assignment_parquet, tmp_path):
    """Test output parquet has correct schema."""
    output_dir = str(tmp_path / "output")
    log_dir = str(tmp_path / "logs")
    os.makedirs(output_dir)
    os.makedirs(log_dir)

    output_file, _, _ = annotate(
        input_file=mock_assignment_parquet,
        output_directory=output_dir,
        germline_database="human",
        log_directory=log_dir,
    )

    # Read output and check columns
    df = pl.read_parquet(output_file)

    # Check key columns exist
    required_columns = [
        "sequence_id",
        "v_gene",
        "j_gene",
        "locus",
        "productive",
    ]
    for col in required_columns:
        assert col in df.columns, f"Missing column: {col}"


def test_annotate_with_debug_logs_succeeded(mock_assignment_parquet, tmp_path):
    """Test debug=True logs succeeded sequences."""
    output_dir = str(tmp_path / "output")
    log_dir = str(tmp_path / "logs")
    os.makedirs(output_dir)
    os.makedirs(log_dir)

    output_file, failed, succeeded = annotate(
        input_file=mock_assignment_parquet,
        output_directory=output_dir,
        germline_database="human",
        log_directory=log_dir,
        debug=True,
    )

    # With debug=True, succeeded log should be returned
    assert succeeded is not None
    assert os.path.exists(succeeded)


def test_annotate_light_chain(mock_light_chain_parquet, tmp_path):
    """Test annotation of light chain (no D gene)."""
    output_dir = str(tmp_path / "output")
    log_dir = str(tmp_path / "logs")
    os.makedirs(output_dir)
    os.makedirs(log_dir)

    output_file, _, _ = annotate(
        input_file=mock_light_chain_parquet,
        output_directory=output_dir,
        germline_database="human",
        log_directory=log_dir,
    )

    df = pl.read_parquet(output_file)

    # Check it's a light chain
    if df.height > 0:
        locus = df["locus"][0]
        assert locus in ["IGK", "IGL"]

        # D gene should be None for light chain
        d_gene = df["d_gene"][0]
        assert d_gene is None


# =============================================
#    annotate_single_sequence() TESTS
# =============================================


def test_annotate_produces_v_gene(annotated_heavy_chain_result):
    """Test full annotation of heavy chain sequence produces V gene."""
    result = annotated_heavy_chain_result

    assert result["v_gene"] is not None
    assert "IGHV" in result["v_gene"]


def test_annotate_produces_j_gene(annotated_heavy_chain_result):
    """Test annotation produces J gene."""
    result = annotated_heavy_chain_result

    assert result["j_gene"] is not None
    assert "IGHJ" in result["j_gene"]


def test_annotate_produces_locus(annotated_heavy_chain_result):
    """Test annotation produces locus."""
    result = annotated_heavy_chain_result

    assert result["locus"] == "IGH"


# =============================================
#       V-GENE PROCESSING TESTS
# =============================================


def test_vgene_alignment_boundaries(annotated_heavy_chain_result):
    """Test V-gene realignment produces correct boundaries."""
    result = annotated_heavy_chain_result

    assert result["v_sequence_start"] is not None
    assert result["v_sequence_end"] is not None
    assert result["v_germline_start"] is not None
    assert result["v_germline_end"] is not None

    # Start should be before end
    assert result["v_sequence_start"] < result["v_sequence_end"]


def test_vgene_sequences_populated(annotated_heavy_chain_result):
    """Test V-gene sequences are populated."""
    result = annotated_heavy_chain_result

    assert result["v_sequence"] is not None
    assert result["v_germline"] is not None
    assert len(result["v_sequence"]) > 0
    assert len(result["v_germline"]) > 0


def test_vgene_gapped_sequences(annotated_heavy_chain_result):
    """Test IMGT-gapped V-gene sequences are generated."""
    result = annotated_heavy_chain_result

    assert result["v_germline_gapped"] is not None
    assert result["v_sequence_gapped"] is not None
    # IMGT-gapped sequences contain dots
    assert "." in result["v_germline_gapped"]


# =============================================
#       JUNCTION/CDR3 TESTS
# =============================================


def test_junction_identification(annotated_heavy_chain_result):
    """Test junction sequence is identified."""
    result = annotated_heavy_chain_result

    # Junction sequence should be identified
    assert result["junction"] is not None
    assert len(result["junction"]) > 0
    # Junction should be nucleotide sequence (divisible by 3 for proper translation)
    assert len(result["junction"]) % 3 == 0


def test_cdr3_extraction(annotated_heavy_chain_result):
    """Test CDR3 is correctly extracted from junction."""
    result = annotated_heavy_chain_result

    assert result["cdr3"] is not None
    assert result["cdr3_aa"] is not None

    # CDR3 is junction without first and last codons
    junction = result["junction"]
    cdr3 = result["cdr3"]
    assert cdr3 == junction[3:-3]


def test_cdr3_length_calculation(annotated_heavy_chain_result):
    """Test CDR3 length is AA length."""
    result = annotated_heavy_chain_result

    cdr3_aa = result["cdr3_aa"]
    cdr3_length = result["cdr3_length"]

    assert cdr3_length == len(cdr3_aa)


def test_junction_aa_present(annotated_heavy_chain_result):
    """Test junction amino acid sequence is present."""
    result = annotated_heavy_chain_result

    assert result["junction_aa"] is not None
    assert len(result["junction_aa"]) > 0


# =============================================
#          REGION EXTRACTION TESTS
# =============================================


def test_all_v_regions_populated(annotated_heavy_chain_result):
    """Test FWR1-3 and CDR1-2 regions are populated."""
    result = annotated_heavy_chain_result

    v_regions = ["fwr1", "cdr1", "fwr2", "cdr2", "fwr3"]
    for region in v_regions:
        assert result[region] is not None, f"V region {region} should not be None"
        assert result[f"{region}_aa"] is not None, f"V region {region}_aa should not be None"


def test_cdr3_region_populated(annotated_heavy_chain_result):
    """Test CDR3 region is populated."""
    result = annotated_heavy_chain_result

    assert result["cdr3"] is not None
    assert result["cdr3_aa"] is not None


def test_fwr4_region_populated(annotated_heavy_chain_result):
    """Test FWR4 region is populated."""
    result = annotated_heavy_chain_result

    assert result["fwr4"] is not None
    assert result["fwr4_aa"] is not None


# =============================================
#          MUTATION ANNOTATION TESTS
# =============================================


def test_v_mutations_annotated(annotated_heavy_chain_result):
    """Test V-gene mutations are identified."""
    result = annotated_heavy_chain_result

    # Mutation count should exist
    assert result["v_mutation_count"] is not None
    assert isinstance(result["v_mutation_count"], int)

    # If there are mutations, the mutation string should be populated
    if result["v_mutation_count"] > 0:
        assert result["v_mutations"] is not None


def test_v_identity_calculation(annotated_heavy_chain_result):
    """Test V-gene identity is correctly calculated."""
    result = annotated_heavy_chain_result

    v_identity = result["v_identity"]
    v_identity_aa = result["v_identity_aa"]

    assert v_identity is not None
    assert v_identity_aa is not None

    # Identity should be between 0 and 1
    assert 0 <= v_identity <= 1
    assert 0 <= v_identity_aa <= 1


# =============================================
#       PRODUCTIVITY ASSESSMENT TESTS
# =============================================


def test_productive_sequence(annotated_heavy_chain_result):
    """Test productive sequence is marked as such."""
    result = annotated_heavy_chain_result

    # Check productivity fields exist
    assert result["productive"] is not None
    assert isinstance(result["productive"], bool)

    # 10E8 is a functional antibody, should be productive
    assert result["productive"] is True


def test_stop_codon_field_exists(annotated_heavy_chain_result):
    """Test stop_codon field exists."""
    result = annotated_heavy_chain_result

    # stop_codon field should exist
    assert "stop_codon" in result.annotations


def test_productivity_issues_field_exists(annotated_heavy_chain_result):
    """Test productivity_issues field exists."""
    result = annotated_heavy_chain_result

    # productivity_issues should be a list (possibly empty)
    issues = result["productivity_issues"]
    assert issues is not None


# =============================================
#           MASK GENERATION TESTS
# =============================================


def test_cdr_mask_generated(annotated_heavy_chain_result):
    """Test CDR mask is generated correctly."""
    result = annotated_heavy_chain_result

    cdr_mask = result["cdr_mask"]
    assert cdr_mask is not None

    # Mask length should match sequence length
    sequence = result["sequence"]
    assert len(cdr_mask) == len(sequence)


def test_cdr_mask_aa_generated(annotated_heavy_chain_result):
    """Test CDR mask AA is generated."""
    result = annotated_heavy_chain_result

    cdr_mask_aa = result["cdr_mask_aa"]
    assert cdr_mask_aa is not None


def test_gene_segment_mask_generated(annotated_heavy_chain_result):
    """Test gene segment mask is generated."""
    result = annotated_heavy_chain_result

    gene_mask = result["gene_segment_mask"]
    assert gene_mask is not None


def test_nongermline_mask_generated(annotated_heavy_chain_result):
    """Test non-germline mask is generated."""
    result = annotated_heavy_chain_result

    ng_mask = result["nongermline_mask"]
    assert ng_mask is not None


# =============================================
#         ALIGNMENT OUTPUT TESTS
# =============================================


def test_sequence_alignment_populated(annotated_heavy_chain_result):
    """Test sequence alignment is populated."""
    result = annotated_heavy_chain_result

    assert result["sequence_alignment"] is not None
    assert result["germline_alignment"] is not None


def test_sequence_aa_alignment_populated(annotated_heavy_chain_result):
    """Test sequence AA alignment is populated."""
    result = annotated_heavy_chain_result

    assert result["sequence_alignment_aa"] is not None
    assert result["germline_alignment_aa"] is not None


# =============================================
#         ASSEMBLED SEQUENCE TESTS
# =============================================


def test_full_sequence_assembled(annotated_heavy_chain_result):
    """Test full VDJ sequence is assembled."""
    result = annotated_heavy_chain_result

    assert result["sequence"] is not None
    assert result["germline"] is not None
    assert len(result["sequence"]) > 0


def test_sequence_aa_assembled(annotated_heavy_chain_result):
    """Test AA sequences are assembled."""
    result = annotated_heavy_chain_result

    assert result["sequence_aa"] is not None
    assert result["germline_aa"] is not None


def test_complete_vdj_field(annotated_heavy_chain_result):
    """Test complete_vdj field is populated."""
    result = annotated_heavy_chain_result

    # Field should exist (may be True or False)
    assert "complete_vdj" in result.annotations
