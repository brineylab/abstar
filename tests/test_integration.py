# Copyright (c) 2025 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

"""
Integration tests for the abstar pipeline.

These tests run the full annotation pipeline on real-world antibody sequence datasets
and validate output correctness. They are slower than unit tests and should be run
selectively during development.
"""

import os

import polars as pl
import pytest
from abutils import Sequence

from abstar.core.abstar import run

TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), "test_data")


# =============================================
#              FIXTURES
# =============================================


@pytest.fixture(scope="module")
def hiv_bnab_hc_results():
    """Run pipeline on HIV bnAb heavy chains (module-scoped for performance)."""
    return run(os.path.join(TEST_DATA_DIR, "test_hiv_bnab_hcs.fasta"))


@pytest.fixture(scope="module")
def hiv_bnab_lc_results():
    """Run pipeline on HIV bnAb light chains (module-scoped for performance)."""
    return run(os.path.join(TEST_DATA_DIR, "test_hiv_bnab_lcs.fasta"))


# =============================================
#      HIV BNAB HEAVY CHAIN TESTS
# =============================================


@pytest.mark.slow
@pytest.mark.integration
class TestHIVBnAbHeavyChains:
    """Integration tests using HIV broadly neutralizing antibody heavy chain sequences."""

    def test_returns_list(self, hiv_bnab_hc_results):
        """Pipeline returns a list for multi-sequence input."""
        assert isinstance(hiv_bnab_hc_results, list)

    def test_result_count(self, hiv_bnab_hc_results):
        """All sequences are processed (some may fail annotation)."""
        assert len(hiv_bnab_hc_results) > 0

    def test_all_results_are_sequences(self, hiv_bnab_hc_results):
        """All results are Sequence objects."""
        assert all(isinstance(s, Sequence) for s in hiv_bnab_hc_results)

    def test_all_have_v_gene(self, hiv_bnab_hc_results):
        """All annotated sequences have V gene assignment."""
        for seq in hiv_bnab_hc_results:
            assert seq["v_gene"] is not None, f"{seq.id} missing v_gene"

    def test_all_are_heavy_chain(self, hiv_bnab_hc_results):
        """All sequences are assigned to IGH locus."""
        for seq in hiv_bnab_hc_results:
            assert seq["locus"] == "IGH", f"{seq.id} has locus {seq['locus']}, expected IGH"

    def test_all_have_j_gene(self, hiv_bnab_hc_results):
        """All annotated sequences have J gene assignment."""
        for seq in hiv_bnab_hc_results:
            assert seq["j_gene"] is not None, f"{seq.id} missing j_gene"

    def test_all_have_cdr3(self, hiv_bnab_hc_results):
        """All annotated sequences have CDR3."""
        for seq in hiv_bnab_hc_results:
            assert seq["cdr3"] is not None, f"{seq.id} missing cdr3"
            assert len(seq["cdr3"]) > 0

    def test_v_genes_are_ighv(self, hiv_bnab_hc_results):
        """All V genes start with IGHV."""
        for seq in hiv_bnab_hc_results:
            assert seq["v_gene"].startswith("IGHV"), (
                f"{seq.id} v_gene={seq['v_gene']} doesn't start with IGHV"
            )

    def test_j_genes_are_ighj(self, hiv_bnab_hc_results):
        """All J genes start with IGHJ."""
        for seq in hiv_bnab_hc_results:
            assert seq["j_gene"].startswith("IGHJ"), (
                f"{seq.id} j_gene={seq['j_gene']} doesn't start with IGHJ"
            )

    def test_identity_in_valid_range(self, hiv_bnab_hc_results):
        """V gene identity is between 0 and 1."""
        for seq in hiv_bnab_hc_results:
            v_id = seq["v_identity"]
            assert v_id is not None, f"{seq.id} missing v_identity"
            assert 0 < v_id <= 1, f"{seq.id} v_identity={v_id} out of range"

    def test_regions_populated(self, hiv_bnab_hc_results):
        """All V-gene regions (FWR1-3, CDR1-2) are populated."""
        for seq in hiv_bnab_hc_results:
            for region in ["fwr1", "cdr1", "fwr2", "cdr2", "fwr3"]:
                assert seq[region] is not None, f"{seq.id} missing {region}"


# =============================================
#      HIV BNAB LIGHT CHAIN TESTS
# =============================================


@pytest.mark.slow
@pytest.mark.integration
class TestHIVBnAbLightChains:
    """Integration tests using HIV broadly neutralizing antibody light chain sequences."""

    def test_returns_list(self, hiv_bnab_lc_results):
        """Pipeline returns a list for multi-sequence input."""
        assert isinstance(hiv_bnab_lc_results, list)

    def test_result_count(self, hiv_bnab_lc_results):
        """All sequences are processed."""
        assert len(hiv_bnab_lc_results) > 0

    def test_all_are_light_chain(self, hiv_bnab_lc_results):
        """All sequences are assigned to IGK or IGL locus."""
        for seq in hiv_bnab_lc_results:
            assert seq["locus"] in {"IGK", "IGL"}, (
                f"{seq.id} has locus {seq['locus']}, expected IGK or IGL"
            )

    def test_no_d_gene(self, hiv_bnab_lc_results):
        """Light chains should not have D gene assignment."""
        for seq in hiv_bnab_lc_results:
            assert seq["d_gene"] is None, f"{seq.id} has d_gene={seq['d_gene']}"

    def test_v_genes_are_igkv_or_iglv(self, hiv_bnab_lc_results):
        """All V genes start with IGKV or IGLV."""
        for seq in hiv_bnab_lc_results:
            v = seq["v_gene"]
            assert v.startswith("IGKV") or v.startswith("IGLV"), (
                f"{seq.id} v_gene={v} not IGKV/IGLV"
            )

    def test_all_have_cdr3(self, hiv_bnab_lc_results):
        """All annotated sequences have CDR3."""
        for seq in hiv_bnab_lc_results:
            assert seq["cdr3"] is not None, f"{seq.id} missing cdr3"


# =============================================
#      DATAFRAME OUTPUT INTEGRATION
# =============================================


@pytest.mark.slow
@pytest.mark.integration
class TestDataFrameOutput:
    """Integration tests for DataFrame output from multi-sequence input."""

    def test_dataframe_from_hc_file(self):
        """Test DataFrame output from heavy chain file."""
        hc_path = os.path.join(TEST_DATA_DIR, "test_hiv_bnab_hcs.fasta")
        result = run(hc_path, as_dataframe=True)
        assert isinstance(result, pl.DataFrame)
        assert result.height > 0
        assert "sequence_id" in result.columns
        assert "v_gene" in result.columns

    def test_dataframe_from_lc_file(self):
        """Test DataFrame output from light chain file."""
        lc_path = os.path.join(TEST_DATA_DIR, "test_hiv_bnab_lcs.fasta")
        result = run(lc_path, as_dataframe=True)
        assert isinstance(result, pl.DataFrame)
        assert result.height > 0


# =============================================
#      FILE OUTPUT INTEGRATION
# =============================================


@pytest.mark.slow
@pytest.mark.integration
class TestFileOutput:
    """Integration tests for file-based output."""

    def test_airr_output_from_hc_file(self, tmp_path):
        """Test AIRR TSV output from heavy chain file."""
        hc_path = os.path.join(TEST_DATA_DIR, "test_hiv_bnab_hcs.fasta")
        project_path = str(tmp_path / "airr_hc")
        run(hc_path, project_path=project_path, output_format="airr")

        airr_dir = os.path.join(project_path, "airr")
        assert os.path.exists(airr_dir)
        tsv_files = [f for f in os.listdir(airr_dir) if f.endswith(".tsv")]
        assert len(tsv_files) >= 1

    def test_parquet_output_from_hc_file(self, tmp_path):
        """Test Parquet output from heavy chain file."""
        hc_path = os.path.join(TEST_DATA_DIR, "test_hiv_bnab_hcs.fasta")
        project_path = str(tmp_path / "parquet_hc")
        run(hc_path, project_path=project_path, output_format="parquet")

        parquet_dir = os.path.join(project_path, "parquet")
        assert os.path.exists(parquet_dir)
        pq_files = [f for f in os.listdir(parquet_dir) if f.endswith(".parquet")]
        assert len(pq_files) >= 1
