# Copyright (c) 2025 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

"""
Tests for the Antibody dataclass.
"""

from abstar.annotation.antibody import Antibody

# =============================================
#          CONSTRUCTION TESTS
# =============================================


def test_antibody_default_construction():
    """Test Antibody can be constructed with all defaults."""
    ab = Antibody()
    assert ab.sequence_id is None
    assert ab.v_gene is None
    assert ab.productive is True
    assert ab.complete_vdj is False
    assert ab.rev_comp is False
    assert ab.stop_codon is False
    assert ab.v_frameshift is False
    assert ab.productivity_issues == []


def test_antibody_with_fields():
    """Test Antibody construction with specific fields."""
    ab = Antibody(
        sequence_id="test_seq",
        v_gene="IGHV1-2*02",
        j_gene="IGHJ4*02",
        d_gene="IGHD3-10*01",
        locus="IGH",
        productive=True,
        cdr3_aa="CARDYW",
        cdr3_length=4,
    )
    assert ab.sequence_id == "test_seq"
    assert ab.v_gene == "IGHV1-2*02"
    assert ab.j_gene == "IGHJ4*02"
    assert ab.d_gene == "IGHD3-10*01"
    assert ab.locus == "IGH"
    assert ab.productive is True
    assert ab.cdr3_aa == "CARDYW"
    assert ab.cdr3_length == 4


def test_antibody_light_chain_no_d_gene():
    """Test Antibody for light chain with no D gene."""
    ab = Antibody(
        sequence_id="lc_test",
        v_gene="IGKV1-39*01",
        j_gene="IGKJ1*01",
        d_gene=None,
        locus="IGK",
    )
    assert ab.d_gene is None
    assert ab.locus == "IGK"


def test_antibody_productivity_issues_default_factory():
    """Test that productivity_issues uses a fresh list per instance."""
    ab1 = Antibody()
    ab2 = Antibody()
    ab1.productivity_issues.append("stop codon")
    assert ab1.productivity_issues == ["stop codon"]
    assert ab2.productivity_issues == []


# =============================================
#          AIRR FIELDS TESTS
# =============================================


def test_airr_fields_populated():
    """Test that airr_fields is populated by __post_init__."""
    ab = Antibody()
    assert hasattr(ab, "airr_fields")
    assert isinstance(ab.airr_fields, list)
    assert len(ab.airr_fields) > 0


def test_airr_fields_contains_core_fields():
    """Test that airr_fields includes core annotation fields."""
    ab = Antibody()
    core_fields = [
        "sequence_id",
        "v_gene",
        "d_gene",
        "j_gene",
        "c_gene",
        "productive",
        "locus",
        "cdr3",
        "cdr3_aa",
        "junction",
        "junction_aa",
        "v_identity",
        "sequence",
        "germline",
    ]
    for field in core_fields:
        assert field in ab.airr_fields, f"{field} should be in airr_fields"


def test_airr_fields_excludes_non_airr_attributes():
    """Test that airr_fields does not include post-init attributes like airr_fields itself."""
    ab = Antibody()
    # airr_fields itself and LoggingMixin attributes should not be in the list
    assert "airr_fields" not in ab.airr_fields


# =============================================
#           TO_DICT TESTS
# =============================================


def test_to_dict_basic():
    """Test basic to_dict conversion."""
    ab = Antibody(sequence_id="test", v_gene="IGHV1-2*02")
    d = ab.to_dict()
    assert isinstance(d, dict)
    assert d["sequence_id"] == "test"
    assert d["v_gene"] == "IGHV1-2*02"
    assert d["j_gene"] is None


def test_to_dict_contains_all_airr_fields():
    """Test that to_dict includes all AIRR fields."""
    ab = Antibody()
    d = ab.to_dict()
    for field in ab.airr_fields:
        assert field in d, f"{field} should be in to_dict output"


def test_to_dict_exclude_single_string():
    """Test to_dict with single string exclude."""
    ab = Antibody(sequence_id="test", v_gene="IGHV1-2*02")
    d = ab.to_dict(exclude="v_gene")
    assert "v_gene" not in d
    assert "sequence_id" in d


def test_to_dict_exclude_iterable():
    """Test to_dict with iterable of excludes."""
    ab = Antibody(sequence_id="test", v_gene="IGHV1-2*02", j_gene="IGHJ4*02")
    d = ab.to_dict(exclude=["v_gene", "j_gene"])
    assert "v_gene" not in d
    assert "j_gene" not in d
    assert "sequence_id" in d


def test_to_dict_include_single_string():
    """Test to_dict with single string include (adds extra fields)."""
    ab = Antibody()
    ab.custom_field = "custom_value"
    d = ab.to_dict(include="custom_field")
    assert "custom_field" in d
    assert d["custom_field"] == "custom_value"


def test_to_dict_include_iterable():
    """Test to_dict with iterable of includes."""
    ab = Antibody()
    ab.extra1 = "val1"
    ab.extra2 = "val2"
    d = ab.to_dict(include=["extra1", "extra2"])
    assert d["extra1"] == "val1"
    assert d["extra2"] == "val2"


def test_to_dict_include_missing_field_returns_none():
    """Test to_dict with include for a field that doesn't exist returns None."""
    ab = Antibody()
    d = ab.to_dict(include="nonexistent")
    assert d["nonexistent"] is None


def test_to_dict_exclude_and_include():
    """Test to_dict with both exclude and include."""
    ab = Antibody(sequence_id="test", v_gene="IGHV1-2*02")
    ab.extra = "extra_value"
    d = ab.to_dict(exclude="v_gene", include="extra")
    assert "v_gene" not in d
    assert "extra" in d
    assert d["extra"] == "extra_value"


# =============================================
#        LOGGING MIXIN INTEGRATION TESTS
# =============================================


def test_antibody_has_logging_methods():
    """Test that Antibody inherits LoggingMixin methods."""
    ab = Antibody()
    assert hasattr(ab, "log")
    assert hasattr(ab, "exception")
    assert hasattr(ab, "format_log")
    assert hasattr(ab, "logs")
    assert hasattr(ab, "exceptions")


def test_antibody_logging():
    """Test that logging works through the Antibody object."""
    ab = Antibody(sequence_id="test")
    ab.log("V gene assignment", "IGHV1-2*02")
    assert len(ab.logs) == 1
    assert "IGHV1-2*02" in ab.logs[0]


def test_antibody_exception_logging():
    """Test that exception logging works through the Antibody object."""
    ab = Antibody(sequence_id="test")
    ab.exception("Failed to assign V gene")
    assert len(ab.exceptions) == 1
    assert "Failed to assign V gene" in ab.exceptions[0]


def test_antibody_format_log():
    """Test format_log output."""
    ab = Antibody(sequence_id="test")
    ab.log("step 1")
    ab.log("step 2")
    log_output = ab.format_log()
    assert "step 1" in log_output
    assert "step 2" in log_output
