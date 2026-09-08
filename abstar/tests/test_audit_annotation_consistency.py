"""Independent oriented-query comparisons must detect equal-length errors too."""
import importlib.util
import json
from pathlib import Path

import polars as pl
import pytest


@pytest.fixture
def audit():
    spec = importlib.util.spec_from_file_location(
        "audit_annotation_consistency",
        Path(__file__).parents[2] / "scripts/audit_annotation_consistency.py",
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def row():
    regions = ("fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4")
    record = dict(
        sequence_id="10E8", annotation_status="annotated", productive=True,
        sequence_oriented="GGGATGAAATGTTTTCCCGGGTGGCCC", frame=1,
        v_sequence_start=3, j_sequence_end=24,
        v_sequence="ATGAAATGT", np1="TTT", d_sequence="CCC", np2="",
        j_sequence="GGGTGG", sequence_alignment="ATGAAA-TGTTTTCCCGGGTGG",
        cdr_mask="000111000222000333000", cdr_mask_aa="0102030",
        gene_segment_mask="VVVVVVVVVNNNDDD" + "J" * 6,
        gene_segment_mask_aa="VVVNDJJ", nongermline_mask="0" * 21,
        nongermline_mask_aa="0" * 7,
    )
    for i, (region, nt, aa) in enumerate(zip(
        regions, ["ATG", "AAA", "TGT", "TTT", "CCC", "GGG", "TGG"], "MKCFPGW",
    )):
        record.update({region: nt, region + "_aa": aa,
                       region + "_start": 3 + 3 * i, region + "_end": 6 + 3 * i})
    return record


def test_oriented_slice_excludes_flanks_and_alignment_gaps(audit, row):
    assert audit.audit_frame(pl.DataFrame([row])).is_empty()


@pytest.mark.parametrize("field,value,check", [
    ("fwr3", "", "regions_nt_content"),
    ("fwr3_aa", "A", "regions_aa_content"),
    ("cdr_mask", "0" * 20, "cdr_mask_length"),
    ("gene_segment_mask_aa", "V" * 6, "gene_segment_mask_aa_length"),
    ("nongermline_mask", None, "nongermline_mask_null"),
    ("sequence_alignment", "A" * 21, "alignment_nt_content"),
    ("sequence_alignment", None, "alignment_nt_content"),
    ("v_sequence", "A" * 9, "segments_nt_content"),
    ("cdr1_start", 5, "cdr1_coordinates"),
    ("frame", 0, "vdj_reference_invalid"),
    ("j_sequence_end", 999, "vdj_reference_invalid"),
])
def test_detects_inconsistency(audit, row, field, value, check):
    row[field] = value
    found = audit.audit_frame(pl.DataFrame([row]))
    assert found.height == 1
    assert check in found.row(0, named=True)["checks"]


@pytest.mark.parametrize("frame,prefix", [(2, "A"), (3, "AC")])
def test_reference_translation_uses_frame_and_ignores_partial_terminal_codon(audit, row, frame, prefix):
    row["sequence_oriented"] = "GGG" + prefix + "ATGAAATGTTTTCCCGGGTGGA" + "CCC"
    row["frame"] = frame
    row["j_sequence_end"] = 25 + len(prefix)
    row["fwr3_aa"] = "A"
    found = audit.audit_frame(pl.DataFrame([row])).row(0, named=True)
    assert found["vdj_aa_length"] == 7
    assert "regions_aa_content" in found["checks"]


def test_missing_columns_are_explicit(audit, row):
    del row["sequence_oriented"]
    with pytest.raises(ValueError, match="sequence_oriented"):
        audit.audit_frame(pl.DataFrame([row]))


def test_report_preserves_duplicate_ids_and_reports_unassigned_separately(audit, row, tmp_path):
    root = tmp_path / "inputs"
    root.mkdir()
    row["fwr3_aa"] = "A"
    unassigned = {k: None for k in row}
    unassigned.update(sequence_id="001", annotation_status="unassigned")
    pl.DataFrame([row, row, unassigned]).write_parquet(root / "sample.parquet")
    destination = tmp_path / "audit.json"
    assert audit.main([str(root), "--report", str(destination)]) == 1
    report = json.loads(destination.read_text())
    assert report["totals"] == {"rows": 3, "annotated": 2, "unassigned": 1, "inconsistent": 2}
    assert [r["sequence_id"] for r in report["candidates"]] == ["10E8", "10E8"]
    assert [r["parquet_row"] for r in report["candidates"]] == [0, 1]
    with pytest.raises(FileExistsError):
        audit.main([str(root), "--report", str(destination)])
    with pytest.raises(ValueError, match="outside"):
        audit.main([str(root), "--report", str(root / "audit.json")])
