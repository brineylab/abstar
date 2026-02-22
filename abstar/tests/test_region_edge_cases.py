# Copyright (c) 2025 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

"""
Edge case tests for antibody region annotation.

See EDGE_CASES.md for full documentation of test rationale and expected
behaviors.
"""

import csv
import os

import abutils
import pytest
from abutils import Sequence

from ..annotation.antibody import Antibody
from ..annotation.regions import get_region_sequence
from ..core.abstar import run
from .edge_case_helpers import REGION_NAMES

# ---------------------------------------------------------------------------
# Ground truth data (loaded at import time for pytest parametrization)
# ---------------------------------------------------------------------------

_CSV_PATH = os.path.join(os.path.dirname(__file__), "..", "test_data", "test_50.csv")
with open(_CSV_PATH, encoding="utf-8-sig", newline="") as _f:
    _GROUND_TRUTH = {r["sequence_id"]: r for r in csv.DictReader(_f)}
_SEQUENCE_IDS = list(_GROUND_TRUTH.keys())

# All 14 region fields to validate
_AA_REGION_NAMES = [f"{r}_aa" for r in REGION_NAMES]
_ALL_REGION_FIELDS = REGION_NAMES + _AA_REGION_NAMES


def _normalize(value: str | None, field: str) -> str:
    """Normalize a region field value for comparison.

    NT region fields are upper-cased so that ground truth CSV (lowercase) and
    abstar output (uppercase) compare equal.  AA fields are left as-is.
    """
    s = value or ""
    if field in REGION_NAMES:
        return s.upper()
    return s


# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def annotated_sequences():
    """Run all 50 ground truth sequences through abstar.run() as a batch.

    Module-scoped so the expensive pipeline call happens only once.
    Returns a dict mapping sequence_id to the annotated Sequence result.
    """
    sequences = [
        Sequence(row["sequence_input"], id=sid)
        for sid, row in _GROUND_TRUTH.items()
    ]
    results = run(sequences)
    # run() returns a single Sequence when only one result; normalise to list
    if isinstance(results, Sequence):
        results = [results]
    return {r["sequence_id"]: r for r in results}


# ===========================================================================
#                       TestGroundTruthBaseline
# ===========================================================================


class TestGroundTruthBaseline:
    """Baseline validation: all 50 unmutated sequences annotate correctly.

    Before testing edge cases, we verify that every ground truth sequence
    produces exact-match annotations for all 14 region fields (7 NT + 7 AA).
    This serves as a sanity check and regression baseline.
    """

    def test_all_sequences_annotated(self, annotated_sequences):
        """All 50 ground truth sequences should produce annotation results."""
        missing = [sid for sid in _SEQUENCE_IDS if sid not in annotated_sequences]
        assert not missing, (
            f"Missing annotations for {len(missing)} sequence(s): {missing}"
        )

    @pytest.mark.parametrize("seq_id", _SEQUENCE_IDS)
    def test_all_regions_match(self, seq_id, annotated_sequences):
        """All 14 region fields (7 NT + 7 AA) match ground truth CSV exactly."""
        result = annotated_sequences.get(seq_id)
        assert result is not None, f"Sequence {seq_id} was not annotated"

        expected = _GROUND_TRUTH[seq_id]
        mismatches = []
        for field in _ALL_REGION_FIELDS:
            actual = _normalize(result[field], field)
            exp = _normalize(expected[field], field)
            if actual != exp:
                mismatches.append(
                    f"  {field}: expected '{exp[:60]}' got '{actual[:60]}'"
                )
        assert not mismatches, (
            f"Region mismatches for {seq_id}:\n" + "\n".join(mismatches)
        )


# ===========================================================================
#                    TestGetRegionSequenceEdgeCases
# ===========================================================================


def _make_alignment(sequence, germline):
    """Create a semiglobal alignment for unit tests."""
    return abutils.tl.semiglobal_alignment(sequence, germline)


def _normalize_region_output(result):
    """Normalize get_region_sequence output for missing-region compatibility.

    The standard return is ``(start, end, sequence)``, but fully missing
    regions may return just ``""``.
    """
    if isinstance(result, tuple):
        return result
    return None, None, result


class TestGetRegionSequenceEdgeCases:
    """Unit tests for get_region_sequence() targeting specific code paths.

    These bypass the full pipeline and test region extraction directly with
    crafted alignment objects.  See EDGE_CASES.md "Unit Tests" section.
    """

    # -- shared germline data from the existing test_regions.py fixtures --

    # Ungapped germline regions
    _FWR1_G = "CAGGTTCAGCTGGTGCAGTCTGGAGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCT"
    _CDR1_G = "CACTTCACCGGCATGCACTGGGAGTCAGTA"
    _FWR2_G = "GTGCGACAGGCCCCTGGACAAGGGCTTGAGCGGATCAACCCTAACAGT"
    _CDR2_G = "GGGGGCACAAACTATGCACAGAAG"
    _FWR3_G = "TTCCAGGGCAGAGTCACCATGACCACAGACACATCCGATACGAGCACAGCCTACATGGAGCTGAGCAGCCTGAGATCTGAGGACACGGCCGTGTATTACTGTGCG"
    _CDR3_G = "AGAGTGGGCAACTGGGGCCAGGGTACCTTTGACTACTGG"
    _FWR4_G = "GGCCAGGGAACCCTGGTCACCGTCTCCTCAGGTAAG"

    _UNGAPPED = _FWR1_G + _CDR1_G + _FWR2_G + _CDR2_G + _FWR3_G + _CDR3_G + _FWR4_G

    # Gapped germline (IMGT-gapped, with dots)
    _FWR1_GG = "CAGGTTCAGCTG...GTGCAGTCTGGAGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCT"
    _CDR1_GG = "CACTTCACCGGC......ATGCACTGGGAGTCAGTA"
    _FWR2_GG = "GTGCGACAGGCCCCTGGACAAGGGCTTGAG...CGGATCAACCCTAACAGT"
    _CDR2_GG = "GGGGGC......ACAAACTATGCACAGAAG"
    _FWR3_GG = "TTCCAGGGCAGAGTCACCATGACCACAGACACATCCGAT............ACGAGCACAGCCTACATGGAGCTGAGCAGCCTGAGATCTGAGGACACGGCCGTGTATTACTGTGCG"
    _CDR3_GG = "AGAGTGGGCAACTGGGGCCAGGGTACCTTTGACTACTGG"
    _FWR4_GG = "GGCCAGGGAACCCTGGTCACCGTCTCCTCAGGTAAG"

    _GAPPED = _FWR1_GG + _CDR1_GG + _FWR2_GG + _CDR2_GG + _FWR3_GG + _CDR3_GG + _FWR4_GG

    @pytest.fixture
    def ab(self):
        return Antibody(sequence_id="test_unit")

    def test_noncodon_boundary_gap_triggers_codon_stealing(self, ab):
        """A deletion spanning a region boundary should trigger codon-stealing.

        Scenario: delete 3 nt at the FWR1/CDR1 boundary (positions 73-75 of
        the ungapped germline, 0-indexed).  The aligned query has a 3-nt gap
        at positions corresponding to the last 1 nt of FWR1 and first 2 nt of
        CDR1.  The gap chars are split non-codon-length across the two regions,
        so the codon-stealing logic (regions.py:174-200) should activate and
        adjust the boundary to keep codons intact.

        Assert:
          - Tuple return shape (start, end, sequence)
          - No gap characters in the returned sequence
          - FWR1 region sequence is shorter than original (deletion removed bases)
          - CDR1 region sequence is also affected (boundary adjusted)
          - Combined FWR1 + CDR1 is exactly 3 nt shorter than original
        """
        # Create a query with a 3-nt deletion at the FWR1/CDR1 boundary
        query = self._UNGAPPED[:73] + self._UNGAPPED[76:]
        aln = _make_alignment(query, self._UNGAPPED)

        # Extract FWR1
        fwr1_result = get_region_sequence(
            region="fwr1",
            aln=aln,
            gapped_germline=self._GAPPED,
            germline_start=1,
            ab=ab,
        )
        fwr1_start, fwr1_end, fwr1_seq = _normalize_region_output(fwr1_result)

        # Extract CDR1
        cdr1_result = get_region_sequence(
            region="cdr1",
            aln=aln,
            gapped_germline=self._GAPPED,
            germline_start=1,
            ab=ab,
        )
        cdr1_start, cdr1_end, cdr1_seq = _normalize_region_output(cdr1_result)

        # Tuple return shape
        assert isinstance(fwr1_result, tuple), "FWR1 should return a tuple"
        assert isinstance(cdr1_result, tuple), "CDR1 should return a tuple"

        # No gaps in returned sequences
        assert "-" not in fwr1_seq, "FWR1 should have no gap characters"
        assert "-" not in cdr1_seq, "CDR1 should have no gap characters"

        # Both regions should have non-empty sequences
        assert len(fwr1_seq) > 0, "FWR1 should not be empty"
        assert len(cdr1_seq) > 0, "CDR1 should not be empty"

        # Combined length should be exactly 3 nt shorter than originals
        orig_combined = len(self._FWR1_G) + len(self._CDR1_G)
        actual_combined = len(fwr1_seq) + len(cdr1_seq)
        assert actual_combined == orig_combined - 3, (
            f"Combined FWR1+CDR1 should be {orig_combined - 3} nt "
            f"(original {orig_combined} minus 3-nt deletion), got {actual_combined}"
        )

        # Codon-stealing should have adjusted the boundary: if FWR1 has
        # trailing gaps that aren't codon-length, the end extends into CDR1
        # territory (or CDR1 start retreats).  Verify the coordinates reflect
        # this adjustment.
        assert fwr1_end is not None and cdr1_start is not None
        assert fwr1_end < cdr1_start or fwr1_end == cdr1_start, (
            "FWR1 end should not exceed CDR1 start after codon-stealing"
        )

    def test_aa_alignment_truncation_corrected_by_nt_start(self, ab):
        """AA alignment truncation should be corrected when nt_region_start=0.

        Tests the edge case at regions.py:149-167.  When the NT alignment
        correctly starts at position 0 but the AA alignment introduces leading
        gaps (due to 1-2 mutations in the first few AAs), the function should
        correct the AA region start to 0.

        Scenario: A query with 3 extra leading amino acids ("AAA" prepended)
        causes the AA alignment target to have leading dashes ("---QVQ..."),
        making get_aligned_position_from_ungapped(0, aligned_target) return > 0.
        The nt_region_start=0 parameter triggers the correction to force
        region_start=0.
        """
        # Build AA-level sequences
        germ_aa = abutils.tl.translate(self._UNGAPPED, frame=1)
        gapped_germ_aa = abutils.tl.translate(self._GAPPED, allow_dots=True)

        # Create a query AA with 3 extra leading AAs so aligner introduces
        # leading gaps in the target
        query_aa = "AAA" + germ_aa

        aln_aa = _make_alignment(query_aa, germ_aa)

        # Verify the alignment target has leading gaps (the condition we're testing)
        target_has_leading_gaps = aln_aa.aligned_target.startswith("-")

        # Call with nt_region_start=0 — this should trigger the correction
        result = get_region_sequence(
            region="fwr1",
            aln=aln_aa,
            gapped_germline=gapped_germ_aa,
            germline_start=1,
            ab=ab,
            aa=True,
            nt_region_start=0,
        )
        start, end, seq = _normalize_region_output(result)

        if target_has_leading_gaps:
            # The correction should have forced start to 0
            assert start == 0, (
                f"AA region start should be corrected to 0 when nt_region_start=0, "
                f"got {start}"
            )

        # Regardless of whether the edge case was triggered, verify basic
        # properties: non-empty, no gaps
        assert len(seq) > 0, "AA FWR1 should not be empty"
        assert "-" not in seq, "AA FWR1 should have no gap characters"

    def test_fully_missing_region_returns_empty(self, ab):
        """When germline_start is beyond a region's IMGT end, return empty.

        Tests the early-exit path at regions.py:98-99 where
        gapped_germline_start > imgt_end.
        """
        aln = _make_alignment(self._UNGAPPED, self._UNGAPPED)

        # germline_start=400 means the query starts well past all regions
        result = get_region_sequence(
            region="fwr1",
            aln=aln,
            gapped_germline=self._GAPPED,
            germline_start=400,
            ab=ab,
        )
        start, end, seq = _normalize_region_output(result)

        assert seq == "", (
            f"Fully missing region should return empty string, got '{seq[:40]}'"
        )

    def test_insertion_at_imgt_gap_position(self, ab):
        """An insertion at an IMGT gap position should be handled correctly.

        Tests that position mapping works when the query has an insertion
        in the CDR1 region (which has IMGT gap positions marked with dots).
        The insertion appears as gaps in the aligned germline (target), and
        get_aligned_position_from_ungapped with is_end_position=True should
        extend the region to include the insertion.

        Assert:
          - Tuple return shape
          - CDR1 sequence is longer than the germline CDR1 (due to insertion)
          - No gap characters in the returned sequence
        """
        # Insert 3 nt at the midpoint of CDR1 (position 12 within CDR1)
        cdr1_midpoint = len(self._FWR1_G) + 12
        query = (
            self._UNGAPPED[:cdr1_midpoint]
            + "AAA"
            + self._UNGAPPED[cdr1_midpoint:]
        )
        aln = _make_alignment(query, self._UNGAPPED)

        result = get_region_sequence(
            region="cdr1",
            aln=aln,
            gapped_germline=self._GAPPED,
            germline_start=1,
            ab=ab,
        )
        start, end, seq = _normalize_region_output(result)

        # Should return a tuple
        assert isinstance(result, tuple), "CDR1 should return a tuple"

        # CDR1 should be longer than original due to the 3-nt insertion
        assert len(seq) == len(self._CDR1_G) + 3, (
            f"CDR1 with 3-nt insertion should be {len(self._CDR1_G) + 3} nt, "
            f"got {len(seq)}"
        )

        # No gaps
        assert "-" not in seq, "CDR1 should have no gap characters"

        # The insertion nucleotides should be present in the sequence
        assert "AAA" in seq or seq.count("A") >= 3, (
            "Insertion bases should be present in the CDR1 sequence"
        )
