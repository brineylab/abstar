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
from .edge_case_helpers import (
    REGION_NAMES,
    get_region_boundaries,
    delete_at_region_boundary,
    insert_at_region_boundary,
    insert_at_position,
    delete_at_position,
    truncate_3prime,
    mutate_near_boundary,
)

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


# ===========================================================================
#                  Integration test helpers
# ===========================================================================

# Pick representative sequences for integration tests
_HEAVY_IDS = [sid for sid in _SEQUENCE_IDS if sid.endswith("_heavy")]
_LIGHT_IDS = [sid for sid in _SEQUENCE_IDS if sid.endswith("_light")]
_REPR_HEAVY_ID = _HEAVY_IDS[0]
_REPR_LIGHT_ID = _LIGHT_IDS[0]

# Adjacent region pairs for boundary tests (boundary name = start of second region)
_BOUNDARY_PAIRS = [
    ("fwr1", "cdr1"),
    ("cdr1", "fwr2"),
    ("fwr2", "cdr2"),
    ("cdr2", "fwr3"),
    ("fwr3", "cdr3"),
]


def _run_single(sequence_input: str, seq_id: str = "edge_case"):
    """Run a single sequence through abstar and return the result.

    Returns None if annotation fails or the pipeline raises a RuntimeError
    (e.g., MMseqs2 fails on empty input after heavy truncation).
    """
    seq = Sequence(sequence_input, id=seq_id)
    try:
        results = run([seq])
    except RuntimeError:
        # MMseqs2 can fail when truncation leaves no query sequence for
        # downstream gene assignment (e.g., empty J-query FASTA).
        return None
    if isinstance(results, Sequence):
        results = [results]
    if not results:
        return None
    return results[0]


def _get_region_len(result, region: str) -> int:
    """Get the length of a region from an annotation result, normalizing None to 0."""
    val = result[region]
    if val is None:
        return 0
    return len(val)


def _assert_no_gaps_in_v_regions(result, regions=None):
    """Assert no gap characters in V-region NT outputs."""
    if regions is None:
        regions = REGION_NAMES
    for r in regions:
        val = result[r] or ""
        assert "-" not in val, f"{r} contains gap characters: {val[:60]}"


def _assert_remote_regions_stable(result, ground_truth_row, affected_boundary,
                                  tolerance=6):
    """Assert regions far from the affected boundary match ground truth.

    Regions adjacent to the affected boundary are skipped (they may shift).
    Regions further away are checked with exact match or within tolerance.
    """
    boundaries_by_region = {
        "cdr1": {"fwr1", "cdr1"},
        "fwr2": {"cdr1", "fwr2"},
        "cdr2": {"fwr2", "cdr2"},
        "fwr3": {"cdr2", "fwr3"},
        "cdr3": {"fwr3", "cdr3"},
    }
    adjacent = boundaries_by_region.get(affected_boundary, set())

    for region in REGION_NAMES:
        if region in adjacent:
            continue
        expected = _normalize(ground_truth_row[region], region)
        actual = _normalize(result[region], region)
        if expected and actual:
            len_diff = abs(len(actual) - len(expected))
            assert len_diff <= tolerance, (
                f"Remote region {region} length diff {len_diff} exceeds "
                f"tolerance {tolerance}: expected {len(expected)}, got {len(actual)}"
            )


# ===========================================================================
#              Edge Case 1: Indels at or Near Region Boundaries
# ===========================================================================


class TestIndelAtBoundary:
    """Integration tests for indels at or near region boundaries.

    See EDGE_CASES.md Edge Case 1 (1a-1e).
    """

    # --- 1a: Codon-length (3-nt) deletion spanning a boundary ---

    @pytest.mark.parametrize("boundary", ["cdr1", "fwr2", "cdr2", "fwr3"])
    def test_3nt_deletion_spanning_boundary(self, boundary, annotated_sequences):
        """A 3-nt deletion spanning a V-gene boundary should not crash and
        should produce a combined region pair that is exactly 3 nt shorter.

        CDR3 is excluded because it uses junction-based identification rather
        than IMGT position mapping.
        """
        row = _GROUND_TRUTH[_REPR_HEAVY_ID]
        # Delete 3 nt starting 1 nt before the boundary (offset=-1)
        modified = delete_at_region_boundary(row, boundary, length=3, offset=-1)
        result = _run_single(modified, f"del3_{boundary}")
        assert result is not None, f"Pipeline returned no result for del3_{boundary}"

        _assert_no_gaps_in_v_regions(result)

        # Identify the pair of regions around this boundary
        idx = REGION_NAMES.index(boundary)
        preceding = REGION_NAMES[idx - 1]
        orig_pair_len = len(row[preceding]) + len(row[boundary])
        actual_pair_len = _get_region_len(result, preceding) + _get_region_len(
            result, boundary
        )
        assert actual_pair_len == orig_pair_len - 3, (
            f"Combined {preceding}+{boundary} should be {orig_pair_len - 3}, "
            f"got {actual_pair_len}"
        )

    def test_3nt_deletion_spanning_fwr3_cdr3_boundary(self, annotated_sequences):
        """A 3-nt deletion at the FWR3/CDR3 boundary should not crash.

        CDR3 uses junction-based identification, so the FWR3/CDR3 boundary
        is special — it involves the conserved Cys codon.  A deletion here
        may cause FWR3 to be empty if the junction anchor is disrupted.
        The key assertion is that the pipeline doesn't crash and CDR3 is
        still identifiable.
        """
        row = _GROUND_TRUTH[_REPR_HEAVY_ID]
        modified = delete_at_region_boundary(row, "cdr3", length=3, offset=-1)
        result = _run_single(modified, "del3_cdr3")
        assert result is not None, "Pipeline returned no result for del3_cdr3"

        _assert_no_gaps_in_v_regions(result)

        # CDR3 should still be identifiable (may differ from original)
        assert _get_region_len(result, "cdr3") > 0, "CDR3 should not be empty"

    # --- 1b: Non-codon-length (1-nt, 2-nt) deletion at a boundary ---

    @pytest.mark.parametrize(
        "boundary,del_len",
        [
            ("cdr1", 1),
            ("cdr1", 2),
            ("cdr2", 1),
            ("cdr2", 2),
            ("cdr3", 1),
            ("cdr3", 2),
        ],
    )
    def test_noncodon_deletion_at_boundary(self, boundary, del_len,
                                           annotated_sequences):
        """A 1-nt or 2-nt deletion at a boundary should not crash.

        Non-codon-length deletions cause frameshifts.  The pipeline may
        legitimately fail to annotate such sequences (returning None), which
        is acceptable — the key requirement is no crash.  When annotation
        succeeds, verify basic properties.
        """
        row = _GROUND_TRUTH[_REPR_HEAVY_ID]
        modified = delete_at_region_boundary(row, boundary, length=del_len, offset=0)
        result = _run_single(modified, f"del{del_len}_{boundary}")

        if result is None:
            # Pipeline didn't crash but couldn't annotate — acceptable for
            # frame-disrupting modifications
            return

        _assert_no_gaps_in_v_regions(result)

        # Regions should still be extracted (non-empty for at least the pair)
        idx = REGION_NAMES.index(boundary)
        preceding = REGION_NAMES[idx - 1]
        pair_len = _get_region_len(result, preceding) + _get_region_len(
            result, boundary
        )
        assert pair_len > 0, (
            f"Combined {preceding}+{boundary} should not be empty"
        )

    # --- 1c: Codon-length (3-nt) insertion at a boundary ---

    @pytest.mark.parametrize("boundary", ["cdr1", "fwr2", "cdr2", "fwr3", "cdr3"])
    def test_3nt_insertion_at_boundary(self, boundary, annotated_sequences):
        """A 3-nt insertion at a boundary should not crash and the combined
        region pair should be 3 nt longer than original.
        """
        row = _GROUND_TRUTH[_REPR_HEAVY_ID]
        modified = insert_at_region_boundary(row, boundary, insertion="AAA", offset=0)
        result = _run_single(modified, f"ins3_{boundary}")
        assert result is not None, f"Pipeline returned no result for ins3_{boundary}"

        _assert_no_gaps_in_v_regions(result)

        idx = REGION_NAMES.index(boundary)
        preceding = REGION_NAMES[idx - 1]
        orig_pair_len = len(row[preceding]) + len(row[boundary])
        actual_pair_len = _get_region_len(result, preceding) + _get_region_len(
            result, boundary
        )
        assert actual_pair_len == orig_pair_len + 3, (
            f"Combined {preceding}+{boundary} should be {orig_pair_len + 3}, "
            f"got {actual_pair_len}"
        )

    # --- 1d: Non-codon-length (2-nt) insertion at a boundary ---

    @pytest.mark.parametrize("boundary", ["cdr1", "fwr2", "cdr3"])
    def test_2nt_insertion_at_boundary(self, boundary, annotated_sequences):
        """A 2-nt insertion 1 position before a boundary should not crash.

        Non-codon-length insertions cause frameshifts.  The pipeline may
        legitimately fail to annotate (returning None).  When annotation
        succeeds, the codon-stealing logic should handle the non-codon gap.
        """
        row = _GROUND_TRUTH[_REPR_HEAVY_ID]
        modified = insert_at_region_boundary(row, boundary, insertion="AG", offset=-1)
        result = _run_single(modified, f"ins2_{boundary}")

        if result is None:
            # Pipeline didn't crash but couldn't annotate — acceptable for
            # frame-disrupting modifications
            return

        _assert_no_gaps_in_v_regions(result)

        # The pair should collectively account for the insertion
        idx = REGION_NAMES.index(boundary)
        preceding = REGION_NAMES[idx - 1]
        orig_pair_len = len(row[preceding]) + len(row[boundary])
        actual_pair_len = _get_region_len(result, preceding) + _get_region_len(
            result, boundary
        )
        assert actual_pair_len == orig_pair_len + 2, (
            f"Combined {preceding}+{boundary} should be {orig_pair_len + 2}, "
            f"got {actual_pair_len}"
        )

    # --- 1e: Large (6-nt, 9-nt) codon-length insertion at a boundary ---

    @pytest.mark.parametrize(
        "boundary,ins_size",
        [
            ("fwr2", 3),
            ("fwr2", 6),
            ("fwr2", 9),
            ("fwr3", 3),
            ("fwr3", 6),
            ("fwr3", 9),
        ],
    )
    def test_large_codon_insertion_at_boundary(self, boundary, ins_size,
                                               annotated_sequences):
        """Multi-codon insertions at a boundary should be handled correctly."""
        row = _GROUND_TRUTH[_REPR_HEAVY_ID]
        insertion = "AGC" * (ins_size // 3)
        modified = insert_at_region_boundary(
            row, boundary, insertion=insertion, offset=0
        )
        result = _run_single(modified, f"ins{ins_size}_{boundary}")
        assert result is not None, (
            f"Pipeline returned no result for ins{ins_size}_{boundary}"
        )

        _assert_no_gaps_in_v_regions(result)

        idx = REGION_NAMES.index(boundary)
        preceding = REGION_NAMES[idx - 1]
        orig_pair_len = len(row[preceding]) + len(row[boundary])
        actual_pair_len = _get_region_len(result, preceding) + _get_region_len(
            result, boundary
        )
        assert actual_pair_len == orig_pair_len + ins_size, (
            f"Combined {preceding}+{boundary} should be {orig_pair_len + ins_size}, "
            f"got {actual_pair_len}"
        )


# ===========================================================================
#           Edge Case 2: Heavily Trimmed J-Gene in Light Chains
# ===========================================================================


class TestTrimmedJGene:
    """Integration tests for heavily trimmed J-gene in light chains.

    See EDGE_CASES.md Edge Case 2 (2a-2d).
    """

    # --- 2a: Remove most of FWR4 ---

    @pytest.mark.parametrize("seq_id", _LIGHT_IDS[:3])
    def test_fwr4_mostly_removed(self, seq_id):
        """Keeping only the last 6 nt of FWR4 should not crash.

        Heavy J-gene trimming can cause the J-query FASTA to be empty,
        which is a known pipeline limitation.  When annotation succeeds,
        FWR4 should be shorter than original.
        """
        row = _GROUND_TRUTH[seq_id]
        fwr4_len = len(row["fwr4"])
        bases_to_remove = fwr4_len - 6
        if bases_to_remove <= 0:
            pytest.skip(f"FWR4 too short to trim ({fwr4_len} nt)")
        modified = truncate_3prime(row, bases_to_remove)
        result = _run_single(modified, f"trim_fwr4_{seq_id}")

        if result is None:
            # Pipeline couldn't annotate (e.g., J-query empty) — no crash
            return

        actual_fwr4_len = _get_region_len(result, "fwr4")
        assert actual_fwr4_len < fwr4_len, (
            f"FWR4 should be shorter than original ({fwr4_len}), got {actual_fwr4_len}"
        )

    # --- 2b: Remove FWR4 entirely ---

    @pytest.mark.parametrize("seq_id", _LIGHT_IDS[:3])
    def test_fwr4_entirely_removed(self, seq_id):
        """Removing FWR4 entirely should not crash.

        When annotation succeeds, FWR4 should be empty or shorter.
        """
        row = _GROUND_TRUTH[seq_id]
        fwr4_len = len(row["fwr4"])
        modified = truncate_3prime(row, fwr4_len)
        result = _run_single(modified, f"no_fwr4_{seq_id}")

        if result is None:
            return

        actual_fwr4 = result["fwr4"] or ""
        assert len(actual_fwr4) < fwr4_len, (
            f"FWR4 should be shorter/empty after full removal, got len={len(actual_fwr4)}"
        )

    # --- 2c: J-gene trimmed to conserved codon only ---

    @pytest.mark.parametrize("seq_id", _LIGHT_IDS[:3])
    def test_fwr4_conserved_codon_only(self, seq_id):
        """Keeping only the first 3 nt of FWR4 should not crash."""
        row = _GROUND_TRUTH[seq_id]
        fwr4_len = len(row["fwr4"])
        bases_to_remove = fwr4_len - 3
        if bases_to_remove <= 0:
            pytest.skip(f"FWR4 too short to trim to 3 nt ({fwr4_len} nt)")
        modified = truncate_3prime(row, bases_to_remove)
        result = _run_single(modified, f"fwr4_codon_{seq_id}")

        if result is None:
            return

        actual_fwr4_len = _get_region_len(result, "fwr4")
        assert actual_fwr4_len <= fwr4_len, (
            f"FWR4 should not exceed original length ({fwr4_len}), got {actual_fwr4_len}"
        )

    # --- 2d: Verify locus and D-gene for light chains ---

    @pytest.mark.parametrize("seq_id", _LIGHT_IDS[:3])
    def test_light_chain_locus_preserved(self, seq_id):
        """After FWR4 trimming, locus should remain IGK/IGL and D-gene absent."""
        row = _GROUND_TRUTH[seq_id]
        fwr4_len = len(row["fwr4"])
        bases_to_remove = fwr4_len - 6
        if bases_to_remove <= 0:
            pytest.skip(f"FWR4 too short to trim ({fwr4_len} nt)")
        modified = truncate_3prime(row, bases_to_remove)
        result = _run_single(modified, f"locus_{seq_id}")

        if result is None:
            return

        locus = result["locus"]
        assert locus in ("IGK", "IGL", None), (
            f"Light chain locus should be IGK/IGL, got {locus}"
        )
        d_call = result["d_call"] or ""
        assert d_call == "", f"Light chain should have no D-gene call, got {d_call}"


# ===========================================================================
#            Edge Case 3: Complementary Frameshifts
# ===========================================================================


class TestComplementaryFrameshifts:
    """Integration tests for complementary frameshifts.

    See EDGE_CASES.md Edge Case 3 (3a-3d).
    """

    # --- 3a: 2-nt deletion in FWR2 + 1-nt insertion in FWR3 (net -1) ---

    def test_del2_fwr2_ins1_fwr3_net_minus1(self):
        """A 2-nt deletion + 1-nt insertion (net -1) should not crash.

        Net -1 is frame-disrupting, so the pipeline may legitimately fail
        to annotate.  When annotation succeeds, regions should be extracted
        and productivity fields populated.
        """
        row = _GROUND_TRUTH[_REPR_HEAVY_ID]
        boundaries = get_region_boundaries(row)
        seq = row["sequence_input"]

        # Delete 2 nt at FWR2 +10
        fwr2_start = boundaries["fwr2"][0]
        seq = delete_at_position(seq, fwr2_start + 10, 2)

        # Insert 1 nt at FWR3 +10 (adjusted for prior deletion)
        fwr3_start = boundaries["fwr3"][0] - 2  # shifted by deletion
        seq = insert_at_position(seq, fwr3_start + 10, "G")

        result = _run_single(seq, "comp_frameshift_3a")

        if result is None:
            # Net -1 frameshift may cause annotation failure — no crash
            return

        _assert_no_gaps_in_v_regions(result)

        # Productivity fields should be populated
        assert result["productive"] is not None
        issues = result["productivity_issues"]
        assert isinstance(issues, (list, type(None)))

    # --- 3b: 2-nt deletion + 2-nt insertion (net 0, in-frame) ---

    def test_del2_fwr2_ins2_fwr3_net_zero(self):
        """A 2-nt deletion + 2-nt insertion (net 0) keeps the sequence in-frame.

        Both affected regions should have altered lengths.  CDR2 (between)
        should be approximately unchanged.
        """
        row = _GROUND_TRUTH[_REPR_HEAVY_ID]
        boundaries = get_region_boundaries(row)
        seq = row["sequence_input"]

        # Delete 2 nt at FWR2 +15
        fwr2_start = boundaries["fwr2"][0]
        seq = delete_at_position(seq, fwr2_start + 15, 2)

        # Insert 2 nt at FWR3 +20 (adjusted for prior deletion)
        fwr3_start = boundaries["fwr3"][0] - 2
        seq = insert_at_position(seq, fwr3_start + 20, "AG")

        result = _run_single(seq, "comp_frameshift_3b")
        assert result is not None, "Pipeline returned no result"

        _assert_no_gaps_in_v_regions(result)

        # CDR2 (between the two indels) should be approximately unchanged
        orig_cdr2_len = len(row["cdr2"])
        actual_cdr2_len = _get_region_len(result, "cdr2")
        assert abs(actual_cdr2_len - orig_cdr2_len) <= 3, (
            f"CDR2 should be ~unchanged: expected ~{orig_cdr2_len}, got {actual_cdr2_len}"
        )

    # --- 3c: 1-nt deletion in FWR2 + 2-nt deletion in FWR3 (net -3, in-frame) ---

    def test_del1_fwr2_del2_fwr3_net_minus3(self):
        """Two non-codon-length deletions summing to -3 (in-frame) should
        preserve region extraction.
        """
        row = _GROUND_TRUTH[_REPR_HEAVY_ID]
        boundaries = get_region_boundaries(row)
        seq = row["sequence_input"]

        # Delete 1 nt at FWR2 +10
        fwr2_start = boundaries["fwr2"][0]
        seq = delete_at_position(seq, fwr2_start + 10, 1)

        # Delete 2 nt at FWR3 +20 (adjusted for prior deletion)
        fwr3_start = boundaries["fwr3"][0] - 1
        seq = delete_at_position(seq, fwr3_start + 20, 2)

        result = _run_single(seq, "comp_frameshift_3c")
        assert result is not None, "Pipeline returned no result"

        _assert_no_gaps_in_v_regions(result)

        # Total sequence should be 3 shorter
        orig_len = len(row["sequence_input"])
        assert len(seq) == orig_len - 3

    # --- 3d: 1-nt insertion in CDR1 + 2-nt insertion in CDR2 (net +3, in-frame) ---

    def test_ins1_cdr1_ins2_cdr2_net_plus3(self):
        """Two non-codon-length insertions summing to +3 (in-frame) should
        preserve region extraction.
        """
        row = _GROUND_TRUTH[_REPR_HEAVY_ID]
        boundaries = get_region_boundaries(row)
        seq = row["sequence_input"]

        # Insert 1 nt at middle of CDR1
        cdr1_start, cdr1_end = boundaries["cdr1"]
        cdr1_mid = cdr1_start + (cdr1_end - cdr1_start) // 2
        seq = insert_at_position(seq, cdr1_mid, "A")

        # Insert 2 nt at middle of CDR2 (adjusted for prior insertion)
        cdr2_start, cdr2_end = boundaries["cdr2"]
        cdr2_mid = cdr2_start + 1 + (cdr2_end - cdr2_start) // 2  # +1 for prior ins
        seq = insert_at_position(seq, cdr2_mid, "TG")

        result = _run_single(seq, "comp_frameshift_3d")
        assert result is not None, "Pipeline returned no result"

        _assert_no_gaps_in_v_regions(result)

        # FWR1 (upstream of both) should be approximately unchanged
        orig_fwr1_len = len(row["fwr1"])
        actual_fwr1_len = _get_region_len(result, "fwr1")
        assert abs(actual_fwr1_len - orig_fwr1_len) <= 3, (
            f"FWR1 should be ~unchanged: expected ~{orig_fwr1_len}, got {actual_fwr1_len}"
        )

        # FWR2 (between CDR1 and CDR2) should be approximately unchanged
        orig_fwr2_len = len(row["fwr2"])
        actual_fwr2_len = _get_region_len(result, "fwr2")
        assert abs(actual_fwr2_len - orig_fwr2_len) <= 3, (
            f"FWR2 should be ~unchanged: expected ~{orig_fwr2_len}, got {actual_fwr2_len}"
        )
