# Copyright (c) 2025 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

"""
Sequence modification utility functions for edge case testing.

These helpers take ground truth rows from test_50.csv and produce modified
sequences with predictable changes for testing the annotation pipeline.

All mutation helpers are deterministic: they use either fixed substitution
maps or a caller-supplied RNG seed.  Base substitutions avoid identity
(``A->A``) and accidental stop codons unless the test explicitly targets
stop codons.
"""

import random

# Ordered list of region names matching the CSV columns
REGION_NAMES = ["fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4"]

# Deterministic base-substitution map.  First alternative is always used by
# ``mutate_position`` when no explicit replacement is given, making single-
# base mutations fully deterministic without an RNG.
_BASE_ALTERNATIVES = {
    "A": ["T", "C", "G"],
    "T": ["C", "G", "A"],
    "C": ["G", "A", "T"],
    "G": ["A", "T", "C"],
    "a": ["t", "c", "g"],
    "t": ["c", "g", "a"],
    "c": ["g", "a", "t"],
    "g": ["a", "t", "c"],
    "N": ["A", "T", "C"],
    "n": ["a", "t", "c"],
}


# ---------------------------------------------------------------------------
# Region boundary mapping
# ---------------------------------------------------------------------------


def get_region_boundaries(row: dict) -> dict:
    """
    Compute (start, end) positions for each region within the full input sequence.

    Returns dict mapping region name to ``(start, end)`` tuple where *start*
    is inclusive and *end* is exclusive (suitable for Python slicing).

    Example for a sequence with FWR1=75nt, CDR1=24nt, FWR2=51nt, ...::

        {"fwr1": (0, 75), "cdr1": (75, 99), "fwr2": (99, 150), ...}
    """
    boundaries = {}
    pos = 0
    for region in REGION_NAMES:
        length = len(row[region])
        boundaries[region] = (pos, pos + length)
        pos += length
    return boundaries


# ---------------------------------------------------------------------------
# Low-level sequence modification
# ---------------------------------------------------------------------------


def insert_at_position(sequence: str, position: int, insertion: str) -> str:
    """Insert *insertion* into *sequence* at *position*.

    Characters at *position* and beyond are shifted right.
    """
    return sequence[:position] + insertion + sequence[position:]


def delete_at_position(sequence: str, position: int, length: int) -> str:
    """Delete *length* characters from *sequence* starting at *position*."""
    return sequence[:position] + sequence[position + length :]


# ---------------------------------------------------------------------------
# Boundary-relative modifications
# ---------------------------------------------------------------------------


def insert_at_region_boundary(
    row: dict, boundary: str, insertion: str, offset: int = 0
) -> str:
    """Insert *insertion* near the start of the *boundary* region.

    Parameters
    ----------
    row : dict
        Ground truth row from test_50.csv.
    boundary : str
        Region name whose start position defines the boundary
        (e.g. ``"cdr1"`` for the FWR1/CDR1 boundary).
    insertion : str
        Nucleotides to insert.
    offset : int
        Offset from the boundary start.  ``0`` inserts exactly at the
        boundary; negative values insert into the preceding region.

    Returns
    -------
    str
        Modified full sequence.
    """
    boundaries = get_region_boundaries(row)
    pos = boundaries[boundary][0] + offset
    return insert_at_position(row["sequence_input"], pos, insertion)


def delete_at_region_boundary(
    row: dict, boundary: str, length: int, offset: int = 0
) -> str:
    """Delete *length* nucleotides near the start of the *boundary* region.

    Parameters
    ----------
    row : dict
        Ground truth row from test_50.csv.
    boundary : str
        Region name whose start position defines the boundary.
    length : int
        Number of nucleotides to delete.
    offset : int
        Offset from the boundary start.  ``0`` deletes starting exactly at the
        boundary; negative values start the deletion in the preceding region.

    Returns
    -------
    str
        Modified full sequence.
    """
    boundaries = get_region_boundaries(row)
    pos = boundaries[boundary][0] + offset
    return delete_at_position(row["sequence_input"], pos, length)


# ---------------------------------------------------------------------------
# Truncation
# ---------------------------------------------------------------------------


def truncate_5prime(row: dict, num_bases: int) -> str:
    """Remove the first *num_bases* nucleotides from the input sequence."""
    return row["sequence_input"][num_bases:]


def truncate_3prime(row: dict, num_bases: int) -> str:
    """Remove the last *num_bases* nucleotides from the input sequence."""
    if num_bases <= 0:
        return row["sequence_input"]
    return row["sequence_input"][:-num_bases]


def truncate_to_region(row: dict, start_region: str) -> str:
    """Truncate the sequence so it starts at *start_region*.

    All regions before *start_region* are removed.
    """
    boundaries = get_region_boundaries(row)
    start = boundaries[start_region][0]
    return row["sequence_input"][start:]


# ---------------------------------------------------------------------------
# Point mutations
# ---------------------------------------------------------------------------


def mutate_position(
    sequence: str, position: int, replacement_base: str | None = None
) -> str:
    """Replace the base at *position* with *replacement_base*.

    If *replacement_base* is ``None``, the first non-identity alternative from
    ``_BASE_ALTERNATIVES`` is used (fully deterministic, no RNG).

    This is a low-level helper and does **not** check for stop codons.  Use
    :func:`mutate_near_boundary` when stop-codon avoidance is required.
    """
    original = sequence[position]
    if replacement_base is None:
        alts = _BASE_ALTERNATIVES.get(original, ["A", "T", "C"])
        replacement_base = alts[0]
    return sequence[:position] + replacement_base + sequence[position + 1 :]


def _would_create_stop(sequence: str, position: int, replacement: str) -> bool:
    """Check if substituting *replacement* at *position* creates a stop codon.

    Assumes the reading frame starts at position 0 of *sequence* (true for all
    ground truth sequences where ``sequence_input == fwr1+cdr1+...+fwr4``).
    """
    codon_start = (position // 3) * 3
    codon_end = codon_start + 3
    if codon_end > len(sequence):
        return False  # partial codon at sequence end — cannot form a stop
    codon = list(sequence[codon_start:codon_end])
    codon[position - codon_start] = replacement
    return "".join(codon).upper() in {"TAA", "TAG", "TGA"}


def mutate_near_boundary(
    row: dict,
    boundary: str,
    num_mutations: int,
    window: int = 6,
    seed: int = 0,
) -> str:
    """Apply *num_mutations* point mutations in a window around *boundary*.

    Mutations are distributed across a ``2 * window`` nucleotide window
    centred on the start of *boundary*.  Uses a seeded RNG for deterministic
    position selection.  Avoids identity substitutions and stop codons.

    Parameters
    ----------
    row : dict
        Ground truth row from test_50.csv.
    boundary : str
        Region name whose start defines the centre of the mutation window.
    num_mutations : int
        Number of point mutations to introduce.
    window : int
        Half-width of the mutation window (default 6 → 12-nt window).
    seed : int
        RNG seed for reproducibility.

    Returns
    -------
    str
        Modified full sequence.
    """
    boundaries = get_region_boundaries(row)
    centre = boundaries[boundary][0]
    seq = row["sequence_input"]

    # Candidate positions clamped to sequence bounds
    start = max(0, centre - window)
    end = min(len(seq), centre + window)
    candidates = list(range(start, end))

    rng = random.Random(seed)
    positions = rng.sample(candidates, min(num_mutations, len(candidates)))
    positions.sort()  # process left-to-right for deterministic codon context

    seq_list = list(seq)
    for pos in positions:
        original = seq_list[pos]
        alts = _BASE_ALTERNATIVES.get(original, ["A", "T", "C"])
        for alt in alts:
            if not _would_create_stop("".join(seq_list), pos, alt):
                seq_list[pos] = alt
                break
        # If every alternative would create a stop (extremely unlikely), the
        # position is left unmutated — this keeps the function safe.

    return "".join(seq_list)


# ---------------------------------------------------------------------------
# CDR3 length modifications
# ---------------------------------------------------------------------------


def extend_cdr3(row: dict, extension: str) -> str:
    """Insert *extension* at the midpoint of CDR3.

    Returns the full modified sequence.
    """
    boundaries = get_region_boundaries(row)
    cdr3_start, cdr3_end = boundaries["cdr3"]
    midpoint = cdr3_start + (cdr3_end - cdr3_start) // 2
    return insert_at_position(row["sequence_input"], midpoint, extension)


def introduce_stop_codon(
    row: dict, region: str, codon_index: int = 1, stop_codon: str = "taa"
) -> str:
    """Replace the codon at *codon_index* within *region* with a stop codon.

    Parameters
    ----------
    row : dict
        Ground truth row from test_50.csv.
    region : str
        Region in which to introduce the stop codon.
    codon_index : int
        0-based codon index within the region (default 1, i.e. second codon).
    stop_codon : str
        Stop codon to use (default ``"taa"``).

    Returns
    -------
    str
        Modified full sequence.
    """
    boundaries = get_region_boundaries(row)
    region_start = boundaries[region][0]
    pos = region_start + codon_index * 3
    seq = row["sequence_input"]
    return seq[:pos] + stop_codon + seq[pos + 3 :]


def shorten_cdr3(row: dict, bases_to_remove: int) -> str:
    """Remove *bases_to_remove* nucleotides from the centre of CDR3.

    Keeps equal portions at the 5' and 3' ends of CDR3.  When
    *bases_to_remove* equals or exceeds the CDR3 length, CDR3 is emptied.

    Returns the full modified sequence.
    """
    boundaries = get_region_boundaries(row)
    cdr3_start, cdr3_end = boundaries["cdr3"]
    cdr3_len = cdr3_end - cdr3_start

    if bases_to_remove >= cdr3_len:
        # Remove entire CDR3
        return row["sequence_input"][:cdr3_start] + row["sequence_input"][cdr3_end:]

    # Keep equal flanks, remove from centre
    keep = cdr3_len - bases_to_remove
    left_keep = keep // 2
    right_keep = keep - left_keep

    del_start = cdr3_start + left_keep
    del_end = cdr3_end - right_keep

    return row["sequence_input"][:del_start] + row["sequence_input"][del_end:]
