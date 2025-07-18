# Copyright (c) 2024 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT


__all__ = [
    "get_gapped_position_from_raw",
    "get_raw_position_from_gapped",
    "get_raw_position_from_aligned",
    "get_gapped_sequence",
]


def get_gapped_position_from_raw(
    position: int,
    gapped_germline: str,
) -> int:
    """
    Get the IMGT-gapped germline position from a raw (ungapped) germline position.

    .. note::
        This function uses 1-based indexing, since that's what IMGT uses for position numbering.

    Parameters
    ----------
    position : int
        The raw (ungapped) position to convert to a gapped position.

    gapped_germline : str
        The gapped germline sequence to convert the raw position to a gapped position.

    Returns
    -------
    int
        The IMGT-gapped position.

    """
    raw = 0
    gapped = 0
    for res in gapped_germline:
        if res == ".":
            gapped += 1
        else:
            raw += 1
            gapped += 1
        if raw == position:
            break
    return gapped


def get_raw_position_from_gapped(
    position: int,
    gapped_germline: str,
    sequence_start: int = 0,
    germline_start: int = 0,
) -> int:
    """
    Get the raw (ungapped) germline position from an IMGT-gapped germline position.

    .. note::
        This function uses 1-based indexing, since that's what IMGT uses for position numbering.

    Parameters
    ----------
    position : int
        The IMGT-gapped position to convert to a raw (ungapped) position.

    gapped_germline : str
        The gapped germline sequence to convert the gapped position to a raw position.

    sequence_start : int, default 0
        The start position of the query sequence in the alignment with germline.

    germline_start : int, default 0
        The start position of the germline sequence in the alignment with the query sequence.

    Returns
    -------
    int
        The raw (ungapped) position.

    """
    raw = len(gapped_germline[:position].replace(".", ""))
    raw -= germline_start
    raw -= sequence_start
    return max([0, raw])


def get_position_from_aligned_reference(
    position: int,
    aligned_sequence: str,
    aligned_reference: str,
) -> int:
    """
    Get the query position from an aligned reference.
    """
    query_position = 0
    reference_position = 0

    for s, r in zip(aligned_sequence, aligned_reference):
        # this has to come first, since we might have an alignment
        # where we want position 0 of the reference, but there are
        # leading gaps for which we need to increment the query first
        if r == "-":
            query_position += 1
        elif reference_position == position:
            return query_position
        elif s == "-":
            reference_position += 1
        else:
            query_position += 1
            reference_position += 1


def get_raw_position_from_aligned(
    position: int,
    aligned_sequence: str,
    aligned_reference: str,
) -> int:
    """
    Gets the aligned position in a sequence, given an aligned reference.

    Parameters
    ----------
    position : int
        The aligned position to convert to a raw (ungapped) position.

    aligned_sequence : str
        The aligned query sequence from which to convert an aligned position to a raw position.

    aligned_reference : str
        The aligned reference sequence.

    Returns
    -------
    int
        The raw (ungapped) position.
    """
    aligned_position = 0
    raw_position = 0
    for s, r in zip(aligned_sequence, aligned_reference):
        if r == "-":
            aligned_position += 1
        else:
            aligned_position += 1
            raw_position += 1
        if raw_position == position:
            return aligned_position


def get_gapped_sequence(
    aligned_sequence: str,
    aligned_germline: str,
    gapped_germline: str,
    germline_start: int,
) -> str:
    """
    Get the gapped sequence from an aligned sequence and aligned germline.

    Parameters
    ----------
    aligned_sequence : str
        The aligned sequence to convert to a gapped sequence.

    aligned_germline : str
        The aligned germline sequence.

    gapped_germline : str
        The gapped germline sequence, which will be used to determine where gaps should be inserted.

    germline_start : int
        The start position of the germline sequence in the alignment with the query sequence.

    Returns
    -------
    str
        The gapped sequence.

    """
    gapped_seq = ""
    seq_pos = 0
    germ_pos = 0
    gap_pos = germline_start

    while gap_pos < len(gapped_germline):
        seq = aligned_sequence[seq_pos]
        germ = aligned_germline[germ_pos]
        gap = gapped_germline[gap_pos]
        if gap == "." and germ != "-":
            gapped_seq += gap
            gap_pos += 1
        else:
            gapped_seq += seq
            seq_pos += 1
            germ_pos += 1
            gap_pos += 1
        if seq_pos >= len(aligned_sequence):
            break

    return gapped_seq
