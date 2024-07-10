# Copyright (c) 2024 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT


__all__ = [
    "get_gapped_position_from_raw",
    "get_raw_position_from_gapped",
    "get_raw_position_from_aligned",
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
