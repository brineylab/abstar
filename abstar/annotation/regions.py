# Copyright (c) 2024 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT


from .antibody import Antibody
from .positions import (
    get_gapped_position_from_raw,
    get_raw_position_from_aligned,
    get_raw_position_from_gapped,
)


def get_region_sequence(
    region: str,
    aligned_sequence: str,
    aligned_germline: str,
    gapped_germline: str,
    germline_start: int,
    ab: Antibody,
    aa: bool = False,
) -> str:
    """
    Get the sequence of a region from the aligned sequence.

    Parameters
    ----------
    region : str
        The region to get the sequence of.

    aligned_sequence : str
        The aligned sequence to get the region from.

    aligned_germline : str
        The aligned germline sequence to get the region from.

    gapped_germline : str
        The gapped germline sequence to get the region from.

    germline_start : int
        The start position of the germline sequence.

    ab : Antibody
        The ``Antibody`` object. Used only for logging. No ``Antibody`` paremeters
        are required or updated.

    aa : bool, default False
        Whether to get the region in amino acids or nucleotides.

    Returns
    -------
    str
        The sequence of the region.

    """
    # IMGT start/end positions for the region
    if aa:
        imgt_start = IMGT_REGION_START_POSITIONS_AA[region]
        imgt_end = IMGT_REGION_END_POSITIONS_AA[region]
    else:
        imgt_start = IMGT_REGION_START_POSITIONS_NT[region]
        imgt_end = IMGT_REGION_END_POSITIONS_NT[region]
    ab.log(f"{region.upper()} IMGT START:", imgt_start)
    ab.log(f"{region.upper()} IMGT END:", imgt_end)

    # if the input sequence is sufficiently truncated at
    # the 5' end that the region is entirely missing...
    gapped_germline_start = get_gapped_position_from_raw(
        position=germline_start + 1,  # needs to be 1-indexed
        gapped_germline=gapped_germline,
    )
    if gapped_germline_start > imgt_end:
        return ""

    # positions in the aligned germline
    aligned_start = get_raw_position_from_gapped(
        position=imgt_start,
        gapped_germline=gapped_germline,
        germline_start=germline_start,
    )
    aligned_end = get_raw_position_from_gapped(
        position=imgt_end,
        gapped_germline=gapped_germline,
        germline_start=germline_start,
    )
    ab.log(f"{region.upper()} ALIGNED START:", aligned_start)
    ab.log(f"{region.upper()} ALIGNED END:", aligned_end)

    # positions in the unaligned sequence
    start = get_raw_position_from_aligned(
        # position=aligned_start + 1,  # needs to be 1-indexed
        position=aligned_start, # no need to do +1 as this is already 1-indexed (cf below)
        aligned_sequence=aligned_sequence,
        aligned_reference=aligned_germline,
    )
    end = get_raw_position_from_aligned(
        # position=aligned_end + 1,  # needs to be 1-indexed
        position=aligned_end, # no need to do +1 as this is already 1-indexed (cf below)
        aligned_sequence=aligned_sequence,
        aligned_reference=aligned_germline,
    )
    ab.log(f"{region.upper()} RAW START:", start - 1)  # actual slicing start
    ab.log(f"{region.upper()} RAW END:", end)

    # collect the sequence
    region_sequence = aligned_sequence[start - 1 : end].replace("-", "")
    return region_sequence


# NOTE: these are the the actual IMGT end positions which are 1-indexed and
# not suitable for slicing in Python. To use these end points in a slice, subtract 1
# from the start position and use the end position unchanged, like so:
# ```
# start = IMGT_REGION_START_POSITIONS_NT["fwr1"] - 1
# end = IMGT_REGION_END_POSITIONS_NT["fwr1"]
# fwr1 = seq[start:end]
# ```
IMGT_REGION_END_POSITIONS_AA = {
    "fwr1": 26,
    "cdr1": 38,
    "fwr2": 55,
    "cdr2": 65,
    "fwr3": 104,
    "cdr3": 117,
    "fwr4": 129,
}

IMGT_REGION_END_POSITIONS_NT = {
    "fwr1": 78,
    "cdr1": 114,
    "fwr2": 165,
    "cdr2": 195,
    "fwr3": 312,
    "cdr3": 351,
    "fwr4": 387,
}

IMGT_REGION_START_POSITIONS_AA = {
    "fwr1": 1,
    "cdr1": 27,
    "fwr2": 39,
    "cdr2": 56,
    "fwr3": 66,
    "cdr3": 105,
    "fwr4": 118,
}

IMGT_REGION_START_POSITIONS_NT = {
    "fwr1": 1,
    "cdr1": 79,
    "fwr2": 115,
    "cdr2": 166,
    "fwr3": 196,
    "cdr3": 313,
    "fwr4": 352,
}
