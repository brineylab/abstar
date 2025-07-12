# Copyright (c) 2024 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT


import abutils

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
        position=aligned_start + 1,  # needs to be 1-indexed
        aligned_sequence=aligned_sequence,
        aligned_reference=aligned_germline,
    )
    end = get_raw_position_from_aligned(
        position=aligned_end + 1,  # needs to be 1-indexed
        aligned_sequence=aligned_sequence,
        aligned_reference=aligned_germline,
    )
    if start is not None:
        ab.log(f"{region.upper()} RAW START:", start - 1)  # actual slicing start
    else:
        ab.log(f"{region.upper()} RAW START: not found")
    if end is not None:
        ab.log(f"{region.upper()} RAW END:", end)
    else:
        ab.log(f"{region.upper()} RAW END: not found")

    # collect the sequence
    if start is not None and end is not None:
        region_sequence = aligned_sequence[start - 1 : end].replace("-", "")
    else:
        region_sequence = ""
    return region_sequence


def identify_cdr3_regions(ab: Antibody) -> Antibody:
    """
    Identify the CDR3 regions of an antibody.

    Parameters
    ----------
    ab : Antibody
        The ``Antibody`` object. Used only for logging. No ``Antibody`` paremeters
        are required or updated.

    Returns
    -------
    Antibody
        The ``Antibody`` object with the CDR3 regions identified.

        The following ``Antibody`` properties are updated:

        - ``cdr3_v``: the V gene region of the CDR3
        - ``cdr3_v_aa``: the V gene region of the CDR3 in amino acids
        - ``cdr3_n1``: the N1 region of the CDR3
        - ``cdr3_n1_aa``: the N1 region of the CDR3 in amino acids
        - ``cdr3_d``: the D region of the CDR3 (only for IGH/TRA/TRD chains with a D-gene call)
        - ``cdr3_d_aa``: the D region of the CDR3 in amino acids (only for IGH/TRA/TRD chains with a D-gene call)
        - ``cdr3_n2``: the N2 region of the CDR3 (only for IGH/TRA/TRD chains with a D-gene call)
        - ``cdr3_n2_aa``: the N2 region of the CDR3 in amino acids (only for IGH/TRA/TRD chains with a D-gene call)
        - ``cdr3_j``: the J gene region of the CDR3
        - ``cdr3_j_aa``: the J gene region of the CDR3 in amino acids

    """
    cdr3_start = ab.sequence.find(ab.fwr3) + len(ab.fwr3)
    cdr3_end = cdr3_start + len(ab.cdr3)

    # V gene region of the CDR3
    cdr3_v_length = len(ab.v_sequence) - cdr3_start
    adjusted_cdr3_v_length = cdr3_v_length - (cdr3_v_length % 3)
    cdr3_v_start = cdr3_start
    cdr3_v_end = cdr3_v_start + adjusted_cdr3_v_length
    ab.cdr3_v = ab.sequence[cdr3_v_start:cdr3_v_end]
    ab.cdr3_v_aa = abutils.tl.translate(ab.cdr3_v)

    # J gene region of the CDR3
    cdr3_j_start = ab.sequence.find(ab.j_sequence)
    cdr3_j_end = cdr3_end
    j_frame = (cdr3_j_start - cdr3_start) % 3
    j_trunc_5 = (-j_frame) % 3  # "wrap-around" modulo for trimming the front of the J
    adjusted_cdr3_j_start = cdr3_j_start + j_trunc_5
    ab.cdr3_j = ab.sequence[adjusted_cdr3_j_start:cdr3_j_end]
    ab.cdr3_j_aa = abutils.tl.translate(ab.cdr3_j)

    # chains with a D-gene call
    if ab.d_call is not None:
        cdr3_d_start = ab.sequence.find(ab.d_sequence)
        cdr3_d_end = cdr3_d_start + len(ab.d_sequence)
        d_start_frame = (cdr3_d_start - cdr3_start) % 3
        d_trunc_5 = (-d_start_frame) % 3  # "wrap-around" modulo
        d_trunc_3 = (cdr3_d_end - cdr3_start) % 3
        adjusted_cdr3_d_start = cdr3_d_start + d_trunc_5
        adjusted_cdr3_d_end = cdr3_d_end - d_trunc_3
        ab.cdr3_d = ab.sequence[adjusted_cdr3_d_start:adjusted_cdr3_d_end]
        ab.cdr3_d_aa = abutils.tl.translate(ab.cdr3_d)

        cdr3_n1_start = cdr3_v_end
        cdr3_n1_end = adjusted_cdr3_d_start
        ab.cdr3_n1 = ab.sequence[cdr3_n1_start:cdr3_n1_end]
        ab.cdr3_n1_aa = abutils.tl.translate(ab.cdr3_n1)

        cdr3_n2_start = adjusted_cdr3_d_end
        cdr3_n2_end = adjusted_cdr3_j_start
        ab.cdr3_n2 = ab.sequence[cdr3_n2_start:cdr3_n2_end]
        ab.cdr3_n2_aa = abutils.tl.translate(ab.cdr3_n2)

    # chains without a D-gene call
    else:
        cdr3_n1_start = cdr3_v_end
        cdr3_n1_end = adjusted_cdr3_j_start
        ab.cdr3_n1 = ab.sequence[cdr3_n1_start:cdr3_n1_end]
        ab.cdr3_n1_aa = abutils.tl.translate(ab.cdr3_n1)


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
