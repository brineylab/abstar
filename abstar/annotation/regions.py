# Copyright (c) 2024 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

from abutils.tools.alignment import PairwiseAlignment

from .antibody import Antibody
from .positions import (
    get_gapped_position_from_raw,
    get_position_from_aligned_reference,
    get_raw_position_from_gapped,
)


def get_region_sequence(
    region: str,
    aln: PairwiseAlignment,
    gapped_germline: str,
    germline_start: int,
    ab: Antibody,
    aa: bool = False,
) -> str:
    """
    Get the sequence of an antibody region (e.g., CDR1).

    Briefly, the germline version of the region of interest is parsed from the
    gapped germline and used to identify the sequence coordinates in the query
    sequence using an alignment of the query sequence to the full germline gene
    (for example, ``v_sg``, which is a semi-global alignment of the oriented input
    sequence and the full germline V gene).

    While this process is a little convoluted, but this is by design, as more
    straight-forward approaches have edge cases that make them less robust:

      - Using a local germline gene alignment (for example, ``v_loc``) can cause
        problems when extensive SHM causes sufficient mismatches at the 5' end
        such that the end of a region (most commonly FWR3) isn't present in the
        local alignment. This causes region parsing to fail.

      - Parsing the ungapped germline region from the ungapped sequence and then
        aligning just the germline region to the oriented input sequence can fail
        when there is an indel near either end of the region, combined with SHM
        in the dangling region on the terminal side of the indel. The gap penalty
        may cause the aligner to incorrectly annotate the indel (essentially, by
        aligning without gaps, the indel is "lost").

    Parameters
    ----------
    region : str
        The region to get the sequence of.

    aln : PairwiseAlignment
        Semi-global alignment of the query sequence to the full (ungapped) germline

    gapped_germline : str
        The gapped germline sequence.

    germline_start : int
        The germline position at which the query sequence begins.

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

    # positions in the ungapped germline
    ungapped_germline_start = get_raw_position_from_gapped(
        position=imgt_start,
        gapped_germline=gapped_germline,
        germline_start=germline_start,
    )
    ungapped_germline_end = get_raw_position_from_gapped(
        position=imgt_end,
        gapped_germline=gapped_germline,
        germline_start=germline_start,
    )
    ab.log(f"{region.upper()} UNGAPPED GERMLINE START:", ungapped_germline_start)
    ab.log(f"{region.upper()} UNGAPPED GERMLINE END:", ungapped_germline_end)

    # get the region start/end using the sequence aligned to germline
    # and the ungapped germline start/end positions
    region_start = get_position_from_aligned_reference(
        position=ungapped_germline_start,
        aligned_sequence=aln.aligned_query,
        aligned_reference=aln.aligned_target,
    )
    region_end = get_position_from_aligned_reference(
        position=ungapped_germline_end,
        aligned_sequence=aln.aligned_query,
        aligned_reference=aln.aligned_target,
    )
    ab.log(f"{region.upper()} REGION START:", region_start)
    ab.log(f"{region.upper()} REGION END:", region_end)

    # NOTE: alignment numbering is inclusive, so we +1 the end position for Python slicing
    region_sequence = aln.query[region_start : region_end + 1]

    return region_sequence


# def get_region_sequence(
#     region: str,
#     aligned_sequence: str,
#     aligned_germline: str,
#     gapped_germline: str,
#     germline_start: int,
#     ab: Antibody,
#     aa: bool = False,
# ) -> str:
#     """
#     Get the sequence of a region from the aligned sequence.

#     Parameters
#     ----------
#     region : str
#         The region to get the sequence of.

#     aligned_sequence : str
#         The aligned sequence to get the region from.

#     aligned_germline : str
#         The aligned germline sequence to get the region from.

#     gapped_germline : str
#         The gapped germline sequence to get the region from.

#     germline_start : int
#         The start position of the germline sequence.

#     ab : Antibody
#         The ``Antibody`` object. Used only for logging. No ``Antibody`` paremeters
#         are required or updated.

#     aa : bool, default False
#         Whether to get the region in amino acids or nucleotides.

#     Returns
#     -------
#     str
#         The sequence of the region.

#     """
#     # IMGT start/end positions for the region
#     if aa:
#         imgt_start = IMGT_REGION_START_POSITIONS_AA[region]
#         imgt_end = IMGT_REGION_END_POSITIONS_AA[region]
#     else:
#         imgt_start = IMGT_REGION_START_POSITIONS_NT[region]
#         imgt_end = IMGT_REGION_END_POSITIONS_NT[region]
#     ab.log(f"{region.upper()} IMGT START:", imgt_start)
#     ab.log(f"{region.upper()} IMGT END:", imgt_end)

#     # if the input sequence is sufficiently truncated at
#     # the 5' end that the region is entirely missing...
#     gapped_germline_start = get_gapped_position_from_raw(
#         position=germline_start + 1,  # needs to be 1-indexed
#         gapped_germline=gapped_germline,
#     )
#     if gapped_germline_start > imgt_end:
#         return ""

#     # positions in the aligned germline
#     aligned_start = get_raw_position_from_gapped(
#         position=imgt_start,
#         gapped_germline=gapped_germline,
#         germline_start=germline_start,
#     )
#     aligned_end = get_raw_position_from_gapped(
#         position=imgt_end,
#         gapped_germline=gapped_germline,
#         germline_start=germline_start,
#     )
#     ab.log(f"{region.upper()} ALIGNED START:", aligned_start)
#     ab.log(f"{region.upper()} ALIGNED END:", aligned_end)

#     # positions in the unaligned sequence
#     start = get_raw_position_from_aligned(
#         position=aligned_start + 1,  # needs to be 1-indexed
#         aligned_sequence=aligned_sequence,
#         aligned_reference=aligned_germline,
#     )
#     end = get_raw_position_from_aligned(
#         position=aligned_end + 1,  # needs to be 1-indexed
#         aligned_sequence=aligned_sequence,
#         aligned_reference=aligned_germline,
#     )
#     if start is not None:
#         ab.log(f"{region.upper()} RAW START:", start - 1)  # actual slicing start
#     else:
#         ab.log(f"{region.upper()} RAW START: not found")
#     if end is not None:
#         ab.log(f"{region.upper()} RAW END:", end)
#     else:
#         ab.log(f"{region.upper()} RAW END: not found")

#     # collect the sequence
#     if start is not None and end is not None:
#         region_sequence = aligned_sequence[start - 1 : end].replace("-", "")
#     else:
#         region_sequence = ""
#     return region_sequence


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
