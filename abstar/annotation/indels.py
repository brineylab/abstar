# Copyright (c) 2024 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT


import re

from .antibody import Antibody
from .positions import get_gapped_position_from_raw

__all__ = [
    "annotate_insertions",
    "annotate_deletions",
]


def annotate_insertions(
    aligned_sequence: str,
    aligned_germline: str,
    gapped_germline: str,
    germline_start: int,
) -> Antibody:
    """
    Annotates non-templated insertions in the V gene region. The notation abstar uses for insertions is:

        ``"123:1>G!"`` (if the insertion causes a frameshift)
        ``"132:3>GGG"`` (if the insertion is in frame)

    In the case of multiple insertions, the insertions are separated by ``"|"`` characters, like so:

        ``"123:1>G!|132:3>GGG"``

    .. note::
        the position in abstar notation is the IMGT position after which the insertion occurs,
        so in the above example, the insertion ``"123:1>G!"`` would occur between positions 123 and 124

    Parameters
    ----------
    aligned_sequence : str
        The aligned query sequence for which insertions will be annotated.

    aligned_germline : str
        The aligned germline sequence.

    gapped_germline : str
        The IMGT-gapped germline sequence.

    germline_start : int
        The start position of the germline sequence in the alignment.

    ab : Antibody
        The ``Antibody`` object.

        The following ``Antibody`` properties are updated:

            - ``v_insertions``: a ``"|"``-separated string of non-templated insertions in the V gene region

        This function does not require any ``Antibody`` properties to be set.

    Returns
    -------
    Antibody
        The updated ``Antibody`` object.

    """
    # annotate insertions
    insertions = []
    for i in re.finditer("-+", aligned_germline):
        raw_start = i.start()
        raw_end = i.end()
        length = raw_end - raw_start
        seq = aligned_sequence[raw_start:raw_end]

        # get the aligned germline position of the insertion start
        start = len(aligned_germline[:raw_start].replace("-", "")) + germline_start
        imgt_start = get_gapped_position_from_raw(start, gapped_germline)

        # format the insertion
        oof = "!" if length % 3 else ""
        insertions.append(f"{imgt_start}:{length}>{seq}{oof}")
    return "|".join(insertions)


def annotate_deletions(
    aligned_sequence: str,
    aligned_germline: str,
    gapped_germline: str,
    germline_start: int,
) -> Antibody:
    """
    Annotates non-templated deletions in the V gene region. The notation abstar uses for deletions is:

        ``123:1>G!`` (if only single nucleotide in length)
        ``123-124:2>GG!`` (if multi-nucleotide and causing a frameshift)
        ``123-125:3>GGG`` (if multi-nucleotide and in frame)

    In the case of multiple deletions, the deletions are separated by ``"|"`` characters, like so:

        ``"123:1>G!|132-134:3>GGG"``

    .. note::
        The position range in abstar notation denotes the IMGT position of the start and end of the deletion,
        inclusive. This conforms to IMGT's convention for deletion position numbering.

    Parameters
    ----------
    aligned_sequence : str
        The aligned query sequence for which insertions will be annotated.

    aligned_germline : str
        The aligned germline sequence.

    gapped_germline : str
        The IMGT-gapped germline sequence.

    germline_start : int
        The start position of the germline sequence in the alignment.

    Returns
    -------
    list[str]
        A list of deletions in abstar notation.

    """
    # annotate deletions
    deletions = []
    for i in re.finditer("-+", aligned_sequence):
        raw_start = i.start()
        raw_end = i.end()
        length = raw_end - raw_start
        seq = aligned_germline[raw_start:raw_end]

        # get the aligned germline position of the deletion start
        #
        # unlike insertions, we want the start position of the actual deletion,
        # which is why we have the +1 on the start (for insertions, we want the
        # position immediately preceeding the insertion)
        start = len(aligned_germline[:raw_start].replace("-", "")) + germline_start + 1
        imgt_start = get_gapped_position_from_raw(start, gapped_germline)
        imgt_end = imgt_start + length

        # format the deletion
        oof = "!" if length % 3 else ""
        if length > 1:
            pos_range = f"{imgt_start}-{imgt_end - 1}"  # end is inclusive
        else:
            pos_range = f"{imgt_start}"
        deletions.append(f"{pos_range}:{length}>{seq}{oof}")
    return "|".join(deletions)


# # ---------------------------------------
# #       VARIABLE REGION INDELS
# # ---------------------------------------


# def annotate_v_insertions(
#     aligned_sequence: str,
#     aligned_germline: str,
#     gapped_germline: str,
#     germline_start: int,
#     ab: Antibody,
# ) -> Antibody:
#     """
#     Annotates non-templated insertions in the V gene region. The notation abstar uses for insertions is:

#         ``"123:1>G!"`` (if the insertion causes a frameshift)
#         ``"132:3>GGG"`` (if the insertion is in frame)

#     In the case of multiple insertions, the insertions are separated by ``"|"`` characters, like so:

#         ``"123:1>G!|132:3>GGG"``

#     .. note::
#         the position in abstar notation is the IMGT position after which the insertion occurs,
#         so in the above example, the insertion ``"123:1>G!"`` would occur between positions 123 and 124

#     Parameters
#     ----------
#     aligned_sequence : str
#         The aligned query sequence for which insertions will be annotated.

#     aligned_germline : str
#         The aligned germline sequence.

#     gapped_germline : str
#         The IMGT-gapped germline sequence.

#     germline_start : int
#         The start position of the germline sequence in the alignment.

#     ab : Antibody
#         The ``Antibody`` object.

#         The following ``Antibody`` properties are updated:

#             - ``v_insertions``: a ``"|"``-separated string of non-templated insertions in the V gene region

#         This function does not require any ``Antibody`` properties to be set.

#     Returns
#     -------
#     Antibody
#         The updated ``Antibody`` object.

#     """
#     # annotate insertions
#     insertions = annotate_insertions(
#         aligned_sequence=aligned_sequence,
#         aligned_germline=aligned_germline,
#         gapped_germline=gapped_germline,
#         germline_start=germline_start,
#     )
#     ab.v_insertions = "|".join(insertions)
#     return ab


# def annotate_v_deletions(
#     aligned_sequence: str,
#     aligned_germline: str,
#     gapped_germline: str,
#     germline_start: int,
#     ab: Antibody,
# ) -> Antibody:
#     """
#     Annotates non-templated deletions in the V gene region. The notation abstar uses for deletions is:

#         ``123:1>G!`` (if only single nucleotide in length)
#         ``123-124:2>GG!`` (if multi-nucleotide and causing a frameshift)
#         ``123-125:3>GGG`` (if multi-nucleotide and in frame)

#     In the case of multiple deletions, the deletions are separated by ``"|"`` characters, like so:

#         ``"123:1>G!|132-134:3>GGG"``

#     .. note::
#         The position range in abstar notation denotes the IMGT position of the start and end of the deletion,
#         inclusive. This conforms to IMGT's convention for deletion position numbering.

#     Parameters
#     ----------
#     aligned_sequence : str
#         The aligned query sequence for which insertions will be annotated.

#     aligned_germline : str
#         The aligned germline sequence.

#     gapped_germline : str
#         The IMGT-gapped germline sequence.

#     germline_start : int
#         The start position of the germline sequence in the alignment.

#     ab : Antibody
#         The ``Antibody`` object.

#         The following ``Antibody`` properties are updated:

#             - ``v_insertions``: a ``"|"``-separated string of non-templated insertions in the V gene region

#         This function does not require any ``Antibody`` properties to be set.

#     Returns
#     -------
#     Antibody
#         The updated ``Antibody`` object.

#     """
#     # annotate deletions
#     deletions = annotate_deletions(
#         aligned_sequence=aligned_sequence,
#         aligned_germline=aligned_germline,
#         gapped_germline=gapped_germline,
#         germline_start=germline_start,
#     )
#     ab.v_deletions = "|".join(deletions)
#     return ab


# # ---------------------------------------
# #       CONSTANT REGION INDELS
# # ---------------------------------------


# def annotate_c_insertions(
#     aligned_sequence: str,
#     aligned_germline: str,
#     gapped_germline: str,
#     germline_start: int,
#     ab: Antibody,
# ) -> Antibody:
#     """
#     Annotates non-templated insertions in the V gene region. The notation abstar uses for insertions is:

#         ``"123:1>G!"`` (if the insertion causes a frameshift)
#         ``"132:3>GGG"`` (if the insertion is in frame)

#     In the case of multiple insertions, the insertions are separated by ``"|"`` characters, like so:

#         ``"123:1>G!|132:3>GGG"``

#     .. note::
#         the position in abstar notation is the IMGT position after which the insertion occurs,
#         so in the above example, the insertion ``"123:1>G!"`` would occur between positions 123 and 124

#     Parameters
#     ----------
#     aligned_sequence : str
#         The aligned query sequence for which insertions will be annotated.

#     aligned_germline : str
#         The aligned germline sequence.

#     gapped_germline : str
#         The IMGT-gapped germline sequence.

#     germline_start : int
#         The start position of the germline sequence in the alignment.

#     ab : Antibody
#         The ``Antibody`` object.

#         The following ``Antibody`` properties are updated:

#             - ``v_insertions``: a ``"|"``-separated string of non-templated insertions in the V gene region

#         This function does not require any ``Antibody`` properties to be set.

#     Returns
#     -------
#     Antibody
#         The updated ``Antibody`` object.

#     """
#     # annotate insertions
#     insertions = annotate_insertions(
#         aligned_sequence=aligned_sequence,
#         aligned_germline=aligned_germline,
#         gapped_germline=gapped_germline,
#         germline_start=germline_start,
#     )
#     ab.c_insertions = "|".join(insertions)
#     return ab


# def annotate_c_deletions(
#     aligned_sequence: str,
#     aligned_germline: str,
#     gapped_germline: str,
#     germline_start: int,
#     ab: Antibody,
# ) -> Antibody:
#     """
#     Annotates non-templated deletions in the V gene region. The notation abstar uses for deletions is:

#         ``123:1>G!`` (if only single nucleotide in length)
#         ``123-124:2>GG!`` (if multi-nucleotide and causing a frameshift)
#         ``123-125:3>GGG`` (if multi-nucleotide and in frame)

#     In the case of multiple deletions, the deletions are separated by ``"|"`` characters, like so:

#         ``"123:1>G!|132-134:3>GGG"``

#     .. note::
#         The position range in abstar notation denotes the IMGT position of the start and end of the deletion,
#         inclusive. This conforms to IMGT's convention for deletion position numbering.

#     Parameters
#     ----------
#     aligned_sequence : str
#         The aligned query sequence for which insertions will be annotated.

#     aligned_germline : str
#         The aligned germline sequence.

#     gapped_germline : str
#         The IMGT-gapped germline sequence.

#     germline_start : int
#         The start position of the germline sequence in the alignment.

#     ab : Antibody
#         The ``Antibody`` object.

#         The following ``Antibody`` properties are updated:

#             - ``v_insertions``: a ``"|"``-separated string of non-templated insertions in the V gene region

#         This function does not require any ``Antibody`` properties to be set.

#     Returns
#     -------
#     Antibody
#         The updated ``Antibody`` object.

#     """
#     # annotate deletions
#     deletions = annotate_deletions(
#         aligned_sequence=aligned_sequence,
#         aligned_germline=aligned_germline,
#         gapped_germline=gapped_germline,
#         germline_start=germline_start,
#     )
#     ab.c_deletions = "|".join(deletions)
#     return ab
