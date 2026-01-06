# Copyright (c) 2025 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT


from .antibody import Antibody


def generate_cdr_mask(
    ab: Antibody,
    aa: bool = False,
    as_string: bool = True,
) -> str | list:
    """

    Create a mask for the CDR regions of an antibody.
    Regions are numbered 0-3, with 0 being any FWR, 1 being the CDR1, 2 being the CDR2, 3 being the CDR3, and 4 being the FWR4.
    Regions are defined according to the `IMGT numbering scheme`_ .

    Parameters:
    ----------
    ab : Antibody
        The antibody to generate a mask for.
    aa : bool, optional
        Whether to generate a mask for the amino acid sequence.
    as_string : bool, optional
        Whether to return the mask as a string. If False, the mask will be returned as a list.

    Returns:
    --------
    str | list: The mask for the CDR regions of the antibody.

    .. _IMGT numbering scheme: https://www.imgt.org/IMGTindex/numbering.php

    """
    # compute the mask
    if aa:
        cdr_mask = _generate_cdr_mask_aa(ab)
    else:
        cdr_mask = _generate_cdr_mask_nt(ab)

    # return
    if as_string:
        return "".join([str(m) for m in cdr_mask])
    return cdr_mask


def generate_gene_segment_mask(
    ab: Antibody,
    aa: bool = False,
    as_string: bool = True,
) -> str | list:
    """
    Create a mask for the gene segments of an antibody.
    """
    if aa:
        segment_mask = _generate_gene_segment_mask_aa(ab)
    else:
        segment_mask = _generate_gene_segment_mask_nt(ab)

    # return
    if as_string:
        return "".join([str(m) for m in segment_mask])
    return segment_mask


def generate_nongermline_mask(
    ab: Antibody,
    aa: bool = False,
    as_string: bool = True,
) -> str | list:
    """
    Create a mask for the non-germline regions of an antibody.

    Parameters
    ----------
    ab : Antibody
        The antibody to generate a mask for.
    aa : bool, optional
        Whether to generate a mask for the amino acid sequence.
    as_string : bool, optional
        Whether to return the mask as a string. If False, the mask will be returned as a list.

    Returns
    -------
    str | list: The mask for the non-germline regions of the antibody.
    """
    # need to use the aligned sequence/germline sequences
    # because indels would cause problems on non-aligned sequences
    if aa:
        segment_mask = ab.gene_segment_mask_aa
        sequence = ab.sequence_alignment_aa
        germline = ab.germline_alignment_aa
    else:
        segment_mask = ab.gene_segment_mask
        sequence = ab.sequence_alignment
        germline = ab.germline_alignment

    nongermline_mask = []
    mask_idx = 0
    for i in range(len(sequence)):
        s = sequence[i]
        g = germline[i]
        m = segment_mask[mask_idx]
        # if there's a deletion, don't increment the segment mask index
        # or add to the nongermline mask
        if s == "-":
            continue
        # N-addition regions are (by definition) non-germline
        if m == "N":
            nongermline_mask.append(1)
        # mutations are non-germline (insertions would be caught here too)
        elif s != g:
            nongermline_mask.append(1)
        # everything else is germline
        else:
            nongermline_mask.append(0)
        mask_idx += 1
        if mask_idx >= len(segment_mask):
            break

    # return
    if as_string:
        return "".join([str(m) for m in nongermline_mask])
    return nongermline_mask


def _generate_cdr_mask_nt(ab: Antibody) -> list:
    """
    Create a mask for the CDR regions of an antibody nucleotide sequence.

    Parameters
    ----------
    ab : Antibody
        The antibody to generate a mask for.

    Returns
    -------
    list: The mask for the CDR regions of the antibody.
    """
    cdr_mask = []
    cdr_mask.extend([0] * len(ab.fwr1))
    cdr_mask.extend([1] * len(ab.cdr1))
    cdr_mask.extend([0] * len(ab.fwr2))
    cdr_mask.extend([2] * len(ab.cdr2))
    cdr_mask.extend([0] * len(ab.fwr3))
    cdr_mask.extend([3] * len(ab.cdr3))
    cdr_mask.extend([0] * len(ab.fwr4))
    return cdr_mask


def _generate_cdr_mask_aa(ab: Antibody) -> list:
    """
    Create a mask for the CDR regions of an antibody amino acid sequence.

    Parameters
    ----------
    ab : Antibody
        The antibody to generate a mask for.

    Returns
    -------
    list: The mask for the CDR regions of the antibody.
    """
    cdr_mask = []
    cdr_mask.extend([0] * len(ab.fwr1_aa))
    cdr_mask.extend([1] * len(ab.cdr1_aa))
    cdr_mask.extend([0] * len(ab.fwr2_aa))
    cdr_mask.extend([2] * len(ab.cdr2_aa))
    cdr_mask.extend([0] * len(ab.fwr3_aa))
    cdr_mask.extend([3] * len(ab.cdr3_aa))
    cdr_mask.extend([0] * len(ab.fwr4_aa))
    return cdr_mask


def _generate_gene_segment_mask_nt(ab: Antibody) -> list:
    """
    Create a mask for the gene segments of an antibody nucleotide sequence.
    """
    segment_mask = []
    cdr3_start = ab.sequence.find(ab.cdr3)
    # there's an edge case in which the V gene is truncated so much that doesn't even contribute the full FWR3 region
    # to fix, we check the lengths of the V-gene region and the start position of the CDR3, and use the minimum of the two, then fill the rest of the FWR3 region with Ns
    if len(ab.v_sequence) < cdr3_start:
        segment_mask.extend(["V"] * len(ab.v_sequence))
        segment_mask.extend(["N"] * (cdr3_start - len(ab.v_sequence)))
    else:
        segment_mask.extend(["V"] * cdr3_start)
    segment_mask.extend(["V"] * len(ab.cdr3_v))
    segment_mask.extend(["N"] * len(ab.cdr3_n1))
    if ab.d_call is not None:
        segment_mask.extend(["D"] * len(ab.cdr3_d))
        segment_mask.extend(["N"] * len(ab.cdr3_n2))
    segment_mask.extend(["J"] * len(ab.cdr3_j))
    # there's an edge case in which the J gene is truncated so much that the FWR4 isn't entirely contained in the J gene
    # this usually happens because of overlap between the V and J gene assignments in light chains, resulting in part of the first FWR4 codon being assigned to the V gene
    # to fix, we extend the mask by the smaller of the J gene and FWR4 lengths
    segment_mask.extend(["J"] * min(len(ab.j_sequence), len(ab.fwr4)))
    return segment_mask


def _generate_gene_segment_mask_aa(ab: Antibody) -> list:
    """
    Create a mask for the gene segments of an antibody amino acid sequence.
    """
    segment_mask = []
    cdr3_start = ab.sequence_aa.find(ab.cdr3_aa)
    # there's an edge case in which the V gene is truncated so much that doesn't even contribute the full FWR3 region
    # to fix, we check the lengths of the V-gene region and the start position of the CDR3, and use the minimum of the two, then fill the rest of the FWR3 region with Ns
    if len(ab.v_sequence_aa) < cdr3_start:
        segment_mask.extend(["V"] * len(ab.v_sequence_aa))
        segment_mask.extend(["N"] * (cdr3_start - len(ab.v_sequence_aa)))
    else:
        segment_mask.extend(["V"] * cdr3_start)
    segment_mask.extend(["V"] * len(ab.cdr3_v_aa))
    segment_mask.extend(["N"] * len(ab.cdr3_n1_aa))
    if ab.d_call is not None:
        segment_mask.extend(["D"] * len(ab.cdr3_d_aa))
        segment_mask.extend(["N"] * len(ab.cdr3_n2_aa))
    segment_mask.extend(["J"] * len(ab.cdr3_j_aa))
    # there's an edge case in which the J gene is truncated so much that the FWR4 isn't entirely contained in the J gene
    # this usually happens because of overlap between the V and J gene assignments in light chains, resulting in part of the first FWR4 codon being assigned to the V gene
    # to fix, we check to see whether the existing segment_mask is longer than the antibody sequence (up to the end of the CDR3)
    # if so, we reduce the length of the J-gene portion of the segment_mask accordingly
    if len(segment_mask) > (cdr3_start + len(ab.cdr3_aa)):
        j_mask_length = len(ab.fwr4_aa) - (
            len(segment_mask) - (cdr3_start + len(ab.cdr3_aa))
        )
        segment_mask.extend(["J"] * j_mask_length)
    else:
        segment_mask.extend(["J"] * len(ab.fwr4_aa))
    return segment_mask
