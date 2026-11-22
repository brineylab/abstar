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
    if len(sequence) != len(germline):
        raise ValueError("Aligned sequence and germline must have equal lengths")
    ungapped_sequence_length = len(sequence.replace("-", ""))
    if len(segment_mask) != ungapped_sequence_length:
        raise ValueError(
            "Gene-segment mask length must equal the ungapped sequence length"
        )
    for s, g in zip(sequence, germline):
        # if there's a deletion, don't increment the segment mask index
        # or add to the nongermline mask
        if s == "-":
            continue
        m = segment_mask[mask_idx]
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
    """Label actual ungapped V/NP1/(D/NP2)/J spans in assembly order."""
    segment_mask = list(
        "V" * len(ab.v_sequence)
        + "N" * len(ab.np1)
        + "D" * len(ab.d_sequence or "")
        + "N" * len(ab.np2 or "")
        + "J" * len(ab.j_sequence)
    )
    if len(segment_mask) != len(ab.sequence.replace("-", "")):
        raise ValueError(
            "Gene-segment mask length must equal the ungapped assembled sequence length"
        )
    return segment_mask


def _generate_gene_segment_mask_aa(ab: Antibody) -> list:
    """Label complete query codons; mixed segment contributions are N.

    ``ab.frame`` is one-based relative to the assembled nucleotide query.
    Leading bases before that frame and terminal partial codons do not translate.
    """
    nt_mask = _generate_gene_segment_mask_nt(ab)
    segment_mask = []
    for start in range(ab.frame - 1, len(nt_mask) - 2, 3):
        codon = nt_mask[start : start + 3]
        segment_mask.append(codon[0] if len(set(codon)) == 1 else "N")
    if len(segment_mask) != len(ab.sequence_aa.replace("-", "")):
        raise ValueError(
            "Gene-segment mask length must equal the ungapped assembled amino acid sequence length"
        )
    return segment_mask
