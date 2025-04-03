# Copyright (c) 2024 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

from typing import Optional

from .antibody import Antibody
from .positions import get_gapped_position_from_raw


def annotate_mutations(
    aligned_sequence: str,
    aligned_germline: str,
    gapped_germline: str,
    germline_start: int,
    ab: Optional[Antibody] = None,
    j_gene: Optional[bool] = False,
    debug: bool = False,
) -> Antibody:
    """
    Annotates mutations in an antibody germline gene region. The notation abstar uses for mutations is:

        ``"123:A>T"`` (if the mutation is a single nucleotide change)

    .. note::
        The mutation position in abstar notation is the IMGT position of the mutation,
        which is 1-indexed

    Parameters
    ----------
    aligned_sequence : str
        The aligned sequence of the query antibody region.

    aligned_germline : str
        The aligned sequence of the reference germline gene.

    gapped_germline : str
        The IMGT-gapped sequence of the reference germline gene.

    germline_start : int
        The starting position of the aligned germline sequence relative to the start of the
        complete germline gene.

    ab : Antibody, optional
        If provided, the ``Antibody`` object is only used for logging purposes.

    debug : bool, optional
        If ``True``, the mutation annotation process is logged to the ``Antibody`` object's
        logger.

    Returns
    -------
    mutations : list[str]
        List of mutations in abstar notation.

    """
    mutations = []
    raw_position = germline_start
    if ab is not None:
        ab.log("ALIGNED SEQUENCE:", aligned_sequence)
        ab.log("ALIGNED GERMLINE:", aligned_germline)
        ab.log("GAPPED GERMLINE:", gapped_germline)
        ab.log("GERMLINE START:", germline_start)
        ab.log("GERM      SEQ       RAW       IMGT      MUTATION  ")
        # ab.log("--------------------------------------------------")
    for s, g in zip(aligned_sequence, aligned_germline):
        log_str = f"{g}{' ' * (10 - len(g))}{s}{' ' * (10 - len(s))}"
        if s == g:
            # sequence and germline match -- no mutation and increment the position
            raw_position += 1
            log_str += f"{raw_position}{' ' * (10 - len(str(raw_position)))}"
            if debug:
                imgt_position = get_gapped_position_from_raw(
                    raw_position, gapped_germline
                )
                log_str += f"{imgt_position}{' ' * (10 - len(str(imgt_position)))}"
        elif g == "-":
            # sequence has insertion -- no mutation and don't increment the position
            log_str += f"{raw_position}{' ' * (10 - len(str(raw_position)))}"
            if debug:
                imgt_position = get_gapped_position_from_raw(
                    raw_position, gapped_germline
                )
                log_str += f"{imgt_position}{' ' * (10 - len(str(imgt_position)))}"
        elif s == "-":
            # sequence has a deletion -- no mutation but increment the position
            raw_position += 1
        else:
            # yay, we found a mutation!
            raw_position += 1
            if j_gene:
                mutation = f"{raw_position}:{g}>{s}"
            else:
                imgt_position = get_gapped_position_from_raw(raw_position, gapped_germline)
                mutation = f"{imgt_position}:{g}>{s}"
            mutations.append(mutation)
            log_str += f"{raw_position}{' ' * (10 - len(str(raw_position)))}"
            log_str += f"{imgt_position}{' ' * (10 - len(str(imgt_position)))}"
            log_str += f"{mutation}"
        if ab is not None:
            ab.log(log_str)
    return mutations


def annotate_v_mutations(
    aligned_sequence: str,
    aligned_germline: str,
    gapped_germline: str,
    germline_start: int,
    is_aa: bool,
    ab: Antibody,
    debug: bool = False,
) -> Antibody:
    """
    Annotates mutations in an antibody germline gene region. The notation abstar uses for mutations is:

        ``"123:A>T"`` (if the mutation is a single nucleotide change)

    .. note::
        The mutation position in abstar notation is the IMGT position of the mutation,
        which is 1-indexed

    Parameters
    ----------
    aligned_sequence : str
        The aligned sequence of the query antibody region.

    aligned_germline : str
        The aligned sequence of the reference germline gene.

    gapped_germline : str
        The IMGT-gapped sequence of the reference germline gene.

    germline_start : int
        The starting position of the aligned germline sequence relative to the start of the
        complete germline gene.

    is_aa : bool
        If ``True``, all input sequences contain amino acids (not nucleotides).

    ab : Antibody
        Antibody object to update with annotation information. The following ``Antibody``
        properties are updated:

        - ``v_mutations`` (if ``is_aa=False``)
        - ``v_mutation_count`` (if ``is_aa=False``)

        - ``v_mutations_aa`` (if ``is_aa=True``)
        - ``v_mutation_count_aa`` (if ``is_aa=True``)

    debug : bool, optional
        If ``True``, the mutation annotation process is logged to the ``Antibody`` object's
        logger.

    Returns
    -------
    ab : Antibody
        Antibody object with updated annotation information.

    """
    mutations = annotate_mutations(
        aligned_sequence=aligned_sequence,
        aligned_germline=aligned_germline,
        gapped_germline=gapped_germline,
        germline_start=germline_start,
        ab=ab,
        debug=debug,
    )
    if is_aa:
        ab.v_mutations_aa = "|".join(mutations)
        ab.v_mutation_count_aa = len(mutations)
    else:
        ab.v_mutations = "|".join(mutations)
        ab.v_mutation_count = len(mutations)
    return ab


def annotate_j_mutations(
	aligned_sequence: str,
    aligned_germline: str,
    gapped_germline: str,
    j_start_position: int,
    is_aa: bool,
    ab: Antibody,
    debug: bool = False,
) -> Antibody:
    """
    Docstring to be completed
    """
    mutations = annotate_mutations(
        aligned_sequence=aligned_sequence,
        aligned_germline=aligned_germline,
        gapped_germline=gapped_germline,
        germline_start=j_start_position,
        ab=ab,
        j_gene=True,
        debug=debug
    )
    if is_aa:
        ab.j_mutations_aa = "|".join(mutations)
        ab.j_mutation_count_aa = len(mutations)
    else:
        ab.j_mutations = "|".join(mutations)
        ab.j_mutation_count = len(mutations)
    return ab


def annotate_c_mutations(
    aligned_sequence: str,
    aligned_germline: str,
    gapped_germline: str,
    germline_start: int,
    is_aa: bool,
    ab: Antibody,
    debug: bool = False,
) -> Antibody:
    """
    Annotates mutations in an antibody germline gene region. The notation abstar uses for mutations is:

        ``"123:A>T"`` (if the mutation is a single nucleotide change)

    .. note::
        The mutation position in abstar notation is the IMGT position of the mutation,
        which is 1-indexed

    Parameters
    ----------
    aligned_sequence : str
        The aligned sequence of the query antibody region.

    aligned_germline : str
        The aligned sequence of the reference germline gene.

    gapped_germline : str
        The IMGT-gapped sequence of the reference germline gene.

    germline_start : int
        The starting position of the aligned germline sequence relative to the start of the
        complete germline gene.

    is_aa : bool
        If ``True``, all input sequences contain amino acids (not nucleotides).

    ab : Antibody
        Antibody object to update with annotation information. The following ``Antibody``
        properties are updated:

        - ``c_mutations`` (if ``is_aa=False``)
        - ``c_mutation_count`` (if ``is_aa=False``)

        - ``c_mutations_aa`` (if ``is_aa=True``)
        - ``c_mutation_count_aa`` (if ``is_aa=True``)

    debug : bool, optional
        If ``True``, the mutation annotation process is logged to the ``Antibody`` object's
        logger.

    Returns
    -------
    ab : Antibody
        Antibody object with updated annotation information.

    """
    mutations = annotate_mutations(
        aligned_sequence=aligned_sequence,
        aligned_germline=aligned_germline,
        gapped_germline=gapped_germline,
        germline_start=germline_start,
        ab=ab,
        debug=debug,
    )
    if is_aa:
        ab.c_mutations_aa = "|".join(mutations)
        ab.c_mutation_count_aa = len(mutations)
    else:
        ab.c_mutations = "|".join(mutations)
        ab.c_mutation_count = len(mutations)
    return ab
