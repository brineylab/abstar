# Copyright (c) 2024 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT


from ..core.antibody import Antibody
from .positions import get_gapped_position_from_raw


def annotate_mutations(
    aligned_sequence: str,
    aligned_germline: str,
    gapped_germline: str,
    germline_start: int,
    is_aa: bool,
    ab: Antibody,
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

    Returns
    -------
    ab : Antibody
        Antibody object with updated annotation information.

    """
    mutations = []
    raw_position = germline_start
    for s, g in zip(aligned_sequence, aligned_germline):
        if s == g:
            # sequence and germline match -- no mutation and increment the position
            raw_position += 1
        elif g == "-":
            # sequence has insertion -- no mutation and don't increment the position
            continue
        elif s == "-":
            # sequence has a deletion -- no mutation but increment the position
            raw_position += 1
        else:
            # yay, we found a mutation!
            raw_position += 1
            imgt_position = get_gapped_position_from_raw(raw_position, gapped_germline)
            mutation = f"{imgt_position}:{g}>{s}"
            mutations.append(mutation)
    if is_aa:
        ab.v_mutations_aa = "|".join(mutations)
        ab.v_mutation_count_aa = len(mutations)
    else:
        ab.v_mutations = "|".join(mutations)
        ab.v_mutation_count = len(mutations)

    return ab
