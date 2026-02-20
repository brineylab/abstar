# Copyright (c) 2024 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT


import os
from typing import Tuple

import abutils
from abutils import Sequence
from natsort import natsorted

from .antibody import Antibody

__all__ = [
    "get_germline_database_path",
    "get_germline",
    # "realign_germline",
    # "reassign_dgene",
    # "process_vgene_alignment",
    # "process_jgene_alignment",
    # "process_dgene_alignment",
]


# ------------------------------
#      GERMLINE GENES
# ------------------------------


def get_germline_database_path(germdb_name: str, receptor: str = "bcr") -> str:
    """
    Get the path to a germline database. The addon directory (which contains user-built
    databases) is checked first, and if a user-built database is not found, the built-in
    database is used if present.

    Parameters
    ----------
    germdb_name : str
        The name of the germline database. Typically the species (like "human"), but can also be a
        custom database (like "humouse") that merges multiple species or is in some other way not
        species-specific.

        .. note::
            All germline database names are lowercase, and identifying germline databases
            by name is case-insensitive.

    receptor : str, default: "bcr"
        The receptor type. Options are "bcr" and "tcr".

    Returns
    -------
    str
        Path to the germline database directory. In the germline directory will be the
        following subdirectories:

          - ``ungapped``: contains ungapped germline sequences in FASTA format (``"v.fasta"``, ``"d.fasta"``, etc.)
          - ``gapped``: contains IMGT-gapped germline sequences in FASTA format (``"v.fasta"``, ``"d.fasta"``, etc.)
          - ``mmseqs``: contains pre-compiled MMseqs2 databases

    Raises
    ------
    FileNotFoundError
        If the germline database is not found.

    ValueError
        If the receptor type is not one of "bcr" or "tcr".

    """
    germdb_name = germdb_name.lower()
    receptor = receptor.lower()
    if receptor not in ["bcr", "tcr"]:
        raise ValueError(f"Receptor type {receptor} not supported")
    # check the addon directory first
    addon_dir = os.path.expanduser(f"~/.abstar/germline_dbs/{receptor}")
    if os.path.isdir(addon_dir):
        if germdb_name.lower() in [os.path.basename(d[0]) for d in os.walk(addon_dir)]:
            return os.path.join(addon_dir, germdb_name)
    # if a user-generated DB isn't found, use the built-in DB
    core_dir = os.path.dirname(os.path.abspath(__file__))  # "abstar/core" directory
    abstar_dir = os.path.dirname(core_dir)  # "abstar" directory
    germdb_path = os.path.join(abstar_dir, f"germline_dbs/{receptor}/{germdb_name}")
    if not os.path.exists(germdb_path):
        raise FileNotFoundError(
            f"Germline database {germdb_name} for receptor {receptor} not found"
        )
    return germdb_path


def get_germline(
    germline_gene: str,
    germdb_name: str,
    receptor: str = "bcr",
    imgt_gapped: bool = False,
    exact_match: bool = False,
    force_constant: bool = False,
    truncate_species: bool = True,
) -> list | Sequence:
    """
    Get the germline sequence for a given germline gene.

    Parameters
    ----------
    germline_gene : str
        The germline gene to get the sequence for. Can be a full germline gene name including allele ('IGHV1-2*02'),
        or a truncated version corresponding to just the gene ('IGHV1-2') or the family ('IGHV1').

        .. warning::
            Querying with gene names that are a substrings of another gene name (for example, 'IGHV1-2' is a substring
            of 'IGHV1-24') will match and return both 'IGHV1-2' and 'IGHV1-24' genes. To prevent this, you can add an
            asterisk to the end of the query gene ('IGHV1-2*') which will cause the query to match only alleles of the
            'IGHV1-2' gene.

        .. note::
            Germline genes are organized into separate files by segment ("v", "d", "j", and "c"). Segments are parsed as
            the fourth character of the germline gene name. For example, 'IGHV1-2*02' is segment 'v'.

            Isotype queries (ex. IGHG1*01 or IGHM*02) are all inferred as segment type 'c'.

    germdb_name : str
        The name of the germline database.

    receptor : str, default: "bcr"
        The receptor type. Options are "bcr" and "tcr".

    imgt_gapped : bool, default: False
        Whether to use the IMGT gapped or ungapped germline sequences.

    force_constant : bool, default: False
        Whether to force the query to match a constant germline gene. This is necessary because both IgD and diversity (D)
        gene names are formatted as IGHD, making it ambiguous whether the query is for IgD or diversity. By default,
        a supplied germline name of the format IGHD is assumed to be diversity.

    truncate_species : bool, default: True
        Whether to truncate the species name from the germline gene name (if it is present). Truncation
        happens before matching, so if `truncate_species` is ``False`` and `exact_match` is ``True``, the query
        must contain the species name if the germline gene name is to be found in the database.

    Returns
    -------
    list or Sequence
        If multiple matches are found (for example, multiple alleles of the query 'IGHV1-2'), a list of Sequence objects is returned.
        If only one match is found (for example, when querying 'IGHV1-2*02'), a single Sequence object is returned.

    Raises
    ------
    FileNotFoundError
        If the germline database file is not found.

    ValueError
        If the germline gene segment is not one of 'v', 'd', 'j', or 'c'. The germline gene segment is inferred as the
        fourth character of the germline gene name. For example, 'IGHV1-2*02' is segment 'v'. Isotype queries (ex.
        IGHG1*01 or IGHM*02) are all inferred as segment type 'c'.

    ValueError
        If no germline genes matching the query are found in the database.

    """
    germdb_path = get_germline_database_path(germdb_name, receptor)
    # parse the segment type
    segment = germline_gene[3].lower()
    if segment in ["a", "g", "e", "m"] or force_constant:
        segment = "c"
    if segment not in ["v", "d", "j", "c"]:
        raise ValueError(
            f"The segment type must be one of 'v', 'd', 'j', or an isotype ('a', 'd', 'e', 'g', 'm'). Your input was {germline_gene}, which corresponds to segment type {segment}."
        )
    # decide on gapped vs ungapped
    if imgt_gapped:
        germdb_file = os.path.join(germdb_path, f"imgt_gapped/{segment}.fasta")
    else:
        germdb_file = os.path.join(germdb_path, f"ungapped/{segment}.fasta")
    if not os.path.exists(germdb_file):
        raise FileNotFoundError(
            f"Germline database {germdb_name} for receptor {receptor} not found"
        )
    # retrieve germline sequences
    all_germs = abutils.io.read_fasta(germdb_file)
    if truncate_species:
        for g in all_germs:
            g.id = g.id.split("__")[0]
    # do the matching
    if exact_match:
        germs = [g for g in all_germs if germline_gene == g.id]
    else:
        germs = [g for g in all_germs if germline_gene in g.id]
    if not germs:
        raise ValueError(
            f"No substring matches to {germline_gene} were found in {germdb_file}"
        )
    if len(germs) == 1:
        return germs[0]
    else:
        return natsorted(germs, key=lambda x: x.id)


# ------------------------------
#     GERMLINE ALILGNMENT
# ------------------------------


def realign_germline(
    sequence: str,
    germline_name: str,
    germdb_name: str,
    imgt_gapped: bool = False,
    semiglobal_aln_params: dict | None = None,
    local_aln_params: dict | None = None,
    skip_semiglobal: bool = False,
    skip_local: bool = False,
    truncate_query: int | None = None,
    truncate_target: int | None = None,
    force_constant: bool = False,
) -> Tuple[abutils.tl.PairwiseAlignment | None]:
    """
    Performs an (optionally) two-step realignment of an assigned germline gene to a query sequence.

      * The first step is a semi-global alignment, which ensures that one or a handful of mismatches
        at the 5' or 3' end of the sequence do not result in incorrect truncation of a sequence that in fact
        should extend to the 5' or 3' end of the germline gene.

      * The second step is a local alignment, which fine-tunes the alignment with optimal alignment
        parameters and allowing truncation at the 3' (for V genes) or 5' (for J genes) end.

    Parameters
    ----------
    sequence : str
        The query sequence to align.

    germline_name : str
        The germline gene to align to.

    imgt_gapped : bool, default: False
        Whether to use the IMGT gapped or ungapped germline sequences. If ``True``, the IMGT gapped
        sequences are used.

    semiglobal_aln_params : dict, default: None
        The parameters to use for the semi-global alignment. Should be a ``dict`` containing one or
        more of the following keys:

        - ``match``: The match score (should be >= 0).
        - ``mismatch``: The mismatch score (must be <= 0).
        - ``gap_open``: The gap open penalty (must be <= 0).
        - ``gap_extend``: The gap extend penalty (must be <= 0).

    local_aln_params : dict, default: None
        The parameters to use for the local alignment. Should be a ``dict`` containing one or
        more of the following keys:

        - ``match``: The match score (should be >= 0).
        - ``mismatch``: The mismatch score (must be <= 0).
        - ``gap_open``: The gap open penalty (must be <= 0).
        - ``gap_extend``: The gap extend penalty (must be <= 0).

    skip_semiglobal : bool, default: False
        Whether to skip the semi-global alignment.

    skip_local : bool, default: False
        Whether to skip the local alignment.

    truncate_query : int, default: None
        The number of bases to truncate from the 3' end of the query sequence.

    truncate_target : int, default: None
        The number of bases to truncate from the 3' end of the target sequence.

    force_constant : bool, default: False
        Whether to force the query to match a constant germline gene. This is necessary because both IgD and diversity (D)
        gene names are formatted as IGHD, making it ambiguous whether the query is for IgD or diversity. By default,
        a supplied germline name of the format IGHD is assumed to be diversity.

    Returns
    -------
    Tuple[Optional[abutils.tl.PairwiseAlignment], Optional[abutils.tl.PairwiseAlignment]]
        The semi-global and local alignments. If ``skip_semiglobal`` is ``True``, the first element of the tuple is ``None``.
        If ``skip_local`` is ``True``, the second element of the tuple is ``None``.

    """
    # initial semi-global alignment
    germ = get_germline(
        germline_name,
        germdb_name,
        imgt_gapped=imgt_gapped,
        exact_match=True,
        force_constant=force_constant,
        truncate_species=False,
    )
    if truncate_query is not None:
        sequence = sequence[: -int(truncate_query)]
    if truncate_target is not None:
        germ = germ[: -int(truncate_target)]
    if skip_semiglobal:
        sg = None
    else:
        semiglobal_aln_params = (
            semiglobal_aln_params if semiglobal_aln_params is not None else {}
        )
        sg = abutils.tl.semiglobal_alignment(sequence, germ, **semiglobal_aln_params)
        sequence = sg.query[sg.query_begin : sg.query_end + 1]
        germ = sg.target[sg.target_begin : sg.target_end + 1]
    # secondary (fine-tuning) local alignment
    if skip_local:
        loc = None
    else:
        local_aln_params = local_aln_params if local_aln_params is not None else {}
        loc = abutils.tl.local_alignment(sequence, germ, **local_aln_params)
    return sg, loc


def reassign_dgene(
    sequence: str, germdb_name: str, aln_params: dict | None = None
) -> abutils.tl.PairwiseAlignment:
    """
    Assigns a D gene to a query sequence using local alignment.

    Parameters
    ----------
    sequence : str
        The query sequence to align.

    germdb_name : str
        The name of the germline database.

    aln_params : dict, default: None
        The parameters to use for the local alignment.


    Returns
    -------
    abutils.tl.PairwiseAlignment
        The highest scoring local alignment. If no suitable alignment is found, ``None`` is returned.

    """
    # fetch all d-genes
    germs = get_germline("IGHD", germdb_name, truncate_species=False)
    aln_params = aln_params if aln_params is not None else {}
    # align the query sequence to all d-genes
    alns = abutils.tl.local_alignment(sequence, targets=germs, **aln_params)
    if len(alns) > 0:
        # return the highest scoring alignment
        return sorted(alns, key=lambda x: x.score, reverse=True)[0]
    return


# ------------------------------
#     ALIGNMENT PROCESSING
# ------------------------------


def process_vgene_alignment(
    semiglobal_aln: abutils.tl.PairwiseAlignment,
    local_aln: abutils.tl.PairwiseAlignment,
    ab: Antibody,
) -> Antibody:
    """
    Processes a V gene alignment and updates ``Antibody`` annotations accordingly.

    .. note:
        all start/end positions are 0-indexed and end postions are
        exclusive, which aligns with Python slicing. This means that
        ``sequence[start : end]`` will work as expected


    Parameters
    ----------
    semiglobal_aln : abutils.tl.PairwiseAlignment
        Semi-global alignment of the query sequence and a V gene.

    local_aln : abutils.tl.PairwiseAlignment
        Local alignment of the query sequence and a V gene.

    ab : Antibody
        Antibody object to update with annotation information.

        The following ``Antibody`` properties are updated:

        - ``v_sequence``: the V gene region of the query sequence
        - ``v_sequence_aa``: the V gene region of the query sequence, translated into amino acids
        - ``v_score``: the score of the local alignment
        - ``v_sequence_start``: start position of the V gene in the query sequence, in positions corresponding to the ``sequence_oriented`` property
        - ``v_sequence_end``: end position of the V gene in the query sequence, in positions corresponding to the ``sequence_oriented`` property
        - ``v_germline``: the V gene region of the assigned germline sequence
        - ``v_germline_aa``: the V gene region of the assigned germline sequence, translated into amino acids
        - ``v_germline_start``: start position of the V gene in the assigned germline sequence
        - ``v_germline_end``: end position of the V gene in the assigned germline sequence
        - ``frame``: reading frame of the aligned query sequence

        Does not require any ``Antibody`` properties to be set.

    Returns
    -------
    Antibody
        Antibody object with updated annotations.

    """
    ab.v_score = local_aln.score
    # sequence/germline start position
    if semiglobal_aln.query_begin < semiglobal_aln.target_begin:
        # the query sequence doesn't extend to the beginning of the germline
        # so we'll use the start position of the local alignment
        ab.v_sequence_start = semiglobal_aln.query_begin + local_aln.query_begin
        ab.v_germline_start = semiglobal_aln.target_begin + local_aln.target_begin
    else:
        # if the query extends 5' at least to the beginning of the germline,
        # we set the start position at the beginning of the germline
        # to avoid an edge case where a mutation at or near the start of the
        # query sequence causes the start position to be incorrect
        ab.v_sequence_start = semiglobal_aln.query_begin
        ab.v_germline_start = semiglobal_aln.target_begin
    # sequence/germline stop position
    # the +1 on v_sequence_end lets us slice nicely going forward
    # so sequence_oriented[v_sequence_start : v_sequence_end] will give the v-gene region
    ab.v_sequence_end = semiglobal_aln.query_begin + local_aln.query_end + 1
    ab.v_germline_end = semiglobal_aln.target_begin + local_aln.target_end + 1
    # v-region sequence and germline
    ab.v_sequence = semiglobal_aln.query[ab.v_sequence_start : ab.v_sequence_end]
    ab.v_germline = semiglobal_aln.target[ab.v_germline_start : ab.v_germline_end]
    # frame -- 1-indexed, which works with abutils.tl.translate()
    ab.v_frame = (3 - (ab.v_germline_start % 3)) % 3 + 1
    ab.frame = ab.v_frame
    # AA sequence and germline
    ab.v_sequence_aa = abutils.tl.translate(ab.v_sequence, frame=ab.frame)
    ab.v_germline_aa = abutils.tl.translate(ab.v_germline, frame=ab.frame)
    return ab


def process_jgene_alignment(
    oriented_input: str,
    v_sequence_end: int,
    semiglobal_aln: abutils.tl.PairwiseAlignment,
    local_aln: abutils.tl.PairwiseAlignment,
    ab: Antibody,
) -> Antibody:
    """
    Processes a J gene alignment and updates ``Antibody`` annotations accordingly.

    .. note:
        all start/end positions are 0-indexed and end postions are
        exclusive, which aligns with Python slicing. This means that
        ``sequence[start : end]`` will work as expected

    Parameters:
    -----------
    oriented_input: str
        The oriented input sequence.

    v_sequence_end: int
        The end position of the V gene in the `oriented_input` sequence.

    semiglobal_aln: abutils.tl.PairwiseAlignment
        Semiglobal alignment of a query sequence to the J gene germline.

    local_aln: abutils.tl.PairwiseAlignment
        Local alignment of a query sequence to the J gene germline.

    ab: Antibody
        Antibody object to update with annotation information.

        The following ``Antibody`` properties are updated:

        - ``j_sequence``: the J gene region of the query sequence
        - ``j_score``: the score of the local alignment
        - ``j_sequence_start``: start position of the J gene in the query sequence, in positions corresponding to the ``sequence_oriented`` property
        - ``j_sequence_end``: end position of the J gene in the query sequence, in positions corresponding to the ``sequence_oriented`` property
        - ``j_germline``: the J gene region of the assigned germline sequence
        - ``j_germline_start``: start position of the J gene in the assigned germline sequence
        - ``j_germline_end``: end position of the J gene in the assigned germline sequence

        Does not require any ``Antibody`` properties to be set.

    """
    ab.j_score = local_aln.score

    ab.log("SEMIGLOBAL ALIGNMENT QUERY START:", semiglobal_aln.query_begin)
    ab.log("SEMIGLOBAL ALIGNMENT QUERY END:", semiglobal_aln.query_end)
    ab.log("LOCAL ALIGNMENT QUERY START:", local_aln.query_begin)
    ab.log("LOCAL ALIGNMENT QUERY END:", local_aln.query_end)
    ab.log("SEMIGLOBAL ALIGNMENT GERMLINE START:", semiglobal_aln.target_begin)
    ab.log("SEMIGLOBAL ALIGNMENT GERMLINE END:", semiglobal_aln.target_end)
    ab.log("LOCAL ALIGNMENT GERMLINE START:", local_aln.target_begin)
    ab.log("LOCAL ALIGNMENT GERMLINE END:", local_aln.target_end)
    ab.log("LOCAL ALIGNMENT QUERY LENGTH:", len(local_aln.query))
    ab.log("LOCAL ALIGNMENT GERMLINE LENGTH:", len(local_aln.target))

    # parse the J sequence and germline start position
    # sequence start is relative to the oriented input sequence
    # germline start is relative to the full (ungapped) germline gene sequence
    # AIRR-C wants 1-based indexing, but 0-based is better so we'll do that instead
    ab.j_sequence_start = (
        v_sequence_end + semiglobal_aln.query_begin + local_aln.query_begin
    )
    ab.j_germline_start = semiglobal_aln.target_begin + local_aln.target_begin

    # like V genes, we need to check whether the local alignment was incorrectly truncated
    # (on the 5' end for J genes) due to mutations at or near the end of the J gene
    # if the local alignment ends before the end of the germline J gene but the full query
    # sequence extends beyond the end of the full germline gene, we extend the alignment
    # to the end of the germline gene
    residual_germ = len(local_aln.target) - (local_aln.target_end + 1)
    residual_seq = len(local_aln.query) - (local_aln.query_end + 1)
    if residual_germ >= 1 and residual_seq >= residual_germ:
        plural = "S" if residual_germ > 1 else ""
        ab.log(
            f"USING SEMIGLOBAL ALIGNMENT END POSITION BECAUSE LOCAL ALIGNMENT WAS TRUNCATED BY {residual_germ} NUCLEOTIDE{plural}"
        )
        ab.j_sequence_end = v_sequence_end + semiglobal_aln.query_end + 1
        ab.j_germline_end = ab.j_germline_start + semiglobal_aln.target_end + 1
    else:
        ab.j_sequence_end = (
            ab.j_sequence_start + (local_aln.query_end - local_aln.query_begin) + 1
        )
        ab.j_germline_end = (
            ab.j_germline_start + (local_aln.target_end - local_aln.target_begin) + 1
        )

    # j-region sequence and germline
    ab.j_sequence = oriented_input[ab.j_sequence_start : ab.j_sequence_end]
    ab.j_germline = semiglobal_aln.target[ab.j_germline_start : ab.j_germline_end]
    return ab


# def process_jgene_alignment(
#     oriented_input: str,
#     v_sequence_end: int,
#     semiglobal_aln: abutils.tl.PairwiseAlignment,
#     local_aln: abutils.tl.PairwiseAlignment,
#     ab: Antibody,
# ) -> Antibody:
#     """
#     Processes a J gene alignment and updates ``Antibody`` annotations accordingly.

#     .. note:
#         all start/end positions are 0-indexed and end postions are
#         exclusive, which aligns with Python slicing. This means that
#         ``sequence[start : end]`` will work as expected

#     Parameters:
#     -----------
#     oriented_input: str
#         The oriented input sequence.

#     v_sequence_end: int
#         The end position of the V gene in the `oriented_input` sequence.

#     semiglobal_aln: abutils.tl.PairwiseAlignment
#         Semiglobal alignment of a query sequence to the J gene germline.

#     local_aln: abutils.tl.PairwiseAlignment
#         Local alignment of a query sequence to the J gene germline.

#     ab: Antibody
#         Antibody object to update with annotation information.

#         The following ``Antibody`` properties are updated:

#         - ``j_sequence``: the J gene region of the query sequence
#         - ``j_score``: the score of the local alignment
#         - ``j_sequence_start``: start position of the J gene in the query sequence, in positions corresponding to the ``sequence_oriented`` property
#         - ``j_sequence_end``: end position of the J gene in the query sequence, in positions corresponding to the ``sequence_oriented`` property
#         - ``j_germline``: the J gene region of the assigned germline sequence
#         - ``j_germline_start``: start position of the J gene in the assigned germline sequence
#         - ``j_germline_end``: end position of the J gene in the assigned germline sequence

#         Does not require any ``Antibody`` properties to be set.

#     """
#     ab.j_score = local_aln.score
#     # AIRR-C wants 1-based indexing, but 0-based is better so we'll do that instead
#     ab.j_sequence_start = (
#         v_sequence_end + semiglobal_aln.query_begin + local_aln.query_begin
#     )
#     ab.j_sequence_end = (
#         ab.j_sequence_start + (local_aln.query_end - local_aln.query_begin) + 1
#     )
#     # germline start/stop positions
#     ab.j_germline_start = semiglobal_aln.target_begin + local_aln.target_begin
#     ab.j_germline_end = (
#         ab.j_germline_start + (local_aln.target_end - local_aln.target_begin) + 1
#     )
#     # j-region sequence and germline
#     ab.j_sequence = oriented_input[ab.j_sequence_start : ab.j_sequence_end]
#     ab.j_germline = semiglobal_aln.target[ab.j_germline_start : ab.j_germline_end]
#     return ab


def process_dgene_alignment(
    oriented_input: str,
    v_sequence_end: int,
    local_aln: abutils.tl.PairwiseAlignment,
    ab: Antibody,
) -> Antibody:
    """
    Processes a D gene alignment and updates ``Antibody`` annotations accordingly.

    .. note:
        all start/end positions are 0-indexed and end postions are
        exclusive, which aligns with Python slicing. This means that
        ``sequence[start : end]`` will work as expected

    Parameters:
    -----------
    oriented_input: str
        The oriented input sequence.

    v_sequence_end: int
        The end position of the V gene in the `oriented_input` sequence.

    semiglobal_aln: abutils.tl.PairwiseAlignment
        Semiglobal alignment of a query sequence to the J gene germline.

    local_aln: abutils.tl.PairwiseAlignment
        Local alignment of a query sequence to the J gene germline.

    ab: Antibody
        Antibody object to update with annotation information.

        The following ``Antibody`` properties are updated:

        - ``d_call``: the D gene call (only updated if the D gene was found during a secondary search, not by the original MMseqs search)
        - ``d_gene``: the D gene name (only updated if the D gene was found during a secondary search, not by the original MMseqs search)
        - ``d_sequence``: the D gene region of the query sequence
        - ``d_score``: the score of the local alignment
        - ``d_sequence_start``: start position of the D gene in the query sequence, in positions corresponding to the ``sequence_oriented`` property
        - ``d_sequence_end``: end position of the D gene in the query sequence, in positions corresponding to the ``sequence_oriented`` property
        - ``d_germline``: the D gene region of the assigned germline sequence
        - ``d_germline_start``: start position of the D gene in the assigned germline sequence
        - ``d_germline_end``: end position of the D gene in the assigned germline sequence
        - ``d_frame``: reading frame of the D gene

        Does not require any ``Antibody`` properties to be set.

    """
    # if the D gene was found during a secondary search, not by the original MMseqs search,
    # we need to update the D gene call and name (since it hasn't already been set)
    if ab.d_call is None:
        ab.d_call = local_aln.target.id
        ab.d_gene = local_aln.target.id.split("*")[0]
    ab.d_score = local_aln.score
    # AIRR-C wants 1-based indexing, but 0-based is better so we'll do that instead
    ab.d_sequence_start = v_sequence_end + local_aln.query_begin
    d_length = local_aln.query_end - local_aln.query_begin
    ab.d_sequence_end = ab.d_sequence_start + d_length + 1
    # germline start/stop positions
    ab.d_germline_start = local_aln.target_begin
    d_germline_length = local_aln.query_end - local_aln.query_begin
    ab.d_germline_end = ab.d_germline_start + d_germline_length + 1
    # d-region sequence and germline
    ab.d_sequence = oriented_input[ab.d_sequence_start : ab.d_sequence_end]
    ab.d_germline = local_aln.target[ab.d_germline_start : ab.d_germline_end]
    # d frame
    ab.d_frame = (3 - (ab.d_germline_start % 3)) % 3 + 1
    return ab


def process_cgene_alignment(
    oriented_input: str,
    j_sequence_end: int,
    semiglobal_aln: abutils.tl.PairwiseAlignment,
    local_aln: abutils.tl.PairwiseAlignment,
    ab: Antibody,
) -> Antibody:
    """
    Processes a constant region alignment and updates ``Antibody`` annotations accordingly.

    .. note:
        all start/end positions are 0-indexed and end postions are
        exclusive, which aligns with Python slicing. This means that
        ``sequence[start : end]`` will work as expected

    Parameters:
    -----------
    oriented_input: str
        The oriented input sequence.

    j_sequence_end: int
        The end position of the J gene in the `oriented_input` sequence.

    semiglobal_aln: abutils.tl.PairwiseAlignment
        Semiglobal alignment of a query sequence to the constant region germline.

    local_aln: abutils.tl.PairwiseAlignment
        Local alignment of a query sequence to the constant region germline.

    ab: Antibody
        Antibody object to update with annotation information.

        The following ``Antibody`` properties are updated:

        - ``c_sequence``: the constant region of the query sequence
        - ``c_score``: the score of the local alignment
        - ``c_sequence_start``: start position of the constant region in the query sequence, in positions corresponding to the ``sequence_oriented`` property
        - ``c_sequence_end``: end position of the constant region in the query sequence, in positions corresponding to the ``sequence_oriented`` property
        - ``c_germline``: the constant region of the assigned germline sequence
        - ``c_germline_start``: start position of the constant region in the assigned germline sequence
        - ``c_germline_end``: end position of the constant region in the assigned germline sequence

        Does not require any ``Antibody`` properties to be set.

    """
    ab.c_score = local_aln.score
    # sequence/germline start position
    if semiglobal_aln.query_begin < semiglobal_aln.target_begin:
        # the query sequence doesn't extend to the beginning of the germline
        # so we'll use the start position of the local alignment
        ab.c_sequence_start = (
            j_sequence_end + semiglobal_aln.query_begin + local_aln.query_begin
        )
        ab.c_germline_start = semiglobal_aln.target_begin + local_aln.target_begin
    else:
        # if the query extends 5' at least to the beginning of the germline,
        # we set the start position at the beginning of the germline
        # to avoid an edge case where a mutation at or near the start of the
        # query sequence causes the start position to be incorrect
        ab.c_sequence_start = j_sequence_end + semiglobal_aln.query_begin
        ab.c_germline_start = semiglobal_aln.target_begin
    # sequence/germline stop positions
    c_length = local_aln.query_end - local_aln.query_begin
    ab.c_sequence_end = ab.c_sequence_start + c_length + 1
    ab.c_germline_end = semiglobal_aln.target_begin + local_aln.target_end + 1
    # j-region sequence and germline
    ab.c_sequence = oriented_input[ab.c_sequence_start : ab.c_sequence_end]
    ab.c_germline = semiglobal_aln.target[ab.c_germline_start : ab.c_germline_end]
    # frame -- 1-indexed, which works with abutils.tl.translate()
    ab.c_frame = (3 - (ab.c_germline_start % 3)) % 3 + 1
    # AA sequence and germline
    ab.c_sequence_aa = abutils.tl.translate(ab.c_sequence, frame=ab.c_frame)
    ab.c_germline_aa = abutils.tl.translate(ab.c_germline, frame=ab.c_frame)
    return ab
