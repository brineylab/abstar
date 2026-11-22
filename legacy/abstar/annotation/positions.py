# Copyright (c) 2024 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT


from .antibody import Antibody


def alignment_columns_for_span(
    aligned_query: str, aligned_reference: str,
    query_origin: int, reference_origin: int,
    query_start: int, query_end: int, reference_start: int, reference_end: int,
) -> tuple[int, int]:
    """Map half-open query/reference spans to their shared alignment columns.

    Origins address the first ungapped residue represented by the trace. Each
    column advances only the coordinate whose residue is present; specifying
    both spaces makes boundaries adjacent to insertions/deletions unambiguous.
    A span that cuts outside the trace or cannot be represented by its columns
    raises ValueError rather than combining evidence from another alignment.
    """
    if len(aligned_query) != len(aligned_reference):
        raise ValueError('Alignment rows must have equal lengths')
    if query_end < query_start or reference_end < reference_start:
        raise ValueError('Alignment spans must be ordered')
    query, reference = query_origin, reference_origin
    boundaries = {(query, reference): 0}
    for column, (q, r) in enumerate(zip(aligned_query, aligned_reference), 1):
        if q == r == '-':
            raise ValueError('Alignment cannot contain double-gap columns')
        query += q != '-'
        reference += r != '-'
        boundaries[query, reference] = column
    try:
        start = boundaries[query_start, reference_start]
        end = boundaries[query_end, reference_end]
    except KeyError as error:
        raise ValueError('Requested span is not represented by the alignment') from error
    if end < start:
        raise ValueError('Alignment columns must be ordered')
    return start, end


__all__ = [
    "alignment_columns_for_span",
    "get_gapped_position_from_raw",
    "get_raw_position_from_gapped",
    "get_raw_position_from_aligned",
    "get_gapped_sequence",
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
        The raw (ungapped) position to convert to a gapped position. Zero
        denotes the boundary before the first germline base and maps to zero.

    gapped_germline : str
        The gapped germline sequence to convert the raw position to a gapped position.

    Returns
    -------
    int
        The IMGT-gapped position.

    """
    if position == 0:
        return 0
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


def get_aligned_position_from_ungapped(
    position: int,
    aligned_sequence: str,
    ab: Antibody | None = None,
    is_end_position: bool = False,
) -> int:
    """
    Get the aligned position from an ungapped position and an aligned reference sequence.

    Parameters
    ----------
    position : int
        The ungapped position to convert to an aligned position.

    aligned_sequence : str
        The aligned sequence to convert the ungapped position to an aligned position.

    ab : Antibody | None, optional
        The ``Antibody`` object. Used only for logging. No ``Antibody`` paremeters
        are required or updated. If not provided, no logging is done.

    is_end_position : bool, default False
        Whether the position is the end position of the region. If True, the function will search for indels immediately following the end position
        and will include them in the region position range.

    Returns
    -------
    int
        The aligned position.
    """
    ungapped_position = 0
    for aligned_position, s in enumerate(aligned_sequence):
        if s == "-":
            # ungapped_position += 1  # this might need to be removed if a semiglobal alignment is used (which may have leading gaps)
            continue
        if ungapped_position == position:
            if is_end_position:
                # edge case where an insertion happens between regions, which would cause both regions to ignore it --
                # the preceding region will end at the position before the insertion, and the subsequent region will skip it as a leading gap
                # our solution is to include insertions (gaps in the aligned_sequence, since this is the germline) in the preceding region
                if (
                    len(aligned_sequence) <= aligned_position + 1
                    or aligned_sequence[aligned_position + 1] != "-"
                ):
                    # no insertion, so we can return the aligned position
                    return aligned_position
                # former_aligned_position = aligned_position
                while aligned_sequence[aligned_position + 1] == "-":
                    aligned_position += 1
            return aligned_position
        else:
            ungapped_position += 1


def get_ungapped_position_from_aligned(
    position: int,
    aligned_sequence: str,
) -> int:
    """
    Get the ungapped position from an aligned position and an aligned sequence.
    """
    ungapped_position = 0
    for aligned_position, s in enumerate(aligned_sequence):
        if s == "-":
            continue
        if aligned_position == position:
            return ungapped_position
        ungapped_position += 1


def get_position_from_aligned_reference(
    position: int,
    aligned_sequence: str,
    aligned_reference: str,
) -> int:
    """
    Get the query position from an aligned reference.
    """
    query_position = 0
    reference_position = 0

    for s, r in zip(aligned_sequence, aligned_reference):
        # this has to come first, since we might have an alignment
        # where we want position 0 of the reference, but there are
        # leading gaps for which we need to increment the query first
        if r == "-":
            query_position += 1
        elif reference_position == position:
            return query_position
        elif s == "-":
            reference_position += 1
        else:
            query_position += 1
            reference_position += 1


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


def get_gapped_sequence(
    aligned_sequence: str,
    aligned_germline: str,
    gapped_germline: str,
    germline_start: int,
) -> str:
    """
    Get the gapped sequence from an aligned sequence and aligned germline.

    Parameters
    ----------
    aligned_sequence : str
        The aligned sequence to convert to a gapped sequence.

    aligned_germline : str
        The aligned germline sequence.

    gapped_germline : str
        The gapped germline sequence, which will be used to determine where gaps should be inserted.

    germline_start : int
        Zero-based start in the ungapped reference. At zero, leading IMGT dots
        are retained; otherwise output begins at the first retained residue.
        Internal dots are emitted before their following reference residue.

    Returns
    -------
    str
        The gapped sequence.

    """
    if len(aligned_sequence) != len(aligned_germline):
        raise ValueError('Alignment rows must have equal lengths')
    residues = [i for i, base in enumerate(gapped_germline) if base != '.']
    if not 0 <= germline_start <= len(residues):
        raise ValueError('Germline start lies outside the reference')
    if germline_start == 0:
        template_position = 0
    elif germline_start == len(residues):
        template_position = len(gapped_germline)
    else:
        template_position = residues[germline_start]
    output = []
    for query, reference in zip(aligned_sequence, aligned_germline):
        if reference != '-':
            while template_position < len(gapped_germline) and gapped_germline[template_position] == '.':
                output.append('.')
                template_position += 1
            if (template_position == len(gapped_germline)
                    or gapped_germline[template_position] != reference):
                raise ValueError('Alignment reference does not match its IMGT template')
            template_position += 1
        # Insertion columns, including terminal insertions, consume no
        # reference residue and must never move the template cursor.
        output.append(query)
    return ''.join(output)
