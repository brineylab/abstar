"""Conservative FWR3 anchor recovery in raw oriented-query coordinates."""

from typing import NamedTuple


class JunctionAnchorError(ValueError):
    """The fallback cannot establish a unique, complete junction anchor."""


class JunctionAnchor(NamedTuple):
    start: int
    end: int
    score: int


def _best_projection(candidates, append=False, position=None):
    """Keep two distinct optimal projections: enough to prove ambiguity.

    Upstream tracebacks with the same anchor projection are equivalent here.
    Distinct partial projections cannot merge when more coordinates are appended.
    """
    best = max(score for score, _ in candidates)
    projections = set()
    for score, paths in candidates:
        if score != best:
            continue
        for path in paths:
            projections.add(path + (position,) if append else path)
            if len(projections) == 2:
                return best, frozenset(projections)
    return best, frozenset(projections)


def recover_junction_anchor(
    query: str, reference: str, *, query_origin: int,
    match: int = 3, mismatch: int = -2, gap_open: int = -35, gap_extend: int = -1,
) -> JunctionAnchor:
    """Fit FWR3 through codon 104 to an upstream-anchored raw query window.

    ``reference`` is ungapped FWR3 ending at the complete anchor codon. ``query``
    starts at its independently mapped upstream boundary and stops before J.
    Both leading ends are fixed; unused query suffix is free. Affine gaps cost
    ``gap_open + (length - 1) * gap_extend``, matching annotation alignment scores.
    Returned coordinates are zero-based, half-open in the oriented input query.

    Scores and anchor projections are propagated through all optimal tracebacks;
    two different optimal projections, gaps within the codon, and insufficient
    sequence support fail explicitly. No C motif or productivity preference is
    used. A mutated/ambiguous codon retains its observed bases and is assessed
    for productivity by the ordinary downstream code.
    """
    if len(reference) < 6 or not query or query_origin < 0:
        raise JunctionAnchorError("junction anchor recovery lacks FWR3 sequence")
    if any(c in "-." for c in query + reference):
        raise JunctionAnchorError("junction anchor recovery requires ungapped sequences")
    query, reference = query.upper(), reference.upper()
    width = len(query) + 1
    impossible = (float("-inf"), frozenset())
    # States: paired bases, deletion in query, insertion in query. Only the
    # preceding row is needed; signatures contain coordinates of the last codon.
    previous = [[impossible] * width for _ in range(3)]
    previous[0][0] = (0, frozenset([()]))
    for j in range(1, width):
        previous[2][j] = (gap_open + (j - 1) * gap_extend, frozenset([()]))
    for i, reference_base in enumerate(reference, 1):
        current = [[impossible] * width for _ in range(3)]
        in_anchor = i > len(reference) - 3
        for j in range(width):
            current[1][j] = _best_projection(
                [(previous[state][j][0] + (gap_extend if state == 1 else gap_open),
                  previous[state][j][1]) for state in range(3)],
                append=in_anchor, position=None,
            )
            if j == 0:
                continue
            substitution = match if query[j - 1] == reference_base and reference_base in "ACGT" else mismatch
            current[0][j] = _best_projection(
                [(previous[state][j - 1][0] + substitution, previous[state][j - 1][1])
                 for state in range(3)],
                append=in_anchor, position=query_origin + j - 1,
            )
            current[2][j] = _best_projection(
                [(current[state][j - 1][0] + (gap_extend if state == 2 else gap_open),
                  current[state][j - 1][1]) for state in range(3)],
            )
        previous = current
    score, projections = _best_projection([cell for row in previous for cell in row])
    if score <= 0 or not projections:
        raise JunctionAnchorError("junction anchor recovery lacks positive FWR3 alignment support")
    if len(projections) != 1:
        raise JunctionAnchorError("junction anchor recovery has competing optimal boundaries")
    positions = next(iter(projections))
    if len(positions) != 3 or None in positions or positions != tuple(range(positions[0], positions[0] + 3)):
        raise JunctionAnchorError("junction anchor recovery has a deleted, interrupted, or truncated codon")
    return JunctionAnchor(positions[0], positions[-1] + 1, int(score))


def recover_fwr3_anchor(ab, retained_query: str, retained_reference: str, **alignment_params):
    """Locate the upstream FWR3 boundary in retained V evidence, then fit raw bases."""
    if len(ab.v_germline_gapped) < 312 or any(c in ".-" for c in ab.v_germline_gapped[309:312]):
        raise JunctionAnchorError("V reference does not contain the complete IMGT anchor")
    reference_start = len(ab.v_germline_gapped[:195].replace(".", ""))
    reference_anchor = len(ab.v_germline_gapped[:309].replace(".", ""))
    query_position = ab.v_sequence_start
    reference_position = ab.v_germline_start
    query_origin = None
    retained_anchor = []
    for query_base, reference_base in zip(retained_query, retained_reference):
        if reference_base != "-" and reference_position == reference_start:
            if query_base != "-":
                query_origin = query_position
        if reference_base != "-" and query_base != "-" and reference_anchor <= reference_position < reference_anchor + 3:
            retained_anchor.append((reference_position - reference_anchor, query_position))
        query_position += query_base != "-"
        reference_position += reference_base != "-"
    if query_origin is None or query_origin >= ab.j_sequence_start:
        raise JunctionAnchorError("junction anchor recovery lacks a mapped upstream FWR3 boundary")
    anchor = recover_junction_anchor(
        ab.sequence_oriented[query_origin:ab.j_sequence_start],
        ab.v_germline_gapped[195:312].replace(".", ""),
        query_origin=query_origin, **alignment_params,
    )
    if any(anchor.start + offset != position for offset, position in retained_anchor):
        raise JunctionAnchorError("junction anchor recovery conflicts with retained V anchor evidence")
    return anchor
