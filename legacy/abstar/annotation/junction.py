"""Conservative query-space recovery of V regions and junction anchors."""

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


def _recover_query_positions(
    query: str, reference: str, *, query_origin: int, projection_positions: tuple[int, ...],
    fixed_positions: dict[int, int | None] | None = None,
    error_prefix: str = "junction anchor recovery",
    match: int = 3, mismatch: int = -2, gap_open: int = -35, gap_extend: int = -1,
) -> tuple[int, tuple[int | None, ...]]:
    """Fit a raw query, retaining joint projections across every optimal path.

    Fixed positions preserve retained reference-to-query base mappings; a None
    value requires the existing deletion. Unused query suffix is free.
    """
    fixed_positions = fixed_positions or {}
    if len(reference) < 6 or not query or query_origin < 0:
        raise JunctionAnchorError(f"{error_prefix} lacks FWR3 sequence")
    if any(c in "-." for c in query + reference):
        raise JunctionAnchorError(f"{error_prefix} requires ungapped sequences")
    query, reference = query.upper(), reference.upper()
    width = len(query) + 1
    impossible = (float("-inf"), frozenset())
    # States: paired bases, deletion in query, insertion in query. Only the
    # preceding row is needed; signatures contain requested boundary/anchor coordinates.
    previous = [[impossible] * width for _ in range(3)]
    previous[0][0] = (0, frozenset([()]))
    for j in range(1, width):
        previous[2][j] = (gap_open + (j - 1) * gap_extend, frozenset([()]))
    for i, reference_base in enumerate(reference, 1):
        current = [[impossible] * width for _ in range(3)]
        in_anchor = i - 1 in projection_positions
        for j in range(width):
            current[1][j] = _best_projection(
                [(previous[state][j][0] + (gap_extend if state == 1 else gap_open),
                  previous[state][j][1]) for state in range(3)],
                append=in_anchor, position=None,
            )
            if i - 1 in fixed_positions and fixed_positions[i - 1] is not None:
                current[1][j] = impossible
            if j == 0:
                continue
            substitution = match if query[j - 1] == reference_base and reference_base in "ACGT" else mismatch
            current[0][j] = _best_projection(
                [(previous[state][j - 1][0] + substitution, previous[state][j - 1][1])
                 for state in range(3)],
                append=in_anchor, position=query_origin + j - 1,
            )
            if i - 1 in fixed_positions and fixed_positions[i - 1] != query_origin + j - 1:
                current[0][j] = impossible
            current[2][j] = _best_projection(
                [(current[state][j - 1][0] + (gap_extend if state == 2 else gap_open),
                  current[state][j - 1][1]) for state in range(3)],
            )
        previous = current
    score, projections = _best_projection([cell for row in previous for cell in row])
    if score <= 0 or not projections:
        raise JunctionAnchorError(f"{error_prefix} lacks positive FWR3 alignment support")
    if len(projections) != 1:
        detail = f" (examples: {sorted(projections, key=repr)})" if error_prefix == "region boundary recovery" else ""
        raise JunctionAnchorError(f"{error_prefix} has competing optimal boundaries{detail}")
    return int(score), next(iter(projections))


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
    score, positions = _recover_query_positions(
        query, reference, query_origin=query_origin,
        projection_positions=tuple(range(len(reference) - 3, len(reference))),
        match=match, mismatch=mismatch, gap_open=gap_open, gap_extend=gap_extend,
    )
    if len(positions) != 3 or None in positions or positions != tuple(range(positions[0], positions[0] + 3)):
        raise JunctionAnchorError("junction anchor recovery has a deleted, interrupted, or truncated codon")
    return JunctionAnchor(positions[0], positions[-1] + 1, score)


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


class RegionBoundaryRecovery(NamedTuple):
    """Oriented-query intervals and anchor, independent of assignment evidence."""

    intervals: dict[str, tuple[int, int]]
    anchor: JunctionAnchor


def recover_v_region_boundaries(ab, retained_query: str, retained_reference: str, **alignment_params):
    """Recover only missing/contradictory mappings beyond a short retained V hit.

    Start at the last mapped region boundary, fit through IMGT codon 104 before
    J, and preserve every retained reference/query base mapping in that window.
    All requested region starts and all three anchor bases must have a single
    joint projection across optimal alignments. Motifs and productivity do not
    influence selection. Return None when ordinary region mapping is supported.
    """
    from .regions import IMGT_REGION_START_POSITIONS_NT

    regions = ("fwr1", "cdr1", "fwr2", "cdr2", "fwr3")
    gapped = ab.v_germline_gapped
    starts = {r: len(gapped[:IMGT_REGION_START_POSITIONS_NT[r] - 1].replace(".", "")) for r in regions}
    anchor_end = len(gapped[:312].replace(".", ""))
    missing_start = ab.v_germline_start <= starts["fwr3"] and ab.v_germline_end <= starts["fwr3"]
    conflicting_end = (ab.junction_start is not None and ab.v_germline_end < anchor_end
                       and ab.v_sequence_end > ab.junction_start + 3)
    if not (missing_start or conflicting_end):
        return None
    if len(gapped) < 312 or any(c in ".-" for c in gapped[309:312]):
        raise JunctionAnchorError("region boundary recovery requires the complete IMGT anchor")

    # Keys are ungapped reference positions; values are oriented-query bases.
    retained = {}
    rp, qp = ab.v_germline_start, ab.v_sequence_start
    for qbase, rbase in zip(retained_query, retained_reference):
        if rbase != "-":
            retained[rp] = qp if qbase != "-" else None
        rp += rbase != "-"
        qp += qbase != "-"
    origins = [r for r in regions if retained.get(starts[r]) is not None]
    if not origins:
        raise JunctionAnchorError("region boundary recovery lacks a mapped upstream region boundary")
    first = regions.index(origins[-1])
    ref_origin = starts[regions[first]]
    query_origin = retained[ref_origin]
    if not (ab.v_sequence_start <= query_origin < ab.j_sequence_start):
        raise JunctionAnchorError("region boundary recovery has an invalid upstream/J interval")
    reference = gapped[:312].replace(".", "")[ref_origin:]
    following = regions[first + 1:]
    cuts = tuple(starts[r] - ref_origin for r in following)
    anchor_positions = tuple(range(len(reference) - 3, len(reference)))
    if len(set(cuts + anchor_positions)) != len(cuts + anchor_positions):
        raise JunctionAnchorError("region boundary recovery has overlapping reference boundaries")
    score, positions = _recover_query_positions(
        ab.sequence_oriented[query_origin:ab.j_sequence_start], reference,
        query_origin=query_origin, projection_positions=cuts + anchor_positions,
        fixed_positions={p - ref_origin: q for p, q in retained.items() if ref_origin <= p < anchor_end},
        error_prefix="region boundary recovery", **alignment_params,
    )
    if None in positions:
        raise JunctionAnchorError("region boundary recovery has a deleted boundary or anchor base")
    anchor_bases = positions[-3:]
    if anchor_bases != tuple(range(anchor_bases[0], anchor_bases[0] + 3)):
        raise JunctionAnchorError("region boundary recovery has an interrupted anchor codon")
    boundaries = (query_origin,) + positions[:-3] + (anchor_bases[-1] + 1,)
    if any(left >= right for left, right in zip(boundaries, boundaries[1:])):
        raise JunctionAnchorError("region boundary recovery has non-increasing region boundaries")
    intervals = dict(zip(regions[first:], zip(boundaries, boundaries[1:])))
    return RegionBoundaryRecovery(intervals, JunctionAnchor(anchor_bases[0], anchor_bases[-1] + 1, score))
