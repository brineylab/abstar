# Copyright (c) 2024 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT


from .antibody import Antibody


VALID_NUCLEOTIDES = frozenset("ACGT")
JUNCTION_MOTIFS = {
    "IGH": "W",
    "IGK": "F",
    "IGL": "F",
    "TRA": "F",
    "TRB": "F",
    "TRD": "F",
    "TRG": "F",
}


def junction_is_in_frame(
    junction_start: int, v_sequence_start: int, frame: int,
) -> bool:
    """Compare oriented-query junction coordinates with the V-region frame.

    Both starts are zero-based positions in the oriented input; ``frame`` is
    one-based relative to the trimmed V region, not the original input.
    """
    return (junction_start - v_sequence_start - (frame - 1)) % 3 == 0


def assess_productivity(ab: Antibody) -> Antibody:
    """
    Checks whether an Antibody is productive and annotates any
    productivity issues.

    Parameters
    ----------

    ab : Antibody
        Antibody object to update with annotation information. The following ``Antibody``
        properties are updated:

        - ``productive``
        - ``productivity_issues``

        The following ``Antibody`` properties must be populated:

        - ``sequence``
        - ``sequence_aa``
        - ``v_call``
        - ``j_call``
        - ``junction_aa``
    """
    if isinstance(ab.productivity_issues, str):
        ab.productivity_issues = [
            issue for issue in ab.productivity_issues.split("|") if issue
        ]

    def add_issue(issue: str) -> None:
        if issue not in ab.productivity_issues:
            ab.productivity_issues.append(issue)

    # Scan for issues. Missing values are explicit failures rather than being
    # treated as evidence of productivity.
    if not ab.sequence_aa:
        add_issue("missing translated sequence")
    elif "*" in ab.sequence_aa:
        add_issue("stop codon(s)")
        ab.stop_codon = True

    if not ab.sequence or any(
        nucleotide not in VALID_NUCLEOTIDES for nucleotide in ab.sequence.upper()
    ):
        add_issue("ambiguous nucleotide(s)")

    v_locus = ab.v_call[:3].upper() if ab.v_call else None
    j_locus = ab.j_call[:3].upper() if ab.j_call else None
    locus = ab.locus.upper() if ab.locus else v_locus
    if v_locus is None or j_locus is None:
        add_issue("missing V/J gene call")
    elif v_locus != j_locus:
        add_issue(
            f"V/J locus mismatch ({ab.v_call} and {ab.j_call})"
        )

    junction_in_frame = True
    if not ab.junction_aa or len(ab.junction_aa) < 2:
        add_issue("missing or truncated junction")
        junction_in_frame = False
    else:
        if ab.junction_aa[0] != "C":
            add_issue("junction does not start with conserved C")
        expected_motif = JUNCTION_MOTIFS.get(locus)
        if expected_motif is None:
            add_issue(f"unsupported locus for junction motif ({locus})")
        elif ab.junction_aa[-1] != expected_motif:
            add_issue(f"junction does not end with conserved {expected_motif}")

    if ab.junction is not None:
        if len(ab.junction) < 6:
            add_issue("missing or truncated junction")
            junction_in_frame = False
        if len(ab.junction) % 3:
            add_issue("junction length is not a multiple of 3")
            junction_in_frame = False
        if any(
            nucleotide not in VALID_NUCLEOTIDES
            for nucleotide in ab.junction.upper()
        ):
            add_issue("ambiguous nucleotide(s) in junction")

    if ab.frame is not None:
        if ab.frame not in (1, 2, 3):
            add_issue(f"invalid reading frame ({ab.frame})")
            junction_in_frame = False
        elif getattr(ab, "junction_start", None) is not None and not junction_is_in_frame(
            ab.junction_start, ab.v_sequence_start, ab.frame
        ):
            add_issue("V/J junction is out of frame")
            junction_in_frame = False

    ab.vj_in_frame = junction_in_frame

    # flag sequences with issues as non-productive
    ab.productive = not ab.productivity_issues
    ab.productivity_issues = "|".join(ab.productivity_issues)

    return ab
