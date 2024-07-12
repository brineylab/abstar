# Copyright (c) 2024 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT


from .antibody import Antibody


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
    # scan for issues
    if "*" in ab.sequence_aa:
        ab.productivity_issues.append("stop codon(s)")
    if "N" in ab.sequence:
        ab.productivity_issues.append("ambiguous nucleotide(s)")
    if ab.v_call[:3] != ab.j_call[:3]:
        ab.productivity_issues.append(
            f"V/J locus mismatch ({ab.v_call} and {ab.j_call})"
        )
    if ab.junction_aa:
        if ab.junction_aa[0] != "C":
            ab.productivity_issues.append("junction does not start with conserved C")
        if ab.junction_aa[-1] not in ["W", "F"]:
            ab.productivity_issues.append("junction does not end with conserved W/F")

    # flag sequences with issues as non-productive
    if ab.productivity_issues:
        ab.productive = False
    ab.productivity_issues = "|".join(ab.productivity_issues)

    return ab
