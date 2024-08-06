#!/usr/bin/env python
# filename: cigar.py

#
# Copyright (c) 2020 Bryan Briney
# License: The MIT license (http://opensource.org/licenses/MIT)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software
# and associated documentation files (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge, publish, distribute,
# sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or
# substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
# BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#


def make_cigar(
    query_start: int,
    germline_start: int,
    aligned_query: str,
    aligned_germline: str,
) -> str:
    """
    Makes a CIGAR string from a realigned germline segment.

    Parameters
    ----------
    query_start : int
        The position in the query sequence where the alignment starts.

    germline_start : int
        The position in the germline sequence where the alignment starts.

    aligned_query : str
        The realigned query sequence.

    aligned_germline : str
        The realigned germline sequence.

    Returns
    -------
    str
        The CIGAR string for the realigned germline segment.

    """
    cigar = ""
    if query_start > 0:
        cigar += "{}S".format(query_start)
    if germline_start > 0:
        cigar += "{}N".format(germline_start)
    cigar += make_alignment_cigar(
        aligned_query,
        aligned_germline,
    )
    return cigar


# def make_cigar(germline_segment):
#     cigar = ""
#     if germline_segment.query_start > 0:
#         cigar += "{}S".format(germline_segment.query_start)
#     if germline_segment.germline_start > 0:
#         cigar += "{}N".format(germline_segment.germline_start)
#     cigar += make_alignment_cigar(
#         germline_segment.realignment.aligned_query,
#         germline_segment.realignment.aligned_target,
#     )
#     return cigar


def get_cigar_code(q: str, t: str) -> str:
    """
    Returns the CIGAR code for a given pair of query and target residues.

    Parameters
    ----------
    q : str
        The query residue.

    t : str
        The target residue.

    Returns
    -------
    str
        The CIGAR code for the given pair of residues.
    """
    if q == "-":
        return "D"
    if t == "-":
        return "I"
    # return "M"
    if q == t:
        return "="
    return "X"


# def get_cigar_code(q, t):
#     if q == "-":
#         return "D"
#     if t == "-":
#         return "I"
#     return "M"


def make_alignment_cigar(query: str, target: str) -> str:
    """
    Makes a CIGAR string from a pair of aligned sequences.

    Parameters
    ----------
    query : str
        The query sequence.

    target : str
        The target sequence.

    Returns
    -------
    str
        The CIGAR string for the aligned sequences.

    """
    prev = None
    count = 0
    cigar = ""
    for q, t in zip(query, target):
        curr = get_cigar_code(q, t)
        if prev is None:
            prev = curr
            count = 1
        elif curr == prev:
            count += 1
        else:
            cigar += "{}{}".format(count, prev)
            prev = curr
            count = 1
    cigar += "{}{}".format(count, prev)
    return cigar


# def make_alignment_cigar(query, target):
#     prev = get_cigar_code(query[0], target[0])
#     count = 1
#     cigar = ""
#     for q, t in zip(query[1:], target[1:]):
#         curr = get_cigar_code(q, t)
#         if prev is None:
#             prev = curr
#         elif curr == prev:
#             count += 1
#         else:
#             count += 1
#             cigar += "{}{}".format(count, prev)
#             prev = curr
#             count = 0
#     cigar += "{}{}".format(count + 1, prev)
#     return cigar
