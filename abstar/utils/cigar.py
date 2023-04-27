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


def make_cigar(germline_segment):
    cigar = ""
    if germline_segment.query_start > 0:
        cigar += "{}S".format(germline_segment.query_start)
    if germline_segment.germline_start > 0:
        cigar += "{}N".format(germline_segment.germline_start)
    cigar += make_alignment_cigar(
        germline_segment.realignment.aligned_query,
        germline_segment.realignment.aligned_target,
    )
    return cigar


def get_cigar_code(q, t):
    if q == "-":
        return "D"
    if t == "-":
        return "I"
    return "M"


def make_alignment_cigar(query, target):
    prev = get_cigar_code(query[0], target[0])
    count = 1
    cigar = ""
    for q, t in zip(query[1:], target[1:]):
        curr = get_cigar_code(q, t)
        if prev is None:
            prev = curr
        elif curr == prev:
            count += 1
        else:
            count += 1
            cigar += "{}{}".format(count, prev)
            prev = curr
            count = 0
    cigar += "{}{}".format(count + 1, prev)
    return cigar
