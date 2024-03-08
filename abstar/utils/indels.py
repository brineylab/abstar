#!/usr/bin/env python
# filename: indels.py

#
# Copyright (c) 2016 Bryan Briney
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


# from __future__ import absolute_import, division, print_function, unicode_literals

import re
import traceback
from typing import Iterable, Optional

# from abutils import Sequence
# from abutils.utils import log
# from abutils.utils.decorators import lazy_property
# from ..core.antibody import Antibody
# from ..core.germline import Germline


class Indel(object):
    """Base class for insertions and deletions"""

    def __init__(self, indel):
        super(Indel, self).__init__()
        self.raw = indel
        self.length = indel.get("len", None)
        self.raw_position = indel.get("pos", None)
        self.sequence = indel.get("seq", None)
        self.fixed = indel.get("fixed", None)
        self.in_frame = self._get_frame()
        # IMGT numbering gets annotated later by the parent Germline object
        # once IMGT numbering for the germline segment is computed.
        self.imgt_position = None
        self.imgt_codon = None

    def __contains__(self, item):
        return item in self.raw.keys()

    def __getitem__(self, key):
        return self.raw.get(key, None)

    def _get_frame(self):
        if "in frame" in self.raw:
            if self.raw["in frame"] in ["yes", "no"]:
                return self.raw["in frame"]
            elif self.raw["in frame"] is True:
                return "yes"
            elif self.raw["in frame"] is False:
                return "no"
        return None


class Insertion(Indel):
    """
    Class for storing and formatting insertions.
    """

    def __init__(self, indel):
        super(Insertion, self).__init__(indel)
        self.type = "insertion"

    @property
    def imgt_formatted(self):
        """
        Returns the insertion in IMGT format.

        IMGT format for insertions is:
            ``123^124insGGG`` (if in frame)
            ``123^124insG#`` (if out of frame)
        Where ``123`` and ``124`` are the germline positions flanking the insertion,
        and ``GGG`` (or ``G``) is the inserted sequence.

        Returns:
        --------
            str: IMGT formatted insertion.

        """
        return "{}^{}>ins^{}{}".format(
            self.imgt_position,
            self.imgt_position + 1,
            self.sequence.lower(),
            "" if self.in_frame else "#",
        )

    @property
    def abstar_formatted(self):
        """
        Returns the insertion in abstar format.

        AbStar format for insertions is:
            ``123:3>GGG`` (if in frame)
            ``123:1!>G`` (if out of frame)
        Where ``123`` is the IMGT position immediately preceding the insertion, ``3`` is
        the length of the insertion, and ``GGG`` is the inserted sequence.

        Returns:
        --------
            str: AbStar formatted insertion.
        """
        return "{}:{}{}>{}".format(
            self.imgt_position,
            self.length,
            "" if self.in_frame else "!",
            self.sequence,
        )

    @property
    def json_formatted(self):
        """
        Returns the insertion in JSON format.

        JSON format for insertions is:
            ``{
                "in_frame": "yes",
                "length": 3,
                "sequence": "GGG",
                "position": "123",
                "codon": "41"
            }``
        Where ``in_frame`` is either "yes" or "no", ``length`` is the length of the insertion,
        ``sequence`` is the inserted sequence, ``position`` is the IMGT position immediately
        preceding the insertion, and ``codon`` is the IMGT codon number immediately preceding
        the insertion.

        Returns:
        --------
            dict: JSON formatted insertion.

        """
        j = {
            "in_frame": "yes" if self.in_frame else "no",
            "length": self.length,
            "sequence": self.sequence,
            "position": str(self.imgt_position),
            "codon": str(self.imgt_codon),
        }
        return j


class Deletion(Indel):
    """
    Class for storing and formatting deletions.
    """

    def __init__(self, indel):
        super(Deletion, self).__init__(indel)
        self.type = "deletion"

    @property
    def imgt_formatted(self):
        """
        Returns the deletion in IMGT format.

        IMGT format for deletions is:
            ``G123>del#`` (if a single nucleotide)
            ``G123-G124>del(2nt)#`` (if multi-nucleotide and out of frame)
            ``G123-G125>del(3nt)`` (if in frame)

        Returns:
        --------
            str: IMGT formatted deletion.
        """
        if self.length == 1:
            return "{}{}>del#".format(self.sequence.lower(), self.imgt_position)
        else:
            return "{}{}-{}{}>del({}nt){}".format(
                self.sequence[0].lower(),
                self.imgt_position,
                self.sequence[-1].lower(),
                self.imgt_position + self.length - 1,
                self.length,
                "" if self.in_frame else "#",
            )

    @property
    def abstar_formatted(self):
        """
        Returns the deletion in abstar format.

        abstar format for deletions is:
            ``123:1!>G`` (if a single nucleotide)
            ``123-124:2!>GG`` (if multi-nucleotide and out of frame)
            ``123-125:3>GGG`` (if in frame)

        Returns:
        --------
            str: abstar formatted deletion.

        """
        if self.length == 1:
            return "{}:{}{}>{}".format(
                self.imgt_position,
                self.length,
                "" if self.in_frame else "!",
                self.sequence,
            )
        else:
            return "{}-{}:{}{}>{}".format(
                self.imgt_position,
                self.imgt_position + self.length - 1,
                self.length,
                "" if self.in_frame else "!",
                self.sequence,
            )

    @property
    def json_formatted(self):
        """
        Returns the deletion in JSON format.

        JSON format for deletions is:
            ``{
                "in_frame": "yes",
                "length": 3,
                "sequence": "GGG",
                "position": "123",
                "codon": "41"
            }``
        Where ``in_frame`` is either "yes" or "no", ``length`` is the length of the deletion,
        ``sequence`` is the deleted sequence, ``position`` is the IMGT position immediately

        Returns:
        --------
            dict: JSON formatted deletion.

        """
        j = {
            "in_frame": "yes" if self.in_frame else "no",
            "length": self.length,
            "sequence": self.sequence,
            "position": str(self.imgt_position),
            "codon": str(self.imgt_codon),
        }
        return j


# ------------
#  Insertions
# ------------


def find_insertions(
    antibody: "Antibody", segment: "Germline"
) -> Optional[Iterable[Insertion]]:
    """
    Identifies and annotates/fixes insertions. Any insertions shorter than a full codon
    will be reverted to germline. All others will be annotated.

    .. note::
        If an Exception is raised during insertion annotation, the exception will be
        recorded in the segment's exception log but will not be raised.

    Parameters
    ----------
    antibody : Antibody
        The Antibody object being annotated.

    segment : Germline
        The Germline object to be annotated.

    Returns
    -------
    list[Insertion]
        A list of ``Insertion`` objects (one for each insertion), or ``None`` if there were
        no annotated insertions.

    """
    try:
        insertions = []
        o = 0
        for i in re.finditer("-+", segment.germline_alignment):
            s = i.start() - o
            e = i.end() - o
            ins_len = e - s
            ins_sequence = segment.query_alignment[s:e]
            if ins_len % 3 == 0 or ins_len > 3:
                insertions.append(
                    _annotate_insertion(s + segment.query_start, ins_len, ins_sequence)
                )
            else:
                insertions.append(
                    _annotate_insertion(
                        s + segment.query_start,
                        ins_len,
                        ins_sequence,
                        fixed=True,
                    )
                )
                _fix_frameshift_insertion(antibody, segment, s, e)
                o += ins_len
        return insertions
    except Exception:
        segment.exception("FIND INSERTIONS ERROR", traceback.format_exc())


def _annotate_insertion(
    start: int, length: int, sequence: str, fixed: bool = False
) -> Insertion:
    """
    Annotates codon-length (non-frameshift) insertions.

    Parameters
    ----------
    start : int
        The starting postion of the insertion.

    length : int
        The length of the insertion.

    sequence : str
        The inserted sequence.

    fixed : bool
        Whether the insertion was fixed (reverted to germline).

    Returns
    -------
    Insertion
        An ``Insertion`` object containing the insertion start position, the insertion length,
        and the inserted sequence.

    """
    in_frame = "yes" if length % 3 == 0 else "no"
    return Insertion(
        {
            "pos": start,
            "len": length,
            "seq": sequence,
            "in frame": in_frame,
            "fixed": fixed,
        }
    )


def _fix_frameshift_insertion(
    antibody: "Antibody", segment: "Germline", s: int, e: int
) -> None:
    """
    Fixes (removes) frameshift insertions.

    Insertion fixing entails several steps:
        1. Remove the frameshift indel from the alignment (query, germline, and midline).
        2. Remove the frameshift indel from the oriented input.
        3. Adjust the alignment end position to account for the removed indel.

    Parameters
    ----------
    antibody : Antibody
        The Antibody object being annotated.

    segment : Germline
        The Germline object to be annotated.

    s : int
        The starting postion of the insertion.

    e : int
        The ending position of the insertion.

    Returns
    -------
    None
        The input ``Germline`` and ``Antibody`` objects are modified in place.

    """
    # remove the frameshift indel from the alignment
    segment.query_alignment = segment.query_alignment[:s] + segment.query_alignment[e:]
    segment.germline_alignment = (
        segment.germline_alignment[:s] + segment.germline_alignment[e:]
    )
    segment.alignment_midline = (
        segment.alignment_midline[:s] + segment.alignment_midline[e:]
    )

    # remove the frameshift indel from the oriented input
    oi_firsthalf = antibody.oriented_input.sequence[: s + segment.query_start]
    oi_secondhalf = antibody.oriented_input.sequence[e + segment.query_start :]
    antibody.oriented_input.sequence = oi_firsthalf + oi_secondhalf

    # adjust the alignment end position
    segment.query_end -= e - s


# ------------
#  Deletions
# ------------


def find_deletions(
    antibody: "Antibody", segment: "Germline"
) -> Optional[Iterable[Deletion]]:
    """
    Identifies and annotates/fixes deletions. Deletions shorter than a full codon will be
    reverted to germline. All others will be annotated.

    .. note::
        If an Exception is raised during deletion annotation, the exception will be
        recorded in the segment's exception log but will not be raised.

    Parameters
    ----------
    antibody : Antibody
        The Antibody object being annotated.

    segment : Germline
        The Germline object to be annotated.

    Returns
    -------
    list[Deletion]
        A list of ``Deletion`` objects (one for each deletion), or ``None`` if there were
        no annotated deletions.

    """
    try:
        deletions = []
        o = 0
        for i in re.finditer("-+", segment.query_alignment):
            s = i.start() - o
            e = i.end() - o
            del_len = e - s
            del_sequence = segment.germline_alignment[s:e]
            if del_len % 3 == 0 or del_len > 3:
                deletions.append(
                    _annotate_deletion(s + segment.query_start, del_len, del_sequence)
                )
            else:
                deletions.append(
                    _annotate_deletion(
                        s + segment.query_start,
                        del_len,
                        del_sequence,
                        fixed=True,
                    )
                )
                _fix_frameshift_deletion(antibody, segment, s, e)
                o += del_len
        return deletions if deletions else []
    except Exception:
        segment.exception("FIND DELETIONS ERROR", traceback.format_exc())


def _annotate_deletion(start, length, sequence, fixed=False):
    """
    Annotates codon-length (non-frameshift) deletions.

    Parameters
    ----------
    start : int
        The starting postion of the deletion.

    length : int
        The length of the deletion.

    sequence : str
        The deleted sequence.

    fixed : bool
        Whether the deletion was fixed (reverted to germline).

    Returns
    -------
    Deletion
        A ``Deletion`` object containing the deletion start position, the deletion length,
        and the deleted sequence.

    """
    in_frame = "yes" if length % 3 == 0 else "no"
    return Deletion(
        {
            "pos": start,
            "len": length,
            "seq": sequence,
            "in frame": in_frame,
            "fixed": fixed,
        }
    )


def _fix_frameshift_deletion(
    antibody: "Antibody", segment: "Germline", s: int, e: int
) -> None:
    """
    Fixes (removes) frameshift deletions.

    Deletion fixing entails several steps:
        1. Remove the frameshift indel from the alignment (query, germline, and midline).
        2. Remove the frameshift indel from the oriented input.
        3. Adjust the alignment end position to account for the reverted indel.

    Parameters
    ----------
    antibody : Antibody
        The Antibody object being annotated.

    segment : Germline
        The Germline object to be annotated.

    s : int
        The starting postion of the deletion.

    e : int
        The ending position of the deletion.

    Returns
    -------
    None
        The input ``Germline`` and ``Antibody`` objects are modified in place.

    """
    # remove the frameshift indel from the alignment
    segment.query_alignment = (
        segment.query_alignment[:s]
        + segment.germline_alignment[s:e]
        + segment.query_alignment[e:]
    )
    segment.alignment_midline = (
        segment.alignment_midline[:s]
        + "".join(["|"] * (e - s))
        + segment.alignment_midline[e:]
    )

    # fix the frameshift deletion in the oriented input
    oi_firsthalf = antibody.oriented_input.sequence[: s + segment.query_start]
    middle = segment.germline_alignment[s:e]
    oi_secondhalf = antibody.oriented_input.sequence[s + segment.query_start :]
    antibody.oriented_input.sequence = oi_firsthalf + middle + oi_secondhalf

    # adjust the alignment end position
    segment.query_end += e - s
