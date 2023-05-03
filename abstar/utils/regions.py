#!/usr/bin/python
# filename: regions.py

#
# Copyright (c) 2015 Bryan Briney
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


from __future__ import absolute_import, division, print_function, unicode_literals

import os
import re
import sys
import traceback

from Bio.Seq import Seq

from abutils.utils import log


# Both IMGT_REGION_END_POSITIONS_AA and IMGT_REGION_END_POSITIONS_NT use the actual
# IMGT end positions, thus they're not suitable for slicing (because Python's uses
# zero-based numbering). To use these end points in a slice, you need to subtract 1
# from the start position (since Python's slicing is exclusive, the actual IMGT
# numbering works fine for the end position of the slice).
IMGT_REGION_END_POSITIONS_AA = {
    "FR1": 26,
    "CDR1": 38,
    "FR2": 55,
    "CDR2": 65,
    "FR3": 104,
    "CDR3": 117,
    "FR4": 129,
}

IMGT_REGION_END_POSITIONS_NT = {
    "FR1": 78,
    "CDR1": 114,
    "FR2": 165,
    "CDR2": 195,
    "FR3": 312,
    "CDR3": 351,
    "FR4": 387,
}

IMGT_REGION_START_POSITIONS_AA = {
    "FR1": 1,
    "CDR1": 27,
    "FR2": 39,
    "CDR2": 56,
    "FR3": 66,
    "CDR3": 105,
    "FR4": 118,
}

IMGT_REGION_START_POSITIONS_NT = {
    "FR1": 1,
    "CDR1": 79,
    "FR2": 115,
    "CDR2": 166,
    "FR3": 196,
    "CDR3": 313,
    "FR4": 352,
}


# def regions(antibody):
#     # V-gene regions
#     antibody.v.regions = VariableRegions(antibody)
#     # J-gene regions
#     antibody.j.regions = JoiningRegions(antibody)


def get_variable_regions(antibody):
    return VariableRegions(antibody)


def get_joining_regions(antibody):
    return JoiningRegions(antibody)


class BaseRegions(object):
    """Base class for Variable and Joining gene regions"""

    def __init__(self, antibody):
        super(BaseRegions, self).__init__()
        self.antibody = antibody

        self._raw_nt_positions = None
        self._nt_seqs = None
        self._aa_seqs = None

    @property
    def raw_nt_positions(self):
        """
        A `dict`, with region names as keys. Values are `list`s containing the start/end points
        (in raw numbering, corresponding to the oriented_input). If raw start/end points
        couldn't be identified, they are `None`. So if an entire region isn't present (either due
        to deletion or query truncation), the value will be `[None, None]`.

        Note that the end points are designed for slicing, so the end point is actually 1 position
        beyond the actual end point. This allows `ab.oriented_input[start:end]` to work. This also
        means that the end point of FR1 should be identical to the start point of CDR1.
        """
        if self._raw_nt_positions is None:
            pos = {}
            for region in self.region_names:
                start = self._raw_region_start_position_nt(region)
                end = self._raw_region_end_position_nt(region)
                pos[region] = [start, end]
            self._raw_nt_positions = pos
        return self._raw_nt_positions

    @property
    def nt_seqs(self):
        if self._nt_seqs is None:
            nt_seqs = {}
            for region in self.region_names:
                start, end = self.raw_nt_positions[region]
                if all([start is None, end is None]):
                    nt_seqs[region] = ""
                    continue
                if start is not None and end is None:
                    err = "{} contains a start position but not an end position".format(
                        region
                    )
                    self.antibody.log("REGION NT SEQUENCE IDENTIFICATION ERROR:", err)
                elif start > end:
                    err = "{} contains a start position {} that is greater than the end position{}".format(
                        region, start, end
                    )
                    self.antibody.log("REGION NT SEQUENCE IDENTIFICATION ERROR:", err)
                else:
                    # if start is None, that means the region begins before the query sequence does. But we
                    # don't necessarily want to start the region sequence at the beginning of the query,
                    # we want to start the region sequence at the point that the query aligns to the germline.
                    if start is None:
                        start = self.segment.query_start
                    nt_seq = self.antibody.oriented_input[start:end]
                    nt_seqs[region] = nt_seq
            self._nt_seqs = nt_seqs
        return self._nt_seqs

    @property
    def aa_seqs(self):
        if self._aa_seqs is None:
            to_translate = {}
            aa_seqs = {}
            first_region = self._get_first_region()
            for i, region in enumerate(self.region_names):
                # in the case of region-spanning indels, we need to process multiple regions
                # simultaneously. The following prevents us from reprocessing a region
                # that contains the second half of a region-spanning deletion.
                if region in list(to_translate.keys()):
                    continue
                # if we're looking at the first region in the query sequence, we need to make
                # sure that we're in the correct reading frame.
                if region == first_region and self.segment.gene_type == "V":
                    region_nt = self.nt_seqs[region][self.antibody.v_rf_offset :]
                else:
                    region_nt = self.nt_seqs[region]
                # check to make sure the region_nt sequence is a multiple of three. If it's not,
                # there are two possibilities:
                # 1) a region-spanning deletion that deletes, for example, 1nt from this region and 2nt
                #    from the following region. In this case, both this region and the next would be out
                #    of frame and we need to fix that.
                # 2) a frameshift indel within this region that's longer than a single codon (since we
                #    don't automatically fix those).
                if len(region_nt) % 3 != 0:
                    # if we're looking at the last region, we obviously can't look to the next region
                    # for the other portion of a region-spanning deletion, so we'll just
                    # try to translate the region as is (as it's either truncated at the 3' end or
                    # we have a frameshift indel that's longer than a single codon which we don't want
                    # to automatically fix).
                    if i == len(self.region_names) - 1:
                        to_translate[region] = self.nt_seqs[region]
                    else:
                        # find the next region that has sequence data. A little convoluted, but this
                        # accounts for the rare case of a large region-spanning deletion that actually
                        # deletes an entire downstream region (for example, a deletion that starts at the end
                        # of FR2, deletes the entire CDR2, and ends near the start of FR3).
                        next_regions = self.region_names[i + 1 :]
                        reg_iter = iter(next_regions)
                        region2 = next(reg_iter)
                        while len(self.nt_seqs[region2].strip("-")) == 0:
                            to_translate[region2] = ""
                            region2 = next(reg_iter)
                        # if the current region and the next region with sequence data combine to form an in-frame
                        # sequence, we can fix the region-spanning deletion so that both regions are translated
                        # correctly.
                        if len(region_nt + self.nt_seqs[region2]) % 3 == 0:
                            region_seq, region2_seq = self._fix_region_spanning_indel(
                                region_nt, self.nt_seqs[region2]
                            )
                            to_translate[region] = region_seq
                            to_translate[region2] = region2_seq
                        # if the combined regions are still out of frame, the current region likely contains a
                        # frameshift indel that shouldn't be automatically fixed (it's longer than a single codon).
                        else:
                            to_translate[region] = region_nt
                else:
                    to_translate[region] = region_nt
            # translate the sequences for each region
            for region, nt_seq in to_translate.items():
                if len(nt_seq) % 3:
                    nt_seq = nt_seq[: -(len(nt_seq) % 3)]
                aa_seq = str(Seq(nt_seq).translate())
                aa_seqs[region] = aa_seq
            self._aa_seqs = aa_seqs
        return self._aa_seqs

    def _raw_region_start_position_nt(self, region):
        # get the raw start position of the region
        imgt_start = IMGT_REGION_START_POSITIONS_NT[region]
        # if there isn't a direct equivalent to the IMGT region start position,
        # that could be because the qyuery sequence is truncated (and the start of the region
        # isn't present in the query sequence) or it could be because a deletion
        # removed the region start position.
        # First, let's check to see if there's any query sequence prior to the start of
        # the region.
        if imgt_start > max(
            [p for p in self.segment._imgt_position_from_raw.values() if p is not None]
        ):
            return None
        # If there's sequence prior to the start of the region but the region start position
        # isn't present, we must have a deletion that removed the region start position. In
        # this case, we just scan forward in the sequence to find the earliest non-deleted
        # position in the region and use that as the region start position
        while self.segment.get_raw_position_from_imgt(imgt_start) is None:
            # we don't want to keep going past the end of the region so only keep looking
            # for the end of the deletion until we reach the start of the next region
            if imgt_start >= IMGT_REGION_END_POSITIONS_NT[region]:
                break
            imgt_start += 1
        return self.segment.get_raw_position_from_imgt(imgt_start)

    def _raw_region_end_position_nt(self, region):
        imgt_end = IMGT_REGION_END_POSITIONS_NT[region]
        # if the end position of the region just happens to be part of a deletion,
        # or if the query sequence doesn't contain the region, then there isn't a
        # corresponding raw sequence position.
        # First, check to see if the query sequence begins after the end of the region
        if imgt_end < min(
            [p for p in self.segment._imgt_position_from_raw.values() if p is not None]
        ):
            return None
        # Now we need to try to find the start of a deletion that includes the region
        # end position
        while self.segment.get_raw_position_from_imgt(imgt_end) is None:
            imgt_end -= 1
            # we don't want to move all the way back into a different region,
            # so only keep looking for the start of the deletion until we
            # reach the start of the region
            if imgt_end < IMGT_REGION_START_POSITIONS_NT[region]:
                break
        return self.segment.get_raw_position_from_imgt(imgt_end) + 1

    def _get_first_region(self):
        for region in self.region_names:
            if self.raw_nt_positions[region][1] is not None:
                return region

    def _fix_region_spanning_indel(self, r1_nt, r2_nt):
        # if the first region contains the majority of the split codon,
        # we add the complete codon to the first region
        if len(r1_nt) % 3 == 2:
            return r1_nt + r2_nt[0], r2_nt[1:]
        # if the second region contains the majority of the split codon,
        # we add the complete codon to the second region
        elif len(r1_nt) % 3 == 1:
            return r1_nt[:-1], r1_nt[-1] + r2_nt
        return r1_nt, r2_nt


class VariableRegions(BaseRegions):
    """docstring for VariableRegions"""

    def __init__(self, antibody):
        super(VariableRegions, self).__init__(antibody)
        self.segment = antibody.v
        self.region_names = ["FR1", "CDR1", "FR2", "CDR2", "FR3"]


class JoiningRegions(BaseRegions):
    """docstring for JoiningRegions"""

    def __init__(self, antibody):
        super(JoiningRegions, self).__init__(antibody)
        self.segment = antibody.j
        self.region_names = ["FR4"]

