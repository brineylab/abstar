#!/usr/bin/python
# filename: junction.py

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


import os
import traceback

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

from abtools import log
from abtools.alignment import global_alignment, local_alignment


def get_junction(vdj, germ=False):
    global logger
    logger = log.get_logger(__name__)
    try:
        return Junction(vdj, germ)
    except Exception, err:
        logger.debug('JUNCTION ERROR: {id}'.format(id=vdj.id))
        logger.debug(traceback.format_exc())


class Junction(object):
    """docstring for Junction"""
    def __init__(self, vdj, germ):
        self.is_germ = germ
        self.in_frame = True
        self._frame_offset = 0
        self.junction_nt_start_pos = self._find_junction_nt_start(vdj)
        self.junction_nt_end_pos = self._find_junction_nt_end(vdj)
        self.junction_nt = self._get_junction_nt_sequence(vdj)
        if len(self.junction_nt) % 3 > 0:
            self.in_frame = False
            self._frame_offset = 3 - len(self.junction_nt) % 3
            n_end = vdj.n2_end if vdj.d is not None else vdj.n_end
            n_end = n_end - self.junction_nt_start_pos
            ns = 'N' * self._frame_offset
            self.junction_nt = self.junction_nt[:n_end] + ns + self.junction_nt[n_end:]
        self.junction_aa = str(Seq(self.junction_nt.replace('-', ''), generic_dna).translate())
        self.cdr3_aa = self.junction_aa[1:-1]
        self.cdr3_nt = self.junction_nt[3:-3]
        if vdj.d:
            self.n1_nt = vdj.vdj_nt[vdj.n1_start:vdj.n1_end]
            self.d_nt = vdj.vdj_nt[vdj.d_start:vdj.d_end]
            self.n2_nt = vdj.vdj_nt[vdj.n2_start:vdj.n2_end]
            self.n_nt = None
            self.d_start_position_nt = self._get_d_start_position_nt(vdj)
            self.d_end_position_nt = self._get_d_end_position(vdj)
            self.d_dist_from_cdr3_start_nt = self.d_start_position_nt
            self.d_dist_from_cdr3_end_nt = self._get_d_dist_from_cdr3_end_nt(vdj)
        else:
            self.n1_nt = None
            self.d_nt = None
            self.n2_nt = None
            self.n_nt = vdj.vdj_nt[vdj.n_start:vdj.n_end]


    def _find_junction_nt_start(self, vdj):
        fr3 = vdj.v.regions.nt_seqs['FR3'][:-3]
        aln = local_alignment(fr3, vdj.vdj_nt)
        return aln.target_end + 1


    def _find_junction_nt_end(self, vdj):
        fr4 = vdj.j.regions.nt_seqs['FR4'][3:]
        aln = local_alignment(fr4, vdj.vdj_nt)
        return aln.target_begin


    def _get_junction_nt_sequence(self, vdj):
        start = self.junction_nt_start_pos
        end = self.junction_nt_end_pos
        if self.is_germ:
            return vdj.vdj_germ_nt[start:end]
        return vdj.vdj_nt[start:end]


    def _get_d_start_position_nt(self, vdj):
        a = local_alignment(self.d_nt, self.cdr3_nt,
                            gap_open_penalty=22, gap_extend_penalty=1)
        d_start = a.target_begin
        return d_start


    def _get_d_end_position(self, vdj):
        return self.d_start_position_nt + len(vdj.d.sequence)


    def _get_d_dist_from_cdr3_end_nt(self, vdj):
        return len(self.cdr3_nt) - (self.d_dist_from_cdr3_start_nt + len(self.d_nt))
