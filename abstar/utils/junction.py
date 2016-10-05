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
from abtools.utils.codons import codons


def get_junction(antibody):
    antibody.log('')
    antibody.log('JUNCTION')
    antibody.log('--------')
    global logger
    logger = log.get_logger(__name__)
    try:
        return Junction(vdj, germ)
    except Exception, err:
        logger.debug('JUNCTION ERROR: {id}'.format(id=vdj.id))
        logger.debug(traceback.format_exc())




class Junction(object):
    """docstring for Junction"""
    def __init__(self, antibody):
        super(Junction, self).__init__()
        self.in_frame = True
        # identify junction sequence
        self.junction_nt_start = self._find_junction_nt_start(antibody)
        self.junction_nt_end = self._find_junction_nt_end(antibody)
        self.junction_nt = antibody.oriented_input[self.junction_nt_start:self.junction_nt_end]
        antibody.log('JUNCTION NT:', self.junction_nt, len(self.junction_nt))
        if len(self.junction_nt) % 3 != 0:
            self.in_frame = False
            antibody.log('WARNING: Junction out of frame!')
        self.junction_aa = str(Seq(self.junction_nt.replace('-', ''), generic_dna).translate())
        antibody.log('JUNCTION AA:', self.junction_aa, len(self.junction_aa))
        # identify CDR3 sequence
        self.cdr3_aa = self.junction_aa[1:-1]
        self.cdr3_nt = self.junction_nt[3:-3]
        # parse N-addition region(s) and D-gene location (if heavy chain)
        if antibody.d is not None:
            n1_start = antibody.v.query_end + 1
            n1_end = antibody.d.query_start
            n2_start = antibody.d.query_end + 1
            n2_end = antibody.j.query_start
            self.v_nt = antibody.oriented_input[self.junction_nt_start:antibody.v.query_end + 1]
            self.n1_nt = antibody.oriented_input[n1_start:n1_end]
            self.d_nt = antibody.oriented_input[n1_end:n2_start]
            self.n2_nt = antibody.oriented_input[n2_start:n2_end]
            self.j_nt = antibody.oriented_input[n2_end:self.junction_nt_end]
            self.d_dist_from_cdr3_start_nt = antibody.d.query_start - self.junction_nt_start - 3
            self.d_dist_from_cdr3_end_nt = self.junction_nt_end - (antibody.d.query_end + 1) - 3
            self.n_nt = None
            antibody.log('V_NT:', self.v_nt)
            antibody.log('N1_NT:', self.n1_nt)
            antibody.log('D_NT:', self.d_nt)
            antibody.log('N2_NT:', self.n2_nt)
            antibody.log('J_NT:', self.j_nt)
            antibody.log('D-GENE DISTANCE FROM CDR3 START:', self.d_dist_from_cdr3_start_nt)
            antibody.log('D-GENE DISTANCE FROM CDR3 END:', self.d_dist_from_cdr3_end_nt)
        else:
            n_start = antibody.v.query_end + 1
            n_end = antibody.j.query_start
            self.v_nt = antibody.oriented_input[self.junction_nt_start:antibody.v.query_end + 1]
            self.n_nt = antibody.oriented_input[n_start:n_end]
            self.j_nt = antibody.oriented_input[n_end:self.junction_nt_end]
            self.n1_nt = None
            self.d_nt = None
            self.n1_nt = None
            self.d_dist_from_cdr3_start_nt = None
            self.d_dist_from_cdr3_end_nt = None
            antibody.log('V_NT:', self.v_nt)
            antibody.log('N_NT:', self.n_nt)
            antibody.log('J_NT:', self.j_nt)


    def _find_junction_nt_start(self, antibody):
        start_codons = ['TGT', 'TGC']
        # position 310 is the start of codon 104, the conserved 2nd-Cys
        junc_start = antibody.v.get_raw_position_from_imgt(310)
        junc_start_codon = antibody.oriented_input[junc_start:junc_start + 3]
        antibody.log('JUNC START:', junc_start_codon, codons[junc_start_codon], junc_start)
        # if the identified junction start isn't normal, use the fallback (and more
        # computationally intensive) method for finding the junction start.
        if junc_start_codon not in start_codons:
            antibody.log('WARNING: Using fallback method to find junction start!')
            return self._fallback_find_junc_nt_start(antibody)
        return junc_start


    def _fallback_find_junc_nt_start(self, antibody):
        # get the FR3 nt sequence of the IMGT gapped germline
        germ_fr3_sequence = antibody.v.imgt_germline.gapped_nt_sequence[196:311].replace('.', '')
        # find the start of the junction (immediately after the end of FR3)
        aln = local_alignment(antibody.oriented_input, germ_fr3_sequence)
        fr3_end = aln.query_end + (len(germ_fr3_sequence) - aln.target_end)
        junc_start_codon = antibody.oriented_input[fr3_end:fr3_end + 3]
        antibody.log('JUNC START:', junc_start_codon, codons[junc_start_codon], fr3_end)
        return fr3_end


    def _find_junction_nt_end(self, antibody):
        if antibody.chain == 'heavy':
            end_codons = ['TGG']
        else:
            end_codons = ['TTT', 'TTC']
        # when calculating the location of the conserved W/F, need to compensate
        # for sequences that don't contain the full J-gene
        joffset = len(antibody.j.raw_germline) - (antibody.j.germline_end + 1)
        # find the end of the J gene, then back up 12 codons to get to
        # the conserved W/F (start of FR4)
        junc_end = antibody.j.query_end - (33 - joffset) + 1
        junc_end_codon = antibody.oriented_input[junc_end:junc_end + 3]
        antibody.log('JUNC END:', junc_end_codon, codons[junc_end_codon], junc_end)
        # if the identified junction end isn't normal, use the fallback (and more
        # computationally intensive) method for finding the junction end.
        if junc_end_codon not in end_codons:
            antibody.log('WARNING: Using fallback method to find junction end!')
            return self._fallback_find_junc_nt_end(antibody)
        return junc_end + 3


    def _fallback_find_junc_nt_end(self, antibody):
        # need to find the start of FR4 in the IMGT germline sequence
        end_res = 'W' if antibody.chain == 'heavy' else 'F'
        for i, res in enumerate(antibody.j.imgt_germline.ungapped_aa_sequence):
            if res == end_res and end_res not in antibody.j.imgt_germline.ungapped_aa_sequence[i + 1:]:
                fr4_nt_start_pos = (antibody.j.imgt_germline.coding_start - 1) + (i * 3)
                break
        germ_fr4_sequence = antibody.j.imgt_germline.gapped_nt_sequence[fr4_nt_start_pos:]
        # find the end of the junction (end of the first codon of FR4)
        aln = local_alignment(antibody.oriented_input, germ_fr4_sequence)
        fr4_start = aln.query_begin - aln.target_begin
        junc_end_codon = antibody.oriented_input[fr4_start:fr4_start + 3]
        antibody.log('JUNC END:', junc_end_codon, codons[junc_end_codon], fr4_start)
        return fr4_start + 3





# def get_junction(vdj, germ=False):
#     global logger
#     logger = log.get_logger(__name__)
#     try:
#         return Junction(vdj, germ)
#     except Exception, err:
#         logger.debug('JUNCTION ERROR: {id}'.format(id=vdj.id))
#         logger.debug(traceback.format_exc())


# class Junction(object):
#     """docstring for Junction"""
#     def __init__(self, vdj, germ):
#         self.is_germ = germ
#         self.in_frame = True
#         self._frame_offset = 0
#         self.junction_nt_start_pos = self._find_junction_nt_start(vdj)
#         self.junction_nt_end_pos = self._find_junction_nt_end(vdj)
#         if any([self.junction_nt_start_pos is None, self.junction_nt_end_pos is None]):
#             logger.debug('JUNCTION IDENTIFICATION ERROR: {}'.format(vdj.id))
#             logger.debug('RAW INPUT: {}'.format(vdj.raw_input))
#             return None
#         self.junction_nt = self._get_junction_nt_sequence(vdj)
#         if len(self.junction_nt) % 3 > 0:
#             logger.debug('OUT OF FRAME JUNCTION: {}'.format(vdj.id))
#             logger.debug('JUNCTION SEQUENCE: {}'.format(self.junction_nt))
#             logger.debug('VDJ_NT: {}'.format(vdj.vdj_nt))
#             logger.debug('FR3 SEQUENCE: {}'.format(vdj.v.regions.nt_seqs['FR3']))
#             logger.debug('FR4 SEQUENCE: {}'.format(vdj.j.regions.nt_seqs['FR4']))
#             logger.debug('JUNCTION START: {}'.format(self.junction_nt_start_pos))
#             logger.debug('JUNCTION END: {}'.format(self.junction_nt_end_pos))
#             self.in_frame = False
#             self._frame_offset = 3 - len(self.junction_nt) % 3
#             n_end = vdj.n2_end if vdj.d is not None else vdj.n_end
#             n_end = n_end - self.junction_nt_start_pos
#             ns = 'N' * self._frame_offset
#             self.junction_nt = self.junction_nt[:n_end] + ns + self.junction_nt[n_end:]
#         self.junction_aa = str(Seq(self.junction_nt.replace('-', ''), generic_dna).translate())
#         self.cdr3_aa = self.junction_aa[1:-1]
#         self.cdr3_nt = self.junction_nt[3:-3]
#         if vdj.d:
#             self.n1_nt = vdj.vdj_nt[vdj.n1_start:vdj.n1_end]
#             self.d_nt = vdj.vdj_nt[vdj.d_start:vdj.d_end]
#             self.n2_nt = vdj.vdj_nt[vdj.n2_start:vdj.n2_end]
#             self.n_nt = None
#             self.d_start_position_nt = self._get_d_start_position_nt(vdj)
#             self.d_end_position_nt = self._get_d_end_position(vdj)
#             self.d_dist_from_cdr3_start_nt = self.d_start_position_nt
#             self.d_dist_from_cdr3_end_nt = self._get_d_dist_from_cdr3_end_nt(vdj)
#         else:
#             self.n1_nt = None
#             self.d_nt = None
#             self.n2_nt = None
#             self.n_nt = vdj.vdj_nt[vdj.n_start:vdj.n_end]


#     def _find_junction_nt_start(self, vdj):
#         fr3 = vdj.v.regions.nt_seqs['FR3'][:-3]
#         aln = local_alignment(fr3, vdj.vdj_nt)
#         if aln:
#             return aln.target_end + 1


#     def _find_junction_nt_end(self, vdj):
#         fr4 = vdj.j.regions.nt_seqs['FR4'][3:]
#         aln = local_alignment(fr4, vdj.vdj_nt)
#         if aln:
#             return aln.target_begin


#     def _get_junction_nt_sequence(self, vdj):
#         start = self.junction_nt_start_pos
#         end = self.junction_nt_end_pos
#         if self.is_germ:
#             return vdj.vdj_germ_nt[start:end]
#         return vdj.vdj_nt[start:end]


#     def _get_d_start_position_nt(self, vdj):
#         a = local_alignment(self.d_nt, self.cdr3_nt,
#                             gap_open_penalty=22, gap_extend_penalty=1)
#         d_start = a.target_begin
#         return d_start


#     def _get_d_end_position(self, vdj):
#         return self.d_start_position_nt + len(vdj.d.sequence)


#     def _get_d_dist_from_cdr3_end_nt(self, vdj):
#         return len(self.cdr3_nt) - (self.d_dist_from_cdr3_start_nt + len(self.d_nt))
