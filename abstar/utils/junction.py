#!/usr/bin/python
# filename: junction.py

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


import math
import os
import traceback

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

from abtools import log
from abtools.alignment import global_alignment, local_alignment
from abtools.utils.codons import codon_lookup as codons


def get_junction(antibody):
    antibody.log('')
    antibody.log('JUNCTION')
    antibody.log('--------')
    try:
        return Junction(antibody)
    except Exception:
        antibody.exception('JUNCTION IDENTIFICATION ERROR', traceback.format_exc())
        raise


class Junction(object):
    """docstring for Junction"""
    def __init__(self, antibody):
        super(Junction, self).__init__()
        self.in_frame = True
        self.fallback_5prime = False
        self.fallback_3prime = False

        # identify junction sequence
        self.junction_nt_start = self._find_junction_nt_start(antibody)
        self.junction_nt_end = self._find_junction_nt_end(antibody)
        self.junction_nt = antibody.oriented_input[self.junction_nt_start:self.junction_nt_end]
        antibody.log('JUNCTION NT:', self.junction_nt, len(self.junction_nt))

        # check junction frame and retry junction identification with the slower
        # (but potentially more accurate) alignment method if it isn't in frame
        if len(self.junction_nt) % 3 != 0:
            antibody.log('WARNING: Junction out of frame!')
            if any([self.fallback_5prime is False, self.fallback_3prime is False]):
                # only try the 5' alignment if it wasn't already used
                if not self.fallback_5prime:
                    antibody.log('Attempting to identify junction start with fallback alignment method')
                    self.junction_nt_start = self._fallback_find_junc_nt_start(antibody)
                # only try the 3' alignment if it wasn't already used
                if not self.fallback_3prime:
                    antibody.log('Attempting to identify junction end with fallback alignment method')
                    self.junction_nt_end = self._fallback_find_junc_nt_end(antibody)
                self.junction_nt = antibody.oriented_input[self.junction_nt_start:self.junction_nt_end]
                antibody.log('JUNCTION NT:', self.junction_nt, len(self.junction_nt))
            if len(self.junction_nt) % 3 != 0:
                self.in_frame = False
        self.junction_aa = str(Seq(self.junction_nt.replace('-', ''), generic_dna).translate())
        antibody.log('JUNCTION AA:', self.junction_aa, len(self.junction_aa))

        # identify CDR3 sequence
        self.cdr3_aa = self.junction_aa[1:-1]
        self.cdr3_nt = self.junction_nt[3:-3]

        # parse N-addition regions
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

        # calculate IMGT numbering for the junction and CDR3
        self.cdr3_imgt_nt_numbering = self._calculate_cdr3_imgt_nt_numbering()
        self.cdr3_imgt_aa_numbering = self._calculate_cdr3_imgt_aa_numbering()
        self.junction_imgt_nt_numbering = ['310', '311', '312'] + self.cdr3_imgt_nt_numbering + ['352', '353', '354']
        self.junction_imgt_aa_numbering = ['104'] + self.cdr3_imgt_aa_numbering + ['118']

        # if the J-gene is long enough to extend into the variable-length portion of the
        # junction numbering, need to adjust the IMGT numbering for the J-gene
        if min(antibody.j.imgt_aa_positions) < 112 or min(antibody.j.imgt_nt_positions) < 334:
            antibody.log('ADJUSTING J-GENE IMGT NUMBERING')
            antibody.log('OLD AA:', ', '.join([str(p) for p in antibody.j.imgt_aa_positions]))
            antibody.log('OLD NT:', ', '.join([str(p) for p in antibody.j.imgt_nt_positions]))
            self._adjust_jgene_imgt_numbering(antibody)
            antibody.log('NEW AA:', ', '.join([str(p) for p in antibody.j.imgt_aa_positions]))
            antibody.log('NEW NT:', ', '.join([str(p) for p in antibody.j.imgt_nt_positions]))
        else:
            antibody.j.correct_imgt_aa_position_from_imgt = {p: p for p in antibody.j.imgt_aa_positions}
            antibody.j.correct_imgt_nt_position_from_imgt = {p: p for p in antibody.j.imgt_nt_positions}


    def _find_junction_nt_start(self, antibody):
        start_codons = ['TGT', 'TGC']
        # position 310 is the start of codon 104, the conserved 2nd-Cys
        junc_start = antibody.v.get_raw_position_from_imgt(310)
        if junc_start is None:
            log_str = 'WARNING: The full 2nd-Cys codon does not appear to be present in the V-gene alignment, '
            log_str += 'likely due to extensive trimming during recombinatinon. '
            log_str += 'Using fallback alignment method to find junction start.'
            antibody.log(log_str)
            return self._fallback_find_junc_nt_start(antibody)
        junc_start_codon = antibody.oriented_input[junc_start:junc_start + 3]
        antibody.log('JUNC START:', junc_start_codon, codons[junc_start_codon], junc_start)
        # if the identified junction start isn't normal, use the fallback (and slower)
        # alignment method for finding the junction start.
        if junc_start_codon not in start_codons:
            antibody.log('WARNING: Using fallback method to find junction start!')
            return self._fallback_find_junc_nt_start(antibody)
        return junc_start


    def _fallback_find_junc_nt_start(self, antibody):
        self.fallback_5prime = True
        # get the FR3 nt sequence of the IMGT gapped germline
        germ_fr3_sequence = antibody.v.imgt_germline.gapped_nt_sequence[196:309].replace('.', '')
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

        # if the identified junction end isn't normal, use the fallback (and slower)
        # alignment method for finding the junction end.
        if junc_end_codon not in end_codons:
            antibody.log('WARNING: Did not identify conserved position 118. Using fallback method to find junction end.')
            return self._fallback_find_junc_nt_end(antibody)
        return junc_end + 3


    def _fallback_find_junc_nt_end(self, antibody):
        self.fallback_3prime = True

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


    def _calculate_cdr3_imgt_nt_numbering(self):
        # if the length of the CDR3 is precisely 13, we don't need to adjust
        # the numbering at all (because the default IMGT numbering allows for exactly
        # 13 AA in the CDR3).
        if len(self.cdr3_aa) == 13:
            return list(range(313, 352))

        # If the CDR3 is shorter than 13 AA, we need to remove position numbers,
        # starting with positions to the left of position 333 (codon 111).
        elif len(self.cdr3_aa) < 13:
            gap = 13 - len(self.cdr3_aa)
            ltrim = int(math.ceil(float(gap) / 2)) * 3
            rtrim = int(math.floor(float(gap) / 2)) * 3
            lnums = list(range(313, 334 - ltrim))
            rnums = list(range(334 + rtrim, 352))
            return lnums + rnums

        # If the CDR3 is longer than 13 AA, we need to add position numbers (with decimal),
        # starting with positions to the right (position 334, or codon 112) and alternating
        # between left (position 333, codon 111) and right.
        # Nucleotides added to 334 will decrement from left to right, nucleotides added to 333
        # will increment from left to right. For example:
        # 331, 332, 333, 331.1, 331.2, 333.3, 334.3, 334.2, 334.1, 334, 335, ...
        # We add the nucleotide numbers in codon-length blocks, so the first codon will be
        # nucleotides 334.3, 334.2, and 334.1.
        elif len(self.cdr3_aa) > 13:
            bonus = len(self.cdr3_aa) - 13
            lbonus = int(math.floor(float(bonus) / 2)) * 3
            rbonus = int(math.ceil(float(bonus) / 2)) * 3
            lbonus_nums = ['333.{}'.format(i) for i in range(1, lbonus + 1)]
            rbonus_nums = ['334.{}'.format(i) for i in reversed(range(1, rbonus + 1))]
            return list(range(313, 334)) + lbonus_nums + rbonus_nums + list(range(334, 352))


    def _calculate_cdr3_imgt_aa_numbering(self):
        # Edge case for extremely short (likely non-productive) CDR3s
        if len(self.cdr3_aa) <= 2:
            return [105, 117][:len(self.cdr3_aa)]

        # if the length of the CDR3 is precisely 13, we don't need to adjust
        # the numbering at all (because the default IMGT numbering allows for exactly
        # 13 AA in the CDR3).
        if len(self.cdr3_aa) == 13:
            return list(range(105, 118))

        # If the CDR3 is shorter than 13 AA, we need to remove position numbers,
        # starting with positions to the left of codon 111 and then alternating left
        # and right.
        elif len(self.cdr3_aa) < 13:
            gap = 13 - len(self.cdr3_aa)
            ltrim = int(math.ceil(float(gap) / 2))
            rtrim = int(math.floor(float(gap) / 2))
            lnums = list(range(105, 112 - ltrim))
            rnums = list(range(112 + rtrim, 118))
            return lnums + rnums

        # If the CDR3 is longer than 13 AA, we need to add position numbers (with decimal),
        # starting with positions to the right (codon 112) and alternating
        # between left (codon 111) and right.
        # AAs added to 112 will decrement from left to right, AAs added to 111
        # will increment from left to right. For example:
        # 109, 110, 111, 111.1, 111.2, 112.3, 112.2, 112.1, 112, 113, ...
        elif len(self.cdr3_aa) > 13:
            bonus = len(self.cdr3_aa) - 13
            lbonus = int(math.floor(float(bonus) / 2))
            rbonus = int(math.ceil(float(bonus) / 2))
            lbonus_nums = ['111.{}'.format(i) for i in range(1, lbonus + 1)]
            rbonus_nums = ['112.{}'.format(i) for i in reversed(range(1, rbonus + 1))]
            return list(range(105, 112)) + lbonus_nums + rbonus_nums + list(range(112, 118))


    def _adjust_jgene_imgt_numbering(self, antibody):
        # adjust AA numbering
        aa_positions_to_replace = [p for p in antibody.j.imgt_aa_positions if p < 118]
        dot_positions = self.cdr3_imgt_aa_numbering[-len(aa_positions_to_replace):]
        new_positions = dot_positions + [p for p in antibody.j.imgt_aa_positions if p >= 118]
        antibody.j.correct_imgt_aa_position_from_imgt = {o: n for o, n in zip(antibody.j.imgt_aa_positions, new_positions)}
        antibody.j.imgt_aa_positions = new_positions
        # adjust NT numbering
        nt_positions_to_replace = [p for p in antibody.j.imgt_nt_positions if p < 352]
        dot_positions = self.cdr3_imgt_nt_numbering[-len(nt_positions_to_replace):]
        new_positions = dot_positions + [p for p in antibody.j.imgt_nt_positions if p >= 352]
        antibody.j.correct_imgt_nt_position_from_imgt = {o: n for o, n in zip(antibody.j.imgt_nt_positions, new_positions)}
        antibody.j.imgt_nt_positions = new_positions
