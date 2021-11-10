#!/usr/bin/env python
# filename: antibody.py

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


from __future__ import absolute_import, division, print_function, unicode_literals

import traceback

from Bio.Seq import Seq

from abutils.utils.alignment import global_alignment, local_alignment

from .germline import get_imgt_germlines
from ..utils import isotype, junction, mutations, productivity, regions
from ..utils.mixins import LoggingMixin


class Antibody(LoggingMixin):
    """
    docstring for Antibody
    """
    def __init__(self, vdj, germ_db):
        super(Antibody, self).__init__()
        LoggingMixin.__init__(self)
        self.id = vdj.sequence.id
        self.raw_input = vdj.sequence
        self.oriented_input = vdj.oriented
        self.v = vdj.v
        self.j = vdj.j
        self.d = vdj.d
        self.chain = vdj.v.chain
        self.germ_db = germ_db
        self.species = self.v.species if self.v is not None else germ_db
        # initialize the log
        self.initialize_log()
        # property vars
        self._strand = None
        self._imgt_position_from_aligned = {}
        self._aligned_position_from_imgt = {}


    @property
    def strand(self):
        if self._strand is None:
            if self.v is None:
                self._strand = None
            else:
                self._strand = self.v.strand
        return self._strand


    @property
    def has_exceptions(self):
        if len(self.exceptions) > 0:
            return True
        segs = [s for s in [self.v, self.d, self.j] if s is not None]
        if any([len(seg.exceptions) > 0 for seg in segs]):
            return True
        return False


    def initialize_log(self):
        log = []
        log.append('=' * len(self.id))
        log.append(self.id)
        log.append('=' * len(self.id))
        log.append('')
        log.append(self.raw_input.fasta)
        log.append('')
        log.append('RAW INPUT: {}'.format(self.raw_input.sequence))
        log.append('ORIENTED INPUT: {}'.format(self.oriented_input.sequence))
        log.append('CHAIN: {}'.format(self.chain))
        log.append('')
        if self.v is not None:
            log.append('V-GENE: {}'.format(self.v.full))
        if self.d is not None:
            log.append('D-GENE: {}'.format(self.d.full))
        if self.j is not None:
            log.append('J-GENE: {}'.format(self.j.full))
        self._log = log



    def annotate(self, uid):
        '''
        Realigns the V(D)J germline genes using optimized SSW alignment
        parameters, followed by more detailed annotation (junction
        identification, ambig correction, parsing of mutations and
        indels, etc).
        '''
        try:
            self._parse_uid(uid)
            self._realign_germlines()
            self._get_junction()
            self._assemble_vdj_sequence()
            self._identify_regions()
            self._mutations()
            self._isotype()
            self._productivity()
        except:
            self.exception('ANNOTATION', traceback.format_exc())


# ---------------------
# GERMLINE REALIGNMENT
# ---------------------

    def _realign_germlines(self):
        # realign V and compute gapped IMGT alignment for the V gene
        self.log('')
        self.log('V-GENE REALIGNMENT')
        self.log('------------------')
        self.v.realign_germline(self)
        self.log('RAW QUERY SEQUENCE:', self.v.raw_query)
        self.log('RAW QUERY LENGTH:', len(self.v.raw_query))
        self.log('RAW GERMLNE LENGTH:', len(self.v.raw_germline))
        self.log('QUERY START:', self.v.query_start)
        self.log('QUERY END:', self.v.query_end)
        self.log('GERMLINE START:', self.v.germline_start)
        self.log('GERMLINE END:', self.v.germline_end)
        self.log('')
        self.log('   QUERY:', self.v.query_alignment)
        self.log('         ', self.v.alignment_midline)
        self.log('GERMLINE:', self.v.germline_alignment)
        self.log('')
        self.log('V-GENE GAPPED IMGT ALIGNMENT')
        self.log('----------------------------')
        self.v.gapped_imgt_realignment()
        self.log('   QUERY:', self.v.imgt_gapped_alignment.aligned_query)
        self.log('         ', self.v.imgt_gapped_alignment.alignment_midline)
        self.log('GERMLINE:', self.v.imgt_gapped_alignment.aligned_target)
        self.log('')
        self.log('IMGT ALIGNMENT START:', self.v.imgt_gapped_alignment.target_begin)
        self.log('IMGT GERMLINE CODING START:', self.v.imgt_germline.coding_start)
        self.log('V-GENE ALIGNMENT READING FRAME:', self.v.alignment_reading_frame)
        self.log('V-GENE CODING REGION:', self.v.coding_region)
        self.log('V-GENE AA SEQUENCE:', self.v.aa_sequence)

        # realign J
        jstart = self.v.query_end + 1
        self.log('')
        self.log('J-GENE REALIGNMENT')
        self.log('------------------')
        self.j.realign_germline(self, query_start=jstart)
        self.log('RAW QUERY SEQUENCE:', self.j.raw_query)
        self.log('RAW QUERY LENGTH:', len(self.j.raw_query))
        self.log('RAW GERMLNE LENGTH:', len(self.j.raw_germline))
        self.log('QUERY START:', self.j.query_start)
        self.log('QUERY END:', self.j.query_end)
        self.log('GERMLINE START:', self.j.germline_start)
        self.log('GERMLINE END:', self.j.germline_end)
        self.log('')
        self.log('   QUERY:', self.j.query_alignment)
        self.log('         ', self.j.alignment_midline)
        self.log('GERMLINE:', self.j.germline_alignment)
        self.log('')
        self.log('J-GENE GAPPED IMGT ALIGNMENT')
        self.log('----------------------------')
        self.j.gapped_imgt_realignment()
        self.log('   QUERY:', self.j.imgt_gapped_alignment.aligned_query)
        self.log('         ', self.j.imgt_gapped_alignment.alignment_midline)
        self.log('GERMLINE:', self.j.imgt_gapped_alignment.aligned_target)
        self.log('')
        self.log('IMGT ALIGNMENT START:', self.j.imgt_gapped_alignment.target_begin)
        self.log('IMGT GERMLINE CODING START:', self.j.imgt_germline.coding_start)
        self.log('J-GENE ALIGNMENT READING FRAME:', self.j.alignment_reading_frame)
        self.log('J-GENE CODING REGION:', self.j.coding_region)
        self.log('J-GENE AA SEQUENCE:', self.j.aa_sequence)

        # realign D (if needed)
        if all([self.chain in ['heavy', 'beta', 'delta'], self.d is not None]):
        # if self.d is not None:
        # if any([self.chain.lower() in ['heavy', 'beta', 'delta'], self.d is not None]):
            try:
                dstart = self.v.query_end + 1
                dend = self.j.query_start
                self.log('')
                self.log('D-GENE REALIGNMENT')
                self.log('------------------')
                self.d.realign_germline(self, query_start=dstart, query_end=dend)
                self.log('RAW QUERY LENGTH:', len(self.d.raw_query))
                self.log('RAW GERMLNE LENGTH:', len(self.d.raw_germline))
                self.log('QUERY START:', self.d.query_start)
                self.log('QUERY END:', self.d.query_end)
                self.log('GERMLINE START:', self.d.germline_start)
                self.log('GERMLINE END:', self.d.germline_end)
                self.log('  QUERY: ', self.d.query_alignment)
                self.log('         ', self.d.alignment_midline)
                self.log('GERMLINE:', self.d.germline_alignment)
            except:
                self.exception('D-GENE REALIGNMENT ERROR', traceback.format_exc())
                self.d = None
        else:
            self.log('')
            self.log('SKIPPING D-GENE REALIGHMENT')
            if self.chain not in ['heavy', 'beta', 'delta']:
                self.log(f'CHAIN ({self.chain}) DOES NOT RECOMBINE D GENES')
            elif self.d is None:
                self.log("D-GENE WAS NOT IDENTIFIED")
            else:
                self.log('UNKNOWN REASON FOR SKIPPING D-GENE REALIGNMENT')



# -------------
# VDJ ASSEMBLY
# -------------

    def _assemble_vdj_sequence(self):
        self.log('')
        self.log('VDJ ASSEMBLY')
        self.log('------------')
        self.gapped_vdj_nt = self._gapped_vdj_nt()
        self.log('GAPPED VDJ NT SEQUENCE:', self.gapped_vdj_nt)
        self.vdj_nt = self.gapped_vdj_nt.replace('-', '')
        self.log('VDJ NT SEQUENCE:', self.vdj_nt)
        self.vdj_aa = self._vdj_aa()
        self.log('VDJ AA SEQUENCE:', self.vdj_aa)
        self.gapped_vdj_germ_nt = self._gapped_vdj_germ_nt()
        self.log('GAPPED GERMLINE VDJ NT SEQUENCE:', self.gapped_vdj_germ_nt)
        self.vdj_germ_nt = self.gapped_vdj_germ_nt.replace('-', '')
        self.log('GERMLINE VDJ NT SEQUENCE:', self.vdj_germ_nt)
        self.vdj_germ_aa = self._vdj_germ_aa()
        self.log('GERMLINE VDJ AA SEQUENCE:', self.vdj_germ_aa)


    def _gapped_vdj_nt(self):
        'Returns the gapped nucleotide sequence of the VDJ region.'
        try:
            gapped_vdj_nt = self.v.query_alignment + \
                self.oriented_input[self.v.query_end + 1:self.j.query_start] + \
                self.j.query_alignment
            return gapped_vdj_nt
        except:
            self.exception('GAPPED VDJ NT ASSEMBLY', traceback.format_exc())


    def _vdj_aa(self):
        'Returns the amino acid sequence of the VDJ region.'
        # self.v_rf_offset = (len(self.oriented_input[self.v.query_start:self.junction.junction_nt_start]) % 3)
        # A correct offset would be expressed by the following statement but because all of the
        # functions operate on the sequence provided by the alignment, which is defines the start of the match
        # the only offset that matters is the offset between the full germline and the portion that matched.
        #self.v_rf_offset = self.oriented_input.sequence.find(self.vdj_nt) % 3
        if self.v.germline_start > 0:
            self.v_rf_offset = 3 - self.v.germline_start
        else:
            self.v_rf_offset = 0
        self.coding_start = self.v.query_start + self.v_rf_offset
        self.coding_end = self.j.query_end - (len(self.oriented_input[self.coding_start:self.j.query_end])) % 3
        self.coding_region = self.oriented_input[self.coding_start:self.coding_end + 1]
        translated_seq = Seq(self.coding_region).translate()
        self.log('READING FRAME OFFSET:', self.v_rf_offset)
        self.log('CODING START:', self.coding_start)
        self.log('CODING END:', self.coding_end)
        self.log('CODING REGION:', self.coding_region)
        return str(translated_seq)


    def _gapped_vdj_germ_nt(self):
        'Returns the gapped germline nucleotide sequence of the VDJ region.'
        try:
            if self.d:
                germ_junction = self.junction.n1_nt + self.d.germline_alignment + self.junction.n2_nt
            else:
                germ_junction = self.junction.n_nt
            germ_vdj_nt = self.v.germline_alignment + germ_junction + self.j.germline_alignment
            return germ_vdj_nt
        except:
            self.exception('GAPPED GERMLINE VDJ NT ASSEMBLY', traceback.format_exc())


    def _vdj_germ_aa(self):
        'Returns the germline amino acid sequence of the VDJ region.'
        trim = len(self.vdj_germ_nt) - (len(self.vdj_germ_nt[self.v_rf_offset:]) % 3)
        translated_seq = Seq(self.vdj_germ_nt[self.v_rf_offset:trim]).translate()
        return str(translated_seq)


    def _parse_uid(self, uid):
        if uid >= 0:
            self.uid = self.raw_input[:uid]
        else:
            self.uid = self.raw_input[uid:]


    def _get_junction(self):
        self.junction = junction.get_junction(self)


    def _identify_regions(self):
        '''
        Identifies and annotates variable/joining gene regions.
        '''
        # from abstar.utils import regions
        self.v.regions = regions.get_variable_regions(self)
        self.j.regions = regions.get_joining_regions(self)
        self.log('')
        self.log('VDJ REGIONS')
        self.log('-----------')
        self.log('FR1 NT sequence:', self.v.regions.nt_seqs['FR1'])
        self.log('CDR1 NT sequence:', self.v.regions.nt_seqs['CDR1'])
        self.log('FR2 NT sequence:', self.v.regions.nt_seqs['FR2'])
        self.log('CDR2 NT sequence:', self.v.regions.nt_seqs['CDR2'])
        self.log('CDR3 NT sequence:', self.junction.cdr3_nt)
        self.log('FR3 NT sequence:', self.v.regions.nt_seqs['FR3'])
        self.log('FR4 NT sequence:', self.j.regions.nt_seqs['FR4'])
        self.log('REGION POSITIONS:', self.v.regions.raw_nt_positions)
        self.log('FR1 AA sequence:', self.v.regions.aa_seqs['FR1'])
        self.log('CDR1 AA sequence:', self.v.regions.aa_seqs['CDR1'])
        self.log('FR2 AA sequence:', self.v.regions.aa_seqs['FR2'])
        self.log('CDR2 AA sequence:', self.v.regions.aa_seqs['CDR2'])
        self.log('FR3 AA sequence:', self.v.regions.aa_seqs['FR3'])
        self.log('CDR3 AA sequence:', self.junction.cdr3_aa)
        self.log('FR4 AA sequence:', self.j.regions.aa_seqs['FR4'])
        self.log('')


    def _mutations(self):
        '''
        Identifies and annotates nucleotide mutations.
        '''
        self.nt_mutations = mutations.nt_mutations(self)
        self.aa_mutations = mutations.aa_mutations(self)


    def _isotype(self):
        self.isotype = isotype.get_isotype(self)


    def _productivity(self):
        self.productivity = productivity.check_productivity(self)


    @staticmethod
    def _realignment_scoring_params(gene):
        '''
        Returns realignment scoring paramaters for a given gene type.

        Args:

            gene (str): the gene type ('V', 'D', or 'J')


        Returns:

            dict: realignment scoring parameters
        '''
        scores = {'V': {'match': 3,
                        'mismatch': -2,
                        'gap_open': -22,
                        'gap_extend': -1},
                  'D': {'match': 3,
                        'mismatch': -2,
                        'gap_open': -22,
                        'gap_extend': -1},
                  'J': {'match': 3,
                        'mismatch': -2,
                        'gap_open': -22,
                        'gap_extend': -1}}
        return scores[gene]
