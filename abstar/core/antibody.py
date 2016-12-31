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


import traceback

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

from abtools.alignment import global_alignment, local_alignment

from ..utils import isotype, junction, mutations, productivity, regions
from ..utils.mixins import LoggingMixin
from .germline import get_imgt_germlines


class Antibody(object, LoggingMixin):
    """
    docstring for Antibody
    """
    def __init__(self, vdj, species):
        super(Antibody, self).__init__()
        LoggingMixin.__init__(self)
        self.id = vdj.sequence.id
        self.raw_input = vdj.sequence
        self.oriented_input = vdj.oriented
        self.v = vdj.v
        self.j = vdj.j
        self.d = vdj.d
        self.chain = vdj.v.chain
        self.species = species.lower()
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


    def initialize_log(self):
        try:
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
        except:
            print('\n')
            print('LOG INITIALIZATION ERROR')
            print(traceback.format_exc())



    def annotate(self, uid):
        '''
        Realigns the V(D)J germline genes using optimized SSW alignment
        parameters, followed by more detailed annotation (junction
        identification, ambig correction, parsing of mutations and
        indels, etc).
        '''
        try:
            print(self.id)
            print('Parsing UIDs...')
            self._parse_uid(uid)
            print('Realigning germlines...')
            self._realign_germlines()
            print('Processing junction...')
            self._get_junction()
            print('Assembling VDJ sequence...')
            self._assemble_vdj_sequence()
            print('Identifying regions...')
            self._identify_regions()
            print('Mutations...')
            self._mutations()
            print('Productivity...')
            self._productivity()
        except:


            print(traceback.format_exc())

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
        self.log('RAW QUERY LENGTH:', len(self.v.raw_query))
        self.log('RAW GERMLNE LENGTH:', len(self.v.raw_germline))
        self.log('QUERY START:', self.v.query_start)
        self.log('QUERY END:', self.v.query_end)
        self.log('GERMLINE START:', self.v.germline_start)
        self.log('GERMLINE END:', self.v.germline_end)
        self.log('')
        self.log('V-GENE GAPPED IMGT ALIGNMENT')
        self.log('----------------------------')
        self.v.gapped_imgt_realignment()
        self.log('  QUERY: ', self.v.imgt_gapped_alignment.aligned_query)
        self.log('         ', self.v.imgt_gapped_alignment.alignment_midline)
        self.log('GERMLINE:', self.v.imgt_gapped_alignment.aligned_target)

        # realign J
        jstart = self.v.query_end + 1
        self.log('')
        self.log('J-GENE REALIGNMENT')
        self.log('------------------')
        self.j.realign_germline(self, query_start=jstart)
        self.log('RAW QUERY LENGTH:', len(self.j.raw_query))
        self.log('RAW GERMLNE LENGTH:', len(self.j.raw_germline))
        self.log('QUERY START:', self.j.query_start)
        self.log('QUERY END:', self.j.query_end)
        self.log('GERMLINE START:', self.j.germline_start)
        self.log('GERMLINE END:', self.j.germline_end)
        self.log('')
        self.log('J-GENE GAPPED IMGT ALIGNMENT')
        self.log('----------------------------')
        self.j.gapped_imgt_realignment()
        self.log('  QUERY: ', self.j.imgt_gapped_alignment.aligned_query)
        self.log('         ', self.j.imgt_gapped_alignment.alignment_midline)
        self.log('GERMLINE:', self.j.imgt_gapped_alignment.aligned_target)

        # realign D (if needed)
        if self.chain == 'heavy':
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
        self.v_rf_offset = (len(self.oriented_input[self.v.query_start:self.junction.junction_nt_start]) % 3)
        self.coding_start = self.v.query_start + self.v_rf_offset
        self.coding_end = self.j.query_end - (len(self.oriented_input[self.coding_start:self.j.query_end])) % 3
        self.coding_region = self.oriented_input[self.coding_start:self.coding_end + 1]
        translated_seq = Seq(self.coding_region, generic_dna).translate()
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
        translated_seq = Seq(self.vdj_germ_nt[self.v_rf_offset:trim], generic_dna).translate()
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
