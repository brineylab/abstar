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


from abstar.utils import junction, mutations, positions, regions
from abstar.vdj.germline import get_imgt_germlines

from abtools.alignment import global_alignment, local_alignment


class Antibody(object):
    """
    docstring for Antibody
    """
    def __init__(self, vdj, species):
        super(Antibody, self).__init__()
        self.id = vdj.sequence.id
        self.raw_input = vdj.sequence
        self.oriented_input = vdj.oriented
        self.v = vdj.v
        self.j = vdj.j
        self.d = vdj.d
        self.chain = vdj.v.chain
        self.species = species.lower()
        self._log = self._initialize_log()
        self._exceptions = []
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


    def annotate(self):
        '''
        Realigns the V(D)J germline genes using optimized SSW alignment
        parameters, followed by more detailed annotation (junction
        identification, ambig correction, parsing of mutations and
        indels, etc).
        '''
        self._realign_germlines()
        self._get_junction()
        self._assemble_vdj_sequence()

        # TODO
        self._imgt_position_numbering()
        self._identify_regions()


    def log(self, *args, **kwargs):
        sep = kwargs.get('sep', ' ')
        lstring = sep.join([str(a) for a in args])
        self._log.append(lstring)


    def exception(self, *args, **kwargs):
        sep = kwargs.get('sep', '\n')
        estring = sep.join([str(a) for a in args])
        self._exceptions.append(estring)


    def format_log(self):
        '''
        Formats the antibody log.

        Log formatting is only performed on sequences that had an
        error during annotation, unless AbStar is run in debug
        mode. In debug mode, all sequences will be logged.

        Returns:
        --------

            str: Formatted log string.
        '''
        self._log += ['', '']
        output = '\n'.join(self._log)
        if self._exceptions:
            output += '\n'
            output += self._format_exceptions
        return output


    def _initialize_log(self):
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
        if self.v is not None:
            # log V-gene stuff
            pass
        if self.j is not None:
            # log J-gene stuff
            pass
        if self.d is not None:
            # log D-gene stuff
            pass
        return log


    def _format_exceptions(self):
        estring = 'EXCEPTIONS'
        estring += '----------'
        estring += ['\n\n'.join(e for e in self._exceptions)]
        return estring


    def _realign_germlines(self):
        # realign V and compute gapped IMGT alignment for the V gene
        self.log('')
        self.log('V-GENE REALIGNMENT')
        self.log('------------------')
        self.v.realign_germline(self.oriented_input)
        self.log('')
        self.log('V-GENE GAPPED IMGT ALIGNMENT')
        self.log('----------------------------')
        self.v.gapped_imgt_realignment()

        # realign J
        jstart = self.v.query_end
        self.log('')
        self.log('J-GENE REALIGNMENT')
        self.log('------------------')
        self.j.realign_germline(self.oriented_input, query_start=jstart)
        self.log('')
        self.log('J-GENE GAPPED IMGT ALIGNMENT')
        self.log('----------------------------')
        self.j.gapped_imgt_realignment()

        # realign D (if needed)
        if self.chain == 'heavy':
            dstart = self.v.query_end
            dend = self.v.query_end + self.j.query_start
            self.d.realign_germline(self.oriented_input, query_start=dstart, query_end=dend)



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
                self.j.input_sequence[:self.j.query_start] + \
                self.j.query_alignment
            return gapped_vdj_nt
        except:
            self.log('VDJ NT ERROR: {}, {}'.format(self.id, self.raw_query))

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
            self.log('VDJ GERM NT ERROR:', self.id, self.raw_query)
            self.log(traceback.format_exc())

    def _vdj_germ_aa(self):
        'Returns the germline amino acid sequence of the VDJ region.'
        trim = len(self.vdj_germ_nt) - (len(self.vdj_germ_nt[self.v_rf_offset:]) % 3)
        translated_seq = Seq(self.vdj_germ_nt[self.v_rf_offset:trim], generic_dna).translate()
        return str(translated_seq)


    def _get_junction(self):
        self.junction = junction.get_junction(self)
        # self.germ_junction = junction.get_junction(self, germ=True)


    def _imgt_position_numbering(self):
        
        pass




    def _imgt_nt_position_numbering(self):
        # V-gene
        v_aln_start = self.v.imgt_position_from_raw(self.v.query_start)  # IMGT position of the first aligned V-gene position
        for i, nt in self.v.query_alignment[:-len(self.junction.v_nt)]:  # don't want to include any of the junction, just through FR3
            imgt_pos = self.v.imgt_position_from_raw(i + v_aln_start)
            self._imgt_nt_position_from_aligned[i] = imgt_pos
            self._aligned_nt_position_from_imgt[imgt_pos] = i

        # J-gene
        junc_start = 310







    def _identify_regions(self):
        '''
        Identifies and annotates variable/joining gene regions.
        '''
        # from abstar.utils import regions
        return regions.regions(self)




    def _nt_mutations(self):
        '''
        Identifies and annotates nucleotide mutations.
        '''
        # from abstar.utils import mutations
        return mutations.nt_mutations(self)


    def _aa_mutations(self):
        '''
        Identifies and annotates amino acid mutations.
        '''
        # from abstar.utils import mutations
        return mutations.aa_mutations(self)


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




