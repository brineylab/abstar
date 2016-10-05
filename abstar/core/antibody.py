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
        self._strand = None


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

        # TODO
        self._assemble_vdj_sequence()
        self._imgt_position_numbering()
        self._identify_regions()


    def log(self, *args, **kwargs):
        sep = kwargs.get('sep', ' ')
        lstring = sep.join([str(a) for a in args])
        self._log.append(lstring)


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
        return '\n'.join(self._log)


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


    # def _realign_germline(self, segment, query_start=None, query_end=None):
    #     '''
    #     Due to restrictions on the available scoring parameters in BLASTn, incorrect truncation
    #     of the v-gene alignment can occur. This function re-aligns the query sequence with
    #     the identified germline variable gene using more appropriate alignment parameters.

    #     Input is a Germline object, corresponding to the gene segment that should be realigned (V, D or J).
    #     '''
    #     germline_seq = self._get_germline_sequence_for_realignment(segment.full)
    #     query = self.oriented_input.sequence[query_start:query_end]
    #     aln_params = self._realignment_scoring_params(segment.gene_type)
    #     alignment = local_alignment(query, germline_seq, **aln_params)
    #     segment.process_realignment(segment, alignment)


    # def _get_germline_sequence_for_realignment(self, germ):
    #     '''
    #     Identifies the appropriate germline variable gene from a database of all
    #     germline variable genes.

    #     Args:

    #         germ (string): Full name of the germline gene using IMGT nomenclature
    #             (ex: 'IGHV1-2*02')

    #     Returns:

    #         Sequence: Requested germline gene, as an AbTools `Sequence` object, or `None`
    #             if the requested germline gene could not be found.
    #     '''
    #     gene_type = germ[3]
    #     mod_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    #     db_file = os.path.join(mod_dir, 'vdj/germline_dbs/fasta/{}_{}.fasta'.format(self.species, gene_type))
    #     for s in SeqIO.parse(open(db_file), 'fasta'):
    #         if s.id == germ:
    #             return Sequence(s)
    #     return None


    # def _process_realignment(self, segment, aln, aln_offset):
    #     aln_offset = aln_offset if aln_offset is not None else 0
    #     segment.raw_query = aln.raw_query
    #     segment.raw_germline = aln.raw_target
    #     segment.query_alignment = aln.aligned_query
    #     segment.germline_alignment = aln.aligned_target
    #     segment.alignment_midline = ''.join(['|' if q == g else ' ' for q, g in zip(aln.aligned_query,
    #                                                                                 aln.aligned_target)])
    #     segment.query_start = aln.query_begin
    #     segment.query_end = aln.query_end
    #     segment.germline_start = aln.target_begin
    #     segment.germline_end = aln.target_end
    #     segment.fix_ambigs()
    #     segment.find_indels()


    # def _gapped_imgt_vgene_realignment(self, v):
    #     imgt_germ = get_imgt_germlines(species=v.species,
    #                                    gene_type=v.gene_type,
    #                                    gene=v.full)
    #     query = v.query_alignment.replace('-', '')
    #     aln_params = self._realignment_scoring_params(v.gene_type)
    #     aln_params['gap_open_penalty'] = 11
    #     v.imgt_gapped_alignment = global_alignment(query, imgt_germ.gapped_nt_sequence, **aln_params)


    # def _imgt_vgene_numbering(self, v):
    #     aln = v.imgt_gapped_alignment
    #     imgt_position_lookup = {}
    #     upos = 0
    #     for i, (g, u) in enumerate(zip(aln.aligned_target, aln.aligned_query)):
    #         if all([u == '-', g == '.']):
    #             continue
    #         imgt_pos = i + 1
    #         imgt_position_lookup[upos] = imgt_pos
    #         upos += 1
    #     v.get_imgt_position_from_raw = imgt_position_lookup
    #     v.get_raw_position_frim_imgt = {v: k for k, v in imgt_position_lookup.items()}


    def _identify_regions(self):
        '''
        Identifies and annotates variable/joining gene regions.
        '''
        # from abstar.utils import regions
        return regions.regions(self)


    def _assemble_vdj_sequence(self):
        self.gapped_vdj_nt = self._gapped_vdj_nt()
        self.vdj_nt = self.gapped_vdj_nt.replace('-', '')
        self.vdj_aa = self._vdj_aa()
        self.gapped_vdj_germ_nt = self._gapped_vdj_germ_nt()
        self.vdj_germ_nt = self.gapped_vdj_germ_nt.replace('-', '')
        self.vdj_germ_aa = self._vdj_germ_aa()

    def _gapped_vdj_nt(self):
        'Returns the gapped nucleotide sequence of the VDJ region.'
        try:
            gapped_vdj_nt = self.v.query_alignment + \
                self.j.input_sequence[:self.j.query_start] + \
                self.j.query_alignment
            return gapped_vdj_nt
        except:
            logger.debug('VDJ NT ERROR: {}, {}'.format(self.id, self.raw_query))

    def _vdj_aa(self):
        'Returns the amino acid sequence of the VDJ region.'
        offset = (self.query_reading_frame * 2) % 3
        trim = len(self.vdj_nt) - (len(self.vdj_nt[offset:]) % 3)
        translated_seq = Seq(self.vdj_nt[offset:trim], generic_dna).translate()
        return str(translated_seq)

    def _gapped_vdj_germ_nt(self):
        'Returns the gapped germline nucleotide sequence of the VDJ region.'
        try:
            if self.d:
                germ_junction = self.junction.n1_nt + \
                    self.d.germline_alignment + \
                    self.junction.n2_nt
            else:
                germ_junction = self.junction.n_nt
            germ_vdj_nt = self.v.germline_alignment + \
                germ_junction + self.j.germline_alignment
            return germ_vdj_nt
        except:
            logger.debug('VDJ GERM NT ERROR: {}, {}'.format(self.id, self.raw_query))
            logger.debug(traceback.format_exc())

    def _vdj_germ_aa(self):
        'Returns the germline amino acid sequence of the VDJ region.'
        offset = (self.query_reading_frame * 2) % 3
        trim = len(self.vdj_germ_nt) - (len(self.vdj_germ_nt[offset:]) % 3)
        translated_seq = Seq(self.vdj_germ_nt[offset:trim], generic_dna).translate()
        return str(translated_seq)


    def _get_junction(self):
        # from abstar.utils import junction
        self.junction = junction.get_junction(self)
        self.germ_junction = junction.get_junction(self, germ=True)


    def _imgt_position_numbering(self):
        pass


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




