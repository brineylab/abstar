#!/usr/bin/python
# filename: mutations.py

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


import math
import re
import traceback

from abtools import log
from abtoolsutils.codons import codons

from . import regions



def nt_mutations(antibody):
    try:
        all_mutations = Mutations()
        for segment in [antibody.v, antibody.j]:
            mutations = Mutations()
            for i, (q, g) in enumerate(zip(segment.query_alignment, segment.germline_alignment)):
                if q != g:
                    raw_pos = i + segment.query_start
                    imgt_pos = segment.get_imgt_position_from_raw(raw_pos)
                    # we don't want to count indels as mutations
                    if any([q == '-', g == '-', imgt_pos is None]):
                        continue
                    imgt_codon = int(math.ceil(imgt_pos / 3.0))
                    if segment.gene_type == 'J':
                        imgt_pos = segment.correct_imgt_nt_position_from_imgt[imgt_pos]
                    mutations.add(Mutation(g, q, raw_pos, imgt_pos, imgt_codon))
            segment.nt_mutations = mutations
            all_mutations.add_many(mutations.mutations)
        antibody.log('')
        antibody.log('NT MUTATIONS')
        antibody.log('------------')
        antibody.log('FR1:', ', '.join([m.abstar_formatted for m in all_mutations.in_region('FR1')]))
        antibody.log('CDR1:', ', '.join([m.abstar_formatted for m in all_mutations.in_region('CDR1')]))
        antibody.log('FR2:', ', '.join([m.abstar_formatted for m in all_mutations.in_region('FR2')]))
        antibody.log('CDR2:', ', '.join([m.abstar_formatted for m in all_mutations.in_region('CDR2')]))
        antibody.log('FR3:', ', '.join([m.abstar_formatted for m in all_mutations.in_region('FR3')]))
        antibody.log('CDR3:', ', '.join([m.abstar_formatted for m in all_mutations.in_region('CDR3')]))
        antibody.log('FR4:', ', '.join([m.abstar_formatted for m in all_mutations.in_region('FR4')]))
        return all_mutations
    except:
        antibody.exception('NT MUTATIONS', traceback.format_exc())
        return Mutations()


def aa_mutations(antibody):
    try:
        all_mutations = Mutations()
        for segment in [antibody.v, antibody.j]:
            mutations = Mutations()
            for i, (q, g) in enumerate(zip(segment.query_alignment, segment.germline_alignment)):
                raw_nt_pos = i + segment.query_start
                imgt_nt_pos = segment.get_imgt_position_from_raw(raw_nt_pos)
                if imgt_nt_pos is None:
                    continue
                if imgt_nt_pos % 3 != 1:
                    continue
                if any([q == '-', g == '-']):
                    continue
                q_codon_end_pos = i + 3
                if q_codon_end_pos >= len(segment.query_alignment):
                    continue
                q_codon = segment.query_alignment[i:q_codon_end_pos]
                while len(q_codon.replace('-', '')) < 3:
                    q_codon_end_pos += 1
                    q_codon = segment.query_alignment[i:q_codon_end_pos]
                    if q_codon_end_pos >= len(segment.query_alignment):
                        break
                if len(q_codon.replace('-', '')) < 3:
                    continue
                g_codon_end_pos = i + 3
                g_codon = segment.germline_alignment[i:g_codon_end_pos]
                while len(g_codon.replace('-', '')) < 3:
                    g_codon_end_pos += 1
                    g_codon = segment.germline_alignment[i:g_codon_end_pos]
                q_aa = codons[q_codon.replace('-', '')]
                g_aa = codons[g_codon.replace('-', '')]
                if q_aa != g_aa:
                    imgt_aa_pos = int(math.ceil(imgt_nt_pos / 3.0))
                    if segment.gene_type == 'J':
                        imgt_aa_pos = segment.correct_imgt_aa_position_from_imgt[imgt_aa_pos]
                    mutations.add(Mutation(g_aa, q_aa, None, imgt_aa_pos, imgt_aa_pos))
            segment.aa_mutations = mutations
            all_mutations.add_many(mutations.mutations)
        antibody.log('')
        antibody.log('AA MUTATIONS')
        antibody.log('------------')
        antibody.log('FR1:', ', '.join([m.abstar_formatted for m in all_mutations.in_region('FR1')]))
        antibody.log('CDR1:', ', '.join([m.abstar_formatted for m in all_mutations.in_region('CDR1')]))
        antibody.log('FR2:', ', '.join([m.abstar_formatted for m in all_mutations.in_region('FR2')]))
        antibody.log('CDR2:', ', '.join([m.abstar_formatted for m in all_mutations.in_region('CDR2')]))
        antibody.log('FR3:', ', '.join([m.abstar_formatted for m in all_mutations.in_region('FR3')]))
        antibody.log('CDR3:', ', '.join([m.abstar_formatted for m in all_mutations.in_region('CDR3')]))
        antibody.log('FR4:', ', '.join([m.abstar_formatted for m in all_mutations.in_region('FR4')]))
        return all_mutations
    except:
        antibody.exception('AA MUTATIONS', traceback.format_exc())


def _get_joining_imgt_mutation_position(codon_num, j):
    pass



class Mutation(object):
    """docstring for Mutation"""
    def __init__(self, was, now, raw_position, imgt_position, imgt_codon):
        super(Mutation, self).__init__()
        self.was = was
        self.now = now
        self.raw_position = raw_position
        self.imgt_position = imgt_position
        self.imgt_codon = imgt_codon


    @property
    def imgt_formatted(self):
        return '{}{}>{}'.format(self.was.lower(),
                                self.imgt_position,
                                self.now.lower())


    @property
    def abstar_formatted(self):
        return '{}:{}>{}'.format(self.imgt_position,
                                 self.was.upper(),
                                 self.now.upper())


    @property
    def json_formatted(self):
        return {'was': self.was,
                'is': self.now,
                'raw_position': self.raw_position,
                'position': self.imgt_position,
                'codon': self.imgt_codon}




class Mutations(object):
    """docstring for Mutations"""
    def __init__(self, mutations=None):
        super(Mutations, self).__init__()
        self.mutations = mutations if mutations is not None else []

    def __iter__(self):
        return (m for m in self.mutations)


    @property
    def count(self):
        return len(self.mutations)


    def add(self, mutation):
        '''
        Add a single mutation.

        Args:

            mutation (Mutation): the Mutation object to be added
        '''
        self.mutations.append(mutation)


    def add_many(self, mutations):
        '''
        Adds multiple mutations.

        Args:

            mutation (list): an iterable of Mutation objects to be added
        '''
        self.mutations += mutations


    def in_region(self, region):
        region_mutations = []
        start = regions.IMGT_REGION_START_POSITIONS_AA[region]
        end = regions.IMGT_REGION_END_POSITIONS_AA[region]
        for mut in self.mutations:
            if all([mut.imgt_codon >= start, mut.imgt_codon <= end]):
                region_mutations.append(mut)
        return region_mutations










# def nt_mutations(blast_result):
#     global logger
#     logger = log.get_logger(__name__)
#     try:
#         return MutationsNT(blast_result)
#     except:
#         logger.debug('NT MUTATIONS ERROR: {}, {}\n'.format(blast_result.id,
#                                                           blast_result.input_sequence))
#         logger.debug('QUERY ALIGNMENT: {}'.format(blast_result.query_alignment))
#         logger.debug('GERMLINE ALIGNMENT: {}'.format(blast_result.germline_alignment))
#         logger.debug(blast_result.regions.raw_positions)
#         logger.debug(blast_result.regions.adjusted_positions)
#         logger.debug(blast_result.regions.nt_seqs)
#         logger.debug(blast_result.regions.germline_nt_seqs)
#         logger.debug(traceback.format_exc())


# def aa_mutations(blast_result):
#     global logger
#     logger = log.get_logger(__name__)
#     try:
#         return MutationsAA(blast_result)
#     except:
#         logger.debug('AA MUTATIONS ERROR: {}, {}\n'.format(blast_result.id,
#                                                           blast_result.input_sequence))
#         logger.debug(traceback.format_exc())



# class MutationsNT(object):
#     '''
#     Structure for identifying and annotating nucleotide mutations.

#     Input is a BlastResult object.
#     '''
#     def __init__(self, blast_result):
#         self.all_mutations = self._find_all_mutations(blast_result)
#         self.mutation_count = self._total_mutation_count(blast_result)
#         self.germline_divergence = self._germline_divergence(blast_result)
#         self.germline_identity = 100. - self.germline_divergence
#         if blast_result.gene_type in ['variable', 'joining']:
#             self.region_mutations = self._find_region_mutations(blast_result)
#             self.region_mutation_count = self._region_mutation_count(blast_result)


#     def _find_all_mutations(self, br):
#         mutations = []
#         q = br.query_alignment
#         g = br.germline_alignment
#         m = br.alignment_midline
#         for i, pos in enumerate(m):
#             if pos == '|':
#                 continue
#             if q[i] == '-' or g[i] == '-':
#                 continue
#             mutations.append({'pos': br.germline_start + i, 'mut': '{}>{}'.format(g[i], q[i])})
#         return mutations


#     def _total_mutation_count(self, br):
#         return len(self.all_mutations)


#     def _find_region_mutations(self, br):
#         mutations = {}
#         start_positions = self._region_start_positions(br)
#         regions = self._region_list(br.gene_type)
#         for r in regions:
#             if br.regions.nt_seqs[r]:
#                 if br.gene_type == 'variable':
#                     start = start_positions[r] if start_positions[r] > br.germline_start else br.germline_start
#                 if br.gene_type == 'joining':
#                     start = br.germline_start + start_positions[r]
#                 mutations[r] = self._mutations_by_region(br, r, start)
#             else:
#                 mutations[r] = []
#         return mutations


#     def _mutations_by_region(self, br, region, start):
#         muts = []
#         q = br.regions.nt_seqs[region]
#         g = br.regions.germline_nt_seqs[region]
#         try:
#             for i, nt in enumerate(q):
#                 if nt == g[i] or nt == '-' or g[i] == '-':
#                     continue
#                 muts.append({'pos': i + start, 'mut': '{}>{}'.format(g[i], q[i])})
#         except IndexError:
#             logger.debug('Query: {}'.format(q))
#             logger.debug('Germline: {}'.format(g))
#             raise
#         return muts


#     def _region_mutation_count(self, br):
#         mutation_counts = {}
#         for r in self._region_list(br.gene_type):
#             if br.regions.nt_seqs[r]:
#                 mutation_counts[r] = len(self.region_mutations[r])
#             else:
#                 mutation_counts[r] = 0
#         return mutation_counts


#     def _region_list(self, gene_type):
#         if gene_type == 'variable':
#             return ['FR1', 'CDR1', 'FR2', 'CDR2', 'FR3']
#         if gene_type == 'joining':
#             return ['FR4', ]


#     def _region_start_positions(self, br):
#         positions = {}
#         if br.gene_type == 'variable':
#             positions['FR1'] = 0
#             for i, region in enumerate(self._region_list(br.gene_type)[1:]):
#                 positions[region] = br.regions.raw_positions[i]
#         if br.gene_type == 'joining':
#             for i, region in enumerate(self._region_list(br.gene_type)):
#                 positions[region] = br.regions.raw_positions[i]
#         return positions


#     def _germline_divergence(self, br):
#         return 100. * self.mutation_count / len(br.query_alignment.replace('-', ''))



# class MutationsAA(object):
#     '''
#     Structure for identifying and annotating amino acid mutations.

#     Input is a BlastResult object.
#     '''
#     def __init__(self, blast_result):
#         self.region_mutations = self._find_region_mutations(blast_result)
#         self.region_mutation_count = self._region_mutation_count(blast_result)
#         self.all_mutations = self._find_all_mutations(blast_result)
#         self.mutation_count = self._total_mutation_count(blast_result)
#         self.germline_divergence = self._germline_divergence(blast_result)
#         self.germline_identity = 100. - self.germline_divergence


#     def _find_region_mutations(self, br):
#         mutations = {}
#         start_positions = self._region_start_positions(br)
#         regions = self._region_list(br.gene_type)
#         for r in regions:
#             if br.regions.aa_seqs[r]:
#                 if br.gene_type == 'variable':
#                     start = start_positions[r] / 3 if start_positions[r] > br.germline_start else br.germline_start / 3
#                 if br.gene_type == 'joining':
#                     start = (br.germline_start + start_positions[r]) / 3
#                 mutations[r] = self._mutations_by_region(br, r, start)
#             else:
#                 mutations[r] = []
#         return mutations


#     def _mutations_by_region(self, br, region, start):
#         muts = []
#         q = br.regions.aa_seqs[region]
#         g = br.regions.germline_aa_seqs[region]
#         for i, aa in enumerate(q):
#             if aa == g[i] or aa == '-' or g[i] == '-':
#                 continue
#             muts.append({'pos': i + start, 'mut': '{}>{}'.format(g[i], q[i])})
#         return muts


#     def _region_mutation_count(self, br):
#         mutation_counts = {}
#         for r in self._region_list(br.gene_type):
#             if br.regions.aa_seqs[r]:
#                 mutation_counts[r] = len(self.region_mutations[r])
#             else:
#                 mutation_counts[r] = 0
#         return mutation_counts


#     def _region_list(self, gene_type):
#         if gene_type == 'variable':
#             return ['FR1', 'CDR1', 'FR2', 'CDR2', 'FR3']
#         if gene_type == 'joining':
#             return ['FR4', ]


#     def _region_start_positions(self, br):
#         positions = {}
#         if br.gene_type == 'variable':
#             positions['FR1'] = 0
#             for i, region in enumerate(self._region_list(br.gene_type)[1:]):
#                 positions[region] = br.regions.raw_positions[i] / 3
#         if br.gene_type == 'joining':
#             for i, region in enumerate(self._region_list(br.gene_type)):
#                 positions[region] = br.regions.raw_positions[i] / 3
#         return positions


#     def _find_all_mutations(self, br):
#         mutations = []
#         for r in self._region_list(br.gene_type):
#             if r in self.region_mutations:
#                 mutations.extend(self.region_mutations[r])
#         return mutations


#     def _total_mutation_count(self, br):
#         return len(self.all_mutations)


#     def _germline_divergence(self, br):
#         aa_length = sum([len(aa.replace('-', '')) for aa in br.regions.aa_seqs.values() if aa])
#         if aa_length:
#             return 100. * self.mutation_count / aa_length
#         return 0.0
