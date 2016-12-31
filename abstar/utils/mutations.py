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


from __future__ import print_function

import math
import re
import traceback

from abtools import log
from abtools.utils.codons import codon_lookup as codons

from .regions import IMGT_REGION_START_POSITIONS_AA, IMGT_REGION_END_POSITIONS_AA



def nt_mutations(antibody):
    try:
        all_mutations = Mutations()
        for segment in [antibody.v, antibody.j]:
            if segment is None:
                continue
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
            segment.nt_identity = 100. - (100. * mutations.count / min(len(segment.query_alignment), len(segment.germline_alignment)))
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


        print(antibody.id, traceback.format_exc(), sep='\n')


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
            segment.aa_identity = 100. - (100. * mutations.count / len(segment.aa_sequence))
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
        start = IMGT_REGION_START_POSITIONS_AA[region]
        end = IMGT_REGION_END_POSITIONS_AA[region]
        for mut in self.mutations:
            if all([mut.imgt_codon >= start, mut.imgt_codon <= end]):
                region_mutations.append(mut)
        return region_mutations


    def in_region_count(self, region):
        return len(self.in_region(region))
