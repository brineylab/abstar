#!/usr/bin/python
# filename: isotype.py

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

import logging
import os
import traceback

from Bio import SeqIO
from Bio.Seq import Seq

from abutils.core.sequence import Sequence
from abutils.utils import log
from abutils.utils.alignment import local_alignment
from abutils.utils.decorators import lazy_property

from ..core.germline import get_germline_database_directory


def get_isotype(antibody):
    try:
        germ_dir = get_germline_database_directory(antibody.germ_db)
        isotype_file = os.path.join(germ_dir, 'isotypes/isotypes.fasta')
        isotype_seqs = [Sequence(s) for s in SeqIO.parse(open(isotype_file, 'r'), 'fasta')]
        return Isotype(antibody, isotype_seqs)
    except:
        antibody.exception('ISOTYPING ERROR', traceback.format_exc())


# def get_isotype(vdj):
#     logger = log.get_logger(__name__)
#     try:
#         mod_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
#         isotype_file = os.path.join(mod_dir, 'ssw/isotypes/{}_isotypes.fasta'.format(vdj.species))
#         isotype_seqs = [Sequence(s) for s in SeqIO.parse(open(isotype_file, 'r'), 'fasta')]
#         isotype_seqs += [Sequence((s.id, s.reverse_complement)) for s in isotype_seqs]
#         return Isotype(vdj, isotype_seqs)
#     except:
#         logger.debug('ISOTYPE ERROR: {}\t{}'.format(vdj.id, vdj.raw_query))
#         logger.debug(traceback.format_exc())



class Isotype(object):
    """docstring for Isotype"""
    def __init__(self, antibody, isotype_seqs):
        super(Isotype, self).__init__()
        self._alignments = self._get_alignments(antibody, isotype_seqs)
        self.alignment = self._alignments[0]


    @lazy_property
    def isotype(self):
        if self.normalized_score < 2:
            return 'unknown'
        return self.alignment.target.id

    @lazy_property
    def score(self):
        return self.alignment.score

    @lazy_property
    def normalized_score(self):
        return float(self.alignment.score) / len(self.alignment)


    def _get_alignments(self, antibody, isotype_seqs):
        query_region = self._get_isotype_query_region(antibody)
        alignments = local_alignment(query_region, targets=isotype_seqs,
            gap_open_penalty=22, gap_extend_penalty=1)
        return sorted(alignments, key=lambda x: x.score, reverse=True)


    def _get_isotype_query_region(self, antibody):
        aln = local_alignment(antibody.vdj_nt, antibody.oriented_input)
        return antibody.oriented_input[aln.target_end:]
