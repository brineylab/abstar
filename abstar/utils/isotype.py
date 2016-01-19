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


import logging
import os
import traceback

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

from abtools import log
from abtools.alignment import local_alignment
from abtools.sequence import Sequence


def get_isotype(vdj):
	logger = log.get_logger(__name__)
	try:
		mod_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
		isotype_file = os.path.join(mod_dir, 'ssw/isotypes/{}_isotypes.fasta'.format(vdj.species))

		isotype_seqs = [Sequence(s) for s in SeqIO.parse(open(isotype_file, 'r'), 'fasta')]
		isotype_seqs += [Sequence((s.id, s.reverse_complement)) for s in isotype_seqs]
		alignments = local_alignment(vdj.raw_query, targets=isotype_seqs,
									 gap_open_penalty=22, gap_extend_penalty=1)
		alignments.sort(key=lambda x: x.score, reverse=True)
		return(alignments[0].target.id)
	except:
		logger.debug('ISOTYPE ERROR: {}\t{}'.format(vdj.id, vdj.raw_query))
		logger.debug(traceback.format_exc())
