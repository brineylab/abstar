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

from abtools.utils import log
from abtools.utils.alignment import local_alignment


def get_isotype(vdj):
	logger = log.get_logger(__name__)
	try:
		mod_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
		isotype_file = os.path.join(mod_dir, 'ssw/isotypes/{}_isotypes.fasta'.format(vdj.species))
		isotypes = [Isotype(i) for i in SeqIO.parse(open(isotype_file, 'r'), 'fasta')]
		for i in isotypes:
			i.align(vdj.raw_query)
		isotypes.sort(key=lambda x: x.top_score, reverse=True)
		return isotypes[0].name
	except:
		logger.debug('ISOTYPE ERROR: {}'.format(vdj.id))
		logger.debug(traceback.format_exc())


class Isotype(object):
	"""
	Holds isotype information for the primer used to
	amplify a single isotype (ex: IgG).

	Input is a single biopython SeqRecord object.
	"""
	def __init__(self, isotype):
		self.name = isotype.id
		self.seq = str(isotype.seq)
		self.rc_seq = self._rev_comp(self.seq)
		self.scores = []
		self.top_score = None


	def align(self, query):
		'''
		Calculates the alignment score (sense and anti-sense) of the
		isotyping primer with the query sequence.
		'''
		query_chunk = query[:len(self.seq) + 50]
		for s in [self.seq, self.rc_seq]:
			self.scores.append(self._score_alignment(query_chunk, s))
		self.top_score = max(self.scores)


	def _score_alignment(self, query, seq):
		# ssw = Aligner(query,
		# 			  match=3,
		# 			  mismatch=2,
		# 		  	  gap_open=22,
		# 		  	  gap_extend=1,
		# 		  	  report_cigar=True)
		# alignment = ssw.align(seq)
		# return alignment.score

		aln = local_alignment(query, seq,
							  gap_open_penalty=22, gap_extend_penalty=1)
		return aln.score

	def _rev_comp(self, seq):
		'''
		Returns the reverse complement of <seq> as a string.
		'''
		s = Seq(self.seq, generic_dna)
		return s.reverse_complement().upper()
