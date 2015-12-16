#!/usr/bin/python
# filename: alignment.py

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


import os
import re

import nwalign as nw
from abstar.ssw.ssw_wrap import Aligner


class NWAlignment(object):
	"""
	Stucture for performing and analyzing a Needleman-Wunch alignment
	using the nwalign package.

	Input is two sequences (query chunk and germline chunk) and the position
	of the query chunk.
	"""
	def __init__(self, query, germ, position=None):
		self.query = query
		self.germ = germ
		self.position = position
		mod_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
		self.matrix = os.path.join(mod_dir, 'nw/blosum62')
		self.alignment = self._nw_alignment()
		self.score = self._score_alignment()


	def _nw_alignment(self):
		return nw.global_align(self.query,
							   self.germ,
							   gap_open=-100,
							   gap_extend=-50,
							   matrix=self.matrix)


	def _score_alignment(self):
		aligned_query = self.alignment[0]
		aligned_germ = self.alignment[1]
		return nw.score_alignment(aligned_query,
								  aligned_germ,
								  gap_open=-100,
								  gap_extend=-50,
								  matrix=self.matrix)


class SSWAlignment(object):
	"""docstring for SSWAlignment"""
	def __init__(self, seq1, seq2):
		self.seq1 = seq1
		self.seq2 = seq2
		self._align()

	def _align(self):
		ssw = Aligner(self.seq1,
					  match=3,
					  mismatch=2,
				  	  gap_open=22,
				  	  gap_extend=1,
				  	  report_cigar=True)
		alignment = ssw.align(self.seq2)
		self._process_alignment(Alignment(self.seq1, self.seq2, alignment))

	def _process_alignment(self, alignment):
		'''
		Processes the result of variable gene realignment and updates BlastResult
		attributes accordingly.

		Input is an Alignment object.
		'''
		self.seq1_alignment, self.seq2_alignment = self._cigar_to_alignment(alignment)
		self.seq1_start = alignment.query_begin
		self.seq1_end = alignment.query_end
		self.seq2_start = alignment.target_begin
		self.seq2_end = alignment.target_end

	def _cigar_to_alignment(self, alignment):
		'''
		Converts the cigar string and aligned region from the query and target
		sequences into a gapped alignment.

		Input is an Alignment object.

		Output is a gapped query sequence and a gapped target sequence.
		'''
		query = ''
		germline = ''
		pattern = re.compile('([MIDNSHPX=])')
		values = pattern.split(alignment.cigar)[:-1]
		pairs = zip(*2*[iter(values)])
		for p in pairs:
			query, germline = self._alignment_chunk_from_cigar_pair(alignment,
																	p,
																	query,
																	germline)
		return query, germline

	def _alignment_chunk_from_cigar_pair(self, alignment, pair, query, germline):
		'''
		Processes a single 'pair' of a cigar string and lengthens the gapped query
		and targed alignment strings.

		Input is an Alignment object, the cigar pair, the partial gapped alignment
		strings for the query and target. If 30M3I250M were to be a cigar string,
		it would have three pairs -- (30, 'M'), (3, 'I') and (250, 'M') --
		which are each formatted as tuples.

		Output is a gapped query alignment and gapped germline alignment that incorporate
		the data from the given cigar string.
		'''
		q = len(query.replace('-', ''))
		g = len(germline.replace('-', ''))
		clength = int(pair[0])
		ctype = pair[1]
		if ctype.upper() == 'S':
			return query, germline
		if ctype.upper() == 'M':
			query += alignment.aligned_query[q:q + clength]
			germline += alignment.aligned_target[g:g + clength]
		elif ctype.upper() == 'D':
			query += alignment.aligned_query[q:q + clength]
			germline += '-' * clength
		elif ctype.upper() == 'I':
			query += '-' * clength
			germline += alignment.aligned_target[g:g + clength]
		return query, germline


class Alignment(object):
	'''
	Data structure to hold the result of a SSW local alignment. Parses the
	raw data from the output of a ssw.Aligner alignment (a PyAlignRes object)
	into something more useful.

	Input is the ID of the top-scoring target, the query sequence, the target
	sequence and the PyAlignRes object.
	'''
	def __init__(self, query_seq, target_seq, alignment):
		super(Alignment, self).__init__()
		self.query_seq = query_seq
		self.target_seq = target_seq
		self.alignment = alignment
		self.cigar = alignment.cigar_string
		self.query_begin = alignment.ref_begin
		self.query_end = alignment.ref_end + 1
		self.target_begin = alignment.query_begin
		self.target_end = alignment.query_end + 1
		self.aligned_query = self._get_aligned_query()
		self.aligned_target = self._get_aligned_target()

	def _get_aligned_query(self):
		'Returns the aligned portion of the query sequence'
		return self.query_seq[self.query_begin:self.query_end]

	def _get_aligned_target(self):
		'Returns the aligned portion of the target sequence'
		return self.target_seq[self.target_begin:self.target_end]
