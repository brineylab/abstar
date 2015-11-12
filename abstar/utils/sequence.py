#!/usr/bin/python
# filename: sequence.py



###########################################################################
#
# Copyright (c) 2014 Bryan Briney.  All rights reserved.
#
# @version: 1.0.0
# @author: Bryan Briney
# @license: MIT (http://opensource.org/licenses/MIT)
#
###########################################################################


from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna


class Sequence(object):
	'''
	Provides basic functions for processing sequences.

	Input is a Biopython SeqRecord object.
	'''
	def __init__(self, seq):
		super(Sequence, self).__init__()
		self.id = seq.id
		self.input = str(seq.seq).upper()
		self.sequence = str(seq.seq).upper()
		self.strand = 'plus'

	def as_fasta(self, start=0, end=None):
		'''
		Returns the sequence, from 'start' to 'end' as a FASTA string. If start
		and end are not provided, returns the entire sequence.
		'''
		if not end:
			end = len(self.sequence)
		return '>{}\n{}'.format(self.id, self.sequence[start:end])

	def region(self, start=0, end=None):
		'''
		Returns the sequence (as a string) from 'start' to 'end'. If start and end
		are not provided, returns the entire sequence.
		'''
		if not end:
			end = len(self.sequence)
		return self.sequence[start:end]

	def reverse_complement(self):
		'''
		Returns a string containing the reverse complement of the sequence.
		'''
		s = Seq(self.sequence, generic_dna)
		rc = s.reverse_complement()
		self.sequence = str(rc)
		self.strand = 'minus'

	def rc(self):
		'''
		Returns a Sequence object of the reverse complement.
		'''
		s = Seq(self.sequence, generic_dna)
		rev_comp = s.reverse_complement()
		return SeqRecord(rev_comp,
						 id=self.id)


