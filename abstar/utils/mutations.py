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


import re
import traceback

from abtools.utils import log


def nt_mutations(blast_result):
	logger = log.get_logger(__name__)
	try:
		return MutationsNT(blast_result)
	except:
		logger.debug('NT MUTATIONS ERROR: {}, {}\n'.format(blast_result.id,
														  blast_result.input_sequence))


def aa_mutations(blast_result):
	logger = log.get_logger(__name__)
	try:
		return MutationsAA(blast_result)
	except:
		logger.debug('AA MUTATIONS ERROR: {}, {}\n'.format(blast_result.id,
														  blast_result.input_sequence))



class MutationsNT(object):
	'''
	Structure for identifying and annotating nucleotide mutations.

	Input is a BlastResult object.
	'''
	def __init__(self, blast_result):
		self.all_mutations = self._find_all_mutations(blast_result)
		self.mutation_count = self._total_mutation_count(blast_result)
		self.germline_divergence = self._germline_divergence(blast_result)
		self.germline_identity = 100. - self.germline_divergence
		if blast_result.gene_type in ['variable', 'joining']:
			self.region_mutations = self._find_region_mutations(blast_result)
			self.region_mutation_count = self._region_mutation_count(blast_result)


	def _find_all_mutations(self, br):
		mutations = []
		q = br.query_alignment
		g = br.germline_alignment
		m = br.alignment_midline
		for i, pos in enumerate(m):
			if pos == '|':
				continue
			if q[i] == '-' or g[i] == '-':
				continue
			mutations.append({'pos': br.germline_start + i, 'mut': '{}>{}'.format(g[i], q[i])})
		return mutations


	def _total_mutation_count(self, br):
		ins_count = len(re.findall('-+', br.germline_alignment))
		del_count = len(re.findall('-+', br.query_alignment))
		return len(self.all_mutations) + ins_count + del_count


	def _find_region_mutations(self, br):
		mutations = {}
		start_positions = self._region_start_positions(br)
		regions = self._region_list(br.gene_type)
		for r in regions:
			if br.regions.nt_seqs[r]:
				if br.gene_type == 'variable':
					start = start_positions[r] if start_positions[r] > br.germline_start else br.germline_start
				if br.gene_type == 'joining':
					start = br.germline_start + start_positions[r]
				mutations[r] = self._mutations_by_region(br, r, start)
			else:
				mutations[r] = []
		return mutations


	def _mutations_by_region(self, br, region, start):
		muts = []
		q = br.regions.nt_seqs[region]
		g = br.regions.germline_nt_seqs[region]
		for i, nt in enumerate(q):
			if nt == g[i] or nt == '-' or g[i] == '-':
				continue
			muts.append({'pos': i + start, 'mut': '{}>{}'.format(g[i], q[i])})
		return muts


	def _region_mutation_count(self, br):
		mutation_counts = {}
		for r in self._region_list(br.gene_type):
			if br.regions.nt_seqs[r]:
				ins_count = len(re.findall('-+', br.regions.germline_nt_seqs[r]))
				del_count = len(re.findall('-+', br.regions.nt_seqs[r]))
				mutation_counts[r] = len(self.region_mutations[r]) + ins_count + del_count
			mutation_counts[r] = 0
		return mutation_counts


	def _region_list(self, gene_type):
		if gene_type == 'variable':
			return ['FR1', 'CDR1', 'FR2', 'CDR2', 'FR3']
		if gene_type == 'joining':
			return ['FR4', ]


	def _region_start_positions(self, br):
		positions = {}
		if br.gene_type == 'variable':
			positions['FR1'] = 0
			for i, region in enumerate(self._region_list(br.gene_type)[1:]):
				positions[region] = br.regions.raw_positions[i]
		if br.gene_type == 'joining':
			for i, region in enumerate(self._region_list(br.gene_type)):
				positions[region] = br.regions.raw_positions[i]
		return positions


	def _germline_divergence(self, br):
		return 100. * self.mutation_count / len(br.query_alignment.replace('-', ''))



class MutationsAA(object):
	'''
	Structure for identifying and annotating amino acid mutations.

	Input is a BlastResult object.
	'''
	def __init__(self, blast_result):
		self.region_mutations = self._find_region_mutations(blast_result)
		self.region_mutation_count = self._region_mutation_count(blast_result)
		self.all_mutations = self._find_all_mutations(blast_result)
		self.mutation_count = self._total_mutation_count(blast_result)
		self.germline_divergence = self._germline_divergence(blast_result)
		self.germline_identity = 100. - self.germline_divergence


	def _find_region_mutations(self, br):
		mutations = {}
		start_positions = self._region_start_positions(br)
		regions = self._region_list(br.gene_type)
		for r in regions:
			if br.regions.aa_seqs[r]:
				if br.gene_type == 'variable':
					start = start_positions[r] / 3 if start_positions[r] > br.germline_start else br.germline_start / 3
				if br.gene_type == 'joining':
					start = (br.germline_start + start_positions[r]) / 3
				mutations[r] = self._mutations_by_region(br, r, start)
			else:
				mutations[r] = []
		return mutations


	def _mutations_by_region(self, br, region, start):
		muts = []
		q = br.regions.aa_seqs[region]
		g = br.regions.germline_aa_seqs[region]
		for i, aa in enumerate(q):
			if aa == g[i] or aa == '-' or g[i] == '-':
				continue
			muts.append({'pos': i + start, 'mut': '{}>{}'.format(g[i], q[i])})
		return muts


	def _region_mutation_count(self, br):
		mutation_counts = {}
		for r in self._region_list(br.gene_type):
			if br.regions.aa_seqs[r]:
				ins_count = len(re.findall('-+', br.regions.germline_nt_seqs[r]))
				del_count = len(re.findall('-+', br.regions.nt_seqs[r]))
				mutation_counts[r] = len(self.region_mutations[r]) + ins_count + del_count
			else:
				mutation_counts[r] = 0
		return mutation_counts


	def _region_list(self, gene_type):
		if gene_type == 'variable':
			return ['FR1', 'CDR1', 'FR2', 'CDR2', 'FR3']
		if gene_type == 'joining':
			return ['FR4', ]


	def _region_start_positions(self, br):
		positions = {}
		if br.gene_type == 'variable':
			positions['FR1'] = 0
			for i, region in enumerate(self._region_list(br.gene_type)[1:]):
				positions[region] = br.regions.raw_positions[i] / 3
		if br.gene_type == 'joining':
			for i, region in enumerate(self._region_list(br.gene_type)):
				positions[region] = br.regions.raw_positions[i] / 3
		return positions


	def _find_all_mutations(self, br):
		mutations = []
		for r in self._region_list(br.gene_type):
			if r in self.region_mutations:
				mutations.extend(self.region_mutations[r])
		return mutations


	def _total_mutation_count(self, br):
		ins_count = len(re.findall('-+', br.germline_alignment))
		del_count = len(re.findall('-+', br.query_alignment))
		return len(self.all_mutations) + ins_count + del_count


	def _germline_divergence(self, br):
		aa_length = sum([len(aa.replace('-', '')) for aa in br.regions.aa_seqs.values() if aa])
		if aa_length:
			return 100. * self.mutation_count / aa_length
		return 0.0
