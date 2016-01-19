#!/usr/bin/python
# filename: codons.py

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
import traceback

from abtools import log


def parse_codons(vdj, gapped):
	logger = log.get_logger(__name__)
	try:
		codons = Codons(vdj, gapped)
		return codons
	except:
		logger.debug('PARSE CODONS ERROR: {}, {}'.format(vdj.id, vdj.raw_input))
		logger.debug(traceback.format_exc())


class Codons(object):
	"""docstring for Codons"""
	def __init__(self, vdj, gapped=False):
		super(Codons, self).__init__()
		self.gapped = gapped
		self.vdj_codons = self._get_vdj_codons(vdj)
		self.vdj_germ_codons = self._get_vdj_codons(vdj)
		self.vdj_codon_regions = self._get_vdj_codon_regions(vdj)
		self.v_codons = self._get_v_codons(vdj)
		self.v_germ_codons = self._get_v_germ_codons(vdj)
		self.j_codons = self._get_j_codons(vdj)
		self.j_germ_codons = self._get_j_germ_codons(vdj)


	def _get_vdj_codons(self, vdj):
		offset = (vdj.query_reading_frame * 2) % 3
		if self.gapped:
			trunc_vdj = vdj.gapped_vdj_nt[offset:]
		else:
			trunc_vdj = vdj.vdj_nt[offset:]
		codons = [trunc_vdj[i:i + 3] for i in range(0, len(trunc_vdj), 3)]
		if len(codons[-1]) != 3:
			codons = codons[:-1]
		return codons

	def _get_vdj_codons(self, vdj):
		offset = (vdj.query_reading_frame * 2) % 3
		if self.gapped:
			trunc_vdj = vdj.gapped_vdj_germ_nt[offset:]
		else:
			trunc_vdj = vdj.vdj_germ_nt[offset:]
		codons = [trunc_vdj[i:i + 3] for i in range(0, len(trunc_vdj), 3)]
		if len(codons[-1]) != 3:
			codons = codons[:-1]
		return codons

	def _get_vdj_codon_regions(self, vdj):
		offset = (vdj.query_reading_frame * 2) % 3
		if self.gapped:
			trunc_region = vdj.gapped_vdj_region_string[offset:]
		else:
			trunc_region = vdj.vdj_region_string[offset:]
		codons = [trunc_region[i:i + 3] for i in range(0, len(trunc_region), 3)]
		if len(codons[-1]) != 3:
			codons = codons[:-1]
		return codons

	def _get_v_codons(self, vdj):
		vbr = vdj.v
		offset = (vdj.query_reading_frame * 2) % 3
		trunc_v = vbr.query_alignment[offset:]
		if not self.gapped:
			trunc_v.replace('-', '')
		codons = [trunc_v[i:i + 3] for i in range(0, len(trunc_v), 3)]
		if len(codons[-1]) != 3:
			codons = codons[:-1]
		return codons

	def _get_v_germ_codons(self, vdj):
		vbr = vdj.v
		offset = (vdj.query_reading_frame * 2) % 3
		trunc_v = vbr.germline_alignment[offset:]
		if not self.gapped:
			trunc_v.replace('-', '')
		codons = [trunc_v[i:i + 3] for i in range(0, len(trunc_v), 3)]
		if len(codons[-1]) != 3:
			codons = codons[:-1]
		return codons

	def _get_j_codons(self, vdj):
		jbr = vdj.j
		rf = jbr.germline_start % 3
		trunc_j = jbr.query_alignment[rf:]
		codons = [trunc_j[i:i + 3] for i in range(0, len(trunc_j), 3)]
		if len(codons[-1]) != 3:
			codons = codons[:-1]
		return codons

	def _get_j_germ_codons(self, vdj):
		jbr = vdj.j
		rf = jbr.germline_start % 3
		trunc_j = jbr.germline_alignment[rf:]
		codons = [trunc_j[i:i + 3] for i in range(0, len(trunc_j), 3)]
		if len(codons[-1]) != 3:
			codons = codons[:-1]
		return codons
