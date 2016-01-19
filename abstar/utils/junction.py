#!/usr/bin/python
# filename: junction.py

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
import traceback

from abtools import log
from abtools.alignment import global_alignment, local_alignment


def get_junction(vdj, germ=False):
	logger = log.get_logger(__name__)
	try:
		return Junction(vdj, germ)
	except Exception, err:
		logger.debug('JUNCTION ERROR: {id}'.format(id=vdj.id))
		logger.debug(traceback.format_exc())


class Junction(object):
	"""docstring for Junction"""
	def __init__(self, vdj, germ):
		self.is_germ = germ
		self.junction_aa_start_pos = self._find_junction_aa_start(vdj)
		self.junction_aa_end_pos = self._find_junction_aa_end(vdj)
		self.junction_aa = self._get_junction_aa_sequence(vdj)
		self.junction_nt_start_pos = self._find_junction_nt_start(vdj)
		self.junction_nt_end_pos = self._find_junction_nt_end(vdj)
		self.junction_nt = self._get_junction_nt_sequence(vdj)
		self.cdr3_aa = self.junction_aa[1:-1]
		self.cdr3_nt = self.junction_nt[3:-3]
		if vdj.d:
			self.n1_nt = vdj.vdj_nt[vdj.n1_start:vdj.n1_end]
			self.d_nt = vdj.vdj_nt[vdj.d_start:vdj.d_end]
			self.n2_nt = vdj.vdj_nt[vdj.n2_start:vdj.n2_end]
			self.n_nt = None
			self.d_start_position_nt = self._get_d_start_position_nt(vdj)
			self.d_end_position_nt = self._get_d_end_position(vdj)
			self.d_dist_from_cdr3_start_nt = self.d_start_position_nt
			self.d_dist_from_cdr3_end_nt = self._get_d_dist_from_cdr3_end_nt(vdj)
		else:
			self.n1_nt = None
			self.d_nt = None
			self.n2_nt = None
			self.n_nt = vdj.vdj_nt[vdj.n_start:vdj.n_end]


	def _find_junction_aa_start(self, vdj):
		fr3_start = sum(vdj.v.regions.aa_lengths.values()) - vdj.v.regions.aa_lengths['FR3']
		alignments = []
		germ = self._get_v_chunk(vdj)
		for i in range(fr3_start, len(vdj.vdj_aa) - 4):
			chunk = vdj.vdj_aa[i:i + 6]
			alignments.append(global_alignment((str(i), chunk), ('germ', germ),
											   aa=True, matrix='blosum62',
											   gap_open=-100, gap_extend=-50))
		alignments.sort(key=lambda x: x.score, reverse=True)
		top = alignments[0]
		if vdj.vdj_aa[int(top.query.id) + 4] != 'C':
			start_pos = self._verify_junction_aa_start(vdj, alignments)
		else:
			start_pos = int(top.query.id) + 4
		return start_pos


	def _verify_junction_aa_start(self, vdj, alignments):
		if alignments[1].score / alignments[0].score >= 0.75:
			if int(alignments[1].query.id) < int(alignments[0].query.id) and vdj.vdj_aa[int(alignments[1].query.id) + 4] == 'C':
				return int(alignments[1].query.id) + 4
		return int(alignments[0].query.id) + 4


	def _find_junction_aa_end(self, vdj):
		alignments = []
		germ = self._get_j_chunk(vdj)
		for i in range(self.junction_aa_start_pos, len(vdj.vdj_aa) - 4):
			chunk = (i, vdj.vdj_aa[i:i + 5])
			alignments.append(global_alignment(chunk, germ, aa=True, matrix='blosum62'))
		alignments.sort(key=lambda x: x.score, reverse=True)
		return int(alignments[0].query.id) + 1


	def _get_junction_aa_sequence(self, vdj):
		start = self.junction_aa_start_pos
		end = self.junction_aa_end_pos
		if self.is_germ:
			return vdj.vdj_germ_aa[start:end]
		return vdj.vdj_aa[start:end]


	def _find_junction_nt_start(self, vdj):
		return self.junction_aa_start_pos * 3 + (vdj.query_reading_frame * 2) % 3


	def _find_junction_nt_end(self, vdj):
		return self.junction_nt_start_pos + len(self.junction_aa) * 3


	def _get_junction_nt_sequence(self, vdj):
		start = self.junction_nt_start_pos
		end = self.junction_nt_end_pos
		if self.is_germ:
			return vdj.vdj_germ_nt[start:end]
		return vdj.vdj_nt[start:end]


	def _get_d_start_position_nt(self, vdj):
		a = local_alignment(self.d_nt, self.cdr3_nt,
							gap_open_penalty=22, gap_extend_penalty=1)
		d_start = a.target_begin
		return d_start


	def _get_d_end_position(self, vdj):
		return self.d_start_position_nt + len(vdj.d.sequence)


	def _get_d_dist_from_cdr3_end_nt(self, vdj):
		return len(self.cdr3_nt) - (self.d_dist_from_cdr3_start_nt + len(self.d_nt))


	def _get_v_chunk(self, vdj):
		trunc_gene = vdj.v.top_germline[:5]
		germs = {
			'human':	{	'IGHV1': 'AVYYCA',
							'IGHV2': 'ATYYCA',
							'IGHV3': 'AVYYCA',
							'IGHV4': 'AVYYCA',
							'IGHV5': 'AMYYCA',
							'IGHV6': 'AVYYCA',
							'IGHV7': 'AVYYCA',
							'IGKV1': 'ATYYCQ',
							'IGKV2': 'GVYYCM',
							'IGKV3': 'AVYYCQ',
							'IGKV4': 'AVYYCQ',
							'IGLV1': 'ADYYCQ',
							'IGLV2': 'ADYYCS',
							'IGLV3': 'ADYYCQ',
							'IGLV4': 'ADYYCQ',
							'IGLV5': 'ADYYCM',
							'IGLV6': 'ADYYCQ',
							'IGLV7': 'AEYYCL',
							'IGLV8': 'SDYYCV',
							'IGLV9': 'SDYHCG'},

			'macaque':	{	'IGHV1': 'AVYYCA',
							'IGHV2': 'ATYYCA',
							'IGHV3': 'AVYYCA',
							'IGHV4': 'AVYYCA',
							'IGHV5': 'ATYYCA',
							'IGHV6': 'AVYYCA',
							'IGHV7': 'AVYYCA',
							'IGKV1': 'ATYYCQ',
							'IGKV2': 'GVYYCM',
							'IGKV3': 'AVYYCQ',
							'IGKV4': 'AVYYCQ',
							'IGKV5': 'AYYFCQ',
							'IGKV6': 'ATYYCQ',
							'IGKV7': 'ADYYCL',
							'IGLV1': 'ADYYCQ',
							'IGLV2': 'ADYYCS',
							'IGLV3': 'ADYYCQ',
							'IGLV4': 'ADYYCQ',
							'IGLV5': 'ADYYCM',
							'IGLV6': 'ADYYCQ',
							'IGLV7': 'AEYYCW',
							'IGLV8': 'SDYYCT',
							'IGLV9': 'SDYHCG'
							},
			
			'mouse':	{	'IGHV1': 'AVYYCA',
							'IGHV2': 'AIYYCA',
							'IGHV3': 'ATYYCA',
							'IGHV4': 'ALYYCA',
							'IGHV5': 'AMYYCA',
							'IGHV6': 'GIYYCT',
							'IGHV7': 'ATYYCA',
							'IGHV8': 'ATYYCA',
							'IGHV9': 'ATYFCA',
							'IGKV1': 'ATYYCQ',
							'IGKV2': 'GVYYCA',
							'IGKV3': 'ATYYCQ',
							'IGKV4': 'ATYYCQ',
							'IGKV5': 'GVYYCQ',
							'IGKV6': 'AVYFCQ',
							'IGKV7': 'THYYCA',
							'IGKV8': 'AVYYCQ',
							'IGKV9': 'ADYYCL',
							'IGLV1': 'AIYFCA',
							'IGLV2': 'AMYFCA',
							'IGLV3': 'AIYICG',
							'IGLV4': 'AIYFCA',
							'IGLV5': 'AIYFCA',
							'IGLV6': 'AIYFCA',
							'IGLV7': 'AIYFCA',
							'IGLV8': 'AIYFCA'},
			
			'rabbit':	{	'IGHV1': 'ATYFCA',
							'IGKV1': 'ATYYCQ',
							'IGLV2': 'AAYYCV',
							'IGLV3': 'ADYYCQ',
							'IGLV4': 'ADYYCG',
							'IGLV5': 'ADYYCA',
							'IGLV6': 'ADYYCQ'},
			
			'b12mouse':	{	'IGHV1': 'AVYYCA',
							'IGHV2': 'AIYYCA',
							'IGHV3': 'ATYYCA',
							'IGHV4': 'ALYYCA',
							'IGHV5': 'AMYYCA',
							'IGHV6': 'GIYYCT',
							'IGHV7': 'ATYYCA',
							'IGHV8': 'ATYYCA',
							'IGHV9': 'ATYFCA',
							'IGKV1': 'ATYYCQ',
							'IGKV2': 'GVYYCA',
							'IGKV3': 'ATYYCQ',
							'IGKV4': 'ATYYCQ',
							'IGKV5': 'GVYYCQ',
							'IGKV6': 'AVYFCQ',
							'IGKV7': 'THYYCA',
							'IGKV8': 'AVYYCQ',
							'IGKV9': 'ADYYCL',
							'IGLV1': 'AIYFCA',
							'IGLV2': 'AMYFCA',
							'IGLV3': 'AIYICG',
							'IGLV4': 'AIYFCA',
							'IGLV5': 'AIYFCA',
							'IGLV6': 'AIYFCA',
							'IGLV7': 'AIYFCA',
							'IGLV8': 'AIYFCA'},
			
			'vrc01mouse':{	'IGHV1': 'AVYYCA',
							'IGHV2': 'AIYYCA',
							'IGHV3': 'ATYYCA',
							'IGHV4': 'ALYYCA',
							'IGHV5': 'AMYYCA',
							'IGHV6': 'GIYYCT',
							'IGHV7': 'ATYYCA',
							'IGHV8': 'ATYYCA',
							'IGHV9': 'ATYFCA',
							'IGKV1': 'ATYYCQ',
							'IGKV2': 'GVYYCA',
							'IGKV3': 'ATYYCQ',
							'IGKV4': 'ATYYCQ',
							'IGKV5': 'GVYYCQ',
							'IGKV6': 'AVYFCQ',
							'IGKV7': 'THYYCA',
							'IGKV8': 'AVYYCQ',
							'IGKV9': 'ADYYCL',
							'IGLV1': 'AIYFCA',
							'IGLV2': 'AMYFCA',
							'IGLV3': 'AIYICG',
							'IGLV4': 'AIYFCA',
							'IGLV5': 'AIYFCA',
							'IGLV6': 'AIYFCA',
							'IGLV7': 'AIYFCA',
							'IGLV8': 'AIYFCA'},

			'9114mouse':{	'IGHV1': 'AVYYCA',
							'IGHV2': 'AIYYCA',
							'IGHV3': 'ATYYCA',
							'IGHV4': 'ALYYCA',
							'IGHV5': 'AMYYCA',
							'IGHV6': 'GIYYCT',
							'IGHV7': 'ATYYCA',
							'IGHV8': 'ATYYCA',
							'IGHV9': 'ATYFCA',
							'IGKV1': 'ATYYCQ',
							'IGKV2': 'GVYYCA',
							'IGKV3': 'ATYYCQ',
							'IGKV4': 'ATYYCQ',
							'IGKV5': 'GVYYCQ',
							'IGKV6': 'AVYFCQ',
							'IGKV7': 'THYYCA',
							'IGKV8': 'AVYYCQ',
							'IGKV9': 'ADYYCL',
							'IGLV1': 'AIYFCA',
							'IGLV2': 'AMYFCA',
							'IGLV3': 'AIYICG',
							'IGLV4': 'AIYFCA',
							'IGLV5': 'AIYFCA',
							'IGLV6': 'AIYFCA',
							'IGLV7': 'AIYFCA',
							'IGLV8': 'AIYFCA'}
			}

		return germs[vdj.species][trunc_gene]


	def _get_j_chunk(self, vdj):
		trunc_gene = vdj.j.top_germline[:5]
		germs = {
			'human':	{	'IGHJ1': 'WGQGT',
							'IGHJ2': 'WGRGT',
							'IGHJ3': 'WGQGT',
							'IGHJ4': 'WGQGT',
							'IGHJ5': 'WGQGT',
							'IGHJ6': 'WGQGT',
							'IGKJ1': 'FGQGT',
							'IGKJ2': 'FGQGT',
							'IGKJ3': 'FGPGT',
							'IGKJ4': 'FGGGT',
							'IGKJ5': 'FGQGT',
							'IGLJ1': 'FGTGT',
							'IGLJ2': 'FGGGT',
							'IGLJ3': 'FGGGT',
							'IGLJ4': 'FGGGT',
							'IGLJ5': 'FGEGT',
							'IGLJ6': 'FGSGT',
							'IGLJ7': 'FGGGT'},

			'macaque':	{	'IGHJ1': 'WGQGA',
							'IGHJ2': 'WGPGT',
							'IGHJ3': 'WGQGL',
							'IGHJ4': 'WGQGV',
							'IGHJ5': 'WGPGV',
							'IGHJ6': 'WGQGV',
							'IGKJ1': 'FGQGT',
							'IGKJ2': 'FGQGT',
							'IGKJ3': 'FGPGT',
							'IGKJ4': 'FGGGT',
							'IGKJ5': 'FGQGT',
							'IGLJ1': 'FGAGT',
							'IGLJ2': 'FGGGT',
							'IGLJ3': 'FGGGT',
							'IGLJ4': 'FCGGT',
							'IGLJ5': 'FGEGT',
							'IGLJ6': 'FGSGT'},

			'mouse':	{	'IGHJ1': 'WGAGT',
							'IGHJ2': 'WGQGT',
							'IGHJ3': 'WGQGT',
							'IGHJ4': 'WGQGT',
							'IGKJ1': 'FGGGT',
							'IGKJ2': 'FGGGT',
							'IGKJ4': 'FGSGT',
							'IGKJ5': 'FGAGT',
							'IGLJ1': 'FGGGT',
							'IGLJ2': 'FGGGT',
							'IGLJ3': 'FGSGT'},

			'rabbit':	{	'IGHJ1': 'WGTGT',
							'IGHJ2': 'WGPGT',
							'IGHJ3': 'WGQGT',
							'IGHJ4': 'WGPGT',
							'IGHJ5': 'WGQGT',
							'IGHJ6': 'WGPGT',
							'IGKJ1': 'FGGGT',
							'IGKJ1': 'FGAGT',
							'IGLJ5': 'FGGGT',
							'IGLJ6': 'FGGGT'},

			'b12mouse':	{	'IGHJ1': 'WGAGT',
							'IGHJ2': 'WGQGT',
							'IGHJ3': 'WGQGT',
							'IGHJ4': 'WGQGT',
							'IGHJ6': 'WGQGT',
							'IGKJ1': 'FGGGT',
							'IGKJ2': 'FGGGT',
							'IGKJ4': 'FGSGT',
							'IGKJ5': 'FGAGT',
							'IGLJ1': 'FGGGT',
							'IGLJ2': 'FGGGT',
							'IGLJ3': 'FGSGT'},

			'vrc01mouse':{	'IGHJ1': 'WGAGT',
							'IGHJ2': 'WGQGT',
							'IGHJ3': 'WGQGT',
							'IGHJ4': 'WGQGT',
							'IGHJ6': 'WGQGT',
							'IGKJ1': 'FGGGT',
							'IGKJ2': 'FGGGT',
							'IGKJ4': 'FGSGT',
							'IGKJ5': 'FGAGT',
							'IGLJ1': 'FGGGT',
							'IGLJ2': 'FGGGT',
							'IGLJ3': 'FGSGT'},

			'9114mouse':{	'IGHJ1': 'WGAGT',
							'IGHJ2': 'WGQGT',
							'IGHJ3': 'WGQGT',
							'IGHJ4': 'WGQGT',
							'IGHJ6': 'WGQGT',
							'IGKJ1': 'FGGGT',
							'IGKJ2': 'FGGGT',
							'IGKJ4': 'FGSGT',
							'IGKJ5': 'FGAGT',
							'IGLJ1': 'FGGGT',
							'IGLJ2': 'FGGGT',
							'IGLJ3': 'FGSGT'}
			}

		return germs[vdj.species][trunc_gene]
