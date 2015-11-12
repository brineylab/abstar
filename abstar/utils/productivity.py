#!/usr/bin/python
# filename: productivity.py

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


def check_productivity(vdj):
	try:
		problems = 0
		problems += stop_codons(vdj)
		problems += ambig_codons(vdj)
		problems += vdj_agreement(vdj)
		problems += conserved_junc_residues(vdj)
		problems += rearrangement(vdj)
		if problems:
			return 'no'
		return 'yes'
	except:
		logging.debug('PRODUCTIVITY CHECK ERROR: {}'.format(vdj.id))


def stop_codons(vdj):
	if '*' in vdj.vdj_aa:
		return 1
	return 0


def ambig_codons(vdj):
	if 'X' in vdj.vdj_aa.upper():
		return 1
	return 0


def vdj_agreement(vdj):
	if vdj.v.top_germline[:3] != vdj.j.top_germline[:3]:
		return 1
	if vdj.d:
		if not vdj.d.top_germline:
			return 0
		if vdj.v.top_germline[:3] != vdj.d.top_germline[:3]:
			return 1
	return 0


def conserved_junc_residues(vdj):
	if not vdj.junction:
		return 1
	start = 'C'
	end = 'W' if vdj.chain == 'heavy' else 'F'
	if vdj.junction.junction_aa[0] != start or vdj.junction.junction_aa[-1] != end:
		return 1
	return 0


def rearrangement(vdj):
	if not vdj.rearrangement:
		return 1
	return 0
