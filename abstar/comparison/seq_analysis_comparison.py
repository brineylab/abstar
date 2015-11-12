#!/usr/bin/env python
# filename: seq_analysis_comparison.py

from __future__ import print_function

import argparse
from pymongo import MongoClient


parser = argparse.ArgumentParser("")
parser.add_argument('-a', '--abblast', dest='abblast_input',required=True, help="The AbBlast input file, tab-delimited. Required.")
parser.add_argument('-g', '--igblast', dest='igblast_input',required=True, help="The IgBLAST input file, tab-delimited. Required.")
parser.add_argument('-i', '--imgt', dest='imgt_input',required=True, help="The IMGT input file, tab-delimited. Required.")
args = parser.parse_args()


class Seq(object):
	"""docstring for Seqs"""
	def __init__(self, seq_id):
		self.seq_id = seq_id

	def add_abblast(self, data):
		self.abblast_v = data[0]
		self.abblast_d = data[1]
		self.abblast_j = data[2]

	def add_igblast(self, data):
		self.igblast_v = data[0]
		self.igblast_d = data[1]
		self.igblast_j = data[2]

	def add_imgt(self, data):
		try:
			self.imgt_v = data[0].split()[1]
		except IndexError:
			self.imgt_v = ''
		try:
			self.imgt_d = data[1].split()[1]
		except IndexError:
			self.imgt_d = ''
		try:
			self.imgt_j = data[2].split()[1]
		except IndexError:
			self.imgt_j = ''


def compare_to_imgt(seqs):
	print('')
	print('===================')
	print('AbBlast versus IMGT')
	print('===================')
	print('')
	# Variable
	print('Variable mismatches')
	print('')
	v_matches = 0
	v_mismatches = 0
	v_missing = 0
	for s in seqs:
		try:
			if s.abblast_v.split('*')[0] == s.imgt_v.split('*')[0]:
				v_matches += 1
			elif s.abblast_v.split('*')[0] == 'IGHV4-b' and s.imgt_v.split('*')[0] == 'IGHV4-38-2':
				v_matches += 1
			elif s.abblast_v.split('*')[0] == 'IGHV2-70' and s.imgt_v.split('*')[0] == 'IGHV2-70D':
				v_matches += 1
			elif s.abblast_v.split('*')[0] == '' or s.imgt_v.split('*')[0] == '':
				v_missing += 1
			else:
				v_mismatches += 1
				# print(s.seq_id, s.abblast_v.split('*')[0], s.imgt_v.split('*')[0])
		except AttributeError:
			v_missing += 1
	# Diversity
	print('Diversity mismatches')
	print('')
	d_matches = 0
	d_mismatches = 0
	d_missing = 0
	for s in seqs:
		try:
			if s.abblast_d.split('*')[0] == s.imgt_d.split('*')[0]:
				d_matches += 1
			else:
				d_mismatches += 1
				# print(s.seq_id, s.abblast_d.split('*')[0], s.imgt_d.split('*')[0])
		except AttributeError:
			d_missing += 1
	# Joining
	print('Joining mismatches')
	print('')
	j_matches = 0
	j_mismatches = 0
	j_missing = 0
	for s in seqs:
		try:
			if s.abblast_j.split('*')[0] == s.imgt_j.split('*')[0]:
				j_matches += 1
			else:
				j_mismatches += 1
				# print(s.seq_id, s.abblast_j.split('*')[0], s.imgt_j.split('*')[0])
		except AttributeError:
			j_missing += 1
	print_results((v_matches, v_mismatches, v_missing, 
				  d_matches, d_mismatches, d_missing, 
				  j_matches, j_mismatches, j_missing))


def compare_to_igblast(seqs):
	print('\n')
	print('======================')
	print('AbBlast versus IgBLAST')
	print('======================')
	print('')
	# Variable
	v_matches = 0
	v_mismatches = 0
	v_missing = 0
	for s in seqs:
		try:
			if s.abblast_v.split('*')[0] == s.igblast_v.split('*')[0]:
				v_matches += 1
			else:
				v_mismatches += 1
		except AttributeError:
			v_missing += 1
	# Diversity
	d_matches = 0
	d_mismatches = 0
	d_missing = 0
	for s in seqs:
		try:
			if s.abblast_d.split('*')[0] == s.igblast_d.split('*')[0]:
				d_matches += 1
			else:
				d_mismatches += 1
		except AttributeError:
			d_missing += 1
	# Joining
	j_matches = 0
	j_mismatches = 0
	j_missing = 0
	for s in seqs:
		try:
			if s.abblast_j.split('*')[0] == s.igblast_j.split('*')[0]:
				j_matches += 1
			else:
				j_mismatches += 1
		except AttributeError:
			j_missing += 1
	print_results((v_matches, v_mismatches, v_missing, 
				  d_matches, d_mismatches, d_missing, 
				  j_matches, j_mismatches, j_missing))

def compare_igblast_to_imgt(seqs):
	print('\n')
	print('===================')
	print('IgBLAST versus IMGT')
	print('===================')
	print('')
	# Variable
	v_matches = 0
	v_mismatches = 0
	v_missing = 0
	for s in seqs:
		try:
			if s.igblast_v.split('*')[0] == s.imgt_v.split('*')[0]:
				v_matches += 1
			elif s.igblast_v.split('*')[0] == 'IGHV4-b' and s.imgt_v.split('*')[0] == 'IGHV4-38-2':
				v_matches += 1
			else:
				v_mismatches += 1
				# print(s.seq_id, s.igblast_v.split('*')[0], s.imgt_v.split('*')[0])
		except AttributeError:
			v_missing += 1
	# Diversity
	d_matches = 0
	d_mismatches = 0
	d_missing = 0
	for s in seqs:
		try:
			if s.igblast_d.split('*')[0] == s.imgt_d.split('*')[0]:
				d_matches += 1
			else:
				d_mismatches += 1
		except AttributeError:
			d_missing += 1
	# Joining
	j_matches = 0
	j_mismatches = 0
	j_missing = 0
	for s in seqs:
		try:
			if s.igblast_j.split('*')[0] == s.imgt_j.split('*')[0]:
				j_matches += 1
			else:
				j_mismatches += 1
		except AttributeError:
			j_missing += 1
	print_results((v_matches, v_mismatches, v_missing, 
				  d_matches, d_mismatches, d_missing, 
				  j_matches, j_mismatches, j_missing))


def print_results(data):
	print('Variable matches: {}'.format(data[0]))
	print('Variable mismatches: {}'.format(data[1]))
	print('Variable missing: {}'.format(data[2]))
	print('Variable accuracy: {}'.format(100. * data[0] / (data[0] + data[1])))
	print('')
	print('Diversity matches: {}'.format(data[3]))
	print('Diversity mismatches: {}'.format(data[4]))
	print('Diversity missing: {}'.format(data[5]))
	print('Diversity accuracy: {}'.format(100. * data[3] / (data[3] + data[4])))
	print('')
	print('Joining matches: {}'.format(data[6]))
	print('Joining mismatches: {}'.format(data[7]))
	print('Joining missing: {}'.format(data[8]))
	print('Joining accuracy: {}'.format(100. * data[6] / (data[6] + data[7])))
	print('')


def parse_input():
	igblast = []
	abblast = []
	imgt = []
	seq_ids = []
	for line in open(args.abblast_input, 'ru'):
		abblast.append(line.split('\t'))
	for line in open(args.igblast_input, 'ru'):
		igblast.append(line.split('\t'))
	imgt_lines = open(args.imgt_input, 'r').read().split('\n')
	print(len(imgt_lines))
	for m in imgt_lines:
		imgt.append(m.split('\t'))
	for a in abblast:
		if a[0] not in seq_ids:
			seq_ids.append(a[0])
	for i in igblast:
		if i[0] not in seq_ids:
			seq_ids.append(i[0])
	for m in imgt:
		if m[0] not in seq_ids:
			seq_ids.append(m[0])
	return seq_ids, abblast, igblast, imgt



def main():
	seq_ids, abblast, igblast, imgt = parse_input()
	seqs = {s : Seq(s) for s in seq_ids}
	for a in abblast:
		seqs[a[0]].add_abblast(a[1:])
	for i in igblast:
		# print(i[0])
		seqs[i[0]].add_igblast(i[1:])
	for m in imgt:
		seqs[m[0]].add_imgt(m[1:])

	compare_to_igblast(seqs.values())
	compare_to_imgt(seqs.values())
	compare_igblast_to_imgt(seqs.values())


if __name__ == '__main__':
	main()

