#!/usr/bin/python
# filename: batch_merge.py

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
import sys
import glob
import shutil
import argparse

import pandaseq

parser = argparse.ArgumentParser("Parses the output of IgBLAST into something suitable for import into a MySQL database")
parser.add_argument('-i', '--in',
					dest='input',
					required=True,
					help="The input directory, containing paired FASTQ files (uncompressed or gzip compressed). Required.")
parser.add_argument('-o', '--out',
					dest='output',
					required=True,
					help="The output directory, will contain merged FASTA files. Required.")
parser.add_argument('-p', '--pandaseq_algo',
					dest="pandaseq_algo",
					default='simple_bayesian',
					choices=['simple_bayesian', 'ea_util', 'flash', 'pear', 'rdp_mle', 'stitch', 'uparse'],
					help="Define merging algorithm to be used by PANDAseq.\
					Options are 'simple_bayesian', 'ea_util', 'flash', 'pear', 'rdp_mle', 'stitch', or 'uparse'.\
					Default is 'simple_bayesian', which is the default PANDAseq algorithm.")
parser.add_argument('-n', '--nextseq',
					dest='nextseq',
					default=False,
					action='store_true',
					help="Use flag if run was performed on a NextSeq sequencer.")
args = parser.parse_args()


def make_direc(d):
	if not os.path.exists(d):
		os.mkdir(d)


def remove_direc(d):
	shutil.rmtree(d)


def list_files(d):
	return sorted([f for f in glob.glob(d + '/*') if os.path.isfile(f)])


def bin_files(files):
	file_bins = {}
	for f in files:
		f_pre = '_'.join(os.path.basename(f).split('_')[:-1])
		if f_pre in file_bins:
			file_bins[f_pre].append(f)
		else:
			file_bins[f_pre] = [f, ]
	return file_bins


def concat(d):
	files = list_files(d)
	file_bins = bin_files(files)
	for bin in file_bins:
		outfile = os.path.join(args.output, '{}.fasta'.format(bin))
		with open(outfile, 'w') as o:
			for f in file_bins[bin]:
				with open(f) as i:
					for line in i:
						o.write(line)


def main():
	make_direc(args.output)
	if args.nextseq:
		temp = os.path.join(args.output, 'temp')
		make_direc(temp)
		o = temp
	else:
		o = args.output
	pandaseq.run(args.input, o, args)
	if args.nextseq:
		print ''
		print 'Concatenating NextSeq lane files for each sample...'
		concat(o)
		remove_direc(o)
		print 'Done.'
		print ''



if __name__ == '__main__':
	main()
