#!/usr/bin/python
# filename: pandaseq.py

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
import subprocess as sp
from multiprocessing import cpu_count


def list_files(d):
    return sorted([f for f in glob.glob(d + '/*') if os.path.isfile(f)])


def pair_files(files, nextseq):
    pairs = {}
    for f in files:
        if nextseq:
            f_prefix = '_'.join(os.path.basename(f).split('_')[:3])
        else:
            f_prefix = '_'.join(os.path.basename(f).split('_')[:2])
        if f_prefix in pairs:
            pairs[f_prefix].append(f)
        else:
            pairs[f_prefix] = [f, ]
    return pairs


def batch_pandaseq(f, r, o, algo):
    cmd = 'pandaseq -f {0} -r {1} -A {2} -d rbfkms -T {3} -w {4}'.format(f, r, algo, cpu_count(), o)
    sp.Popen(cmd, shell=True, stderr=sp.STDOUT, stdout=sp.PIPE).communicate()


def merge_reads(files, output, args, i):
    files.sort()
    f = files[0]
    r = files[1]
    if args.nextseq:
        lane = os.path.basename(f).split('_')[-3]
        sample_id = os.path.basename(f).split('_')[0]
        sample = sample_id + '_' + lane
    else:
        sample = os.path.basename(f).split('_')[0]
    print_sample_info(i, sample)
    o = os.path.join(output, '{}.fasta'.format(sample))
    batch_pandaseq(f, r, o, args.pandaseq_algo)
    # print_sample_end(log)


def print_start_info():
    sys.stdout.write('\n')
    sys.stdout.write('\n')
    sys.stdout.write('========================================\n')
    sys.stdout.write('Merging reads with PANDAseq\n')
    sys.stdout.write('========================================\n')
    sys.stdout.write('\n')


def print_input_info(files):
    sys.stdout.write('The input directory contains {} pair(s) of files to be merged.\n\n'.format(len(files) / 2))


def print_sample_info(i, sample):
    sys.stdout.write('[ {} ]  Processing sample {}\n'.format(str(i), sample))


def print_sample_end():
    sys.stdout.write('Done.\n')


def run(input, output, args, log=False):
    if log:
        import logging
    try:
        nextseq = args.nextseq
    except:
        args.nextseq = False
    print_start_info()
    files = list_files(input)
    print_input_info(files)
    pairs = pair_files(files, args.nextseq)
    for i, pair in enumerate(sorted(pairs.keys())):
        if len(pairs[pair]) == 2:
            merge_reads(pairs[pair], output, args, i)
            if log:
                file1 = os.path.basename(pairs[pair][0])
                file2 = os.path.basename(pairs[pair][1])
                logging.info('Merging {} and {}'.format(file1, file2))
