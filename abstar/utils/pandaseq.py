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

from abtools import log

logger = log.get_logger('basespace')


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


def merge_reads(files, output, algo, nextseq, i):
    files.sort()
    f = files[0]
    r = files[1]
    if nextseq:
        lane = os.path.basename(f).split('_')[-3]
        sample_id = os.path.basename(f).split('_')[0]
        sample = sample_id + '_' + lane
    else:
        sample = os.path.basename(f).split('_')[0]
    print_sample_info(i, sample)
    o = os.path.join(output, '{}.fasta'.format(sample))
    batch_pandaseq(f, r, o, algo)
    return o


def print_start_info():
    logger.info('')
    logger.info('')
    logger.info('========================================')
    logger.info('Merging reads with PANDAseq')
    logger.info('========================================')
    logger.info('')


def print_input_info(files):
    logger.info('The input directory contains {} pair(s) of files to be merged.\n'.format(len(files) / 2))


def print_sample_info(i, sample):
    logger.info('[ {} ]  Processing sample {}'.format(str(i + 1), sample))


def print_sample_end():
    logger.info('Done.')


def run(input, output, algorithm='simple_bayesian', nextseq=False):
    '''
    Merge paired-end FASTQ files with PANDAseq.

    Examples:

        To merge a directory of raw (gzip compressed) files from a MiSeq run::

            merged_files = run('/path/to/input', '/path/to/output')

        Same as above, but using the Pear_ read merging algorithm::

            merged_files = run('/path/to/input', '/path/to/output', algorithm='pear')

        To merge a list of file pairs::

            file_pairs = [(sample1_R1.fastq, sample1_R2.fastq),
                          (sample2_R1.fastq.gz, sample2_R2.fastq.gz),
                          (sample3_R1.fastq, sample3_R2.fastq)]
            merged_files = run(file_pairs, '/path/to/output')

        .. _Pear: http://sco.h-its.org/exelixis/web/software/pear/


    Args:

        input (str, list): Input can be one of three things:

                1. path to a directory of paired FASTQ files
                2. a list of paired FASTQ files
                3. a list of read pairs, with each read pair being a list/tuple
                   containing paths to two paired read files

            Regardless of what input type is provided, paired FASTQ files can be either
            gzip compressed or uncompressed.

            When providing a list of files or a directory of files, it is assumed that
            all files follow Illumina naming conventions. If your file names aren't
            Illumina-like, submit your files as a list of read pairs to ensure that
            the proper pairs of files are merged.

        output (str): Path to an output directory, into which merged FASTQ files will be deposited.
            To determine the filename for the merged file, the R1 file (or the first file in the read
            pair) is split at the first occurance of the '_' character. Therefore, the read pair
            ``['my-sequences_R1.fastq', 'my-sequences_R2.fastq']`` would be merged into ``my-sequences.fasta``.

        algorithm (str): PANDAseq algorithm to be used for merging reads. Choices are: 'simple_bayesian',
            'ea_util', 'flash', 'pear', 'rdp_mle', 'stitch', or 'uparse'. Default is 'simple_bayesian',
            which is the default PANDAseq algorithm.

        nextseq (bool): Set to ``True`` if the sequencing data was generated on a NextSeq. Needed
            because the naming conventions for NextSeq output files differs from MiSeq output.


    Returns:

        list: a list of merged file paths
    '''
    print_start_info()
    if os.path.isdir(input):
        files = list_files(input)
        pairs = pair_files(files, nextseq)
    elif type(input) in [list, tuple]:
        if all([type(i) in [list, tuple] for i in input]) and all([len(i) == 2 for i in input]):
            files = [f for sublist in input for f in sublist]
            pairs = {n: i for n, i in zip(range(len(input)), input)}
        elif all([os.path.isfile(i) for i in input]):
            files = input
            pairs = pair_files(files, nextseq)
        else:
            err = 'ERROR: Invalid input. Input may be one of three things:\n'
            err += '  1. a directory path\n'
            err += '  2. a list of file paths\n'
            err += '  3. a list of file pairs (lists/tuples containing exactly 2 file paths)'
            raise RuntimeError(err)
    else:
        err = 'ERROR: Invalid input. Input may be one of three things:\n'
        err += '  1. a directory path\n'
        err += '  2. a list of file paths\n'
        err += '  3. a list of file pairs (lists/tuples containing exactly 2 file paths)'
        raise RuntimeError(err)
    print_input_info(files)
    merged_files = []
    for i, pair in enumerate(sorted(pairs.keys())):
        if len(pairs[pair]) == 2:
            logger.info('Merging {} and {}'.format(pairs[pair][0], pairs[pair][1]))
            mf = merge_reads(pairs[pair], output, algorithm, nextseq, i)
            merged_files.append(mf)
    return merged_files
