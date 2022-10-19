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


from __future__ import absolute_import, division, print_function, unicode_literals

import os
import sys
import glob
import subprocess as sp
from multiprocessing import cpu_count

from natsort import natsorted

from abutils.utils import log

logger = log.get_logger('basespace')


def list_files(d):
    return sorted([f for f in glob.glob(d + '/*') if os.path.isfile(f)])


def pair_files(files):
    return group_files(files)


def group_files(files, delim_count=3):
    groups = {}
    for f in files:
        f_prefix = '_'.join(os.path.basename(f).split('_')[:-delim_count])
        if f_prefix not in groups:
            groups[f_prefix] = []
        groups[f_prefix].append(f)
    return groups


def concat_lanes(lane_files, merged_file):
    with open(merged_file, 'w') as outfile:
        for fname in lane_files:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)
            outfile.write('\n')
    return merged_file


def pandaseq(f, r, o, algo):
    cmd = 'pandaseq -f {0} -r {1} -A {2} -d rbfkms -T {3} -w {4}'.format(f, r, algo, cpu_count(), o)
    sp.Popen(cmd, shell=True, stderr=sp.STDOUT, stdout=sp.PIPE).communicate()


def merge_reads(sample_name, sample_files, output, algo, i):
    if len(sample_files) > 2:
        subgroups = group_files(sample_files, delim_count=2)
        subgroups = {n: f for n, f in subgroups.items() if len(f) == 2}
        names = natsorted(subgroups.keys())
        file_pairs = [subgroups[n] for n in names]
    else:
        names = [sample_name, ]
        file_pairs = [sample_files, ]
    output_files = []
    for name, files in zip(names, file_pairs):
        files.sort()
        f = files[0]
        r = files[1]
        print_sample_info(i, name)
        o = os.path.join(output, f'{name}.fasta')
        pandaseq(f, r, o, algo)
        output_files.append(o)
    if len(output_files) > 1:
        merged_file = os.path.join(output, f'{sample_name}.fasta')
        o = concat_lanes(output_files, merged_file)
        for of in output_files:
            os.unlink(of)
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
    logger.info('[ {} ]  Merging sample {}'.format(str(i + 1), sample))


def print_sample_end():
    logger.info('Done.')


def run(input, output, algorithm='simple_bayesian', nextseq=False):
    '''
    Merge paired-end FASTQ files with PANDAseq. File name format must match that produced 
    by ``bcl2fastq``:

        ``<sample_name>_<sample_num>_<lane>_<read>_001.fastq.gz``

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
                3. a dict with sample names as keys and a list/tuple
                   containing paths to two paired read files as values

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


    Returns:

        list: a list of merged file paths
    '''
    print_start_info()
    if os.path.isdir(input):
        files = list_files(input)
        groups = group_files(files)
    elif type(input) in [dict, ]:
        if all([type(i) in [list, tuple] for i in input]) and all([len(i) == 2 for i in input]):
            files = [f for sublist in input for f in sublist]
            groups = {n: {'0': i} for n, i in zip(range(len(input)), input)}
        elif type(input) in [list, tuple] and all([os.path.isfile(i) for i in input]):
            files = input
            groups = group_files(files)
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
    for i, sample_name, sample_files in enumerate(natsorted(groups.keys())):
        mf = merge_reads(sample_name, sample_files, output, algorithm, i)
        merged_files.append(mf)
    return merged_files
