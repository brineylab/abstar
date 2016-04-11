#!/usr/bin/python
# filename: preprocess.py

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


from multiprocessing import cpu_count
import os
from subprocess import Popen, PIPE
import sys

from Bio import SeqIO

from utils.pandaseq import pair_files

from abtools.log import get_logger
from abtools.pipeline import list_files, make_dir


logger = get_logger('preprocess')


def quality_trim(input_directory=None, output_directory=None,
        quality_cutoff=15, length_cutoff=50, force_individual=False,
        quality_type='sanger', compress_output=True, file_pairs=None,
        singles_directory=None, nextseq=False, paired_reads=True,
        allow_5prime_trimming=False):
    if input_directory is None and any([file_pairs is None, output_directory is None]):
        err = '\nERROR: Either an input_directory must be provided or '
        err += 'both file_pairs and an output_directory must be provided.\n'
        print(err)
        sys.exit(1)
    if file_pairs:
        files = file_pairs
    else:
        input_directory = os.path.normpath(input_directory)
        if output_directory is None:
            oparent = os.path.dirname(input_directory)
            output_directory = os.path.join(oparent, 'quality_trimmed')
        make_dir(output_directory)
        if paired_reads:
            files = list_files(input_directory)
            file_pairs = pair_files(files, nextseq)
            files = file_pairs.values()
        else:
            files = [[f] for f in list_files(input_directory)]
    for f in files:

        logger.info(f)

        if len(f) == 2:
            paired_end = True
        elif len(f) == 1:
            paired_end = False
        else:
            err = 'ERROR: Each batch of files must contain either 1 (single-end reads) or '
            err += '2 (paired-end reads) files. This batch contains {} files:\n{}'.format(
                len(f), '\n'.join(f))
            err2 = 'If your filenames do not follow default the Illumina naming scheme, '
            err2 += 'you can force the processing of each file individually with the <force_individual> option.\n'
            err2 += 'Note that this can cause some problems for paired-end reads in which one read passes '
            err2 += 'passes trimming and the other does not.\n'
            err2 += 'If you have paired-end reads that do not follow the Illumina naming scheme, '
            err2 += 'you can pass pairs of filenames (a list of lists/tuples) with the <file_pairs> option. '
            err2 += 'If using <file_pairs>, the output directory must also be provided.'
            logger.info(err)
            logger.info(err2)
            continue
        f.sort()
        # set basic sickle cmd options
        sickle = 'sickle pe' if paired_end else 'sickle se'
        sickle += ' -t {}'.format(quality_type)
        sickle += ' -l {}'.format(length_cutoff)
        sickle += ' -q {}'.format(quality_cutoff)
        if compress_output:
            sickle += ' -g'
        if not allow_5prime_trimming:
            sickle += ' -x'
        # compute input/output filenames, add to sickle cmd
        sickle += ' -f {}'.format(f[0])
        o1_basename = os.path.basename(f[0]).rstrip('.gz')
        if compress_output:
            o1_basename += '.gz'
        sickle += ' -o {}'.format(os.path.join(output_directory, o1_basename))
        if paired_end:
            sickle += ' -r {}'.format(f[1])
            o2_basename = os.path.basename(f[1]).rstrip('.gz')
            if compress_output:
                o2_basename += '.gz'
            sickle += ' -p {}'.format(os.path.join(output_directory, o2_basename))
        # compute singles output filename, add to sickle cmd
        if singles_directory is not None:
            sfilename = '{}_{}_singles.fastq'.format(
                o1_basename.rstrip('.gz').rstrip('.fastq').rstrip('.fq'),
                o2_basename.rstrip('.gz').rstrip('.fastq').rstrip('.fq'))
            if compress_output:
                sfilename += '.gz'
            sickle += ' -s {}'.format(os.path.join(singles_directory, sfilename))
        else:
            sickle += ' -s /dev/null'
        # run sickle
        p = Popen(sickle, stdout=PIPE, stderr=PIPE, shell=True)
        stdout, stderr = p.communicate()
        logger.debug(stdout)
        logger.debug(stderr)
    return output_directory


def adapter_trim(input_directory, output_directory=None,
        adapter_5prime=None, adapter_3prime=None,
        adapter_5prime_anchored=None, adapter_3prime_anchored=None,
        adapter_both=None, compress_output=True):
    input_directory = os.path.normpath(input_directory)
    if output_directory is None:
        oparent = os.path.dirname(input_directory)
        output_directory = os.path.join(oparent, 'adapter_trimmed')
    make_dir(output_directory)
    files = list_files(input_directory)
    # parse adapter FASTA files, compile adapter option list
    adapters = []
    opts = ['-g', '-a', '-b']
    adapt_files = [adapter_5prime, adapter_3prime, adapter_both]
    for o, a in zip(opts, adapt_files):
        if a is None:
            continue
        adapts = [str(s.seq) for s in SeqIO.parse(open(a, 'r'), 'fasta')]
        adapters += [' '.join(z) for z in zip([o] * len(adapts), adapts)]
    if adapter_5prime_anchored is not None:
        adapts = ['^{}'.format(str(s.seq)) for s in SeqIO.parse(open(adapter_5prime_anchored, 'r'), 'fasta')]
        adapters += ['-g {}'.format(a) for a in adapts]
    if adapter_3prime_anchored is not None:
        adapts = ['{}$'.format(str(s.seq)) for s in SeqIO.parse(open(adapter_3prime_anchored, 'r'), 'fasta')]
        adapters += ['-a {}'.format(a) for a in adapts]
    # process input files
    for ifile in files:
        oname = os.path.basename(ifile).rstrip('.gz')
        if compress_output:
            oname += '.gz'
        ofile = os.path.join(output_directory, oname)
        # set up cutadapt command
        adapter_string = ' '.join(adapters)
        cutadapt = 'cutadapt -o {} {} {}'.format(ofile, adapter_string, ifile)



        print(cutadapt)



        # run cutadapt
        p = Popen(cutadapt, stdout=PIPE, stderr=PIPE, shell=True)
        stdout, stderr = p.communicate()
        logger.debug(stdout)
        logger.debug(stderr)
    return output_directory


def fastqc(input_directory, output_directory=None, threads=-1):
    '''
    Performs FASTQC analysis on raw NGS data.
    '''
    input_directory = os.path.normpath(input_directory)
    if output_directory is None:
        oparent = os.path.dirname(input_directory)
        output_directory = os.path.join(oparent, 'fastqc_reports')
    make_dir(output_directory)
    files = list_files(input_directory)
    if threads < 0:
        threads = cpu_count()
    fastqc_cmd = 'fastqc --noextract -o={} -t={} {}'.format(output_directory,
        threads, ' '.join(files))
    p = Popen(fastqc_cmd, stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = p.communicate()
    logger.debug(stdout)
    logger.debug(stderr)
    return output_directory
