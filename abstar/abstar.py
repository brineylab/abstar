#!/usr/bin/env python
# filename: abstar.py

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


from __future__ import print_function, unicode_literals

from argparse import ArgumentParser
from glob import glob
from multiprocessing import Pool
from subprocess import Popen, PIPE

# External
from Bio import SeqIO
from abtools import log
from abtools.sequence import Sequence

import gzip
import logging
import os
import re
import sys
import time
import traceback
import warnings


#####################################################################
#
#                             ARGUMENTS
#
#####################################################################


def parse_arguments(print_help=False):
    parser = ArgumentParser("Performs germline assignment and other relevant annotation on antibody sequence data from NGS platforms.")
    parser.add_argument('-p', '--project', dest='data_dir', default=None,
                        help="The data directory, where files will be downloaded (or have previously \
                        been download), temp files will be stored, and output files will be \
                        written. During StarCluster configuration, all ephemeral drives on the master \
                        node will be combined into a single RAID volume and mounted at the data directory. \
                        Not necessary if '-o', '-t' and '-i' are provided. \
                        Default is '/data'.")
    parser.add_argument('-i', '--in', dest='input', default=None,
                        help="The input file or directory, to be split and processed in parallel. \
                        If a directory is given, all files in the directory will be iteratively processed.")
    parser.add_argument('-o', '--out', dest='output', default=None,
                        help="The output directory, into which the JSON-formatted output files will be deposited. \
                        If the directory does not exist, it will be created.")
    parser.add_argument('-l', '--log', dest='log', default=None,
                        help="The log file, to which log info will be written. \
                        Default is <output_directory>/abstar.log if '-o/--out' is specificied. \
                        If '-o/--out' isn't provided, default is <data_directory>/abstar.log")
    parser.add_argument('-t', '--temp', dest='temp', default=None,
                        help="The directory in which temp files will be stored. If the directory doesn't exist, \
                        it will be created. Required.")
    parser.add_argument('--sequences', dest='sequences', default=None,
                        help="Only used when passing sequences directly through the API.")
    parser.add_argument('--use-test-data', dest='use_test_data', default=False, action='store_true',
                        help="If set, and if '-i/--input' isn't provided, AbStar's included test data \
                        will be used as input.")
    parser.add_argument('-k', '--chunksize', dest='chunksize', default=250, type=int,
                        help="Approximate number of sequences in each distributed job. \
                        Defaults to 250. \
                        Set to 0 if you want file splitting to be turned off \
                        Don't change unless you know what you're doing.")
    parser.add_argument('-T', '--output_type', dest="output_type", choices=['json', 'imgt', 'hadoop'], default='json',
                        help="Select the output type. Options are 'json', 'imgt' and 'impala'. \
                        IMGT output mimics the Summary table produced by IMGT High-V/Quest, \
                        to maintain some level of compatibility with existing IMGT-based pipelines. \
                        JSON output is much more detailed. \
                        Hadoop output is columnar and easily converted to binary HDFS-friendly formats \
                        (Parquet, Avro) for use in Impala or other Hadoop query engines (Pig, Hive, Spark). \
                        Defaults to JSON output.")
    parser.add_argument('-m', '--merge', dest="merge", action='store_true', default=False,
                        help="Use if the input files are paired-end FASTQs \
                        (either gzip compressed or uncompressed) from Illumina platforms. \
                        Prior to running the germline assignment pipeline, \
                        paired reads will be merged with PANDAseq.")
    parser.add_argument('-P', '--pandaseq_algo', dest="pandaseq_algo", default='simple_bayesian',
                        choices=['simple_bayesian', 'ea_util', 'flash', 'pear', 'rdp_mle', 'stitch', 'uparse'],
                        help="Define merging algorithm to be used by PANDAseq.\
                        Options are 'simple_bayesian', 'ea_util', 'flash', 'pear', 'rdp_mle', 'stitch', or 'uparse'.\
                        Default is 'simple_bayesian', which is the default PANDAseq algorithm.")
    parser.add_argument('-n', '--nextseq', dest="nextseq", action='store_true', default=False,
                        help="Use if the run was performed on a NextSeq sequencer.")
    parser.add_argument('-u', '--umid', dest="uaid", type=int, default=0,
                        help="Length of the unique antibody identifiers (UAIDs) \
                        used when preparing samples, if used. \
                        Default is unbarcoded (UAID length of 0).")
    parser.add_argument('-I', '--isotype', dest="isotype", action='store_false', default=True,
                        help="If set, the isotype will not be determined for heavy chains.\
                        If not set, isotyping sequences for the appropriate species will be used.")
    parser.add_argument('-b', '--basespace', dest="basespace", default=False, action='store_true',
                        help="Use if files should be downloaded directly from BaseSpace. \
                        Files will be downloaded into the input directory.")
    parser.add_argument('-c', '--cluster', dest="cluster", default=False, action='store_true',
                        help="Use if performing computation on a Celery cluster. \
                        If set, input files will be split into many subfiles and passed \
                        to a Celery queue. If not set, input files will still be split, but \
                        will be distributed to local processors using multiprocessing.")
    parser.add_argument('-S', '--starcluster', dest="starcluster", default=False, action='store_true',
                        help="Use if performing analysis on a StarCluster instance. \
                        If set, the cluster will be configured to NFS share all ephemeral drives \
                        on the master node and Celery workers will be started on all worker nodes. \
                        Configuration only needs to be run once per cluster, so additional runs on\
                        an already-configured cluster should be run without this option.")
    parser.add_argument('-D', '--debug', dest="debug", action='count', default=0,
                        help="If set, logs additional debug information. \
                        Use -DD to print verbose exception information to screen in addition to writing to log.")
    parser.add_argument('-s', '--species', dest='species', default='human',
                        choices=['human', 'macaque', 'mouse', 'rabbit', 'b12mouse', 'vrc01mouse', '9114mouse'])
    parser.add_argument('-z', '--gzip', dest='gzip', default=False, action='store_true',
                        help="Compress the output to a gzipped file")
    parser.add_argument('--pretty', dest='pretty', default=False, action='store_true',
                        help='Pretty format json file')
    parser.add_argument('--no-padding', dest='padding', default=True, action='store_false',
                        help="If passed, will eliminate padding from json file. \
                        Don't use if you don't know what you are doing")
    if print_help:
        parser.print_help()
    else:
        args = parser.parse_args()
        return args


class Args(object):
    def __init__(self, project_dir=None, input=None, output=None, log=None, temp=None,
                 use_test_data=False, sequences=None, chunksize=250, output_type='json',
                 merge=False, pandaseq_algo='simple_bayesian',
                 nextseq=False, umid=0, isotype=False,
                 basespace=False, cluster=False, starcluster=False,
                 debug=False, print_debug=False, species='human', gzip=False):
        super(Args, self).__init__()
        self.sequences = sequences
        self.data_dir = str(project_dir) if project_dir is not None else project_dir
        self.input = str(input) if input is not None else input
        self.output = str(output) if output is not None else output
        self.log = str(log) if log is not None else log
        self.temp = str(temp) if temp is not None else temp
        args.use_test_data = use_test_data
        self.chunksize = int(chunksize)
        self.output_type = str(output_type)
        self.merge = True if basespace else merge
        self.pandaseq_algo = str(pandaseq_algo)
        self.nextseq = nextseq
        self.uaid = int(umid)
        self.isotype = isotype
        self.basespace = basespace
        self.cluster = cluster
        self.starcluster = starcluster
        self.debug = 1 if debug else 0
        self.gzip = gzip
        if print_debug and self.debug > 0:
            self.debug == 2
        self.species = species


def validate_args(args):
    if all([args.use_test_data, args.input is None]):
        mod_dir = os.path.dirname(os.path.abspath(__file__))
        args.input = os.path.join(mod_dir, 'test_data')
        print(args.input)
    if not any([args.data_dir,
                args.sequences,
                all([args.input, args.output, args.temp])]):
        parse_arguments(print_help=True)
        sys.exit(1)
    if args.sequences:
        args.output_type = 'json'
    if all([args.sequences is not None, args.temp is None]):
        args.temp = './tmp'


#####################################################################
#
#                        FILES AND DIRECTORIES
#
#####################################################################


def make_directories(args):
    full_paths = []
    if args.data_dir and not os.path.exists(args.data_dir):
        _make_direc(args.data_dir, args)
    indir = args.input if args.input else os.path.join(args.data_dir, 'input')
    outdir = args.output if args.output else os.path.join(args.data_dir, 'output')
    tempdir = args.temp if args.temp else os.path.join(args.data_dir, 'temp')
    for d in [indir, outdir, tempdir]:
        d = os.path.abspath(d)
        full_paths.append(d)
        _make_direc(d, args)
    return full_paths


def _make_direc(d, args):
    if not os.path.exists(d):
        os.makedirs(d)
    if args.cluster:
        cmd = 'sudo chmod 777 {}'.format(d)
        p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
        stdout, stderr = p.communicate()


def make_merge_dir(args):
    merge_dir = os.path.join(args.data_dir, 'merged')
    _make_direc(merge_dir, args)
    return merge_dir


def setup_logging(args):
    log_dir = args.output if args.output else args.data_dir
    logfile = args.log if args.log else os.path.join(log_dir, 'abstar.log')
    debug = True if args.debug > 0 else False
    print_debug = True if args.debug == 2 else False
    log.setup_logging(logfile, debug=args.debug)
    global logger
    logger = log.get_logger('abstar')


def log_options(input_dir, output_dir, temp_dir, args):
    logger.info('')
    logger.info('')
    logger.info('-' * 25)
    logger.info('ABSTAR')
    logger.info('-' * 25)
    logger.info('')
    logger.info('SPECIES: {}'.format(args.species))
    logger.info('CHUNKSIZE: {}'.format(args.chunksize))
    logger.info('OUTPUT TYPE: {}'.format(args.output_type))
    if args.merge or args.basespace:
        logger.info('PANDASEQ ALGORITHM: {}'.format(args.pandaseq_algo))
    logger.info('UAID: {}'.format(args.uaid))
    logger.info('ISOTYPE: {}'.format('yes' if args.isotype else 'no'))
    logger.info('EXECUTION: {}'.format('cluster' if args.cluster else 'local'))
    logger.info('DEBUG: {}'.format('True' if args.debug > 0 else 'False'))
    logger.debug('INPUT: {}'.format(input_dir))
    logger.debug('OUTPUT: {}'.format(output_dir))
    logger.debug('TEMP: {}'.format(temp_dir))


def list_files(d, log=False):
    if os.path.isdir(d):
        glob_pattern = os.path.join(os.path.expanduser(d), '*')
        # need to fix square brackets, or glob freaks out
        if any(['[' in glob_pattern, ']' in glob_pattern]):
            glob_pattern = re.sub(r'\[', '[[]', glob_pattern)
            glob_pattern = re.sub(r'(?<!\[)\]', '[]]', glob_pattern)
        files = sorted(glob(glob_pattern))
    else:
        files = [d, ]
    if log:
        fnames = [os.path.basename(f) for f in files]
        logger.info('FILES: {}'.format(', '.join(fnames)))
    return files


def concat_output(input_file, jsons, output_dir, args):
    bname = os.path.basename(input_file)
    if '.' in bname:
        split_name = bname.split('.')
        oprefix = '.'.join(split_name[:-1])
    else:
        oprefix = bname
    osuffix = '.json' if args.output_type == 'json' else '.txt'
    oname = oprefix + osuffix
    ofile = os.path.join(output_dir, oname)
    logger.info('Concatenating {} job outputs into a single output file.'.format(len(jsons)))
    logger.info('')
    if args.gzip:
        ohandle = gzip.open(ofile + ".gz", 'wb')
    else:
        ohandle = open(ofile, 'w')

    with ohandle as out_file:
        if args.output_type in ['json', 'hadoop']:
            for json in jsons:
                with open(json) as f:
                    for line in f:
                        out_file.write(line)
                out_file.write('\n')
        if args.output_type == 'imgt':
            for i, json in enumerate(jsons):
                with open(json) as f:
                    for j, line in enumerate(f):
                        if i == 0:
                            out_file.write(line)
                        elif j >= 1:
                            out_file.write(line)
                    out_file.write('\n')
    return ofile


def clear_temp_files(temp_files):
    for f in temp_files:
        os.unlink(f)


#####################################################################
#
#                       INPUT PROCESSING
#
#####################################################################


def download_files(input_dir):
    from utils.basespace import BaseSpace
    bs = BaseSpace()
    logger.info('')
    logger.info('BASESPACE PROJECT NAME: {}'.format(bs.project_name))
    logger.info('BASESPACE PROJECT ID: {}'.format(bs.project_id))
    num_files = bs.download(input_dir)
    logger.info('')
    logger.info('BASESPACE FILE DOWNLOAD COUNT: {}'.format(num_files))


def merge_reads(input_dir, args):
    from utils import pandaseq
    merge_dir = make_merge_dir(args)
    pandaseq.run(input_dir,
                 merge_dir,
                 args.pandaseq_algo,
                 args.nextseq)
    return merge_dir


def format_check(input_list):
    formats = []
    for f in input_list:
        formats.append(_get_format(f))
    return formats


def _get_format(in_file):
    with open(in_file) as f:
        line = f.next()
        while line.strip() == '':
            line = f.next()
        if line.lstrip().startswith('>'):
            return 'fasta'
        elif line.lstrip().startswith('@'):
            return 'fastq'
        else:
            return None


def split_file(f, fmt, temp_dir, args):
    file_counter = 0
    seq_counter = 0
    total_seq_counter = 0
    subfiles = []
    fastas = []
    if '.' in os.path.basename(f):
        out_prefix = '.'.join(os.path.basename(f).split('.')[:-1])
    else:
        out_prefix = os.path.basename(f)
    if args.chunksize != 0:
        try:
            for seq in SeqIO.parse(open(f, 'r'), fmt.lower()):
                fastas.append('>{}\n{}'.format(seq.id, str(seq.seq)))
                seq_counter += 1
                total_seq_counter += 1
                if seq_counter == args.chunksize:
                    out_file = os.path.join(temp_dir, '{}_{}'.format(out_prefix, file_counter))
                    ohandle = open(out_file, 'w')
                    ohandle.write('\n'.join(fastas))
                    ohandle.close()
                    fastas = []
                    seq_counter = 0
                    file_counter += 1
                    subfiles.append(out_file)
        except ValueError:
            print('')
            print('ERROR: invalid file.')
            print('{} is not properly formatted'.format(f))
            print(traceback.format_exception_only)
            clear_temp_files(subfiles)
            sys.exit(1)

        # We don't want our files split
    else:
        for seq in SeqIO.parse(open(f, 'r'), fmt.lower()):
            total_seq_counter += 1
        subfiles.append(f)
        file_counter = 1

    # unless the input file is an exact multiple of args.chunksize,
    # need to write the last few sequences to a split file.
    if seq_counter:
        file_counter += 1
        out_file = os.path.join(temp_dir, '{}_{}'.format(out_prefix, file_counter))
        open(out_file, 'w').write('\n' + '\n'.join(fastas))
        subfiles.append(out_file)
    logger.info('SEQUENCES: {}'.format(total_seq_counter))
    logger.info('JOBS: {}'.format(file_counter))
    return subfiles, total_seq_counter


#####################################################################
#
#                            PRINTING
#
#####################################################################


def print_input_file_info(f, fmt):
    fname = os.path.basename(f)
    logger.info('')
    logger.info('')
    logger.info(fname)
    logger.info('-' * len(fname))
    logger.info('FORMAT: {}'.format(fmt.lower()))


def print_job_stats(total_seqs, good_seqs, start_time, end_time):
    run_time = end_time - start_time
    logger.info('{} sequences contained an identifiable rearrangement'.format(good_seqs))
    logger.info('AbStar completed in {} seconds'.format(run_time))


def update_progress(finished, jobs, failed=None):
    pct = int(100. * finished / jobs)
    ticks = pct / 2
    spaces = 50 - ticks
    if failed:
        prog_bar = '\r({}/{}) |{}{}|  {}% ({}, {})'.format(finished, jobs, '|' * ticks, ' ' * spaces, pct, finished - failed, failed)
    else:
        prog_bar = '\r({}/{}) |{}{}|  {}%'.format(finished, jobs, '|' * ticks, ' ' * spaces, pct)
    sys.stdout.write(prog_bar)
    sys.stdout.flush()


#####################################################################
#
#                             JOBS
#
#####################################################################


def run_jobs(files, output_dir, args):
    sys.stdout.write('\nRunning VDJ...\n')
    if args.cluster:
        return _run_jobs_via_celery(files, output_dir, args)
    elif args.debug or args.chunksize == 0:
        return _run_jobs_singlethreaded(files, output_dir, args)
    else:
        return _run_jobs_via_multiprocessing(files, output_dir, args)


def _run_jobs_singlethreaded(files, output_dir, args):
    from utils.vdj import run as run_vdj
    results = []
    for i, f in enumerate(files):
        try:
            results.append(run_vdj(f, output_dir, args))
            update_progress(i + 1, len(files))
        except:
            logger.debug('FILE-LEVEL EXCEPTION: {}'.format(f))
            logging.debug(traceback.format_exc())
    logger.info('')
    return results


def _run_jobs_via_multiprocessing(files, output_dir, args):
    from utils.vdj import run as run_vdj
    p = Pool(maxtasksperchild=50)
    async_results = []
    for f in files:
        async_results.append((f, p.apply_async(run_vdj, (f,
                                                         output_dir,
                                                         args))))
    monitor_mp_jobs([ar[1] for ar in async_results])
    results = []
    for a in async_results:
        try:
            results.append(a[1].get())
        except:
            logger.debug('FILE-LEVEL EXCEPTION: {}'.format(a[0]))
            if args.debug:
                traceback.print_exc()
            logging.debug(''.join(traceback.format_exc()))
            continue
    p.close()
    p.join()
    return results


def monitor_mp_jobs(results):
    finished = 0
    jobs = len(results)
    while finished < jobs:
        time.sleep(1)
        ready = [ar for ar in results if ar.ready()]
        finished = len(ready)
        update_progress(finished, jobs)
    sys.stdout.write('\n\n')


def _run_jobs_via_celery(files, output_dir, args):
    from utils.vdj import run as run_vdj
    async_results = []
    for f in files:
        async_results.append(run_vdj.delay(f,
                                           output_dir,
                                           args))
    succeeded, failed = monitor_celery_jobs(async_results)
    failed_files = [f for i, f in enumerate(files) if async_results[i].failed()]
    for ff in failed_files:
        logger.debug('FAILED FILE: {}'.format(f))
    return [s.get() for s in succeeded]


def monitor_celery_jobs(results):
    finished = 0
    jobs = len(results)
    while finished < jobs:
        time.sleep(1)
        succeeded = [ar for ar in results if ar.successful()]
        failed = [ar for ar in results if ar.failed()]
        finished = len(succeeded) + len(failed)
        update_progress(finished, jobs, failed=len(failed))
    sys.stdout.write('\n\n')
    return succeeded, failed


def run(*args, **kwargs):
    '''
    Runs AbStar.

    Input sequences can be provided in several different formats:

        1) individual sequences as positional arguments: ``run(seq1, seq2, temp=temp, output=output)``
        2) a list of sequences, as an argument: ``run([seq1, seq2], temp=temp, output=output)``
        3) a single FASTA/Q-formatted input file, passed via ``input``
        4) a directory of FASTA/Q-formatted files, passed via ``input``

    When passing sequences (not FASTA/Q files), the sequences can be in any format recognized
    by ``abtools.sequence.Sequence``, including:

        - a raw nucleotide sequence, as a string (a random sequence ID will be assigned)
        - a list/tuple of the format ``[sequence_id, sequence]``
        - a BioPython SeqRecord object
        - an AbTools Sequence object

    Either sequences, ``project_dir``, or all of ``input``, ``output`` and ``temp`` are required.


    Examples:

        If processing a single sequence, you can pass the raw sequence, as a string::

            import abstar

            result = abstar.run('ATGC')

        or a list/tuple of the format ``[sequence_id, sequence]``::

            result = abstar.run(['seq1', 'ATGC'])

        If you pass just the raw sequence, a random sequence ID will be generated with ``uuid.uuid4()``.
        In either case, when given a single sequence, ``abstar.run()`` will return a single AbTools ``Sequence``
        object. If running multiple sequences, you can either pass each sequence as a positional argument::

            result_list = run(['seq1', 'ATGC'], ['seq2', 'CGTA'])

        or you can pass a list of sequences as the first argument, in this case using sequences parsed from a
        FASTA file using Biopython::

            from Bio import SeqIO

            fasta = open('my_sequences.fasta', 'r')
            seqs = [s for s in SeqIO.parse(fasta, 'fasta')]
            result_list = abstar.run(seqs)

        When given multiple sequences, ``abstar.run()`` will return a list of AbTools ``Sequence`` objects,
        one per input sequence.

        If you'd prefer not to parse the FASTQ/A file into a list (for example, if the input file is
        extremely large), you can pass the input file path directly, along with a temp directory and output
        directory::

            result_files = abstar.run(input='/path/to/my_sequences.fasta',
                                      temp='/path/to/temp',
                                      output='/path/to/output')

        Given a file path, ``abstar.run()`` returns a list of output file paths. In the above case,
        ``result_files`` will be a list containing a single output file path:
        ``/path/to/output/my_sequences.json``.

        If you have a directory containing multiple FASTQ/A files, you can pass the directory path
        using ``input``::

            result_files = abstar.run(input='/path/to/input',
                                      temp='/path/to/temp',
                                      output='/path/to/output')

        As before, ``result_files`` will contain a list of output file paths.

        If your input directory contains paired FASTQ files (gzip compressed or uncompressed)
        that need to be merged prior to processing with AbStar::

            result_files = abstar.run(input='/path/to/input',
                                      temp='/path/to/temp',
                                      output='/path/to/output',
                                      merge=True)

        The paired read files in ``input`` will be merged with PANDAseq prior to processing with AbStar.
        By default, PANDAseq's 'simple bayesian' read merging algorithm is used, although alternate
        algorithms can be selected with ``pandaseq_algo``.

        AbStar also provides an alternate CSV-formatted output type that mimics the `IMGT Summary file`_.
        This option is provided to minimize the effort needed to convert existing
        IMGT-based pipelines to AbStar. Alternate output is only available when passing an input file or
        directory; passing individual sequences or a list of sequences will always return Sequence objects.
        To produce IMGT-formatted output::

            result_files = abstar.run(input='/path/to/input',
                                      temp='/path/to/temp',
                                      output='/path/to/output',
                                      output_type='imgt')

        .. _IMGT Summary file: http://www.imgt.org/IMGT_vquest/share/textes/imgtvquest.html#Esummary


    Args:

        project_dir (str): Path to the project directory. Most useful when directly downloading
            files from BaseSpace, and all subdirectories will be created by AbStar.

        input (str): Path to input directory, containing FASTA/Q files. If performing
            read merging with PANDAseq, paired FASTQ files may be gzip compressed.

        output (str): Path to output directory.

        temp (str): Path to temp directory, where intermediate job files will be stored.

        log (str): Path to log file. If not provided and ``project_dir`` is provided,
            the log will be written to ``/path/to/project_dir/abstar.log``. If output is
            provided, log will be written to ``/path/to/output/abstar.log``.

        species (str): Species of the antibody sequences. Choices are 'human', 'macaque',
            'mouse' and 'rabbit'. Default is 'human'.

        isotype (bool): If True, the isotype will infered by aligning the sequence region
            downstream of the J-gene. If False, the isotype will not be determined.
            Default is True.

        umid (int): Length (in nucleotides) of the Unique Molecular ID used to barcode input RNA.
            A positive integer results in the UMID being parsed from the start of the read (or merged
            read), a negative integer results in parsing from the end of the read. Default is 0,
            which results in no UMID parsing.

        gzip (bool): If True, compresses output files with gzip. Default is False.

        pretty (bool): If True, formats JSON output files to be more human-readable. If False,
            JSON output files contain one record per line. Default is False.

        output_type (str): Options are 'json' or 'imgt'. IMGT output mimics the Summary
            table produced by IMGT High-V/Quest, to maintain a level of compatibility with
            existing IMGT-based pipelines. JSON output is much more detailed. Default is 'json'.

        merge (bool): If True, input must be paired-read FASTA files (gzip compressed or uncompressed)
            which will be merged with PANDAseq prior to processing with AbStar. If ``basespace`` is True,
            ``merge`` is automatically set to True. Default is False.

        pandaseq_algo (str): Define merging algorithm to be used by PANDAseq. Options are
            'simple_bayesian', 'ea_util', 'flash', 'pear', 'rdp_mle', 'stitch', or 'uparse'. Default is
            'simple_bayesian', which is the default PANDAseq algorithm.

        debug (bool): If ``True``, ``abstar.run()`` runs in single-threaded mode, the log is much more verbose,
            and temporary files are not removed. Default is ``False``.


    Returns:

        If the input is a single sequence, ``run`` returns a single AbTools ``Sequence`` object.

        If the input is a list of sequences, ``run`` returns a list of AbTools ``Sequence`` objects.

        If the input is a file or a directory of files, ``run`` returns a list of output files.
    '''

    warnings.filterwarnings("ignore")
    if len(args) == 1:
        # if there's a single arg, need to check if it's a single sequence...
        try:
            sequences = [Sequence(args[0]), ]
        except:
            # ...or a list of sequences
            try:
                sequences = [Sequence(s) for s in args[0]]
            except:
                print('ERROR: invalid format for sequence input:')
                for a in args:
                    print(a)
                sys.exit(1)
    # if multiple args, assume each is a sequence
    elif len(args) > 1:
        try:
            sequences = [Sequence(s) for s in args]
        except:
            print('ERROR: invalid format for sequence input:')
            for a in args:
                print(a)
            sys.exit(1)
    kwargs['sequences'] = sequences
    args = Args(**kwargs)
    validate_args(args)
    global logger
    logger = log.get_logger('abstar')
    output = main(args)
    if args.sequences is not None:
        output = [Sequence(o) for o in output]
        if len(output) == 1:
            return output[0]
    return output


def run_standalone(args):
    output_dir = main(args)


def main(args):
    if args.sequences is not None:
        from utils import output, vdj
        vdjs = vdj.process_sequences(args.sequences, args)
        return output.as_dict(vdjs)
    else:
        input_dir, output_dir, temp_dir = make_directories(args)
        setup_logging(args)
        log_options(input_dir, output_dir, temp_dir, args)
        if args.basespace:
            args.merge = True
            download_files(input_dir)
        if args.merge:
            input_dir = merge_reads(input_dir, args)
        if args.isotype:
            args.isotype = args.species
        input_files = [f for f in list_files(input_dir, log=True) if os.stat(f).st_size > 0]
        output_files = []
        for f, fmt in zip(input_files, format_check(input_files)):
            # skip the non-FASTA/Q files
            if fmt is None:
                continue
            start_time = time.time()
            print_input_file_info(f, fmt)
            subfiles, seq_count = split_file(f, fmt, temp_dir, args)
            run_info = run_jobs(subfiles, temp_dir, args)
            temp_json_files = [r[0] for r in run_info if r is not None]
            processed_seq_counts = [r[1] for r in run_info if r is not None]
            vdj_end_time = time.time()
            output_file = concat_output(f, temp_json_files, output_dir, args)
            output_files.append(output_file)
            if not args.debug:
                clear_temp_files(subfiles + temp_json_files)
            print_job_stats(seq_count, sum(processed_seq_counts), start_time, vdj_end_time)
        return output_files


if __name__ == '__main__':
    warnings.filterwarnings("ignore")
    args = parse_arguments()
    output_dir = main(args)
    sys.stdout.write('\n\n')
