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

import sys
import traceback
import warnings
import logging
import os
import time
import gzip


# External
#import skbio

#####################################################################
#
#                             ARGUMENTS
#
#####################################################################


def parse_arguments(print_help=False):
    parser = ArgumentParser("Performs germline assignment and other relevant annotation on antibody sequence data from NGS platforms.")
    parser.add_argument('-d', '--data', dest='data_dir', default=None,
                        help="The data directory, where files will be downloaded (or have previously \
                        been download), temp files will be stored, and output files will be \
                        written. During StarCluster configuration, all ephemeral drives on the master \
                        node will be combined into a single RAID volume and mounted at the data directory. \
                        Not necessary if '-o', '-t' and '-i' are provided. \
                        Default is '/data'.")
    parser.add_argument('-i', '--in', dest='input', default=None,
                        help="The input file or directory, to be split and processed in parallel. \
                        If a directory is given, all files in the directory will be iteratively processed. \
                        Required.")
    parser.add_argument('-o', '--out', dest='output', default=None,
                        help="The output directory, into which the JSON-formatted output files will be deposited. \
                        If the directory does not exist, it will be created. \
                        Required.")
    parser.add_argument('-l', '--log', dest='log', default=None,
                        help="The log file, to which log info will be written. \
                        Default is <output_directory>/abstar.log if '-o/--out' is specificied. \
                        If '-o/--out' isn't provided, default is <data_directory>/abstar.log")
    parser.add_argument('-t', '--temp', dest='temp', default=None,
                        help="The directory in which temp files will be stored. If the directory doesn't exist, \
                        it will be created. Required.")
    parser.add_argument('-k', '--chunksize', dest='chunksize', default=250, type=int,
                        help="Approximate number of sequences in each distributed job. \
                        Defaults to 250. \
                        Set to 0 if you want file splittint to be turned off \
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
    parser.add_argument('-p', '--pandaseq_algo', dest="pandaseq_algo", default='simple_bayesian',
                        choices=['simple_bayesian', 'ea_util', 'flash', 'pear', 'rdp_mle', 'stitch', 'uparse'],
                        help="Define merging algorithm to be used by PANDAseq.\
                        Options are 'simple_bayesian', 'ea_util', 'flash', 'pear', 'rdp_mle', 'stitch', or 'uparse'.\
                        Default is 'simple_bayesian', which is the default PANDAseq algorithm.")
    parser.add_argument('-n', '--nextseq', dest="nextseq", action='store_true', default=False,
                        help="Use if the run was performed on a NextSeq sequencer.")
    parser.add_argument('-u', '--uaid', dest="uaid", type=int, default=0,
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
    """Holds arguments, mimics argparse's Namespace when running abstar as an imported module"""

    def __init__(self, data_dir=None, input=None, output=None, log=None, temp=None,
                 chunksize=250, output_type='json',
                 merge=False, pandaseq_algo='simple_bayesian',
                 nextseq=False, uaid=0, isotype=False,
                 basespace=False, cluster=False, starcluster=False,
                 debug=False, print_debug=False, species='human', gzip=False):
        super(Args, self).__init__()
        self.data_dir = str(data_dir) if data_dir is not None else data_dir
        self.input = str(input) if input is not None else input
        self.output = str(output) if output is not None else output
        self.log = str(log) if log is not None else log
        self.temp = str(temp) if temp is not None else temp
        self.chunksize = int(chunksize)
        self.output_type = str(output_type)
        self.merge = True if basespace else merge
        self.pandaseq_algo = str(pandaseq_algo)
        self.nextseq = nextseq
        self.uaid = int(uaid)
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
    if not args.data_dir and not all([args.input, args.output, args.temp]):
        parse_arguments(print_help=True)
        sys.exit(1)


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
    log.setup_logging(logfile, debug=args.debug, print_debug=print_debug)
    global logger
    logger = log.get_logger('abstar')


def log_options(args):
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
    if args.debug > 0:
        logger.info('DEBUG LEVEL: {}'.format('print exceptions' if args.debug == 2 else 'log exceptions'))


def list_files(d, log=False):
    if os.path.isdir(d):
        expanded_dir = os.path.expanduser(d)
        files = sorted(glob(expanded_dir + '/*'))
    else:
        files = [d, ]
    if log:
        fnames = [os.path.basename(f) for f in files]
        logger.info('FILES: {}'.format(', '.join(fnames)))
    return files


def concat_output(input_file, temp_dir, output_dir, args):
    bname = os.path.basename(input_file)
    if '.' in bname:
        split_name = bname.split('.')
        oprefix = '.'.join(split_name[:-1])
    else:
        oprefix = bname
    osuffix = '.json' if args.output_type == 'json' else '.txt'
    oname = oprefix + osuffix
    ofile = os.path.join(output_dir, oname)
    #open(ofile, 'w').write('')
    jsons = [f for f in list_files(temp_dir) if os.path.basename(f).startswith(oprefix) and f.endswith(osuffix)]
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


def clear_temp_dir(temp_dir):
    for f in list_files(temp_dir):
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
                 args,
                 log=True)
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


def run(**kwargs):
    '''
    Runs AbStar.

    Either ::data_dir:: or all of ::input::, ::output::, and ::temp::
    are required.

    To parse unique antibody IDs (UAIDs, or molecular barcodes), set ::uaid::
    to the length of the barcode. It is assumed that the barcode is at the start
    if the sequence.

    Default options:

    data_dir=None
    input=None
    output=None
    log=None
    temp=None
    chunksize=250
    output_type='json'
    merge=False
    pandaseq_algo='simple_bayesian'
    next_seq=False
    uaid=0
    isotype=False
    basespace=False
    cluster=False
    starcluster=False
    species='human'
    debug=False
    '''
    warnings.filterwarnings("ignore")
    args = Args(**kwargs)
    validate_args(args)
    global logger
    logger = log.get_logger('abstar')
    output_dir = main(args)
    return list_files(output_dir)


def run_standalone(args):
    output_dir = main(args)


def main(args):
    input_dir, output_dir, temp_dir = make_directories(args)
    setup_logging(args)
    log_options(args)
    if args.basespace:
        args.merge = True
        download_files(input_dir)
    if args.merge:
        input_dir = merge_reads(input_dir, args)
    if args.isotype:
        args.isotype = args.species
    input_files = [f for f in list_files(input_dir, log=True) if os.stat(f).st_size > 0]
    for f, fmt in zip(input_files, format_check(input_files)):
        # skip the non-FASTA/Q files
        if fmt is None:
            continue
        start_time = time.time()
        print_input_file_info(f, fmt)
        subfiles, seq_count = split_file(f, fmt, temp_dir, args)
        job_stats = run_jobs(subfiles, temp_dir, args)
        vdj_end_time = time.time()
        concat_output(f, temp_dir, output_dir, args)
        clear_temp_dir(temp_dir)
        print_job_stats(seq_count, sum(job_stats), start_time, vdj_end_time)
    return output_dir


if __name__ == '__main__':
    warnings.filterwarnings("ignore")
    args = parse_arguments()
    output_dir = main(args)
    sys.stdout.write('\n\n')
