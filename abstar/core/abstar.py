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


from __future__ import absolute_import, division, print_function, unicode_literals

import warnings

warnings.simplefilter(action="ignore", category=FutureWarning)

from argparse import ArgumentParser
import csv
from glob import glob
import gzip
import json
import logging
from multiprocessing import Pool, cpu_count
import os
import pkg_resources
import re
from subprocess import Popen, PIPE
import sys
import tempfile
import time
import traceback
from typing import Iterable, Optional
import warnings
import shutil

from Bio import SeqIO

import dask.dataframe as dd

from abutils.core.sequence import Sequence, read_json, read_csv
from abutils.utils import log
from abstar.utils.parquet_schema import schema

# from abutils.utils.pipeline import list_files

from .antibody import Antibody
from ..assigners.assigner import BaseAssigner
from ..assigners.registry import ASSIGNERS

# from ..utils import output
from ..utils.output import (
    get_abstar_result,
    get_abstar_results,
    get_output,
    write_output,
    get_header,
    get_output_suffix,
    get_output_separator,
    get_parquet_dtypes,
)
from ..utils.queue.celery import celery


if sys.version_info[0] > 2:
    STR_TYPES = [
        str,
    ]
else:
    STR_TYPES = [str, unicode]


from ..version import __version__

# __version__ = pkg_resources.require("abstar")[0].version

# ASSIGNERS = {cls.__name__.lower(): cls for cls in vars()['BaseAssigner'].__subclasses__()}


#####################################################################
#
#                             ARGUMENTS
#
#####################################################################


def create_parser() -> ArgumentParser:
    parser = ArgumentParser(
        prog="abstar",
        description="VDJ assignment and antibody sequence annotation. Scalable from a single sequence to billions of sequences.",
    )
    parser.add_argument(
        "-p",
        "--project",
        dest="project_dir",
        default=None,
        help="The data directory, where files will be downloaded (or have previously \
                        been download), temp files will be stored, and output files will be \
                        written. During StarCluster configuration, all ephemeral drives on the master \
                        node will be combined into a single RAID volume and mounted at the data directory. \
                        Not necessary if '-o', '-t' and '-i' are provided. \
                        Default is '/data'.",
    )
    parser.add_argument(
        "-i",
        "--input",
        dest="input",
        default=None,
        help="The input file or directory, to be split and processed in parallel. \
                        If a directory is given, all files in the directory will be iteratively processed. \
                        Required, unless <project> is provided.",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        default=None,
        help="The output directory, into which the JSON-formatted output files will be deposited. \
                        If the directory does not exist, it will be created. \
                        Required, unless <project> is provided.",
    )
    parser.add_argument(
        "-l",
        "--log",
        dest="log",
        default=None,
        help="The log directory, into which log files will be deposited. \
                        Default is <project_directory>/log if <project> is supplied, otherwise \
                        the parent directory of <output> will be used.",
    )
    parser.add_argument(
        "-t",
        "--temp",
        dest="temp",
        default=None,
        help="The directory in which temp files will be stored. If the directory doesn't exist, \
                        it will be created. Required, unless <project> is provided.",
    )
    parser.add_argument(
        "--sequences",
        dest="sequences",
        default=None,
        help="Only used when passing sequences directly through the API.",
    )
    parser.add_argument(
        "-a",
        "--assigner",
        dest="assigner",
        default="blastn",
        help="VDJ germline assignment method to use. \
                        Options are: {}. Default is blastn".format(
            ", ".join(ASSIGNERS.keys())
        ),
    )
    parser.add_argument(
        "-k",
        "--chunksize",
        dest="chunksize",
        default=500,
        type=int,
        help="Approximate number of sequences in each distributed job. \
                        Defaults to 500. \
                        Set to 0 if you want file splitting to be turned off \
                        Don't change unless you know what you're doing.",
    )
    parser.add_argument(
        "-O",
        "--output-type",
        dest="output_type",
        action="append",
        choices=["json", "imgt", "tabular", "airr"],
        help="Select the output type. Options are 'json', 'imgt', 'csv' and 'airr'. \
                        IMGT output mimics the Summary table produced by IMGT High-V/Quest, \
                        to maintain some level of compatibility with existing IMGT-based pipelines. \
                        JSON output is much more detailed, and is suitable for direct import into MongoDB. \
                        Tabular output contains a subset of the JSON output format in CSV format. \
                        AIRR output is tab-delimited and conforms to AIRR data standards. \
                        Defaults to JSON output.",
    )
    parser.add_argument(
        "--parquet",
        dest="parquet",
        action="store_true",
        default=False,
        help="For tabular output formats, create parquet formatted output in addition to \
                        the tabular delimited output. Default is False.",
    )
    parser.add_argument(
        "-m",
        "--merge",
        dest="merge",
        action="store_true",
        default=False,
        help="Use if the input files are paired-end FASTQs \
                        (either gzip compressed or uncompressed) from Illumina platforms. \
                        Prior to running the germline assignment pipeline, \
                        paired reads will be merged with PANDAseq.",
    )
    parser.add_argument(
        "-P",
        "--pandaseq-algo",
        dest="pandaseq_algo",
        default="simple_bayesian",
        choices=[
            "simple_bayesian",
            "ea_util",
            "flash",
            "pear",
            "rdp_mle",
            "stitch",
            "uparse",
        ],
        help="Define merging algorithm to be used by PANDAseq.\
                        Options are 'simple_bayesian', 'ea_util', 'flash', 'pear', 'rdp_mle', 'stitch', or 'uparse'.\
                        Default is 'simple_bayesian', which is the default PANDAseq algorithm.",
    )
    parser.add_argument(
        "-n",
        "--nextseq",
        dest="nextseq",
        action="store_true",
        default=False,
        help="Use if the run was performed on a NextSeq sequencer.",
    )
    parser.add_argument(
        "-u",
        "--uid",
        dest="uid",
        type=int,
        default=0,
        help="Length of the unique identifier (UID, or molecular barcode) if used. \
                        Default is unbarcoded (UID length of 0).",
    )
    parser.add_argument(
        "-I",
        "--isotype",
        dest="isotype",
        action="store_false",
        default=True,
        help="If set, the isotype will not be determined for heavy chains.\
                        If not set, isotyping sequences for the appropriate species will be used.",
    )
    parser.add_argument(
        "-b",
        "--basespace",
        dest="basespace",
        default=False,
        action="store_true",
        help="Use if files should be downloaded directly from BaseSpace. \
                        Files will be downloaded into the input directory.",
    )
    parser.add_argument(
        "-c",
        "--cluster",
        dest="cluster",
        default=False,
        action="store_true",
        help="Use if performing computation on a Celery cluster. \
                        If set, input files will be split into many subfiles and passed \
                        to a Celery queue. If not set, input files will still be split, but \
                        will be distributed to local processors using multiprocessing.",
    )
    parser.add_argument(
        "-N",
        "--num-cores",
        dest="num_cores",
        default=0,
        type=int,
        help="Number of cores used by abstar. Default is `0`, which uses \
                        all available cores.",
    )
    parser.add_argument(
        "-D",
        "--debug",
        dest="debug",
        action="store_true",
        default=False,
        help="If set, logs information about successfully assigned sequences as well \
                        as unsuccessful sequences. Useful for debugging small test datasets, \
                        as the logging is fairly detailed. \
                        Default is to only log unsuccessful sequences.",
    )
    parser.add_argument(
        "-T",
        "--use-test-data",
        dest="use_test_data",
        action="store_true",
        default=False,
        help="If set, AbStar will run on a small sample of test data (1000 sequences). \
                        These samples are of human origin, but AbStar doesn't force the use of the human \
                        germline database. Of the 1000 sequences, only 999 contain an antibody rearrangement, \
                        so one sequence should fail the germline assignment process. \
                        Default is False. \
                        Providing temp and output directories (or a project directory) is required.",
    )
    parser.add_argument(
        "--germ_db",
        dest="germ_db",
        default="human",
        help="Germline database to use. Built-in \
                        options include 'human', 'mouse', 'macaaque' and 'humouse'. Custom databases can alsu \
                        be built and used. Default is 'human'.",
    )
    parser.add_argument(
        "-r",
        "--receptor",
        default="bcr",
        help="Receptor type. Options are 'bcr' and 'tcr'.",
    )
    parser.add_argument(
        "-z",
        "--gzip",
        dest="gzip",
        default=False,
        action="store_true",
        help="Compress the output to a gzipped file",
    )
    parser.add_argument(
        "--raw",
        dest="raw",
        default=False,
        action="store_true",
        help="Returns raw output (a dict).",
    )
    parser.add_argument(
        "--json-keys",
        dest="json_keys",
        default=None,
        help='If supplied, allows selection of a subset of the normal JSON output. \
                        Should be provided as a list of JSON key names, separated by commas: \
                        "--json-keys seq_id,v_gene,d_gene,j_gene,cdr3_aa". \
                        Only top-level keys are supported (that is, cannot select a single nested element). \
                        Note that if --add-padding is set, padding will be included whether or not the \
                        padding key name is included in --json-keys.',
    )
    parser.add_argument(
        "--pretty",
        dest="pretty",
        default=False,
        action="store_true",
        help="Pretty format json file",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s {version}".format(version=__version__),
    )
    parser.add_argument(
        "--add-padding",
        dest="padding",
        default=False,
        action="store_true",
        help="If passed, will add padding to the json output file. \
                        Really only useful if you're using an old version of MongoDB.",
    )
    parser.add_argument(
        "--quiet",
        dest="verbose",
        default=True,
        action="store_false",
        help="If set, suppresses logging and printing progress to screen",
    )

    return parser


class Args(object):
    def __init__(
        self,
        project_dir=None,
        input=None,
        output=None,
        log=None,
        temp=None,
        sequences=None,
        chunksize=500,
        output_type=["json",],
        assigner="blastn",
        merge=False,
        pandaseq_algo="simple_bayesian",
        use_test_data=False,
        parquet=False,
        nextseq=False,
        uid=0,
        isotype=False,
        pretty=False,
        num_cores=0,
        basespace=False,
        cluster=False,
        padding=False,
        raw=False,
        json_keys=None,
        debug=False,
        germ_db="human",
        receptor="bcr",
        gzip=False,
        verbose=True,
    ):
        super(Args, self).__init__()
        self.sequences = sequences
        self.project_dir = (
            os.path.abspath(project_dir) if project_dir is not None else project_dir
        )
        self.input = os.path.abspath(input) if input is not None else input
        self.output = os.path.abspath(output) if output is not None else output
        self.use_test_data = use_test_data
        self.log = os.path.abspath(log) if log is not None else log
        self.temp = os.path.abspath(temp) if temp is not None else temp
        self.chunksize = int(chunksize)
        self.output_type = [output_type,] if output_type in STR_TYPES else output_type
        self.parquet = parquet
        self.assigner = assigner
        self.merge = True if basespace else merge
        self.pandaseq_algo = str(pandaseq_algo)
        self.nextseq = nextseq
        self.uid = int(uid)
        self.isotype = isotype
        self.basespace = basespace
        self.cluster = cluster
        self.num_cores = num_cores
        self.pretty = pretty
        self.debug = debug
        self.gzip = gzip
        self.raw = raw
        self.json_keys = json_keys
        self.padding = padding
        self.germ_db = germ_db
        self.receptor = receptor.lower()
        self.verbose = verbose


def validate_args(args):
    # make sure all necessary args were provided
    if not any(
        [
            args.project_dir,
            args.sequences,
            all([any([args.input, args.use_test_data]), args.output, args.temp]),
        ]
    ):
        create_parser().print_help()
        sys.exit(1)
        # alter output type if abstar is being run interactively
        # if args.sequences:
        #     if len(args.output_type) > 1:
        #         print('\nWARNING: Multiple output formats are not supported when runing abstar in interactive mode.')
        #         print('\nUsing {} format'.format(args.output_type[0]))
        #         args.output_type = args.output_type[:1]
        args.raw = True

    # find the number of available cores if abstar is being run on all of them
    if args.num_cores == 0:
        args.num_cores = cpu_count()

    # set default temp directory if not provided
    if all([args.sequences is not None, args.temp is None]):
        args.temp = "/tmp"

    # set default output type if not provided
    if args.output_type is None:
        args.output_type = [
            "json",
        ]
    if isinstance(args.output_type, str):
        args.output_type = [
            args.output_type,
        ]

    # process JSON key string if provided
    if args.json_keys is not None:
        args.json_keys = args.json_keys.split(",")

    # check that the supplied receptor is appropriate:
    if args.receptor not in ["bcr", "tcr"]:
        print(f"\nERROR: Supplied receptor type {args.receptor} is not supported.")
        print("Abstar supports the following receptor types: BCR and TCR.")

    # check to ensure a germline database exists for the requested species
    addon_species_dbs = []
    builtin_species_dbs = []
    addon_germline_dbs_dir = os.path.expanduser(
        f"~/.abstar/germline_dbs/{args.receptor}"
    )
    if os.path.isdir(addon_germline_dbs_dir):
        addon_species_dbs += [
            os.path.basename(d[0]) for d in os.walk(addon_germline_dbs_dir)
        ]
    mod_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    builtin_germline_dbs_dir = os.path.join(
        mod_dir, f"assigners/germline_dbs/{args.receptor}"
    )
    if os.path.isdir(builtin_germline_dbs_dir):
        builtin_species_dbs += [
            os.path.basename(d) for d in list_files(builtin_germline_dbs_dir)
        ]
    if not any(
        [
            args.germ_db.lower() in addon_species_dbs,
            args.germ_db.lower() in builtin_species_dbs,
        ]
    ):
        print(
            f"\nERROR: A germline database was not found for the requested species ({args.species})."
        )
        print(
            f"\nBuilt-in {args.receptor.upper()} databases exist for the following species:"
        )
        print(", ".join(builtin_species_dbs))
        if addon_species_dbs:
            print(
                f"\nUser-generated {args.receptor.upper()} databases exist for the following species:"
            )
            print(", ".join(addon_species_dbs))
        else:
            print("No user-generated databases exist on this computer.")
        print(
            "\nCustom germline databases can be created with the build_abstar_germline_db command"
        )
        sys.exit(1)


#####################################################################
#
#                        FILES AND DIRECTORIES
#
#####################################################################


def make_directories(args):
    full_paths = []
    if args.project_dir is not None and not os.path.exists(args.project_dir):
        _make_direc(args.project_dir, args)
    if args.use_test_data:
        indir = None
    else:
        indir = (
            args.input
            if args.input is not None
            else os.path.join(args.project_dir, "input")
        )
    outdir = (
        args.output
        if args.output is not None
        else os.path.join(args.project_dir, "output")
    )
    tempdir = (
        args.temp if args.temp is not None else os.path.join(args.project_dir, "temp")
    )
    log_parent = (
        args.project_dir if args.project_dir is not None else os.path.dirname(outdir)
    )
    logdir = args.log if args.log is not None else os.path.join(log_parent, "log")
    logtempdir = os.path.join(logdir, "temp")
    for d in [indir, outdir, tempdir, logdir, logtempdir]:
        if d is None:
            full_paths.append(None)
            continue
        d = os.path.abspath(d)
        full_paths.append(d)
        _make_direc(d, args.cluster)
    if tempdir is not None:
        for subdir in ["input"] + args.output_type:
            _make_direc(os.path.join(tempdir, subdir), args.cluster)
    if outdir is not None:
        for subdir in args.output_type:
            _make_direc(os.path.join(outdir, subdir), args.cluster)
    return full_paths[:-1]


def _make_direc(d, cluster):
    if not os.path.exists(d):
        os.makedirs(d)
    if cluster:
        cmd = "sudo chmod 777 {}".format(d)
        p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
        stdout, stderr = p.communicate()


def make_merge_dir(args):
    merge_parent = (
        args.project_dir
        if args.project_dir is not None
        else os.path.dirname(args.input)
    )
    merge_dir = os.path.abspath(os.path.join(merge_parent, "merged"))
    _make_direc(merge_dir, args)
    return merge_dir


def setup_logging(log_dir, debug):
    logfile = os.path.join(log_dir, "abstar.log")
    debug = True if debug > 0 else False
    # print_debug = True if debug == 2 else False
    log.setup_logging(logfile, debug=debug)
    global logger
    logger = log.get_logger("abstar")


def log_options(input_dir, output_dir, temp_dir, args):
    logger.info("")
    logger.info("")
    logger.info("-" * 25)
    logger.info("ABSTAR")
    logger.info("-" * 25)
    logger.info("")
    logger.info("GERMLINE DB: {}".format(args.germ_db))
    logger.info("RECEPTOR: {}".format(args.receptor))
    logger.info("CHUNKSIZE: {}".format(args.chunksize))
    logger.info("OUTPUT TYPE: {}".format(", ".join(args.output_type)))
    if args.merge or args.basespace:
        logger.info("PANDASEQ ALGORITHM: {}".format(args.pandaseq_algo))
    logger.info("UID: {}".format(args.uid))
    logger.info("ISOTYPE: {}".format("yes" if args.isotype else "no"))
    logger.info("EXECUTION: {}".format("cluster" if args.cluster else "local"))
    logger.info("DEBUG: {}".format("True" if args.debug else "False"))
    logger.debug("INPUT: {}".format(input_dir))
    logger.debug("OUTPUT: {}".format(output_dir))
    logger.debug("TEMP: {}".format(temp_dir))


def list_files(d, log=False):
    if os.path.isdir(d):
        glob_pattern = os.path.join(os.path.expanduser(d), "*")
        # need to fix square brackets, or glob freaks out
        if any(["[" in glob_pattern, "]" in glob_pattern]):
            glob_pattern = re.sub(r"\[", "[[]", glob_pattern)
            glob_pattern = re.sub(r"(?<!\[)\]", "[]]", glob_pattern)
        files = sorted(glob(glob_pattern))
    else:
        files = [
            d,
        ]
    if log:
        fnames = [os.path.basename(f) for f in files]
        logger.info("")
        logger.info("FILE COUNT: {}".format(len(fnames)))
        logger.info("FILES: {}".format(", ".join(fnames)))
    return files


# def get_output_suffix(output_format):
#     osuffixes = {'json': '.json',
#                  'imgt': '.csv',
#                  'tabular': '.txt'}
#     return osuffixes[output_format.lower()]


def build_output_base(output_types):
    needs_header = ["imgt", "tabular", "airr"]
    output_dict = {}
    for output_type in output_types:
        if output_type in needs_header:
            output_dict[output_type] = [
                get_header(output_type),
            ]
        else:
            output_dict[output_type] = []
    return output_dict


def concat_outputs(input_file, temp_output_file_dicts, output_dir, args):
    # temp_output_file_dicts is a list of dicts with output format as key and output file as value
    ofiles = []
    bname = os.path.basename(input_file)
    if "." in bname:
        split_name = bname.split(".")
        oprefix = ".".join(split_name[:-1])
    else:
        oprefix = bname
    logger.info("")
    for output_type in sorted(args.output_type):
        temp_files = [d[output_type] for d in temp_output_file_dicts]
        osuffix = get_output_suffix(output_type)
        oname = oprefix + osuffix
        output_subdir = os.path.join(output_dir, output_type)
        if not os.path.exists(output_subdir):
            os.makedirs(output_subdir)
        ofile = os.path.join(output_subdir, oname)
        # temp_files = [tf for tf in temp_output_files if tf.endswith(osuffix)]
        logger.info(
            "Concatenating {} {}-formatted job outputs into a single output file".format(
                len(temp_files), output_type.upper()
            )
        )
        if args.gzip:
            ohandle = gzip.open(ofile + ".gz", "wb")
        else:
            ohandle = open(ofile, "wb")

        with ohandle as out_file:
            # JSON-formatted files don't have headers, so we don't worry about it
            if output_type == "json" and not args.parquet:
                for temp_file in temp_files:
                    with open(temp_file, "rb") as f:
                        shutil.copyfileobj(
                            f, out_file, length=16 * 1024 ** 2
                        )  # Increasing buffer size to 16MB for faster transfer
            # For file formats with headers, only keep headers from the first file
            elif output_type in ["imgt", "tabular", "airr"]:
                for i, temp_file in enumerate(temp_files):
                    with open(temp_file, "rb") as f:
                        for j, line in enumerate(f):
                            if i == 0:
                                out_file.write(line)
                            elif j >= 1:
                                out_file.write(line)

        if args.parquet:
            logger.info("Converting concatenated output to parquet format")
            # Make clear the output format from which the parquet file is generated.
            # If the parquet files are generated from json for example, the filename would be f"{oprefix}_from_json".
            pname = f"{oprefix}_from_{output_type}"
            pfile = os.path.join(output_subdir, pname)
            dtypes = get_parquet_dtypes(output_type)

            if output_type == "json":
                df = dd.read_parquet(
                    os.path.dirname(temp_files[0])
                )  # Read in all parquet files in temp dir
                df.repartition(partition_size="100MB").to_parquet(
                    pfile,
                    engine="pyarrow",
                    compression="snappy",
                    write_index=False,
                    schema=schema,
                )
            else:
                df = dd.read_csv(
                    ofile, sep=get_output_separator(output_type), dtype=dtypes
                )
                df.to_parquet(
                    pfile, engine="pyarrow", compression="snappy", write_index=False
                )

        ofiles.append(ofile)
    return ofiles


def concat_logs(input_file, logs, log_dir, log_type):
    bname = os.path.basename(input_file)
    if "." in bname:
        split_name = bname.split(".")
        lprefix = ".".join(split_name[:-1])
    else:
        lprefix = bname
    lfile = os.path.join(log_dir, "{}.{}".format(lprefix, log_type))
    lhandle = open(lfile, "w")
    with lhandle as logfile:
        for log in logs:
            with open(log) as f:
                for line in f:
                    logfile.write(line)
    return lfile


def clear_temp_files(temp_files):
    for f in temp_files:
        os.unlink(f)


#####################################################################
#
#                       INPUT PROCESSING
#
#####################################################################


def download_files(input_dir):
    from ..utils.basespace import BaseSpace

    bs = BaseSpace()
    logger.info("")
    logger.info("BASESPACE PROJECT NAME: {}".format(bs.project_name))
    logger.info("BASESPACE PROJECT ID: {}".format(bs.project_id))
    num_files = bs.download(input_dir)
    logger.info("")
    logger.info("BASESPACE FILE DOWNLOAD COUNT: {}".format(num_files))


def merge_reads(input_dir, args):
    from ..utils import pandaseq

    merge_dir = make_merge_dir(args)
    pandaseq.run(input_dir, merge_dir, args.pandaseq_algo, args.nextseq, args.debug)
    return merge_dir


def format_check(input_list):
    formats = []
    for f in input_list:
        formats.append(_get_format(f))
    return formats


def _get_format(in_file):
    with open(in_file) as f:
        line = f.readline()
        while line.strip() == "":
            line = f.readline()
        if line.lstrip().startswith(">"):
            return "fasta"
        elif line.lstrip().startswith("@"):
            return "fastq"
        else:
            return None


def split_file(f, fmt, temp_dir, args):
    file_counter = 0
    seq_counter = 0
    total_seq_counter = 0
    subfiles = []
    sequences = []
    if "." in os.path.basename(f):
        out_prefix = ".".join(os.path.basename(f).split(".")[:-1])
    else:
        out_prefix = os.path.basename(f)
    if args.chunksize != 0:
        try:
            with open(f, "r") as f_handle:
                for seq in SeqIO.parse(f_handle, fmt.lower()):
                    sequences.append(seq.format(fmt.lower()))
                    seq_counter += 1
                    total_seq_counter += 1
                    if seq_counter == args.chunksize:
                        out_file = os.path.join(
                            temp_dir, "{}_{}".format(out_prefix, file_counter)
                        )
                        ohandle = open(out_file, "w")
                        ohandle.write("\n".join(sequences))
                        ohandle.close()
                        sequences = []
                        seq_counter = 0
                        file_counter += 1
                        subfiles.append(out_file)
        except ValueError:
            print("")
            print("ERROR: invalid file.")
            print("{} is not properly formatted".format(f))
            print(traceback.format_exception_only)
            clear_temp_files(subfiles)
            sys.exit(1)

    # We don't want our files split
    else:
        with open(f, "r") as f_handle:
            for seq in SeqIO.parse(f_handle, fmt.lower()):
                total_seq_counter += 1
            subfiles.append(f)
            file_counter = 1

    # unless the input file is an exact multiple of args.chunksize,
    # need to write the last few sequences to a split file.
    if seq_counter:
        out_file = os.path.join(temp_dir, "{}_{}".format(out_prefix, file_counter))
        open(out_file, "w").write("\n".join(sequences))
        subfiles.append(out_file)
        file_counter += 1
    logger.info("SEQUENCES: {}".format(total_seq_counter))
    logger.info("JOBS: {}".format(file_counter))
    return subfiles, total_seq_counter


#####################################################################
#
#                            PRINTING
#
#####################################################################


def print_input_file_info(f, fmt):
    fname = os.path.basename(f)
    logger.info("")
    logger.info("")
    logger.info(fname)
    logger.info("-" * len(fname))
    logger.info("FORMAT: {}".format(fmt.lower()))


def log_job_stats(total_seqs, good_seq_counts, start_time, end_time):
    run_time = end_time - start_time
    zero_files = sum([c == 0 for c in good_seq_counts])
    if zero_files > 0:
        logger.info(
            "{} files contained no successfully processed sequences".format(zero_files)
        )
    good_seqs = sum(good_seq_counts)
    logger.info("")
    logger.info(
        "{} sequences contained an identifiable rearrangement".format(good_seqs)
    )
    logger.info("AbStar completed in {} seconds".format(run_time))


def print_job_stats(total_seqs, good_seq_counts, start_time, end_time):
    run_time = end_time - start_time
    zero_files = sum([c == 0 for c in good_seq_counts])
    if zero_files > 0:
        logger.info(
            "{} files contained no successfully processed sequences".format(zero_files)
        )
    good_seqs = sum(good_seq_counts)
    print("{} sequences contained an identifiable rearrangement".format(good_seqs))
    print("abstar completed in {} seconds".format(round(run_time, 2)))
    print("")


#####################################################################
#
#                             JOBS
#
#####################################################################


@celery.task
def run_abstar(seq_file, output_dir, log_dir, file_format, arg_dict):
    """
    Wrapper function to multiprocess (or not) the assignment of V, D and J
    germline genes. Also writes the JSON-formatted output to file.

    Input is a a FASTA- or FASTQ-formatted file of antibody sequences, the output directory,
    and a dictionary of runtime args.

    Output is a tuple containing (0) path to the output file, (1) the number of successfully
    annotated antibody sequences, (2) path to the log file for successfully annotated sequences (an
    empty string unless args.debug is True) and (3) path to the log file for unsuccessfully annotated
    sequences (only an eompty string if all sequences in the input file were successful).
    """
    try:
        # Args instances can't be serialized by Celery, so we need to pass them in
        # as a dict and re-convert to an Args object inside the Celery job
        args = Args(**arg_dict)
        # identify output file
        output_filename = os.path.basename(seq_file)
        output_suffixes = [
            get_output_suffix(output_type) for output_type in sorted(args.output_type)
        ]
        output_files = [
            os.path.join(output_dir, output_filename + output_suffix)
            for output_suffix in output_suffixes
        ]
        # setup log files
        annotated_logfile = os.path.join(
            log_dir, "temp/{}.annotated".format(output_filename)
        )
        annotated_loghandle = open(annotated_logfile, "a")
        failed_logfile = os.path.join(log_dir, "temp/{}.failed".format(output_filename))
        failed_loghandle = open(failed_logfile, "a")
        unassigned_logfile = os.path.join(
            log_dir, "temp/{}.unassigned".format(output_filename)
        )
        unassigned_loghandle = open(unassigned_logfile, "a")
        # start assignment
        assigner_class = ASSIGNERS[args.assigner]
        assigner = assigner_class(args.germ_db, args.receptor)
        assigner(seq_file, file_format)  # call the assigner
        # process all of the successfully assigned sequences
        outputs_dict = build_output_base(args.output_type)
        assigned = [Antibody(vdj, args.germ_db) for vdj in assigner.assigned]
        successful = 0
        for ab in assigned:
            try:
                ab.annotate(args.uid)
                result = get_abstar_result(
                    ab,
                    pretty=args.pretty,
                    padding=args.padding,
                    raw=args.raw,
                    keys=args.json_keys,
                )
                # results.add_result(result)
                for i, output_type in enumerate(args.output_type):
                    try:
                        output = get_output(result, output_type)
                    except:
                        ab.exception("OUTPUT CREATION ERROR", traceback.format_exc())
                    if output is not None:
                        outputs_dict[output_type].append(
                            get_output(result, output_type)
                        )
                        # only write debug log data once
                        if i == 0:
                            successful += 1
                            if args.debug:
                                annotated_loghandle.write(ab.format_log())
                    # only write failed log data once
                    elif i == 0:
                        failed_loghandle.write(ab.format_log())
            except:
                ab.exception("ANNOTATION ERROR", traceback.format_exc())
                failed_loghandle.write(ab.format_log())
        # outputs = [outputs_dict[ot] for ot in sorted(args.output_type)]
        # write_output(outputs, output_files)
        output_file_dict = write_output(
            outputs_dict, output_dir, output_filename, args.parquet
        )
        # capture the log for all unsuccessful sequences
        for vdj in assigner.unassigned:
            unassigned_loghandle.write(vdj.format_log())
        # close the log handles
        unassigned_loghandle.close()
        annotated_loghandle.close()
        failed_loghandle.close()
        # return the number of successful assignments
        return (
            output_file_dict,
            successful,
            annotated_logfile,
            failed_logfile,
            unassigned_logfile,
        )
    except:
        logging.debug(traceback.format_exc())

        print(traceback.format_exc())


def process_sequences(sequences, args):
    project_temp_dir = tempfile.TemporaryDirectory(dir=args.temp)
    args.project_dir = project_temp_dir.name
    project_temp_dir.cleanup()
    os.makedirs(args.project_dir)

    args.temp = os.path.join(args.project_dir, "temp")
    args.log = os.path.join(args.project_dir, "log")
    args.output = os.path.join(args.project_dir, "output")

    seq_file = tempfile.NamedTemporaryFile(dir=args.project_dir, delete=False)
    seq_file.close()
    with open(seq_file.name, "w") as f:
        f.write("\n".join([s.fasta for s in sequences]))
    args.input = seq_file.name


def run_jobs(files, output_dir, log_dir, file_format, args):
    if args.sequences is not None:
        if args.verbose:
            sys.stdout.write("\nRunning abstar...\n")
    else:
        if args.verbose:
            sys.stdout.write("\nRunning VDJ...\n")
    if args.cluster:
        return _run_jobs_via_celery(files, output_dir, log_dir, file_format, args)
    elif args.debug or args.chunksize == 0:
        return _run_jobs_singlethreaded(files, output_dir, log_dir, file_format, args)
    else:
        return _run_jobs_via_multiprocessing(
            files, output_dir, log_dir, file_format, args
        )


def _run_jobs_singlethreaded(files, output_dir, log_dir, file_format, args):
    results = []
    update_progress(0, len(files))
    for i, f in enumerate(files):
        try:
            results.append(run_abstar(f, output_dir, log_dir, file_format, vars(args)))
            update_progress(i + 1, len(files))
        except:
            logger.debug("FILE-LEVEL EXCEPTION: {}".format(f))
            logging.debug(traceback.format_exc())
    logger.info("")
    return results


def _run_jobs_via_multiprocessing(files, output_dir, log_dir, file_format, args):
    p = Pool(processes=args.num_cores)
    async_results = []
    if args.verbose:
        update_progress(0, len(files))
    for f in files:
        async_results.append(
            (
                f,
                p.apply_async(
                    run_abstar, (f, output_dir, log_dir, file_format, vars(args))
                ),
            )
        )
    monitor_mp_jobs([ar[1] for ar in async_results], print_progress=args.verbose)
    results = []
    for a in async_results:
        try:
            results.append(a[1].get())
        except:
            logger.debug("FILE-LEVEL EXCEPTION: {}".format(a[0]))
            # if args.debug:
            #     traceback.print_exc()
            logging.debug("".join(traceback.format_exc()))
            continue
    p.close()
    p.join()
    return results


def monitor_mp_jobs(results, print_progress=True):
    finished = 0
    jobs = len(results)
    while finished < jobs:
        time.sleep(1)
        ready = [ar for ar in results if ar.ready()]
        finished = len(ready)
        if print_progress:
            update_progress(finished, jobs)
    if print_progress:
        sys.stdout.write("\n\n")


def _run_jobs_via_celery(files, output_dir, log_dir, file_format, args):
    async_results = []
    for f in files:
        async_results.append(
            run_abstar.delay(f, output_dir, log_dir, file_format, vars(args))
        )
    succeeded, failed = monitor_celery_jobs(async_results)
    failed_files = [f for i, f in enumerate(files) if async_results[i].failed()]
    for ff in failed_files:
        logger.debug("FAILED FILE: {}".format(f))
    return [s.get() for s in succeeded]


def monitor_celery_jobs(results, print_progress=True):
    finished = 0
    jobs = len(results)
    while finished < jobs:
        time.sleep(1)
        succeeded = [ar for ar in results if ar.successful()]
        failed = [ar for ar in results if ar.failed()]
        finished = len(succeeded) + len(failed)
        if print_progress:
            update_progress(finished, jobs, failed=len(failed))
    if print_progress:
        sys.stdout.write("\n\n")
    return succeeded, failed


def update_progress(finished, jobs, failed=None):
    pct = int(100.0 * finished / jobs)
    ticks = int(pct / 2)
    spaces = int(50 - ticks)
    if failed:
        prog_bar = "\r({}/{}) |{}{}|  {}% ({}, {})".format(
            finished, jobs, "|" * ticks, " " * spaces, pct, finished - failed, failed
        )
    else:
        prog_bar = "\r({}/{}) |{}{}|  {}%".format(
            finished, jobs, "|" * ticks, " " * spaces, pct
        )
    sys.stdout.write(prog_bar)
    sys.stdout.flush()


#####################################################################
#
#                             MAIN
#
#####################################################################


def run(*args, **kwargs):
    """
    Runs abstar.

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
        - an abtools Sequence object

    .. Caution:: Supplying a single input sequence in list/tuple format is not supported, as abstar 
                 assumes that each element of an iterable is a separate sequence if an iterable 
                 is the only argument. Either convert the list/tuple to a ``Sequence`` object before 
                 calling ``abstar.run()`` or supply a nested list containing the sequence information, 
                 for example: ``[[sequence_id, sequence], ]``.
                 

    Either sequences, ``project_dir``, or all of ``input``, ``output`` and ``temp`` are required.


    Examples:

        If processing a single sequence, you can pass the raw sequence, as a string::

            import abstar

            result = abstar.run('ATGC')

        or as a ``Sequence`` object::

            sequence = Sequence('ATGC', id='seq1')

            result = abstar.run(sequence)

        If you pass just the raw sequence, a random sequence ID will be generated with ``uuid.uuid4()``.
        In either case, when given a single sequence, ``abstar.run()`` will return a single ``Sequence``
        object. If running multiple sequences, you can either pass each sequence as a positional argument::

            result_list = run(['seq1', 'ATGC'], ['seq2', 'CGTA'])

        or you can pass a list of sequences as the first argument, in this case using sequences parsed from a
        FASTA file using Biopython::

            from Bio import SeqIO

            fasta = open('my_sequences.fasta', 'r')
            seqs = [s for s in SeqIO.parse(fasta, 'fasta')]
            result_list = abstar.run(seqs)

        When given multiple sequences, ``abstar.run()`` will return a list of abtools ``Sequence`` objects,
        one per input sequence.

        If you'd prefer not to parse the FASTQ/A file into a list (for example, if the input file is
        extremely large), you can pass the input file path directly, along with a temp directory and output
        directory::

            result_files = abstar.run(input='/path/to/my_sequences.fasta',
                                      temp='/path/to/temp',
                                      output='/path/to/output')

        Given a file path, ``abstar.run()`` returns a list of output file paths. In the above case,
        ``result_files`` will be a list containing a single output file path:
        ``/path/to/output/json/my_sequences.json``.

        If you have a directory containing multiple FASTQ/A files, you can pass the directory path
        using ``input``::

            result_files = abstar.run(input='/path/to/input',
                                      temp='/path/to/temp',
                                      output='/path/to/output')

        As before, ``result_files`` will contain a list of output file paths.

        If your input directory contains paired FASTQ files (gzip compressed or uncompressed)
        that need to be merged prior to processing with abstar::

            result_files = abstar.run(input='/path/to/input',
                                      temp='/path/to/temp',
                                      output='/path/to/output',
                                      merge=True)

        The paired read files in ``input`` will be merged with PANDAseq prior to processing with abstar.
        By default, PANDAseq's 'simple bayesian' read merging algorithm is used, although alternate
        algorithms can be selected with ``pandaseq_algo``.

        abstar provides several output format options. By default, abstar will produce JSON-formatted 
        output file. abstar's output format options include:

        json 
            abstar's default format, in Javascript Object Notation (JSON) format. This format is the most comprehensive. 
            JSON's nesting and inclusion of programmatic objects (such as lists) make this format extremely 
            flexible and well-suited to adaptive immune receptor sequence data, particularly for cases in which
            addtitional fields may be added in the future (such as clonality-related annotations). 

        airr 
            Tab-delimited format compatible with the Adaptive Immune Receptor Repertoires Community's (AIRR-C)
            `schema guidelines`_. This format contains all required fields, several "optional" fields, and 
            several fields that are not part of the schema but conform to the naming conventions of existing schema
            fields (examples include ``v_mutations`` and ``v_mutations_aa``).

        tabular 
            Comma-delimited format containing a subset of the fields contained in the default JSON output format. 
            This format was originally conceived for extremely large datasets, for which output size and compatibility
            with tabular databases (such as MySQL and Apache Spark) were high priorities.
            
        imgt 
            Comma-delimited format that mimics the `IMGT Summary file`_. This output option is provided 
            to minimize the effort needed to convert existing IMGT-based pipelines to abstar. 
            
            
        Multiple output formats can be produced in a single run of abstar, although this is only available when 
        passing an input file or directory; passing individual sequences or a list of sequences (which returns 
        ``Sequence`` objects) can only return a single output format. To produce AIRR output::

            result_files = abstar.run(input='/path/to/input',
                                      temp='/path/to/temp',
                                      output='/path/to/output',
                                      output_type='airr')
        
        To produce both JSON and AIRR-formatted outputs::

            result_files = abstar.run(input='/path/to/input',
                                      temp='/path/to/temp',
                                      output='/path/to/output',
                                      output_type=['json', 'airr'])

        In interactive mode (providing ``Sequence`` objects rather than an input file or directory), returning
        AIRR-formated data can be accomplished by::

            results = abstar.run(sequences, output_type='airr')

        .. _schema guidelines: https://docs.airr-community.org/en/latest/datarep/rearrangements.html
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

        germ_db (str): Germline database to be used. Choices are 'human', 'macaque',
            'mouse', 'humouse', and 'rabbit'. The 'humouse' database contains all germline genes 
            from human and mouse databaes, and is designed to process data from humanized mouse 
            models expressing one or more human germline genes as well as mouse germline genes. 
            Default is 'human'.

        isotype (bool): If True, the isotype will infered by aligning the sequence region
            downstream of the J-gene. If False, the isotype will not be determined.
            Default is True.

        uid (int): Length (in nucleotides) of the Unique Molecular ID used to barcode input RNA.
            A positive integer results in the UMID being parsed from the start of the read (or merged
            read), a negative integer results in parsing from the end of the read. Default is 0,
            which results in no UMID parsing.

        gzip (bool): If True, compresses output files with gzip. Default is False.

        pretty (bool): If True, formats JSON output files to be more human-readable. If False,
            JSON output files contain one record per line. Default is False.

        output_type (str): Options are 'json' 'airr', 'tabular', or 'imgt'. JSON output is the most 
            detailed. Default is 'json'.

        merge (bool): If True, input must be paired-read FASTA files (gzip compressed or uncompressed)
            which will be merged with PANDAseq prior to processing with AbStar. If ``basespace`` is True,
            ``merge`` is automatically set to True. Default is False.

        pandaseq_algo (str): Define merging algorithm to be used by PANDAseq. Options are
            'simple_bayesian', 'ea_util', 'flash', 'pear', 'rdp_mle', 'stitch', or 'uparse'. Default is
            'simple_bayesian', which is the default PANDAseq algorithm.

        debug (bool): If ``True``, ``abstar.run()`` runs in single-threaded mode, the log is much more verbose,
            and temporary files are not removed. Default is ``False``.

        verbose (bool): If ``True``, progress is logged and printed to screen. If ``False``, logging and 
            progress printing are suppressed. Default is ``True``.


    Returns:

        If the input is a single sequence, ``run()`` returns a single abtools ``Sequence`` object.

        If the input is a list of sequences, ``run()`` returns a list of abtools ``Sequence`` objects.

        If the input is a file or a directory of files, ``run()`` returns a list of output files.
    """

    warnings.filterwarnings("ignore")
    if len(args) == 1:
        try:
            if isinstance(args[0], (list, tuple, set)):
                sequences = [Sequence(s) for s in args[0]]
            else:
                sequences = [
                    Sequence(args[0]),
                ]
        except:
            print("ERROR: invalid format for sequence input:")
            for a in args:
                print(a)
            sys.exit(1)
    elif len(args) > 1:
        try:
            sequences = [Sequence(s) for s in args]
        except:
            print("ERROR: invalid format for sequence input:")
            for a in args:
                print(a)
            sys.exit(1)
    else:
        sequences = None

    kwargs["sequences"] = sequences
    args = Args(**kwargs)
    validate_args(args)
    global logger
    if args.verbose:
        logger = log.get_logger("abstar")
        logger.handles = []
    else:
        logger = logging.getLogger("abstar")

    if args.sequences is not None:
        process_sequences(args.sequences, args)
        if len(args.output_type) > 1:
            args.output_type = args.output_type[:1]

    output_files = main(args)

    if args.sequences is not None:
        ofmt = args.output_type[0]
        ofile = output_files[0]
        if ofmt == "tabular":
            output = read_csv(
                ofile, delimiter=",", id_key="seq_id", sequence_key="vdj_nt"
            )
            # with open(ofile) as f:
            # reader = csv.DictReader(f, delimiter=',')
            # output = [Sequence(r) for r in reader]
        if ofmt == "airr":
            output = read_csv(
                ofile, delimiter="\t", id_key="sequence_id", sequence_key="sequence"
            )
            # with open(ofile) as f:
            #     reader = csv.DictReader(f, delimiter='\t')
            #     output = [Sequence(r, id_key='sequence_id', seq_key='sequence')  for r in reader]
        if ofmt == "json":
            output = read_json(ofile, id_key="seq_id", sequence_key="vdj_nt")
            # with open(ofile) as f:
            #     output = [Sequence(json.loads(line)) for line in f]
        # output = [Sequence(o) for o in output]
        if len(output) == 1:
            return output[0]
        else:
            return output
    else:
        return output_files


def run_standalone(args):
    output_dir = main(args)


def main(args):
    input_dir, output_dir, temp_dir, log_dir = make_directories(args)
    if args.sequences is None:
        setup_logging(log_dir, args.debug)
        log_options(input_dir, output_dir, temp_dir, args)
    if args.use_test_data:
        mod_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        input_files = [
            os.path.join(mod_dir, "test_data/test_1k.fasta"),
        ]
    else:
        if args.basespace:
            args.merge = True
            download_files(input_dir)
        if args.merge:
            input_dir = merge_reads(input_dir, args)
        if args.isotype:
            args.isotype = args.germ_db
        input_files = [
            f for f in list_files(input_dir, log=True) if os.stat(f).st_size > 0
        ]
    output_files = []
    for f, fmt in zip(input_files, format_check(input_files)):
        # skip the non-FASTA/Q files
        if fmt is None:
            continue
        start_time = time.time()
        print_input_file_info(f, fmt)
        input_tempdir = os.path.join(temp_dir, "input")
        subfiles, seq_count = split_file(f, fmt, input_tempdir, args)
        run_info = run_jobs(subfiles, temp_dir, log_dir, fmt, args)
        temp_output_file_dicts = [r[0] for r in run_info if r is not None]
        processed_seq_counts = [r[1] for r in run_info if r is not None]
        annotated_log_files = [r[2] for r in run_info if r is not None]
        failed_log_files = [r[3] for r in run_info if r is not None]
        unassigned_log_files = [r[4] for r in run_info if r is not None]
        vdj_end_time = time.time()
        _output_files = concat_outputs(f, temp_output_file_dicts, output_dir, args)
        unassigned_file = concat_logs(f, unassigned_log_files, log_dir, "unassigned")
        failed_file = concat_logs(f, failed_log_files, log_dir, "failed")
        if args.debug:
            annotated_file = concat_logs(f, annotated_log_files, log_dir, "annotated")
        output_files.extend(_output_files)
        if not args.debug:
            flat_temp_files = [
                f for subdict in temp_output_file_dicts for f in subdict.values()
            ]
            clear_temp_files(
                subfiles
                + flat_temp_files
                + annotated_log_files
                + failed_log_files
                + unassigned_log_files
            )
        if args.sequences is not None:
            if args.verbose:
                print_job_stats(
                    seq_count, processed_seq_counts, start_time, vdj_end_time
                )
        else:
            log_job_stats(seq_count, processed_seq_counts, start_time, vdj_end_time)
    return output_files


def run_main(arg_list: Optional[Iterable[str]] = None):
    warnings.filterwarnings("ignore")
    args = create_parser().parse_args(args=arg_list)
    validate_args(args)
    output_dir = main(args)
    sys.stdout.write("\n\n")
    return output_dir


if __name__ == "__main__":
    run_main()
