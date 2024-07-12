#!/usr/bin/env python
# filename: abstar.py

#
# Copyright (c) 2024 Bryan Briney
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
from typing import Iterable, Optional, Union

import abutils
from abutils import Sequence

from ..preprocess.merging import merge_fastqs

#  TODO: inputs/returns
#  --------------------
#
#  - inputs
#    - single Sequence (or something that can be processed by Sequence)
#    - Iterable of Sequence objects
#    - single FASTA/Q file
#    - directory of FASTA/Q files
#    - directory of FASTA/Q files with subdirectories (samples in subdirectories)
#
#  - returns
#    - Sequence objects
#      - should this happen automatically if an output directory is not provided?
#      - running via the API should allow outputs to be returned as Sequence objects or as files
#    - directory of AIRR/parquet files
#      - always return AIRR/parquet files when running at CLI -- maybe a check to ensure an output directory is provided?
#
#
#  the cleanest approach might be to just have two i/o arguments: 1) sequences and 2) project_path:
#
#  - sequences: Iterable[str, Sequence, Iterable]
#    - if str and is a file: process the file
#    - if str and is a directory: process the files in the directory
#    - if str and not either file or directory: assume str is a sequence, convert to Sequence and process
#    - if Sequence or Iterable[Sequence]: process accordingly
#
#  - project_path: Optional[str]
#    - if provided, str should be a directory path into which tmp, log and output files will be deposited
#    - if not provided:
#      - tmp files will be written to /tmp
#      - logs will not be generated
#      - outputs will be returned as Sequence objects


#  TODO: overall workflow
#  ----------------------
#
#  - set up directory structure
#    - user provides a "project path" and we make tmp, log and output (airr and/or parquet) directories within it
#
#  - preprocessing steps
#    - for now, just sequence merging for paired FASTQ inputsm but maybe primer/adapter trimming for FASTA inputs in the future?
#    - all samples are processed before proceeding to assignment
#    - results get deposited in project/merged
#
#  - VDJC assignment
#    - processes an entire input/merged file in a single job (thanks, MMseqs!)
#    - output is a single parquet file, deposited in project/tmp/assignment
#    - log failed assignments in project/logs/assignment as a single txt file per sample
#
#  - split the assignment output into job files
#    - deposited in project/tmp/assignment
#    - once the splitting is done, should call cleanup() on the assigner to remove tmp files that are no longer needed
#
#  - annotation
#    - annotation jobs run in parallel (multiprocessing)
#    - outputs are deposited in project/tmp/annotation
#    - log failed annotations (and successful annotations, if debug=True) as separate failed/succeeeded files in project/log/annotation
#    - temporary output and log files are collected into separate lists so they can be merged and then deleted.
#
#  - create outputs
#    - concat output files (using polars) and write to single tsv (AIRR) and/or parquet files in project/airr and project/parquet
#    - concatenate log files into single failed/succeeded text files in project/logs/annotation
#    - remove all of the temporary output and log files


def run(
    sequences: Optional[Iterable[str, Sequence, Iterable]] = None,
    project_path: Optional[str] = None,
    germline_database: str = "human",
    output_format: Union[str, Iterable[str]] = "airr",
    merge: bool = False,
    merge_kwargs: Optional[dict] = None,
    chunksize: int = 500,
    temp_path: Optional[str] = None,
    quiet: bool = False,
):
    # set up directories
    log_dir = os.path.join(project_path, "log")
    abutils.io.make_dir(log_dir)

    if isinstance(output_format, str):
        output_format = [output_format]


def run_on_sequences(
    sequences: Iterable[str, abutils.Sequence, Iterable],
    project_path: Optional[str] = None,
    germ_db: str = "human",
    output_format: Union[str, Iterable[str]] = "airr",
    chunksize: int = 500,
    temp_path: Optional[str] = None,
    quiet: bool = False,
):
    pass


def run_on_files(
    files: Iterable[str],
    project_path: Optional[str] = None,
    germ_db: str = "human",
    output_format: Union[str, Iterable[str]] = "airr",
    merge: bool = False,
    merge_kwargs: Optional[dict] = None,
    chunksize: int = 500,
    temp_path: Optional[str] = None,
    log_path: Optional[str] = None,
    quiet: bool = False,
):
    if merge:
        # set up merge directories
        merge_dir = os.path.join(project_path, "merged")
        merge_log_dir = os.path.join(log_path, "merge_fastqs")
        abutils.io.make_dir(merge_dir)
        abutils.io.make_dir(merge_log_dir)
        # merge fastqs
        merge_kwargs = merge_kwargs if merge_kwargs is not None else {}
        files = merge_fastqs(files, merge_dir, merge_log_dir, **merge_kwargs)
