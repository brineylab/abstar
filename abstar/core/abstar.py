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
from Bio.SeqRecord import SeqRecord

from ..preprocess.merging import merge_fastqs


def run(
    sequences: Optional[Iterable[str, abutils.Sequence, SeqRecord, Iterable]] = None,
    project_path: Optional[str] = None,
    germ_db: str = "human",
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
    sequences: Iterable[str, abutils.Sequence, SeqRecord, Iterable],
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
        merge_kwargs = merge_kwargs if merge_kwargs is not None else {}
        merge_dir = os.path.join(project_path, "merged")
        merge_log_dir = os.path.join(log_path, "merge_fastqs")
        abutils.io.make_dir(merge_dir)
        abutils.io.make_dir(merge_log_dir)

        # merge fastqs
        merge_fastqs(files, merge_dir, merge_log_dir, **merge_kwargs)
