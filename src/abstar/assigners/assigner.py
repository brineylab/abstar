# Copyright (c) 2024 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT


from __future__ import annotations

import logging
from pathlib import Path

import abutils

from ..annotation.germline_alignment import get_germline_database_path

__all__ = ["AssignerBase"]


class AssignerBase:
    """
    docstring for AssignerBase
    """

    def __init__(
        self,
        output_directory: str | Path,
        germdb_name: str,
        receptor: str,
        log_directory: str | Path,
        logger: logging.Logger | None = None,
        concise_logging: bool = False,
        debug: bool = False,
    ):
        self.output_directory = Path(output_directory).resolve()
        self.log_directory = Path(log_directory).resolve()
        self.receptor = receptor
        self.germdb_name = germdb_name
        self.debug = debug
        self.logger = logger if logger is not None else abutils.log.null_logger()
        self.concise_logging = concise_logging
        self.germdb_path = get_germline_database_path(germdb_name, receptor)
        self.to_delete = []  # files that should be deleted during cleanup

        abutils.io.make_dir(self.output_directory)

    def __call__(self, sequence_file: str) -> str:
        """
        All classes subclassing ``AssignerBase`` must implement this method.

        The only expected argument is `sequence_file`, which is a single FASTA or FASTQ-formatted
        file containing input sequences. Gzipped files are supported.

        The return value should be the path to a parquet-formatted file containing the following
        columns:
          * ``sequence_id``: a unique identifier for each input sequence.
          * ``sequence``: the input sequence.
          * ``quality``: the quality score for each input sequence.
          * ``rev_comp``: a boolean indicating whether the input sequence is in the reverse_complement orientation
          * ``v_call``: the V-gene call for each input sequence.
          * ``v_support``: the V-gene evalue for each input sequence.
          * ``d_call``: the D-gene call for each input sequence.
          * ``d_support``: the D-gene evalue for each input sequence.
          * ``j_call``: the J-gene call for each input sequence.
          * ``j_support``: the J-gene evalue for each input sequence.

        """
        raise NotImplementedError("Subclasses must implement this method")

    def cleanup(self):
        abutils.io.delete_files(self.to_delete)
        self.to_delete = []
