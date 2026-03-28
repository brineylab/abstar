# Copyright (c) 2024 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

"""Preprocessing utilities for sequencing data.

from __future__ import annotations

Namespace module following scanpy conventions. Provides functions for
preprocessing raw sequencing data before annotation, including paired-end
read merging via fastp.

Usage::

    import abstar

    # merge paired-end FASTQ files
    abstar.pp.merge_fastqs(forward, reverse, merged_output)

    # group paired FASTQ files by sample
    pairs = abstar.pp.group_paired_fastqs("/path/to/fastq_dir")
"""

from .preprocess.merging import *

__all__ = [
    "merge_fastqs",
    "group_paired_fastqs",
]
