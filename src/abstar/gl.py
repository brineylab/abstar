# Copyright (c) 2024 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

"""Germline database access and management.

Namespace module following scanpy conventions. Provides functions for
querying germline gene databases and building custom germline databases.

Usage::

    import abstar

    # get the path to a built-in germline database
    db_path = abstar.gl.get_germline_database_path("human")

    # retrieve a single germline gene sequence
    germline = abstar.gl.get_germline("IGHV1-2*02", "human")

    # build a custom germline database
    abstar.gl.build_germline_database("custom_db", fasta_files=["seqs.fasta"])
"""

from .annotation.germline_alignment import *
from .core.germline_builder import *

__all__ = [
    # from annotation.germline_alignment
    "get_germline_database_path",
    "get_germline",
    "realign_germline",
    "reassign_dgene",
    "process_vgene_alignment",
    "process_jgene_alignment",
    "process_dgene_alignment",
    # from core.germline_builder
    "build_germline_database",
]
