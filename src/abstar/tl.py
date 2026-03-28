# Copyright (c) 2024 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

"""Analysis tools and utilities.

Namespace module following scanpy conventions. Provides tools for
post-annotation analysis, including UMI parsing.

Usage::

    import abstar

    # parse UMIs from annotated sequences
    abstar.tl.parse_umis(sequences, umi_pattern="NNNNNNNN")
"""

from .annotation.umi import *

__all__ = [
    "parse_umis",
]
