#!/usr/bin/env python
# filename: germline.py

#
# Copyright (c) 2016 Bryan Briney
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
from typing import Union

import abutils
from abutils import Sequence
from natsort import natsorted

__all__ = ["get_germdb_path", "get_germline"]


def get_germdb_path(germdb_name: str, receptor: str = "bcr") -> str:
    """
    Get the path to a germline database. The addon directory (which contains user-built
    databases) is checked first, and if a user-built database is not found, the built-in
    database is used if present.

    Parameters
    ----------
    germdb_name : str
        The name of the germline database. Typically the species (like "human"), but can also be a
        custom database (like "humouse") that merges multiple species or is in some other way not
        species-specific.

        .. note::
            All germline database names are lowercase, and identifying germline databases
            by name is case-insensitive.

    receptor : str, default: "bcr"
        The receptor type. Options are "bcr" and "tcr".

    Returns
    -------
    str
        The path to the germline database.

    Raises
    ------
    FileNotFoundError
        If the germline database is not found.

    """
    germdb_name = germdb_name.lower()
    receptor = receptor.lower()
    if receptor not in ["bcr", "tcr"]:
        raise ValueError(f"Receptor type {receptor} not supported")

    # check the addon directory first
    addon_dir = os.path.expanduser(f"~/.abstar/germline_dbs/{receptor}")
    if os.path.isdir(addon_dir):
        if germdb_name.lower() in [os.path.basename(d[0]) for d in os.walk(addon_dir)]:
            return os.path.join(addon_dir, f"{receptor}/{germdb_name}")

    # if a user-generated DB isn't found, use the built-in DB
    core_dir = os.path.dirname(os.path.abspath(__file__))  # "abstar/core" directory
    abstar_dir = os.path.dirname(core_dir)  # "abstar" directory
    germdb_path = os.path.join(abstar_dir, f"germline_dbs/{receptor}/{germdb_name}")
    if not os.path.exists(germdb_path):
        raise FileNotFoundError(
            f"Germline database {germdb_name} for receptor {receptor} not found"
        )
    return germdb_path


def get_germline(
    germline_gene: str,
    germdb_name: str,
    receptor: str = "bcr",
    imgt_gapped: bool = False,
) -> Union[list, Sequence]:
    """
    Get the germline sequence for a given germline gene.

    Parameters
    ----------
    germline_gene : str
        The germline gene to get the sequence for. Can be a full germline gene name including allele ('IGHV1-2*02"), or a truncated
        version corresponding to just the gene ('IGHV1-2') or the family ('IGHV1').

            .. warning:
                Querying with gene names that are a substrings
                of another gene name (for example, 'IGHV1-2' is a substring of 'IGHV1-24') will match and return both
                'IGHV1-2' and 'IGHV1-24' genes. To prevent this, you can add an asterisk to the end of the query gene
                ('IGHV1-2*') which will cause the query to match only alleles of the 'IGHV1-2' gene.

    germdb_name : str
        The name of the germline database.

    receptor : str, default: "bcr"
        The receptor type. Options are "bcr" and "tcr".

    imgt_gapped : bool, default: False
        Whether to use the IMGT gapped or ungapped germline sequences.

    Returns
    -------
    list or Sequence
        If multiple matches are found (for example, multiple alleles of the query 'IGHV1-2'), a list of Sequence objects is returned.
        If only one match is found (for example, when querying 'IGHV1-2*02'), a single Sequence object is returned.

    Raises
    ------
    ValueError
        If the germline gene segment is not one of 'V', 'D', or 'J'. The germline gene segment is inferred as the
        fourth character of the germline gene name. For example, 'IGHV1-2*02' is segment 'V'.

    ValueError
        If the germline gene is not found in the database or if more than one germline
        sequence is found for the same gene.

    """
    germdb_path = get_germdb_path(germdb_name, receptor)
    segment = germline_gene[3].lower()
    if segment not in ["v", "d", "j"]:
        raise ValueError(
            f"The segment type must be one of 'V', 'D', or 'J'. Your input was {germline_gene}, which corresponds to segment type {segment}."
        )

    if imgt_gapped:
        germdb_file = os.path.join(germdb_path, f"imgt_gapped/{segment}.fasta")
    else:
        germdb_file = os.path.join(germdb_path, f"ungapped/{segment}.fasta")
    germs = abutils.io.read_fasta(germdb_file)
    germs = [g for g in germs if germline_gene in g.id]
    if not germs:
        raise ValueError(
            f"Expected at least one germline sequence in {germdb_file}, but found {len(germs)}"
        )
    if len(germs) == 1:
        return germs[0]
    else:
        return natsorted(germs, key=lambda x: x.id)
