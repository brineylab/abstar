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
import subprocess as sp
import sys
from typing import Optional, Union

import abutils
import click
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
        If the germline gene segment is not one of 'V', 'D', 'J', or 'C'. The germline gene segment is inferred as the
        fourth character of the germline gene name. For example, 'IGHV1-2*02' is segment 'V'.

    ValueError
        If the germline gene is not found in the database.

    """
    germdb_path = get_germdb_path(germdb_name, receptor)
    segment = germline_gene[3].lower()
    if segment not in ["v", "d", "j", "c"]:
        raise ValueError(
            f"The segment type must be one of 'V', 'D', 'J', or 'C'. Your input was {germline_gene}, which corresponds to segment type {segment}."
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


# ===============================
#
#   CUSTOM GERMLINE DATABASES
#
# ===============================


@click.command()
@click.option(
    "-v",
    "--variable",
    type=str,
    required=True,
    help="Path to a FASTA-formatted file containing gapped Variable gene sequences. Sequences for both heavy (or alpha and delta) and light (or beta and gamma) chains should be included in a single file.",
)
@click.option(
    "-d",
    "--diversity",
    type=str,
    required=True,
    help="Path to a FASTA-formatted file containing gapped Diversity gene sequences. Sequences for both heavy (or alpha and delta) and light (or beta and gamma) chains should be included in a single file.",
)
@click.option(
    "-j",
    "--joining",
    type=str,
    required=True,
    help="Path to a FASTA-formatted file containing gapped Joining gene sequences. Sequences for both heavy (or alpha and delta) and light (or beta and gamma) chains should be included in a single file.",
)
@click.option(
    "-c",
    "--constant",
    type=str,
    required=False,
    help="Path to a FASTA-formatted file containing gapped Constant gene sequences. Sequences for both heavy (or alpha and delta) and light (or beta and gamma) chains should be included in a single file.",
)
@click.option(
    "-n",
    "--name",
    type=str,
    required=True,
    help="Name of the custom germline database",
)
@click.option(
    "-r",
    "--receptor",
    type=click.Choice(["bcr", "tcr"], case_sensitive=False),
    show_default=True,
    default="bcr",
    help="Receptor type",
)
@click.option(
    "-l",
    "--location",
    type=str,
    default=None,
    help="Location into which the new germline databases will be deposited. This option is provided primarily to test database creation without overwriting current databases of the same name.",
)
@click.option(
    "-m",
    "--manifest",
    type=str,
    show_default=True,
    default=None,
    help="Path to a plain-text file containing information about the germline database. Format is not important, but this is the place for optional inforamtion like the origin of the germline database, the date of download, etc.",
)
@click.option(
    "--verbose/--quiet",
    default=True,
    help="Print verbose output",
)
@click.option(
    "--debug",
    is_flag=True,
    default=False,
    help="Run in debug mode, which will print more verbose information including stdout/stderr from command line tools.",
)
def build_germline_database(
    variable: str,
    diversity: str,
    joining: str,
    name: str,
    receptor: str,
    location: Optional[str] = None,
    manifest: Optional[str] = None,
    constant: Optional[str] = None,
    verbose: bool = True,
    debug: bool = False,
):
    addon_dir = get_addon_directory(location, receptor)
    check_for_existing_db(addon_dir, name)
    make_db_directories(addon_dir, name)
    for segment, input_file in [
        ("variable", variable),
        ("diversity", diversity),
        ("joining", joining),
        ("constant", constant),
    ]:
        if input_file is None:
            continue
        if verbose:
            print_segment_info(segment, input_file)
        gapped_file = make_gapped_db(
            input_file=input_file,
            addon_directory=addon_dir,
            segment=segment[0].lower(),
            dbname=name,
            verbose=verbose,
        )
        ungapped_file = make_ungapped_db(
            input_file=gapped_file,
            addon_directory=addon_dir,
            segment=segment[0].lower(),
            dbname=name,
            verbose=verbose,
        )
        make_mmseqs_db(
            input_file=ungapped_file,
            addon_directory=addon_dir,
            segment=segment[0].lower(),
            dbname=name,
            verbose=verbose,
            debug=debug,
        )

    if manifest is not None:
        if verbose:
            print_manifest_info(manifest)
        transfer_manifest_data(manifest, addon_dir, name)


def get_addon_directory(db_location: Optional[str], receptor: str) -> str:
    """
    Get the path to the addon directory.

    Parameters
    ----------
    db_location : Optional[str]
        The path to the addon directory. If not provided, the default location (~/.abstar/) is used.

    receptor : str
        The receptor type.

    Returns
    -------
    str
        The path to the addon directory.

    """
    if db_location is not None:
        print("\n")
        print(
            "NOTE: You have selected a non-default location for the germline directory."
        )
        string = "abstar only looks in the default location (~/.abstar/) for user-created germline databases, "
        string += "so this database will not be used by abstar. The custom database location option is primarily "
        string += "provided so that users can test the database creation process without overwriting existing databases."
        print(string)
        addon_dir = db_location
    else:
        addon_dir = os.path.expanduser("~/.abstar/germline_dbs")
    addon_dir = os.path.join(addon_dir, receptor.lower())
    abutils.io.make_dir(addon_dir)
    return addon_dir


def check_for_existing_db(addon_dir: str, dbname: str) -> None:
    """
    Check if a germline database already exists in the addon directory.

    Parameters
    ----------
    addon_dir : str
        The path to the addon directory.

    dbname : str
        The name of the germline database.

    """
    dbs = [os.path.basename(d[0]) for d in os.walk(addon_dir)]
    if dbname.lower() in dbs:
        print("\n")
        print("WARNING: A {} germline database already exists.".format(dbname.lower()))
        print("Creating a new database with that name will overwrite the old one.")
        keep_going = input("Do you want to continue? [y/N]: ")
        if keep_going.lower() not in ["y", "yes"]:
            print("")
            print("Aborting germline database creation.")
            print("\n")
            sys.exit()


def make_db_directories(addon_dir: str, dbname: str) -> None:
    """
    Make the main directory for the germline database and the subdirectories for each file type.

    Parameters
    ----------
    addon_dir : str
        The path to the addon directory.

    dbname : str
        The name of the germline database.

    """
    # make the main DB directory
    dbname_dir = os.path.join(addon_dir, dbname.lower())
    if not os.path.isdir(dbname_dir):
        abutils.io.make_dir(dbname_dir)

    # make subdirectories
    db_names = ["imgt_gapped", "ungapped", "mmseqs"]
    for db_name in db_names:
        db_dir = os.path.join(dbname_dir, db_name)
        abutils.io.make_dir(db_dir)


def transfer_manifest_data(manifest: str, addon_directory: str, dbname: str) -> str:
    """
    Transfer manifest data to the new germline database.

    Parameters
    ----------
    manifest : str
        The path to the manifest file.

    addon_directory : str
        The path to the addon directory.

    dbname : str
        The name of the germline database.

    Returns
    -------
    str
        The path to the manifest file.

    """
    # read manifest data
    with open(manifest, "r") as f:
        manifest_data = f.read()

    # write to manifest file
    manifest_file = os.path.join(addon_directory, f"{dbname.lower()}/manifest.txt")
    with open(manifest_file, "w") as f:
        f.write(manifest_data)
    return manifest_file


def make_gapped_db(
    input_file: str,
    addon_directory: str,
    segment: str,
    dbname: str,
    verbose: bool = False,
) -> str:
    """
    Make the IMGT-gapped FASTA file for a given segment.

    Parameters
    ----------
    input_file : str
        The path to the input FASTA file.

    addon_directory : str
        The path to the addon directory.

    segment : str
        The segment type.

    dbname : str
        The name of the germline database.

    verbose : bool, default: False
        Whether to print verbose output.

    Returns
    -------
    str
        The path to the IMGT-gapped database file.

    """
    if verbose:
        print("  - IMGT-gapped FASTA")

    # read input sequences
    seqs = sorted(abutils.io.read_fasta(input_file), key=lambda x: x.id)
    fastas = [f">{s.id}\n{s.seq.upper()}" for s in seqs]

    # write to the output DB file
    output_file = os.path.join(
        addon_directory,
        f"{dbname.lower()}/imgt_gapped/{segment.lower()}.fasta",
    )
    with open(output_file, "w") as f:
        f.write("\n".join(fastas))
    return output_file


def make_ungapped_db(
    input_file: str,
    addon_directory: str,
    segment: str,
    dbname: str,
    verbose: bool = False,
) -> str:
    """
    Make the IMGT-gapped FASTA file for a given segment.

    Parameters
    ----------
    input_file : str
        The path to the input FASTA file.

    addon_directory : str
        The path to the addon directory.

    segment : str
        The segment type.

    dbname : str
        The name of the germline database.

    verbose : bool, default: False
        Whether to print verbose output.

    Returns
    -------
    str
        The path to the ungapped database file.

    """
    if verbose:
        print("  - ungapped FASTA")

    # read input sequences
    seqs = sorted(abutils.io.read_fasta(input_file), key=lambda x: x.id)
    fastas = [f">{s.id}\n{s.seq.upper().replace('.', '')}" for s in seqs]

    # write to the output DB file
    output_file = os.path.join(
        addon_directory,
        f"{dbname.lower()}/ungapped/{segment.lower()}.fasta",
    )
    with open(output_file, "w") as f:
        f.write("\n".join(fastas))
    return output_file


def make_mmseqs_db(
    input_file: str,
    addon_directory: str,
    segment: str,
    dbname: str,
    verbose: bool = False,
    debug: bool = False,
) -> str:
    """
    Make the MMseqs2 database for a given segment.

    Parameters
    ----------
    input_file : str
        The path to the input (ungapped) FASTA file.

    addon_directory : str
        The path to the addon directory.

    segment : str
        The segment type.

    dbname : str
        The name of the germline database.

    verbose : bool, default: False
        Whether to print verbose output.

    debug : bool, default: False
        Whether to print debug output.

    Returns
    -------
    str
        The path to the ungapped database file.

    """
    if verbose:
        print("  - MMseqs2 database")

    # create MMseqs2 database
    output_file = os.path.join(
        addon_directory,
        f"{dbname.lower()}/mmseqs/{segment.lower()}",
    )
    createdb_cmd = f"mmseqs createdb {input_file} {output_file}"
    p = sp.Popen(createdb_cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    if debug:
        print(createdb_cmd)
        print(stdout)
        print(stderr)

    # create MMseqs2 index
    createindex_cmd = f"mmseqs createindex {output_file} /tmp"
    p = sp.Popen(createindex_cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    if debug:
        print(createindex_cmd)
        print(stdout)
        print(stderr)

    return output_file


def print_segment_info(segment: str, input_file: str) -> None:
    """
    Print information about the segment (variable, diversity, joining, constant).
    """
    seqs = abutils.io.read_fasta(input_file)
    seg_string = "  " + segment.upper() + "  "
    print("\n")
    print("-" * len(seg_string))
    print(seg_string)
    print("-" * len(seg_string))
    print(input_file)
    print("input file contains {} sequences".format(len(seqs)))
    print("")
    print("Building germline databases:")


def print_manifest_info(manifest: str) -> None:
    seg_string = "  MANIFEST  "
    print("\n")
    print("-" * len(seg_string))
    print(seg_string)
    print("-" * len(seg_string))
    print(manifest)
    print("")
    print("Transferring manifest data...")
