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


import datetime
import difflib
import json
import os
import shutil
import subprocess as sp
import sys
from concurrent.futures import ThreadPoolExecutor
from typing import Iterable, Optional, Union

# from weakref import ref
import abutils
import parasail
from abutils import Sequence

# import click
# from natsort import natsorted
# from networkx import max_flow_min_cost
from tqdm.auto import tqdm

from ..gl import get_germline, get_germline_database_path
from ..utils import MATRIX_PATH

__all__ = ["build_germline_database"]

# ===============================
#
#   CUSTOM GERMLINE DATABASES
#
# ===============================


#  TODO: inputs/returns
#  --------------------
#
#  - inputs
#    - one or more FASTA or JSON files containing VDJ gene segments
#       - probably easiest to use a separate command-line flag for each file type that can be used multiple times -- "-f" for FASTA and "-j" for JSON
#       - no requirement for homogeneity of gene segments in the input files -- it doesn't have to be one file for V, one for J, etc
#       - combining heterogeneity with multiple JSON/FASTA files should make it much easier to build multi-species databases
#    - separate FASTA file(s) containing constant regions (since D genes and IgD genes both start with "IGHD")
#       - no need to support JSON at this point, because OGRDB doesn't have constant regions (yet)
#       - should be able to use a single flag that can be used multiple times -- "-c"?
#    - should also copy the unprocessed input data into a subdirectory ("raw"?) so that it's linked to the resulting germline database
#


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


def build_germline_database(
    name: str,
    fastas: Optional[Union[str, Iterable[str]]] = None,
    jsons: Optional[Union[str, Iterable[str]]] = None,
    constants: Optional[Union[str, Iterable[str]]] = None,
    receptor: str = "bcr",
    manifest: Optional[str] = None,
    include_species_in_name: bool = True,
    location: Optional[str] = None,
    reference: str = "human",
    verbose: bool = True,
    debug: bool = False,
) -> None:
    """
    Builds a custom germline database.

    Parameters
    ----------
    name : str
        The name of the germline database. This will be how the database can be invoked when running ``abstar``.

        .. warning::
            Although custom databases are stored in a different location than built-in databases (and thus will not actually
            overwrite them on disk), custom databases are given priority over built-in databases when running ``abstar``.
            This means that a custom database named ``human`` will be used instead of the built-in ``human`` database.

    fastas : str | Iterable[str], optional
        The path to the FASTA file(s) containing the VDJ segments, either as a single file or a list of files.

    jsons : str | Iterable[str], optional
        The path to the JSON file(s) containing the VDJ segments, either as a single file or a list of files.

    constants : str | Iterable[str], optional
        The path to the FASTA file containing the constant regions, either as a single file or a list of files.

    receptor : str, default: "bcr"
        The type of receptor to build the database for.

    manifest : str, optional
        The path to the manifest file containing the metadata for the germline database.

    include_species_in_name : bool, default: True
        Whether to include the species name in the database name.

    location : str, optional
        The location of the germline database.

    reference : str, default: "human"
        The reference species to use for adding IMGT gaps. Only used if ungapped VDJ sequences are provided.

    verbose : bool, default: True
        Whether to print verbose output.

    debug : bool, default: False
        Whether to print debug output.

    Returns
    -------
    None
        The function does not return anything.

    """
    gapped_vdjs = []
    gapped_constants = []

    # set up database directory structure
    database_dir = get_database_directory(location, receptor)
    if check_for_existing_db(name, receptor, database_dir):
        confirm_overwrite_existing_db(name)
    database_dir = os.path.join(database_dir, name.lower())
    sub_dirs = make_db_directories(database_dir)
    raw_dir = sub_dirs["raw"]
    gapped_dir = sub_dirs["imgt_gapped"]
    ungapped_dir = sub_dirs["ungapped"]
    mmseqs_dir = sub_dirs["mmseqs"]

    # process FASTA-formatted VDJ segments
    if fastas:
        if verbose:
            print("processing FASTA-formatted VDJ segments")
        if isinstance(fastas, str):
            fastas = [fastas]
        for fasta in fastas:
            if not os.path.isfile(fasta):
                raise FileNotFoundError(f"The file {fasta} does not exist.")
            copy_to_raw(fasta, raw_dir)
            sequences = process_fasta(
                fasta_file=fasta, include_species_in_name=include_species_in_name
            )
            gapped_vdjs.extend(sequences)

    # process JSON-formatted VDJ segments
    if jsons:
        if verbose:
            print("processing JSON-formatted VDJ segments")
        if isinstance(jsons, str):
            jsons = [jsons]
        for _json in jsons:
            if not os.path.isfile(_json):
                raise FileNotFoundError(f"The file {_json} does not exist.")
            copy_to_raw(_json, raw_dir)
            sequences = process_json(
                _json, include_species_in_name=include_species_in_name
            )
            gapped_vdjs.extend(sequences)

    # add IMGT gapps (if they don't already exist)
    gapped_vdjs = add_imgt_gaps(gapped_vdjs, reference=reference)

    # process FASTA-formatted constant regions
    if constants:
        if verbose:
            print("processing FASTA-formatted constant regions")
        if isinstance(constants, str):
            constants = [constants]
        for constant in constants:
            if not os.path.isfile(constant):
                raise FileNotFoundError(f"The file {constant} does not exist.")
            copy_to_raw(constant, raw_dir)
            sequences = process_fasta(
                fasta_file=constant, include_species_in_name=False
            )
            gapped_constants.extend(sequences)

    # IMGT-gapped database
    if verbose:
        print("")
        print("building IMGT-gapped database")
    make_fasta_dbs(
        vdjs=gapped_vdjs,
        constants=gapped_constants,
        database_dir=gapped_dir,
        verbose=verbose,
    )

    # ungapped database
    if verbose:
        print("")
        print("building ungapped database")
    ungapped_vdjs = [
        Sequence(s.sequence.replace(".", ""), id=s.id) for s in gapped_vdjs
    ]
    ungapped_constants = [
        Sequence(s.sequence.replace(".", ""), id=s.id) for s in gapped_constants
    ]
    make_fasta_dbs(
        vdjs=ungapped_vdjs,
        constants=ungapped_constants,
        database_dir=ungapped_dir,
        verbose=verbose,
    )

    # MMseqs database
    if verbose:
        print("")
        print("building MMseqs database")
    make_mmseqs_dbs(
        database_dir=mmseqs_dir,
        ungapped_dir=ungapped_dir,
        verbose=verbose,
        debug=debug,
    )

    # manifest
    if manifest is not None:
        if verbose:
            print("")
            print("transferring manifest data")
        transfer_manifest_data(manifest, database_dir)


def add_imgt_gaps(
    germlines: abutils.Sequence | Iterable[abutils.Sequence], reference: str = "human"
) -> Iterable[abutils.Sequence]:
    """
    Add IMGT gaps to germline sequences.

    Parameters
    ----------
    germlines : abutils.Sequence | Iterable[abutils.Sequence]
        The germline sequences to add IMGT gaps to.

        .. note::
            Gaps are only added to V genes. D and J genes are always returned unchanged.

    reference : str, default: "human"
        The reference species to use for the IMGT gaps.

    Returns
    -------
    Iterable[abutils.Sequence]
        The germline sequences with IMGT gaps added. If any of the V gene sequences already
        have gaps, all sequences are returned unchanged.

    """
    if isinstance(germlines, abutils.Sequence):
        germlines = [germlines]
    # check for gaps -- if they're present, return the gapped sequences
    vgenes = [g for g in germlines if g.id[3] == "V"]
    others = [g for g in germlines if g.id[3] != "V"]
    if any(["." in v.sequence for v in vgenes]):
        return germlines

    # get IMGT-gapped germline V genes
    germs = get_germline("IGHV", germdb_name=reference, imgt_gapped=True)

    # build parasail matrix that more heavily penalizes mismatches to gaps
    matrix_file = os.path.join(MATRIX_PATH, "imgt_gapped.txt")
    matrix = parasail.Matrix(matrix_file)

    # add IMGT gaps
    gapped_vgenes = []
    for ungapped in vgenes:
        # find the best germline match
        alns = abutils.tl.semiglobal_alignment(
            ungapped, targets=germs, matrix=matrix, gap_open=-25
        )
        top_aln = alns[0]
        # build the gapped sequence
        gapped = ""
        for q, t in zip(top_aln.aligned_query, top_aln.aligned_target):
            if q == "-":
                gapped += "."
            else:
                gapped += q
        # remove trailing gaps, but keep any 5' gaps to preserve IMGT position information
        while True:
            if gapped[-1] == ".":
                gapped = gapped[:-1]
            else:
                break
        gapped_vgenes.append(abutils.Sequence(gapped, id=ungapped.id))
    return gapped_vgenes + others


# -------------------------
#   GENERATING CUSTOM DATABASE FROM IgDiscover
# -------------------------


# def build_germdb_from_igdiscover(
#     name: str,
#     igdiscover_output: str,
#     constants: Optional[str] = None,
#     receptor: str = "bcr",
#     species: str = "human",
#     location: Optional[str] = None,
#     verbose: bool = True,
#     debug: bool = False,
# ) -> None:
#     """
#     Builds a custom reference database using the output from IgDiscover

#     """

#     abutils.io.make_dir("/tmp/refs")

#     origin = get_germline_database_path(receptor=receptor, germdb_name=species)
#     shutil.copyfile(os.path.join(origin, "manifest.txt"), "/tmp/refs/manifest.txt")

#     with open("/tmp/refs/manifest.txt", "a") as f:
#         f.write("\n\n")
#         f.write("/!\\ CUSTOMIZED REFERENCE SPECIFIC TO DONOR /!\\\n\n")
#         f.write("Custom database modified with IgDiscover output\n")
#         f.write(f"Donor identifier: {name}\n")
#         f.write(f"Database created on {str(datetime.date.today())}\n")

#     files = abutils.io.list_files(igdiscover_output, extension="fasta")
#     for file_in in files:
#         filename = os.path.basename(file_in).lower()
#         file_out = os.path.join("/tmp/refs/", filename)

#         if "v" in filename:  # We only need to gap the V gene file
#             gapped_sequences = []
#             sequences = abutils.io.read_fasta(file_in)
#             v_refs = abutils.io.read_fasta(
#                 os.path.join(origin, "imgt_gapped", "v.fasta")
#             )
#             v_refs_hyphen = [
#                 Sequence(s.sequence.replace(".", "-"), id=s.id) for s in v_refs
#             ]

#             with ThreadPoolExecutor() as executor:
#                 if verbose:
#                     gapped_sequences = list(
#                         tqdm(
#                             executor.map(
#                                 lambda seq: pairwise_gap_sequence(
#                                     seq, v_refs_hyphen, verbose=verbose
#                                 ),
#                                 sequences,
#                             ),
#                             desc="Gapping new sequences...",
#                             total=len(sequences),
#                         )
#                     )

#                 else:
#                     gapped_sequences = list(
#                         executor.map(
#                             lambda seq: pairwise_gap_sequence(
#                                 seq, v_refs_hyphen, verbose=verbose
#                             ),
#                             sequences,
#                         ),
#                     )

#             with open(file_out, "w") as f:
#                 for s in gapped_sequences:
#                     f.write(s.fasta)
#                     f.write("\n")

#         else:  # D and J gene fils don't need to be gapped
#             shutil.copyfile(file_in, file_out)

#     # Constant genes are never part of the IgDiscover output, hence they need to be imported
#     if constants == None:
#         shutil.copyfile(
#             os.path.join(origin, "imgt_gapped", "c.fasta"), "/tmp/refs/c.fasta"
#         )

#     fastas = abutils.io.list_files("/tmp/refs/", extension="fasta")

#     # Building the germline_database using prepared files
#     build_germline_database(
#         name=name,
#         fastas=fastas,
#         constants=constants,
#         receptor=receptor,
#         manifest="/tmp/refs/manifest.txt",
#         include_species_in_name=False,
#         location=location,
#         verbose=verbose,
#         debug=debug,
#     )

#     if not debug:
#         shutil.rmtree("/tmp/refs")

#     return


# # def gap_sequence(sequence, reference, gaps='.'):
# #     to_align = [sequence, ] + reference
# #     aln = abutils.tools.alignment.mafft(sequences = to_align, mafft_bin='mafft')
# #     gapped_seq = [s for s in aln if s.id == sequence.id][0].sequence.replace("-", gaps)
# #     gapped = Sequence(gapped_seq, id=sequence.id)

# #     return gapped


# def pairwise_gap_sequence(sequence, references, gaps=".", verbose=False):
#     gene_name = sequence.id.split("_")[0]

#     # Getting the correct reference gene from IMGT to gap new sequence
#     try:
#         reference = [s for s in references if s.id == gene_name][0]
#     except:
#         match = difflib.get_close_matches(
#             gene_name, [r.id for r in references], n=1, cutoff=0.6
#         )[0]
#         if verbose:
#             print(
#                 f"No exact match for {gene_name}. Using {match} to perform alignement..."
#             )
#         reference = [s for s in references if s.id == match][0]

#     # Performing the alignment
#     aln = abutils.aln.semiglobal_alignment(
#         query=sequence, target=reference, mismatch=-5
#     )
#     gapped_seq = aln.aligned_query.replace("-", gaps)
#     gapped = Sequence(gapped_seq, id=sequence.id)

#     return gapped


# -------------------------
#   DATABASE DIRECTORIES
# -------------------------


def get_database_directory(receptor: str, db_location: Optional[str]) -> str:
    """
    Get the path to the receptor-level germline database directory.

    Parameters
    ----------
    db_location : Optional[str]
        The path to the addon directory. If not provided, the default location (~/.abstar/) will be used.

    receptor : str
        The receptor type.

    Returns
    -------
    str
        The path to the addon directory.

    """
    if db_location is not None:
        print("")
        print(
            "NOTE: You have selected a non-default location for the germline directory."
        )
        string = "abstar only looks in the default location (~/.abstar/) for user-created germline databases, "
        string += "so this database will not be used by abstar. The custom database location option is primarily "
        string += "provided so that users can test the database creation process without overwriting existing databases.\n"
        print(string)
        database_dir = db_location
    else:
        database_dir = os.path.expanduser("~/.abstar/germline_dbs")
    database_dir = os.path.join(database_dir, receptor.lower())
    abutils.io.make_dir(database_dir)
    return database_dir


def check_for_existing_db(
    name: str, receptor: str, location: Optional[str] = None
) -> bool:
    """
    Check if a germline database already exists in the addon directory.

    Parameters
    ----------
    name : str
        The name of the germline database.

    receptor : str
        The receptor type.

    location : Optional[str], default: None
        The path to a non-standard directory that may contain the germline database.

    Returns
    -------
    bool
        True if the germline database already exists, False otherwise.

    """
    if location is None:
        location = get_database_directory(receptor)
    dbs = [os.path.basename(d[0]) for d in os.walk(location)]
    if name.lower() in dbs:
        return True
    else:
        return False


def confirm_overwrite_existing_db(name: str) -> bool:
    """
    Confirm that the user wants to overwrite an existing germline database.

    Parameters
    ----------
    name : str
        The name of the germline database.

    Returns
    -------
    bool
        True if the user wants to continue, False otherwise.
    """
    print("")
    print(f"WARNING: A {name.lower()} germline database already exists.")
    print("Creating a new database with that name will overwrite the old one.")
    keep_going = input("Do you want to continue? [y/N]: ")
    if keep_going.lower() not in ["y", "yes"]:
        print("")
        print("Aborting germline database creation.")
        print("\n")
        sys.exit()
    else:
        print("")


def make_db_directories(database_dir: str) -> None:
    """
    Make the main directory for the germline database and the subdirectories for each file type.

    Parameters
    ----------
    database_dir : str
        The path to the database directory.

    """
    # make the main DB directory
    abutils.io.make_dir(database_dir)

    # make subdirectories
    subdirs = {}
    subdir_names = ["raw", "imgt_gapped", "ungapped", "mmseqs"]
    for subdir_name in subdir_names:
        subdir = os.path.join(database_dir, subdir_name)
        abutils.io.make_dir(subdir)
        subdirs[subdir_name] = subdir
    return subdirs


def transfer_manifest_data(manifest_file: str, database_dir: str) -> None:
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
    with open(manifest_file, "r") as f:
        manifest_data = f.read()

    # write to manifest file
    manifest_file = os.path.join(database_dir, "manifest.txt")
    with open(manifest_file, "w") as f:
        f.write(manifest_data)


def copy_to_raw(fasta: str, raw_dir: str) -> None:
    """
    Copy a FASTA file to the raw directory.

    Parameters
    ----------
    fasta : str
        The path to the input FASTA file.

    raw_dir : str
        The path to the raw directory.
    """
    shutil.copy(fasta, raw_dir)


# -------------------------
#   PROCESS INPUT FILES
# -------------------------


def process_fasta(
    fasta_file: str, include_species_in_name: bool = True
) -> Iterable[Sequence]:
    """
    Process a FASTA file and return a list of Sequence objects.
    """
    seqs = abutils.io.read_fasta(fasta_file)
    if not include_species_in_name:
        for s in seqs:
            s.id = s.id.split("__")[0]
    return seqs


def process_json(
    json_file: str, include_species_in_name: bool = True
) -> Iterable[Sequence]:
    """
    Process a JSON file and return a list of Sequence objects.
    """
    seqs = []
    with open(json_file, "r") as f:
        jdata = json.load(f)
    for entry in jdata["GermlineSet"][0]["allele_descriptions"]:
        name = entry["label"]
        if entry["sequence_type"] == "V":
            # V-genes are the only ones with IMGT-gapped sequences,
            # and it's in a different location than D/J sequences
            gapped = [
                d
                for d in entry["v_gene_delineations"]
                if d["delineation_scheme"] == "IMGT"
            ][0]["aligned_sequence"]
        else:
            gapped = entry["coding_sequence"]
        species = entry["species"]["label"].lower().replace(" ", "_")
        if include_species_in_name:
            name = f"{name}__{species}"
        seqs.append(Sequence(gapped, id=name))
    return seqs


# -------------------------
#     MAKE DATABASES
# -------------------------


def make_fasta_dbs(
    vdjs: Iterable[Sequence],
    constants: Iterable[Sequence],
    database_dir: str,
    verbose: bool = False,
) -> None:
    """
    Make the IMGT-gapped and ungapped FASTA databases.

    Parameters
    ----------
    vdjs : Iterable[Sequence]
        The VDJ genes.

    constants : Iterable[Sequence]
        The constant regions.

    database_dir : str
        The path to the database directory.

    verbose : bool, default: False
        Whether to print verbose output.
    """
    # VDJ genes
    for segment in ["V", "D", "J"]:
        if verbose:
            if segment == "V":
                print("  V", end="")
            else:
                print(f" | {segment}", end="")
        seqs = [s for s in vdjs if s.id[3] == segment]
        if seqs:
            abutils.io.to_fasta(
                seqs, os.path.join(database_dir, f"{segment.lower()}.fasta")
            )

    # constant regions
    if constants:
        if verbose:
            print(" | C")
        abutils.io.to_fasta(constants, os.path.join(database_dir, "c.fasta"))


def make_mmseqs_dbs(
    database_dir: str,
    ungapped_dir: str,
    verbose: bool = False,
    debug: bool = False,
) -> None:
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

    """
    # VDJ genes
    for segment in ["V", "D", "J"]:
        if verbose:
            if segment == "V":
                print("  V", end="")
            else:
                print(f" | {segment}", end="")
        ungapped_file = os.path.join(ungapped_dir, f"{segment.lower()}.fasta")
        output_file = os.path.join(database_dir, f"{segment.lower()}")
        _make_mmseqs_db(ungapped_file, output_file, debug=debug)

    # constant regions
    segment = "C"
    ungapped_file = os.path.join(ungapped_dir, f"{segment.lower()}.fasta")
    if os.path.exists(ungapped_file):
        if verbose:
            print(f" | {segment}")
        output_file = os.path.join(database_dir, f"{segment.lower()}")
        _make_mmseqs_db(ungapped_file, output_file, debug=debug)
    elif verbose:
        print("")


def _make_mmseqs_db(input_file: str, output_file: str, debug: bool = False) -> None:
    """
    Make an MMseqs2 database.

    Parameters
    ----------
    input_file : str
        The path to the input file.

    output_file : str
        The path to the output file.

    debug : bool, default: False
        Whether to print debug output.
    """
    # create MMseqs2 database
    mmseqs_bin = abutils.bin.get_path("mmseqs")
    createdb_cmd = f"{mmseqs_bin} createdb {input_file} {output_file}"
    p = sp.Popen(createdb_cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    if debug:
        print(createdb_cmd)
        print(stdout)
        print(stderr)

    # create MMseqs2 index
    createindex_cmd = f"{mmseqs_bin} createindex {output_file} /tmp"
    p = sp.Popen(createindex_cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    if debug:
        print(createindex_cmd)
        print(stdout)
        print(stderr)


# def print_segment_info(segment: str, input_file: str) -> None:
#     """
#     Print information about the segment (variable, diversity, joining, constant).
#     """
#     seqs = abutils.io.read_fasta(input_file)
#     seg_string = "  " + segment.upper() + "  "
#     print("\n")
#     print("-" * len(seg_string))
#     print(seg_string)
#     print("-" * len(seg_string))
#     print(input_file)
#     print("input file contains {} sequences".format(len(seqs)))
#     print("")
#     print("Building germline databases:")


# def print_manifest_info(manifest: str) -> None:
#     seg_string = "  MANIFEST  "
#     print("\n")
#     print("-" * len(seg_string))
#     print(seg_string)
#     print("-" * len(seg_string))
#     print(manifest)
#     print("")
#     print("Transferring manifest data...")
