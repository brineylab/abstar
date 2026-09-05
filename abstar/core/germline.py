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
import fcntl
import json
import os
import shlex
import shutil
import stat
import subprocess as sp
import sys
import tempfile
import uuid
import warnings
from contextlib import contextmanager
from concurrent.futures import ThreadPoolExecutor
from typing import Iterable

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


class GermlineBuildExternalToolError(RuntimeError):
    """A custom-database build command failed with inspectable diagnostics."""

    def __init__(self, message, *, command=(), returncode=None, stdout="", stderr=""):
        self.command = tuple(command)
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr
        if self.command:
            message += (
                f"\nCOMMAND: {shlex.join(self.command)}\nEXIT STATUS: {returncode}"
                f"\nSTDOUT:\n{stdout}\nSTDERR:\n{stderr}"
            )
        super().__init__(message)


class GermlinePublicationError(RuntimeError):
    """Publication and rollback both failed, leaving a recoverable backup."""

    def __init__(self, primary_error, rollback_error, backup_path):
        self.primary_error = primary_error
        self.rollback_error = rollback_error
        self.backup_path = backup_path
        super().__init__(
            f"database publication failed: {primary_error}; rollback failed: "
            f"{rollback_error}; existing database retained at {backup_path}"
        )

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
    fastas: str | Iterable[str] | None = None,
    jsons: str | Iterable[str] | None = None,
    constants: str | Iterable[str] | None = None,
    receptor: str = "bcr",
    manifest: str | None = None,
    include_species_in_name: bool = True,
    location: str | None = None,
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
    receptor = receptor.lower()
    if receptor not in {"bcr", "tcr"}:
        raise ValueError("receptor must be 'bcr' or 'tcr'")
    if not name or name in {".", ".."} or os.path.basename(name) != name:
        raise ValueError("database name must be a non-empty path-safe name")

    database_root = get_database_directory(receptor, location)
    normalized_name = name.lower()
    database_dir = os.path.join(database_root, normalized_name)
    initial_destination = inspect_database_destination(database_dir)
    replacing = initial_destination is not None
    if replacing:
        confirm_overwrite_existing_db(name)

    staging_dir = tempfile.mkdtemp(prefix=f".{normalized_name}.staging-", dir=database_root)
    gapped_vdjs = []
    gapped_constants = []
    try:
        sub_dirs = make_db_directories(staging_dir)
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

        # add IMGT gaps (if they don't already exist)
        validate_germlines(gapped_vdjs, gapped_constants, receptor)
        gapped_vdjs = add_imgt_gaps(
            gapped_vdjs, reference=reference, receptor=receptor
        )

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
                    fasta_file=constant, include_species_in_name=include_species_in_name
                )
                gapped_constants.extend(sequences)

        validate_germlines(gapped_vdjs, gapped_constants, receptor)

        # IMGT-gapped database
        if verbose:
            print("")
            print("building IMGT-gapped database")
        make_fasta_dbs(gapped_vdjs, gapped_constants, gapped_dir, verbose)

        # ungapped database
        if verbose:
            print("")
            print("building ungapped database")
        ungapped_vdjs = [
            Sequence(s.sequence.replace(".", ""), id=s.id) for s in gapped_vdjs
        ]
        ungapped_constants = [
            Sequence(s.sequence.replace(".", ""), id=s.id)
            for s in gapped_constants
        ]
        make_fasta_dbs(ungapped_vdjs, ungapped_constants, ungapped_dir, verbose)

        # MMseqs database
        if verbose:
            print("")
            print("building MMseqs database")
        make_mmseqs_dbs(mmseqs_dir, ungapped_dir, verbose, debug)

        # manifest
        if manifest is not None:
            if not os.path.isfile(manifest):
                raise FileNotFoundError(f"The file {manifest} does not exist.")
            transfer_manifest_data(manifest, staging_dir)

        validate_staged_database(
            staging_dir, receptor, require_manifest=manifest is not None
        )
        with database_build_lock(database_root, normalized_name):
            current_destination = inspect_database_destination(database_dir)
            if current_destination != initial_destination:
                raise RuntimeError(
                    f"germline database destination {database_dir} changed during build"
                )
            publish_database(staging_dir, database_dir, replacing)
        staging_dir = None
    finally:
        if staging_dir is not None:
            shutil.rmtree(staging_dir, ignore_errors=True)


def add_imgt_gaps(
    germlines: abutils.Sequence | Iterable[abutils.Sequence],
    reference: str = "human",
    receptor: str = "bcr",
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
        Germline sequences with IMGT gaps added to each ungapped V gene. Already-gapped
        V genes and all D/J genes are returned unchanged.

    """
    if isinstance(germlines, abutils.Sequence):
        germlines = [germlines]
    vgenes = [g for g in germlines if g.id[3] == "V"]
    others = [g for g in germlines if g.id[3] != "V"]

    # build parasail matrix that more heavily penalizes mismatches to gaps
    matrix_file = os.path.join(MATRIX_PATH, "imgt_gapped.txt")
    matrix = parasail.Matrix(matrix_file)

    # add IMGT gaps
    gapped_vgenes = []
    for ungapped in vgenes:
        if "." in ungapped.sequence:
            gapped_vgenes.append(ungapped)
            continue
        germs = get_germline(
            f"{ungapped.id[:3]}V",
            germdb_name=reference,
            receptor=receptor,
            imgt_gapped=True,
        )
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


def get_database_directory(receptor: str, db_location: str | None = None) -> str:
    """
    Get the path to the receptor-level germline database directory.

    Parameters
    ----------
    receptor : str
        The receptor type.

    db_location : Optional[str]
        The path to the addon directory. If not provided, the default location (~/.abstar/) will be used.

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
    receptor = receptor.lower()
    if receptor not in {"bcr", "tcr"}:
        raise ValueError("receptor must be 'bcr' or 'tcr'")
    database_dir = os.path.join(database_dir, receptor)
    abutils.io.make_dir(database_dir)
    return database_dir


def check_for_existing_db(
    name: str, receptor: str, location: str | None = None
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
    return os.path.isdir(os.path.join(location, name.lower()))


def publish_database(staging_dir: str, database_dir: str, replacing: bool) -> None:
    """Atomically expose a completed staged database, preserving the old one."""
    backup_dir = None
    if replacing:
        backup_dir = f"{database_dir}.backup-{uuid.uuid4().hex}"
        os.replace(database_dir, backup_dir)
    try:
        os.replace(staging_dir, database_dir)
    except Exception as primary_error:
        if backup_dir is not None:
            try:
                os.replace(backup_dir, database_dir)
            except Exception as rollback_error:
                raise GermlinePublicationError(
                    primary_error, rollback_error, backup_dir
                ) from primary_error
        raise
    if backup_dir is not None:
        try:
            shutil.rmtree(backup_dir)
        except OSError as error:
            warnings.warn(
                f"database published successfully, but backup cleanup failed: "
                f"{error}; retained backup at {backup_dir}",
                RuntimeWarning,
                stacklevel=2,
            )


def inspect_database_destination(database_dir: str) -> tuple[int, int] | None:
    """Return a safe directory identity without following destination links."""
    try:
        metadata = os.lstat(database_dir)
    except FileNotFoundError:
        return None
    if stat.S_ISLNK(metadata.st_mode):
        raise ValueError(
            f"germline database destination must not be a symlink: {database_dir}"
        )
    if not stat.S_ISDIR(metadata.st_mode):
        raise ValueError(
            f"germline database destination must be a directory: {database_dir}"
        )
    return metadata.st_dev, metadata.st_ino


@contextmanager
def database_build_lock(database_root: str, name: str):
    """Coordinate publication by abstar builders targeting the same database."""
    lock_path = os.path.join(database_root, f".{name}.lock")
    flags = os.O_CREAT | os.O_RDWR
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(lock_path, flags, 0o600)
    try:
        metadata = os.fstat(descriptor)
        if not stat.S_ISREG(metadata.st_mode):
            raise ValueError(f"database build lock must be a regular file: {lock_path}")
        fcntl.flock(descriptor, fcntl.LOCK_EX)
        yield lock_path
    finally:
        try:
            fcntl.flock(descriptor, fcntl.LOCK_UN)
        finally:
            os.close(descriptor)


_MMSEQS_COMPONENT_SUFFIXES = (
    "",
    ".dbtype",
    ".index",
    "_h",
    "_h.dbtype",
    "_h.index",
)


def validate_staged_database(
    database_dir: str, receptor: str, *, require_manifest: bool = False
) -> None:
    """Validate a complete custom database before it becomes discoverable."""
    database_dir = os.path.abspath(database_dir)
    staging_metadata = os.lstat(database_dir)
    if stat.S_ISLNK(staging_metadata.st_mode) or not stat.S_ISDIR(staging_metadata.st_mode):
        raise ValueError("staged germline database root must be a real directory")
    missing = []
    segment_components = {}
    for segment in "vdjc":
        components = [
            os.path.join("imgt_gapped", f"{segment}.fasta"),
            os.path.join("ungapped", f"{segment}.fasta"),
            *[
                os.path.join("mmseqs", f"{segment}{suffix}")
                for suffix in _MMSEQS_COMPONENT_SUFFIXES
            ],
        ]
        segment_components[segment] = components
        present = any(os.path.lexists(os.path.join(database_dir, path)) for path in components)
        if segment in "vj" or present:
            missing.extend(
                path for path in components
                if not os.path.lexists(os.path.join(database_dir, path))
            )
    if require_manifest:
        manifest = os.path.join(database_dir, "manifest.txt")
        if not os.path.lexists(manifest):
            missing.append("manifest.txt")
    if missing:
        raise ValueError(
            "incomplete staged germline database: missing " + ", ".join(missing)
        )

    for segment, components in segment_components.items():
        if not os.path.lexists(os.path.join(database_dir, components[0])):
            continue
        for component in components:
            _validate_staged_file(
                database_dir,
                component,
                require_nonempty=True,
            )
    if require_manifest:
        _validate_staged_file(database_dir, "manifest.txt", require_nonempty=True)

    for segment, components in segment_components.items():
        gapped_path = os.path.join(database_dir, components[0])
        if not os.path.lexists(gapped_path):
            continue
        ungapped_path = os.path.join(database_dir, components[1])
        gapped = list(abutils.io.read_fasta(gapped_path))
        ungapped = list(abutils.io.read_fasta(ungapped_path))
        gapped_ids = [str(sequence.id) for sequence in gapped]
        ungapped_ids = [str(sequence.id) for sequence in ungapped]
        if not gapped:
            raise ValueError(
                f"staged germline database {components[0]} contains no sequences"
            )
        if len(gapped_ids) != len(set(gapped_ids)):
            raise ValueError(
                f"staged germline database {components[0]} contains duplicate IDs"
            )
        if len(ungapped_ids) != len(set(ungapped_ids)):
            raise ValueError(
                f"staged germline database {components[1]} contains duplicate IDs"
            )
        if gapped_ids != ungapped_ids:
            raise ValueError(
                f"staged germline database {segment.upper()} gapped and ungapped IDs differ"
            )
        for gapped_sequence, ungapped_sequence in zip(gapped, ungapped):
            if gapped_sequence.sequence.replace(".", "") != ungapped_sequence.sequence:
                raise ValueError(
                    f"staged germline database {segment.upper()} gapped and ungapped "
                    f"sequences differ for {gapped_sequence.id}"
                )


def _validate_staged_file(
    database_dir: str, relative_path: str, *, require_nonempty: bool
) -> None:
    """Require an owned regular file beneath real staged directories."""
    current = database_dir
    parts = relative_path.split(os.sep)
    for part in parts[:-1]:
        current = os.path.join(current, part)
        metadata = os.lstat(current)
        if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISDIR(metadata.st_mode):
            raise ValueError(
                f"staged germline database component parent must be a real directory: "
                f"{relative_path}"
            )
    path = os.path.join(database_dir, relative_path)
    metadata = os.lstat(path)
    if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISREG(metadata.st_mode):
        raise ValueError(
            f"staged germline database component must be a regular non-symlink file: "
            f"{relative_path}"
        )
    if require_nonempty and metadata.st_size == 0:
        raise ValueError(f"staged germline database component is empty: {relative_path}")


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
    destination = os.path.join(raw_dir, os.path.basename(fasta))
    if os.path.exists(destination):
        raise ValueError(f"duplicate raw input filename: {os.path.basename(fasta)}")
    shutil.copy2(fasta, destination)


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


def validate_germlines(
    vdjs: Iterable[Sequence], constants: Iterable[Sequence], receptor: str
) -> None:
    """Validate identifiers, loci, alphabets, uniqueness, and required segments."""
    vdjs = list(vdjs)
    constants = list(constants)
    allowed_loci = {
        "bcr": {"IGH", "IGK", "IGL"},
        "tcr": {"TRA", "TRB", "TRD", "TRG"},
    }[receptor]
    allowed_d_loci = {"bcr": {"IGH"}, "tcr": {"TRB", "TRD"}}[receptor]
    allowed_bases = set("ACGTRYSWKMBDHVN.")
    seen = set()

    for sequence in [*vdjs, *constants]:
        identifier = str(sequence.id)
        bare_id = identifier.split("__", 1)[0]
        if len(bare_id) < 4 or any(character.isspace() for character in identifier):
            raise ValueError(f"invalid germline identifier: {identifier!r}")
        if identifier in seen:
            raise ValueError(f"duplicate germline identifier: {identifier}")
        seen.add(identifier)
        locus = bare_id[:3].upper()
        if locus not in allowed_loci:
            raise ValueError(
                f"germline {identifier} has locus {locus}, incompatible with {receptor}"
            )
        normalized = sequence.sequence.upper()
        if not normalized or set(normalized) - allowed_bases:
            raise ValueError(f"invalid nucleotide sequence for germline {identifier}")
        sequence.sequence = normalized

    segments_by_locus = {}
    for sequence in vdjs:
        bare_id = str(sequence.id).split("__", 1)[0]
        segment = bare_id[3].upper()
        locus = bare_id[:3].upper()
        if segment not in {"V", "D", "J"}:
            raise ValueError(f"invalid V/D/J segment identifier: {sequence.id}")
        if segment == "D" and locus not in allowed_d_loci:
            raise ValueError(f"locus does not contain D genes: {sequence.id}")
        if segment != "V" and "." in sequence.sequence:
            raise ValueError(f"only V genes may contain IMGT gaps: {sequence.id}")
        segments_by_locus.setdefault(locus, set()).add(segment)

    if not segments_by_locus:
        raise ValueError("custom germline database requires V and J genes")
    for locus, segments in segments_by_locus.items():
        if {"V", "J"} - segments:
            raise ValueError(
                f"custom germline database locus {locus} requires V and J genes"
            )

    allowed_constant_segments = {"bcr": set("CADEGM"), "tcr": {"C"}}[receptor]
    for sequence in constants:
        bare_id = str(sequence.id).split("__", 1)[0]
        if bare_id[3].upper() not in allowed_constant_segments:
            raise ValueError(f"invalid constant-region identifier: {sequence.id}")


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
        ungapped_file = os.path.join(ungapped_dir, f"{segment.lower()}.fasta")
        if not os.path.exists(ungapped_file):
            continue
        if verbose:
            if segment == "V":
                print("  V", end="")
            else:
                print(f" | {segment}", end="")
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
    createdb_cmd = [mmseqs_bin, "createdb", input_file, output_file]
    createdb = _run_mmseqs_command(createdb_cmd)
    if debug:
        print(" ".join(createdb_cmd))
        print(createdb.stdout)
        print(createdb.stderr)

    # create MMseqs2 index
    with tempfile.TemporaryDirectory(
        prefix=f".{os.path.basename(output_file)}.index-",
        dir=os.path.dirname(output_file),
    ) as index_tmp:
        createindex_cmd = [
            mmseqs_bin,
            "createindex",
            output_file,
            index_tmp,
            "--search-type",
            "3",
        ]
        createindex = _run_mmseqs_command(createindex_cmd)
        if debug:
            print(" ".join(createindex_cmd))
            print(createindex.stdout)
            print(createindex.stderr)


def _run_mmseqs_command(command: list[str]) -> sp.CompletedProcess:
    """Run MMseqs and retain its diagnostics when a build step fails."""
    try:
        return sp.run(command, check=True, capture_output=True, text=True)
    except sp.CalledProcessError as error:
        raise GermlineBuildExternalToolError(
            "MMseqs command failed",
            command=command,
            returncode=error.returncode,
            stdout=error.stdout or "",
            stderr=error.stderr or "",
        ) from error
    except OSError as error:
        raise GermlineBuildExternalToolError(
            "MMseqs command could not be started",
            command=command,
            stderr=str(error),
        ) from error


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
