# Copyright (c) 2024 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

import multiprocessing as mp
import os
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Iterable, Optional, Union

import abutils
import click
import polars as pl
from abutils import Sequence
from natsort import natsorted
from tqdm.auto import tqdm

from ..annotation.annotator import annotate
from ..assigners.mmseqs import MMseqs
from ..preprocess.merging import merge_fastqs
from ..utils.callbacks import parse_dict_from_string

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


@click.command()
@click.argument(
    "sequences",
    type=str,
    help="Path to a FASTA/Q file or a directory of FASTA/Q files. Gzip-compressed files are supported.",
)
@click.argument(
    "project_path",
    type=str,
    help="Path to a directory in which tmp, log and output files will be deposited",
)
@click.option(
    "--germline_database",
    type=str,
    show_default=True,
    default="human",
    help="Name of the germline database to be used for assignment/annotation",
)
@click.option(
    "--receptor",
    type=click.Choice(["bcr", "tcr"], case_sensitive=False),
    show_default=True,
    default="bcr",
    help="Name of the receptor to be used for assignment/annotation",
)
@click.option(
    "-O" "--output_format",
    type=click.Choice(["airr", "parquet"], case_sensitive=False),
    multiple=True,
    show_default=True,
    default=["airr"],
    help="Format of the output files",
)
@click.option(
    "--umi_pattern",
    type=str,
    default=None,
    help="Pattern to match for extracting the UMI sequence, or name of a built-in pattern",
)
@click.option(
    "--umi_length",
    type=int,
    default=None,
    help="Length of the UMI sequence to extract. If positive, the UMI will be parsed from the start of the sequence. If negative, the UMI will be parsed from the end of the sequence.",
)
@click.option(
    "--merge",
    is_flag=True,
    default=False,
    help="Whether to merge FASTQ files prior to assignment",
)
@click.option(
    "--merge_kwargs",
    type=str,
    callback=parse_dict_from_string,
    default=None,
    help="Keyword arguments to pass to the merge_fastqs function. Format must be 'key1=val1,key2=val2'",
)
@click.option(
    "--interleaved_fastq",
    is_flag=True,
    default=False,
    help="Whether the input FASTQ files are interleaved",
)
@click.option(
    "--chunksize",
    type=int,
    show_default=True,
    default=500,
    help="Number of sequences to process at a time",
)
@click.option(
    "--n_processes",
    type=int,
    show_default=True,
    default=None,
    help="Number of processes to use for annotation",
)
@click.option(
    "--verbose/--quiet",
    default=True,
    help="Whether to print verbose output",
)
@click.option(
    "--debug",
    is_flag=True,
    default=False,
    help="Whether to run in debug mode, which results in temporary files being retained and additional logging.",
)
def run(
    sequences: Union[str, Sequence, Iterable[Sequence]],
    project_path: Optional[str] = None,
    germline_database: str = "human",
    receptor: str = "bcr",
    output_format: Union[str, Iterable[str]] = "airr",
    umi_pattern: Optional[str] = None,
    umi_length: Optional[int] = None,
    merge: bool = False,
    merge_kwargs: Optional[dict] = None,
    interleaved_fastq: bool = False,
    chunksize: int = 500,
    n_processes: Optional[int] = None,
    verbose: bool = False,
    debug: bool = False,
):
    """
    Annotate antibody or TCR sequences.

    Parameters
    ----------

    sequences : Union[str, Sequence, Iterable[Sequence]]
        The sequences to annotate. Can be one of the following:

          - ``str``: path to a FASTA/Q file, path to a directory of FASTA/Q files, or a single sequence, as a string
          - ``Sequence``: a single ``abutils.Sequence`` object
          - ``Iterable[Sequence]``: an iterable of ``abutils.Sequence`` objects

        .. note::
            If `sequences` is a directory path, files in the directory will be consumed recursively, including
            files in any subfolders.

    project_path : Optional[str] = None
        If provided, the path to a directory in which tmp, log and output files will be deposited. If not provided,
        annotated sequecnes will be returned as ``Sequence`` objects and log/tmp directories will be placed in
        ``"/tmp"``.

        .. warning::
            If ``debug=False`` (the default), the temporary directory will be removed during cleanup. If a
            directory named ``"tmp"`` exists in the provided `project_path`, it will be deleted.

    germline_database : str = "human",
        Name of the germline database to be used for assignment/annotation. Built-in options are
        "human", "mouse", and "macaque" and "humouse".

    receptor : str = "bcr",
        Name of the receptor to be used for assignment/annotation. Options are "bcr" and "tcr".

    output_format : Union[str, Iterable[str]] = "airr",
        Format of the output files. Options are "airr" and "parquet". If more than one output format is
        desired, a list of multiple output formats can be provided.

    umi_pattern : Optional[str], default=None
        Pattern to match for parsing the UMI sequence, or name of a built-in pattern.

    umi_length : Optional[int], default=None
        Length of the UMI sequence.

    merge : bool = False,
        Whether to merge FASTQ files prior to assignment. If ``True``, the ``merge_fastqs`` function will be used to merge
        FASTQ files.

    merge_kwargs : Optional[dict] = None,
        Keyword arguments to pass to the ``merge_fastqs`` function.

    interleaved_fastq : bool = False,
        Whether the input FASTQ files are interleaved. This is common when using data downloaded from the
        SRA_, which reuqires that paired FASTQ files be submitted as a single interleaved file rather than two
        separate files.

    chunksize : int = 500,
        Number of sequences to process at a time.

    n_processes : Optional[int] = None,
        Number of processes to use for annotation. If ``None``, the number of processes will be set to the number of
        available CPU cores.

    verbose : bool = False,
        Whether to print verbose output.

    debug : bool = False,
        If ``True``, the following additional things will happen:
          - successfully annotated sequences will be logged in addition to sequences that errored during annotation
          - all tmp files will be retained
          - stdout and stderr from various third party tools (MMseqs, fastp, etc) will be captured and logged


    .. _SRA: https://www.ncbi.nlm.nih.gov/sra


    """
    # output format
    if isinstance(output_format, str):
        output_format = [output_format]

    # set up log/output/temp directories
    if project_path is not None:
        return_sequences = False
        project_path = os.path.abspath(project_path)
    else:
        return_sequences = True
        sequences_to_return = []
        output_format = ["parquet"]
        project_path = tempfile.TemporaryDirectory(prefix="abstar", dir="/tmp")
    log_dir = os.path.join(project_path, "logs")
    abutils.io.make_dir(log_dir)
    temp_dir = os.path.join(project_path, "tmp")
    abutils.io.make_dir(temp_dir)
    for fmt in output_format:
        abutils.io.make_dir(os.path.join(project_path, fmt))

    # process input sequences
    sequence_files = _process_inputs(sequences, temp_dir)

    # merge FASTQ files
    if merge:
        merge_dir = os.path.join(project_path, "merged")
        sequence_files = merge_fastqs(
            sequence_files, merge_dir, interleaved=interleaved_fastq, **merge_kwargs
        )

    # annotation config
    if n_processes is None:
        n_processes = mp.cpu_count()
    annot_kwargs = {
        "output-directory": temp_dir,
        "germline_database": germline_database,
        "log_directory": log_dir,
        "umi_pattern": umi_pattern,
        "umi_length": umi_length,
        "debug": debug,
    }

    # initialize the assigner
    assigner = MMseqs(
        output_directory=temp_dir,
        log_directory=log_dir,
        germdb_name=germline_database,
        receptor=receptor,
        verbose=verbose,
        debug=debug,
    )

    # annotate sequences
    for sequence_file in natsorted(sequence_files):
        to_delete = []

        # parse sample name
        sample_name = ".".join(
            os.path.basename(sequence_file).rstrip(".gz").split(".")[:-1]
        )
        if verbose:
            print(f"  {sample_name}")
            print("-" * (len(sample_name) + 4))

        # assign VDJC genes
        assign_file = assigner(sequence_file)  # returns a parquet file
        assigner.cleanup()

        # split into annotation jobs
        split_assign_files = abutils.io.split_parquet(
            assign_file, temp_dir, num_rows=chunksize
        )

        # run annotation jobs
        annotated_files = []
        failed_log_files = []
        succeeded_log_files = []
        if verbose:
            progress_bar = tqdm(
                total=len(split_assign_files),
                desc="  - annotating",
            )
        with ProcessPoolExecutor(max_workers=n_processes) as executor:
            futures = [
                executor.submit(annotate, f, **annot_kwargs) for f in split_assign_files
            ]
            for future in as_completed(futures):
                annotated, failed, succeeded = future.result()
                annotated_files.append(annotated)
                failed_log_files.append(failed)
                succeeded_log_files.append(succeeded)
                if verbose:
                    progress_bar.update(1)
        if verbose:
            progress_bar.close()

        # get output Sequences
        if return_sequences:
            sequences_to_return.extend(abutils.io.read_parquet(annotated_files))

        # or assemble output files (including logs)
        else:
            output_df = pl.scan_parquet(annotated_files)
            if "airr" in output_format:
                airr_file = os.path.join(project_path, f"airr/{sample_name}.tsv")
                output_df.sink_csv(airr_file, separator="\t")
            if "parquet" in output_format:
                parquet_file = os.path.join(
                    project_path, f"parquet/{sample_name}.parquet"
                )
                output_df.sink_parquet(parquet_file)
            # assemble logs
            failed_log_file = os.path.join(log_dir, f"{sample_name}.failed")
            _assemble_logs(failed_log_files, failed_log_file)
            if debug:
                # only log succeeded sequences if we're in debug mode
                succeeded_log_file = os.path.join(log_dir, f"{sample_name}.succeeded")
                _assemble_logs(succeeded_log_files, succeeded_log_file)

        # collect files for removal
        # NOTE: failed_log_files aren't here, because they're retained regardless of debug status
        to_delete.append(assign_file)
        to_delete.extend(split_assign_files)
        to_delete.extend(annotated_files)
        to_delete.extend(succeeded_log_files)

        # remove tmp files
        if not debug:
            _delete_files(to_delete)

    if return_sequences:
        if len(sequences_to_return) == 1:
            return sequences_to_return[0]
        return sequences_to_return
    else:
        return


def _process_inputs(
    sequences: Union[str, Sequence, Iterable[Sequence]], temp_dir: str
) -> Iterable[str]:
    """
    Process the various inputs accepted by abstar and return a list of one or more sequence files.

    Parameters
    ----------
    sequences : Union[str, Sequence, Iterable[Sequence]]
        The sequences to process.

    temp_dir : str
        The path to a directory in which tmp files will be deposited.

    Returns
    -------
    sequence_files : Iterable[str]
        A list of one or more sequence files.
    """
    sequence_files = None
    if isinstance(sequences, str):
        if os.path.isfile(sequences):
            sequence_files = [os.path.abspath(sequences)]
        elif os.path.isdir(sequences):
            sequence_files = abutils.io.list_files(sequences, recursive=True)
        else:
            sequences = Sequence(sequences)
    if isinstance(sequences, (Sequence, Iterable[Sequence])):
        if isinstance(sequences, Sequence):
            sequences = [sequences]
        temp_file = tempfile.NamedTemporaryFile(delete=False, dir=temp_dir, mode="w")
        fastas = [seq.fasta for seq in sequences]
        temp_file.write("\n".join(fastas))
        temp_file.close()
        sequence_files = [temp_file.name]
    if sequence_files is None:
        raise ValueError(
            "Invalid input sequences. Must be a path to a file or directory, a single sequence, or an iterable of sequences."
        )
    return sequence_files


def _assemble_logs(log_files: Iterable[str], combined_log_file: str) -> None:
    """
    Assemble log files into a single file.

    Parameters
    ----------
    log_files : Iterable[str]
        A list of log files to be assembled.

    combined_log_file : str
        The path to the file in which the combined log will be written.

    Returns
    -------
    None
    """
    with open(combined_log_file, "w") as f:
        for log_file in log_files:
            if log_file is not None:
                with open(log_file, "r") as log:
                    f.write(log.read())


def _delete_files(files: Iterable[str]) -> None:
    """
    Delete a list of files.
    """
    for f in files:
        if f is not None:
            if os.path.exists(f):
                os.remove(f)
