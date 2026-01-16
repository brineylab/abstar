# Copyright (c) 2024 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

import multiprocessing as mp
import os
import shutil
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from typing import Iterable

import abutils
import polars as pl
from abutils import Sequence
from natsort import natsorted
from tqdm.auto import tqdm

from ..annotation.annotator import annotate
from ..assigners.mmseqs import MMseqs
from ..preprocess.merging import merge_fastqs

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


def run(
    sequences: str | Sequence | Iterable[Sequence],
    project_path: str | None = None,
    germline_database: str = "human",
    receptor: str = "bcr",
    output_format: str | Iterable[str] = "airr",
    umi_pattern: str | None = None,
    umi_length: int | None = None,
    merge: bool = False,
    merge_kwargs: dict | None = None,
    interleaved_fastq: bool = False,
    chunksize: int = 500,
    mmseqs_chunksize: int = 1e6,
    mmseqs_threads: int | None = None,
    n_processes: int | None = None,
    copy_inputs_to_project: bool = False,
    verbose: bool = False,
    concise_logging: bool = False,
    as_dataframe: bool = False,
    started_from_cli: bool = False,
    debug: bool = False,
) -> Iterable[Sequence] | Sequence | None:
    """
    Annotate antibody or TCR sequences.

    \b
    command line arguments:
      SEQUENCES can be a FASTA/Q file or a directory of FASTA/Q files.
      PROJECT_PATH is the path to a directory in which tmp, log and output files will be deposited.

    \f

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
        Number of sequences to process at a time for annotation.

    mmseqs_chunksize : int = 1e6,
        Number of sequences to process at a time for MMseqs2 searches (VDJC assignment).

    n_processes : Optional[int] = None,
        Number of processes to use for annotation. If ``None``, the number of processes will be set to the number of
        available CPU cores.

    copy_inputs_to_project : bool = False,
        Whether to copy input sequences to the project directory.

    as_dataframe : bool = False,
        Whether to return the output as a polars DataFrame. If ``True``, the output will be returned as a polars DataFrame.
        Note that this option is only available when running interactively (via the API) and when the `project_path` is not provided.

    concise_logging : bool = False,
        Whether to use concise logging. If ``True``, the logging is more concise, more suitable for running
        abstar inside some other pipeline where only basic progress information is needed.

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
        project_path = tempfile.TemporaryDirectory(prefix="abstar", dir="/tmp").name
    log_dir = os.path.join(project_path, "logs")
    abutils.io.make_dir(log_dir)
    temp_dir = os.path.join(project_path, "tmp")
    abutils.io.make_dir(temp_dir)
    for fmt in output_format:
        abutils.io.make_dir(os.path.join(project_path, fmt))

    # setup logging
    global logger
    if started_from_cli:
        logger = _setup_logging(
            log_dir,
            add_stream_handler=verbose,
            single_line_handler=True,
            debug=debug,
        )
    elif verbose or concise_logging:
        verbose = True
        logger = abutils.log.NotebookLogger(verbose=verbose, end="")
    else:
        logger = abutils.log.null_logger()

    if started_from_cli:
        _log_run_parameters(
            project_path=project_path,
            germline_database=germline_database,
            receptor=receptor,
            output_format=output_format,
            umi_pattern=umi_pattern,
            umi_length=umi_length,
            merge=merge,
            merge_kwargs=merge_kwargs,
            interleaved_fastq=interleaved_fastq,
            chunksize=chunksize,
            n_processes=n_processes,
            copy_inputs_to_project=copy_inputs_to_project,
            verbose=verbose,
            debug=debug,
        )

    # process input sequences
    sequence_files = _process_inputs(sequences, temp_dir)
    if copy_inputs_to_project:
        _copy_inputs_to_project(sequence_files, project_path)

    # merge FASTQ files
    if merge or interleaved_fastq:
        merge_dir = os.path.join(project_path, "merged")
        abutils.io.make_dir(merge_dir)
        merge_log_dir = os.path.join(log_dir, "merge_fastqs")
        abutils.io.make_dir(merge_log_dir)
        # log merge info
        if not concise_logging:
            logger.info("\n\n\n")
            logger.info("MERGE FASTQS\n")
            logger.info("============\n")
            logger.info(f"merge directory: {merge_dir}\n")
            # merging
        sequence_files = merge_fastqs(
            sequence_files,
            merge_dir,
            interleaved=interleaved_fastq,
            show_progress=verbose,
            log_directory=merge_log_dir,
            **merge_kwargs,
        )

    # print sequence file info
    if started_from_cli:
        _log_sequence_file_info(sequence_files)

    # annotation config
    if n_processes is None:
        n_processes = mp.cpu_count()
    annot_kwargs = {
        "output_directory": temp_dir,
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
        logger=logger,
        concise_logging=concise_logging,
        chunksize=mmseqs_chunksize,
        threads=mmseqs_threads,
        debug=debug,
    )

    # annotate sequences
    for sequence_file in natsorted(sequence_files):
        start_time = datetime.now()
        to_delete = []

        # parse sample name
        sample_name = ".".join(
            os.path.basename(sequence_file).rstrip(".gz").split(".")[:-1]
        )
        if not sample_name:  # if the input was Sequence object(s), not file(s)
            sample_name = "sequences"
        # log sample info
        if started_from_cli:
            logger.info("\n\n")
            logger.info("-" * (len(sample_name) + 4))
            logger.info("\n")
            logger.info(f"  {sample_name}\n")
            logger.info("-" * (len(sample_name) + 4))
            logger.info("\n")

        # assign VDJC genes, the returned assign_file is in parquet format
        assign_file, raw_sequence_count = assigner(sequence_file)
        assigner.cleanup()

        # split into annotation jobs
        split_assign_files = abutils.io.split_parquet(
            assign_file, temp_dir, num_rows=chunksize
        )

        # run annotation jobs
        annotated_files = []
        failed_log_files = []
        succeeded_log_files = []
        # log annotation info
        if concise_logging:
            logger.info("\nsequence annotation: ")
        else:
            logger.info("\n")
            logger.info("sequence annotation:\n")
        if verbose and started_from_cli:
            progress_bar = tqdm(
                total=len(split_assign_files),
                bar_format="{desc:<2.5}{percentage:3.0f}%|{bar:25}{r_bar}",
            )
        elif verbose:
            progress_bar = tqdm(
                total=len(split_assign_files),
            )
        with ProcessPoolExecutor(
            max_workers=n_processes,
            mp_context=mp.get_context("spawn"),
        ) as executor:
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
            if as_dataframe:
                sequence_df = pl.read_parquet(annotated_files)
                sequence_count = sequence_df.height
            else:
                annotated_sequences = abutils.io.read_parquet(annotated_files)
                sequences_to_return.extend(annotated_sequences)
                sequence_count = len(annotated_sequences)

            # log results summary
            duration = datetime.now() - start_time
            _log_results_summary(
                sequence_count=sequence_count,
                sequences_per_second=raw_sequence_count / duration.total_seconds(),
                seconds=duration.total_seconds(),
                concise_logging=concise_logging,
            )

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

            # log results summary
            sequence_count = output_df.select(pl.count()).collect().row(0)[0]
            duration = datetime.now() - start_time
            _log_results_summary(
                sequence_count=sequence_count,
                sequences_per_second=raw_sequence_count / duration.total_seconds(),
                seconds=duration.total_seconds(),
                concise_logging=concise_logging,
            )

        # collect files for removal
        to_delete.append(assign_file)
        to_delete.extend(split_assign_files)
        to_delete.extend(annotated_files)
        to_delete.extend(succeeded_log_files)
        to_delete.extend(failed_log_files)

        # remove tmp files
        if not debug:
            _delete_files(to_delete)

    if return_sequences:
        if as_dataframe:
            return sequence_df
        if len(sequences_to_return) == 1:
            return sequences_to_return[0]
        return sequences_to_return
    else:
        return


def _process_inputs(
    sequences: str | Sequence | Iterable[Sequence],
    temp_dir: str,
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
        # input is a string -- either a file/directory path or a raw sequence string
        if os.path.isfile(sequences):
            sequence_files = [os.path.abspath(sequences)]
        elif os.path.isdir(sequences):
            sequence_files = abutils.io.list_files(
                sequences,
                recursive=True,
                extension=["fasta", "fa", "fastq", "fq", "fasta.gz", "fastq.gz"],
            )
        else:
            sequences = Sequence(sequences)
    if isinstance(sequences, Sequence) or (
        isinstance(sequences, Iterable)
        and all(isinstance(s, Sequence) for s in sequences)
    ):
        # input is a Sequence or an iterable of Sequences
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
    return natsorted(sequence_files)


def _copy_inputs_to_project(sequence_files: Iterable[str], project_path: str) -> None:
    """
    Copy input sequences to the project directory.

    Parameters
    ----------
    sequence_files : Iterable[str]
        A list of sequence files to copy.

    project_path : str
        The path to the project directory.
    """
    inputs_path = os.path.join(project_path, "input")
    abutils.io.make_dir(inputs_path)
    for f in sequence_files:
        shutil.copy(f, inputs_path)


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


# ===============================
#
#       PRINTING/LOGGING
#
# ===============================


def _setup_logging(
    log_dir: str,
    add_stream_handler: bool,
    single_line_handler: bool = False,
    debug: bool = False,
) -> None:
    logfile = os.path.join(log_dir, "abstar.log")
    abutils.log.setup_logging(
        logfile=logfile,
        add_stream_handler=add_stream_handler,
        single_line_handler=single_line_handler,
        debug=debug,
    )
    # global logger
    logger = abutils.log.get_logger(
        name="abstar",
        add_stream_handler=add_stream_handler,
        single_line_handler=single_line_handler,
    )
    return logger


def _log_run_parameters(
    project_path: str,
    germline_database: str,
    receptor: str,
    output_format: str | Iterable[str],
    umi_pattern: str | None,
    umi_length: int | None,
    merge: bool,
    merge_kwargs: dict | None,
    interleaved_fastq: bool,
    chunksize: int,
    n_processes: int | None,
    copy_inputs_to_project: bool,
    verbose: bool,
    debug: bool,
) -> None:
    logger.info("\n")
    # printing the splash line-by-line makes the log files look nicer
    for line in ABSTAR_SPLASH.split("\n"):
        logger.info(line + "\n")
    logger.info("\n")
    logger.info("RUN PARAMETERS\n")
    logger.info("===============\n")
    logger.info(f"PROJECT PATH: {project_path}\n")
    logger.info(f"GERMLINE DATABASE: {germline_database}\n")
    logger.info(f"RECEPTOR: {receptor}\n")
    logger.info(f"OUTPUT FORMAT: {', '.join(output_format)}\n")
    logger.info(f"UMI PATTERN: {umi_pattern}\n")
    logger.info(f"UMI LENGTH: {umi_length}\n")
    logger.info(f"MERGE: {merge}\n")
    logger.info(f"MERGE KWARGS: {merge_kwargs}\n")
    logger.info(f"INTERLEAVED FASTQ: {interleaved_fastq}\n")
    logger.info(f"CHUNKSIZE: {chunksize}\n")
    logger.info(
        f"NUM PROCESSES: {n_processes if n_processes is not None else 'auto'}\n"
    )
    logger.info(f"COPY INPUTS TO PROJECT: {copy_inputs_to_project}\n")
    logger.info(f"VERBOSE: {verbose}\n")
    logger.info(f"DEBUG: {debug}\n")


def _log_sequence_file_info(sequence_files: Iterable[str]) -> None:
    num_files = len(sequence_files)
    plural = "files" if num_files > 1 else "file"
    logger.info("\n\n\n")
    logger.info("INPUT FILES\n")
    logger.info("===========\n")
    logger.info(f"found {num_files} input {plural}:\n")
    if num_files < 6:
        for f in sequence_files:
            logger.info(f"  {os.path.basename(f)}\n")
    else:
        for f in sequence_files[:5]:
            logger.info(f"  {os.path.basename(f)}\n")
        logger.info(f"  ... and {num_files - 5} more\n")


def _log_results_summary(
    sequence_count: int,
    sequences_per_second: int,
    seconds: float,
    concise_logging: bool = False,
) -> None:
    if seconds < 60:
        hours = 0
        minutes = 0
        seconds = seconds
    elif seconds < 3600:
        hours = 0
        minutes = int(seconds / 60)
        seconds = int(seconds % 60)
    else:
        hours = int(seconds / 3600)
        minutes = (seconds % 3600) / 60
        seconds = seconds % 60
    duration_string = f"{hours:02}:{minutes:02}:{seconds:02.2f}"
    if concise_logging:
        logger.info(f"annotated sequences: {sequence_count:,}\n")
    else:
        logger.info("\n")
        logger.info(f"{sequence_count:,} sequences had an identifiable rearrangement\n")
    logger.info(
        f"time elapsed: {duration_string} ({sequences_per_second:,.2f} sequences/sec)\n"
    )


ABSTAR_SPLASH = """
         __         __ 
  ____ _/ /_  _____/ /_____ ______
 / __ `/ __ \/ ___/ __/ __ `/ ___/
/ /_/ / /_/ (__  ) /_/ /_/ / /    
\__,_/_.___/____/\__/\__,_/_/      
"""
