# Copyright (c) 2024 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

from typing import Iterable, Optional, Union

import click
from abutils import Sequence

from ..core.abstar import run as _run
from ..core.germline import build_germline_database as _build_germline_database
from ..utils.callbacks import HiddenClickOption, parse_dict_from_string


@click.group()
def cli():
    pass


# cli.add_command(build_germline_database)


@cli.command()
@click.argument(
    "input_path",
    type=str,
    # help="Path to a FASTA/Q file or a directory of FASTA/Q files. Gzip-compressed files are supported.",
)
@click.argument(
    "project_path",
    type=str,
    # help="Path to a directory in which tmp, log and output files will be deposited",
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
    "-o",
    "--output-format",
    type=click.Choice(["airr", "parquet"], case_sensitive=False),
    multiple=True,
    show_default=True,
    default=["airr"],
    help="Format of the output files",
)
@click.option(
    "--umi-pattern",
    type=str,
    default=None,
    help="Pattern to match for extracting the UMI sequence, or name of a built-in pattern",
)
@click.option(
    "--umi-length",
    type=int,
    default=None,
    help="Length of the UMI sequence to extract. If positive, the UMI will be parsed from the start of the sequence. If negative, the UMI will be parsed from the end of the sequence.",
)
@click.option(
    "-m",
    "--merge",
    is_flag=True,
    default=False,
    help="Whether to merge FASTQ files prior to assignment",
)
@click.option(
    "--merge-kwargs",
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
    "-c",
    "--chunksize",
    type=int,
    show_default=True,
    default=500,
    help="Number of sequences to process at a time",
)
@click.option(
    "--n-processes",
    type=int,
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
@click.option(
    "--copy-inputs-to-project",
    cls=HiddenClickOption,
    is_flag=True,
    default=True,
)
@click.option(
    "--started_from_cli",
    cls=HiddenClickOption,
    is_flag=True,
    default=True,
)
def run(
    input_path: Union[str, Sequence, Iterable[Sequence]],
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
    copy_inputs_to_project: bool = False,
    verbose: bool = False,
    started_from_cli: bool = False,
    debug: bool = False,
) -> Optional[Union[Iterable[Sequence], Sequence]]:
    """
    Annotate antibody or TCR sequences.

    \b
    command line arguments:
      INPUT_PATH can be a FASTA/Q file or a directory of FASTA/Q files.
      PROJECT_PATH is the path to a directory in which tmp, log and output files will be deposited.
    """
    return _run(
        sequences=input_path,
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
        started_from_cli=started_from_cli,
        debug=debug,
    )


@cli.command()
@click.argument(
    "name",
    type=str,
    required=True,
)
@click.option(
    "-f",
    "--fasta",
    "fasta_files",
    type=str,
    multiple=True,
    default=None,
    help="Path to a FASTA-formatted file containing gapped germline gene sequences. Can be used multiple times to provide multiple files.",
)
@click.option(
    "-j",
    "--json",
    "json_files",
    type=str,
    multiple=True,
    default=None,
    help="Path to a JSON-formatted file containing germline gene sequences. Can be used multiple times to provide multiple files.",
)
@click.option(
    "-c",
    "--constant",
    "constant_files",
    type=str,
    multiple=True,
    default=None,
    help="Path to a FASTA-formatted file containing gapped constant region sequences. Can be used multiple times to provide multiple files.",
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
    "--include-species-in-name/--exclude-species-from-name",
    default=True,
    help="Whether to include the species in the name of each sequence, like so: 'IGHV1-2*02__homo_sapiens'. Only applies if input data is JSON-formatted or FASTA-formatted with species names already included in the sequence names.",
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
    name: str,
    fasta_files: Optional[Union[str, Iterable[str]]] = None,
    json_files: Optional[Union[str, Iterable[str]]] = None,
    constant_files: Optional[Union[str, Iterable[str]]] = None,
    receptor: str = "bcr",
    manifest: Optional[str] = None,
    include_species_in_name: bool = True,
    location: Optional[str] = None,
    verbose: bool = True,
    debug: bool = False,
):
    return _build_germline_database(
        name=name,
        fastas=fasta_files,
        jsons=json_files,
        constants=constant_files,
        receptor=receptor,
        manifest=manifest,
        include_species_in_name=include_species_in_name,
        location=location,
        verbose=verbose,
        debug=debug,
    )
