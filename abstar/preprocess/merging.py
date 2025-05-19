# Copyright (c) 2024 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT


import os
import subprocess as sp
import tempfile
from typing import Iterable, Optional, Union

from abutils.bin import get_path as get_binary_path
from abutils.io import (
    concatenate_files,
    delete_files,
    list_files,
    make_dir,
    parse_fastx,
    rename_file,
)
from natsort import natsorted
from tqdm.auto import tqdm

__all__ = ["merge_fastqs", "group_paired_fastqs"]


class FASTQFile:
    """
    Base class for processing paired FASTQ filenames.

    Parameters
    ----------
    file : str
        Path to a FASTQ file.

    """

    def __init__(self, file: str):
        self.file = file
        self.path = os.path.abspath(file)
        self.basename = os.path.basename(self.path)
        self.dir = os.path.dirname(self.path)
        self.filename = self.basename.rstrip(".gz").rstrip(".fastq").rstrip(".fq")

    def __eq__(self, other):
        return self.name == other.name


class IlluminaFile(FASTQFile):
    """
    Class for processing Illumina-formatted FASTQ filenames.
    Also compatible with Element's "legacy" filename structure.

    Parameters
    ----------
    file : str
        Path to a FASTQ file.

    """

    def __init__(self, file: str):
        super().__init__(file)
        self.schema = "illumina"

    @property
    def name(self):
        return "_".join(self.filename.split("_")[:-4])

    @property
    def number(self):
        return self.filename.split(".")[0].split("_")[-1]

    @property
    def read(self):
        return self.filename.split("_")[-2]

    @property
    def lane(self):
        return self.filename.split("_")[-3]

    @property
    def sample(self):
        return self.filename.split("_")[-4]


class ElementFile(FASTQFile):
    """
    Class for processing Element-formatted FASTQ filenames.

    Parameters
    ----------
    file : str
        Path to a FASTQ file.

    """

    def __init__(self, file: str):
        super().__init__(file)
        self.schema = "element"

    @property
    def name(self):
        return "_".join(self.filename.split("_")[:-1])

    @property
    def number(self):
        return ""

    @property
    def read(self):
        return self.filename.split("_")[-1]

    @property
    def lane(self):
        return ""

    @property
    def sample(self):
        return ""


class MergeGroup:
    """
    Group of FASTQ files to be merged.

    Parameters
    ----------
    name : str
        Name of the sample.

    files : Iterable[FASTQFile]
        Iterable of FASTQFile objects.

    """

    def __init__(self, name: str, files: Iterable[FASTQFile]):
        self.name = name
        self.files = files
        self.merged_file = None  # assigned by merge function

    def merge(
        self,
        merged_directory: str,
        log_directory: Optional[str] = None,
        format: str = "fastq",
        algo: str = "fastp",
        binary_path: Optional[str] = None,
        minimum_overlap: int = 30,
        allowed_mismatches: int = 5,
        allowed_mismatch_percent: float = 20.0,
        trim_adapters: bool = True,
        adapter_file: Optional[str] = None,
        quality_trim: bool = True,
        window_size: int = 4,
        quality_cutoff: int = 20,
        merge_args: Optional[str] = None,
        compress_output: bool = False,
        # show_progress: bool = False,
        debug: bool = False,
    ) -> str:
        groups = self._group_by_lane()
        n_groups = len(groups)
        # if show_progress:
        #     print(f"\n{self.name}")
        #     # print("-" * len(self.name))
        #     if n_groups > 1:
        #         groups = tqdm(
        #             groups,
        #             desc="  - merging lanes",
        #             total=n_groups,
        #             unit="lane",
        #             leave=False,
        #         )
        #     else:
        #         print("  - merging paired FASTQs")

        # merge function
        if algo.lower() == "fastp":
            merge_func = merge_fastqs_fastp
        else:
            raise ValueError(f"Invalid merge algorithm: {algo}. Must be 'fastp'.")

        # merge files
        merged_files = []
        compress_suffix = ".gz" if compress_output else ""
        self.merged_file = os.path.join(
            merged_directory, f"{self.name}.{format.lower()}{compress_suffix}"
        )
        # if show_progress:
        #     groups = tqdm(groups, total=n_groups, leave=True)
        for group in groups:
            r1, r2 = natsorted(group, key=lambda x: x.read)
            # add name and compression suffix, if necessary
            name = self.name if n_groups == 1 else f"{self.name}_{r1.lane}"
            merged_path = os.path.join(
                merged_directory, f"{name}.{format.lower()}{compress_suffix}"
            )
            # merge
            merge_func(
                r1.path,
                r2.path,
                merged_path,
                name=name,
                log_directory=log_directory,
                trim_adapters=trim_adapters,
                quality_trim=quality_trim,
                window_size=window_size,
                quality_cutoff=quality_cutoff,
                minimum_overlap=minimum_overlap,
                allowed_mismatches=allowed_mismatches,
                allowed_mismatch_percent=allowed_mismatch_percent,
                adapter_file=adapter_file,
                binary_path=binary_path,
                additional_args=merge_args,
                debug=debug,
            )
            merged_files.append(merged_path)

        # concatenate merged files if dataset has multiple lanes
        if len(merged_files) > 1:
            # if verbose:
            #     print("  - concatenating merged files")
            concatenate_files(merged_files, self.merged_file)
            delete_files(merged_files)
        else:
            rename_file(merged_files[0], self.merged_file)
        return self.merged_file

    def _group_by_lane(self):
        lane_dict = {}
        for f in self.files:
            if f.lane not in lane_dict:
                lane_dict[f.lane] = []
            lane_dict[f.lane].append(f)
        return [
            group for lane, group in natsorted(lane_dict.items(), key=lambda x: x[0])
        ]


def merge_fastqs(
    files: Union[str, Iterable],
    output_directory: str,
    output_format: str = "fastq",
    log_directory: Optional[str] = None,
    schema: str = "illumina",
    algo: str = "fastp",
    binary_path: Optional[str] = None,
    merge_args: Optional[str] = None,
    minimum_overlap: int = 30,
    allowed_mismatches: int = 5,
    allowed_mismatch_percent: float = 20.0,
    trim_adapters: bool = True,
    adapter_file: Optional[str] = None,
    quality_trim: bool = True,
    window_size: int = 4,
    quality_cutoff: int = 20,
    interleaved: bool = False,
    compress_output: bool = False,
    debug: bool = False,
    show_progress: bool = False,
) -> Iterable[str]:
    """
    Merge paired-end fastq files.

    Parameters
    ----------
    files : Union[str, Iterable]
        Path to a directory containing paired-end fastq files, or an iterable of paired-end fastq files.

    output_directory : str
        Path to the directory in which to save the merged files.

    output_format : str, optional
        Output format. Must be 'fastq' or 'fasta'. Default is 'fastq'.

    log_directory : str, optional
        Path to the directory in which to save the fastp reports.
        If not provided, fastp reports will not be saved.

    schema : str, optional
        Schema of the file names. Can be 'illumina' or 'element'. Default is 'illumina'.

    algo : str, optional
        Algorithm to use for merging. Must be 'fastp'. Default is 'fastp'.

    binary_path : str, optional
        Path to the merge algorithm binary. If not provided, the binary packaged with abutils will be used.

    merge_args : str, optional
        Additional arguments (as a string) to pass to the merge function.

    minimum_overlap : int, optional
        Minimum overlap between reads. Default is 30.

    allowed_mismatches : int, optional
        Allowed mismatches between reads. Default is 5.

    allowed_mismatch_percent : float, optional
        Allowed mismatch percentage between reads. Default is 20.0.

    trim_adapters : bool, optional
        If True, trim adapters. Default is True.

    adapter_file : str, optional
        Path to a FASTA file containing adapter sequences. If not provided, the default
        adapter sequences (Illumina TruSeq) will be used.

    quality_trim : bool, optional
        If True, trim low-quality bases. Default is True.

    window_size : int, optional
        Sliding window size for quality trimming. Default is 4.

    quality_cutoff : int, optional
        Mean quality cutoff for quality trimming. Default is 20.

    compress_output : bool, optional
        If True, compress the output file using gzip. Default is False.

    debug : bool, optional
        If True, print debug output. Default is False.

    show_progress : bool, optional
        If True, show a progress bar. Default is False.

    Returns
    -------
    list
        A list of paths to the merged FASTA/Q files.

    """
    # input/output files and directories
    if isinstance(files, str):
        if os.path.isdir(files):
            files = list_files(
                files, extension=[".fastq", ".fq", ".fastq.gz", ".fq.gz"]
            )
        elif os.path.isfile(files):
            files = [files]
        else:
            raise ValueError(
                f"If files are supplied as a string, it must be either a directory or a file. The supplied ({files}) is neither."
            )
    make_dir(output_directory)

    merged_files = []

    if interleaved:
        if show_progress:
            _print_files_to_merge(
                [os.path.basename(f) for f in files], name="interleaved FASTQ file"
            )
            files = tqdm(
                files,
                desc="merging FASTQ files:",
                leave=True,
                bar_format="{desc}{percentage:3.0f}%|{bar:25}{r_bar}",
            )
        for f in files:
            compress_suffix = ".gz" if compress_output else ""
            name = os.path.basename(f).rstrip(".gz").rstrip(".fastq")
            merged_file = os.path.join(
                output_directory, f"{name}.{output_format.lower()}{compress_suffix}"
            )
            # fastp seems to occasionally have problems with multi-line interleaved FASTQ
            # files, so we'll make a temp file with sequences/qualities each on a single line
            tmp = tempfile.NamedTemporaryFile(delete=False, dir=os.path.dirname(f))
            for seq in parse_fastx(f):
                tmp.write(f"{seq.fastq}\n")
            tmp.close()
            merge_fastqs_fastp(
                forward=tmp.name,
                merged=merged_file,
                binary_path=binary_path,
                log_directory=log_directory,
                minimum_overlap=minimum_overlap,
                allowed_mismatches=allowed_mismatches,
                allowed_mismatch_percent=allowed_mismatch_percent,
                trim_adapters=trim_adapters,
                adapter_file=adapter_file,
                quality_trim=quality_trim,
                window_size=window_size,
                quality_cutoff=quality_cutoff,
                additional_args=merge_args,
                interleaved=True,
                debug=debug,
            )
            merged_files.append(merged_file)
            os.unlink(tmp.name)
    else:
        # group files by sample
        file_pairs = group_paired_fastqs(files, schema=schema)
        if show_progress:
            _print_files_to_merge(
                [f.name for f in file_pairs], name="paired-end FASTQ group"
            )
            file_pairs = tqdm(
                file_pairs,
                desc="merging FASTQ files:",
                leave=True,
                bar_format="{desc}{percentage:3.0f}%|{bar:25}{r_bar}",
            )
        # merge files
        for fp in file_pairs:
            merged_file = fp.merge(
                merged_directory=output_directory,
                format=output_format,
                algo=algo,
                binary_path=binary_path,
                log_directory=log_directory,
                minimum_overlap=minimum_overlap,
                allowed_mismatches=allowed_mismatches,
                allowed_mismatch_percent=allowed_mismatch_percent,
                trim_adapters=trim_adapters,
                adapter_file=adapter_file,
                quality_trim=quality_trim,
                window_size=window_size,
                quality_cutoff=quality_cutoff,
                merge_args=merge_args,
                compress_output=compress_output,
                debug=debug,
                # show_progress=show_progress,
            )
            merged_files.append(merged_file)
    return merged_files


def group_paired_fastqs(
    files: Union[str, Iterable], schema: str = "illumina"
) -> Iterable[MergeGroup]:
    """
    Group paired-end fastq files by sample name. If a sample has multiple lanes, the files will be combined.

    Parameters
    ----------
    files : Union[str, Iterable]
        Path to a directory containing paired-end fastq files, or an iterable of paired-end fastq files.

    schema : str, optional
        Schema of the file names. Can be 'illumina' or 'element'. Default is 'illumina'.

    Returns
    -------
    list
        A list of MergeGroup objects, one for each sample. Note that if a sample has multiple lanes,
        its MergeGroup will contain the paired reads for all lanes. Moreover, calling
        `MergeGroup.merge()` will merge the paired reads for all lanes and concatenate the results into
        a single file.

    """
    # process input files
    if isinstance(files, str):
        if not os.path.isdir(files):
            err = f"The supplied file path ({files}) does not exist or is not a directory."
            raise ValueError(err)
        files = list_files(files, extension=[".fastq", ".fq", ".fastq.gz", ".fq.gz"])
    if schema.lower() == "illumina":
        files = [IlluminaFile(f) for f in files]
    elif schema.lower() == "element":
        files = [ElementFile(f) for f in files]
    else:
        raise ValueError(f"Invalid schema: {schema}. Must be 'illumina' or 'element'.")

    # group files by sample
    group_dict = {}
    for f in files:
        if f.name not in group_dict:
            group_dict[f.name] = []
        group_dict[f.name].append(f)
    groups = [MergeGroup(name, group_files) for name, group_files in group_dict.items()]
    return natsorted(groups, key=lambda x: x.name)


def merge_fastqs_vsearch(
    forward: str,
    reverse: Optional[str],
    merged_file: str,
    output_format: str = "fasta",
    binary_path: Optional[str] = None,
    minimum_overlap: int = 30,
    allowed_mismatches: int = 5,
    allowed_mismatch_percent: float = 20.0,
    additional_args: Optional[str] = None,
    interleaved: bool = False,
    debug: bool = False,
) -> str:
    """
    Merge paired-end reads using vsearch.

    Parameters
    ----------
    forward : str
        Path to the forward read file.

    reverse : str
        Path to the reverse read file.

    merged_file : str
        Path to the merged read file.

    output_format : str, optional
        Output format. Must be 'fasta' or 'fastq'. Default is 'fasta'.

    binary_path : str, optional
        Path to a vsearch binary. If not provided, the vsearch binary packaged with abutils will be used.

    additional_args : str, optional
        Additional arguments (as a string) to pass to vsearch.

    debug : bool, optional
        If True, print vsearch stdout and stderr. Default is False.

    Returns
    -------
    str
        Path to the merged read file.

    """
    # validate input files
    if not os.path.isfile(forward):
        err = f"The supplied forward read file path ({forward}) does not exist or is not a file."
        raise ValueError(err)
    if reverse is not None and not os.path.isfile(reverse):
        err = f"The supplied reverse read file path ({reverse}) does not exist or is not a file."
        raise ValueError(err)
    if reverse is None and not interleaved:
        err = "Reverse read file path is required when not using interleaved mode."
        raise ValueError(err)
    # make output directory
    out_dir = os.path.dirname(merged_file)
    make_dir(out_dir)
    # get the vsearch binary
    if binary_path is None:
        binary_path = get_binary_path("vsearch")
    # compile the vsearch command
    cmd = f"{binary_path} --fastq_mergepairs {forward} --reverse {reverse}"
    if output_format.lower() == "fasta":
        cmd += " --fastaout {merged_file}"
    elif output_format.lower() == "fastq":
        cmd += " --fastqout {merged_file}"
    else:
        err = f"Invalid output format: {output_format}. Must be 'fasta' or 'fastq'."
        raise ValueError(err)
    if additional_args is not None:
        cmd += f" {additional_args}"
    # merge reads
    p = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    if debug:
        print(stdout.decode())
        print(stderr.decode())
    if p.returncode != 0:
        err = f"Error merging reads with vsearch: {stderr.decode()}"
        raise ValueError(err)
    return merged_file


def merge_fastqs_fastp(
    forward: str,
    reverse: Optional[str] = None,
    merged: Optional[str] = None,
    binary_path: Optional[str] = None,
    minimum_overlap: int = 30,
    allowed_mismatches: int = 5,
    allowed_mismatch_percent: int = 20,
    correct_overlap_region: bool = False,
    trim_adapters: bool = True,
    adapter_file: Optional[str] = None,
    quality_trim: bool = True,
    window_size: int = 4,
    quality_cutoff: int = 20,
    name: Optional[str] = None,
    log_directory: Optional[str] = None,
    interleaved: bool = False,
    additional_args: Optional[str] = None,
    debug: bool = False,
) -> str:
    """
    Merge paired-end reads using fastp.

    Parameters
    ----------
    forward : str
        Path to the forward read file.

    reverse : str
        Path to the reverse read file. If `interleaved` is ``True``, the merged output file can optionally
        be provided to this parameter. This allows calling the function using positional arguments for the
        input/output files, like so:

        ```python
        merged_file = merge_fastqs_fastp("path/to/forward.fastq", "path/to/merged.fastq", interleaved=True)
        ```

    merged : str
        Path to the merged read file. If the merged file ends in ".gz", the merged file
        will be gzip compressed.

    binary_path : str, optional
        Path to a fastp binary. If not provided, the fastp binary packaged with abutils will be used.

    minimum_overlap : int, optional
        Minimum overlap between reads. Default is 30.

    allowed_mismatches : int, optional
        Allowed mismatches between reads. Default is 5.

    allowed_mismatch_percent : float, optional
        Allowed mismatch percentage between reads. Default is 20.0.

    correct_overlap_region : bool, optional
        If True, correct the overlap region. Default is False.

    trim_adapters : bool, optional
        If True, trim adapters. Default is True.

    adapter_file : str, optional
        Path to a FASTA file containing adapter sequences. If not provided, the default
        adapter sequences (Illumina TruSeq) will be used.

    quality_trim : bool, optional
        If True, trim low-quality bases. Default is True.

    window_size : int, optional
        Sliding window size for quality trimming. Default is 4.

    quality_cutoff : int, optional
        Mean quality cutoff for quality trimming. Default is 20.

    name : str, optional
        Name of the sample. If not provided, the name of the merged file will be used.

    log_directory : str, optional
        Path to the directory in which to save the fastp report.
        If not provided, fastp reports will not be saved.

    additional_args : str, optional
        Additional arguments (as a string) to pass to fastp.

    debug : bool, optional
        If True, print vsearch stdout and stderr. Default is False.

    Returns
    -------
    str
        Path to the merged read file.

    """
    # validate input files
    if not os.path.isfile(forward):
        err = f"The supplied forward read file path ({forward}) does not exist or is not a file."
        raise ValueError(err)
    if reverse is not None and not os.path.isfile(reverse):
        err = f"The supplied reverse read file path ({reverse}) does not exist or is not a file."
        raise ValueError(err)
    if reverse is None and not interleaved:
        err = "Reverse read file path is required when not using interleaved mode."
        raise ValueError(err)
    if merged is None:
        if interleaved and reverse is not None:
            merged = os.path.abspath(reverse)
            reverse = None
        else:
            raise ValueError("You must supply a path for the merged output file.")

    # make output directory
    merged = os.path.abspath(merged)
    out_dir = os.path.dirname(merged)
    make_dir(out_dir)

    # get the fastp binary
    if binary_path is None:
        binary_path = get_binary_path("fastp")

    # compile the fastp command
    cmd = f"{binary_path} -i '{forward}' --merge --merged_out '{merged}'"
    if not interleaved:
        cmd += f" -I '{reverse}'"
    cmd += f" --overlap_len_require {int(minimum_overlap)}"
    cmd += f" --overlap_diff_limit {int(allowed_mismatches)}"
    cmd += f" --overlap_diff_percent_limit {int(allowed_mismatch_percent)}"
    if correct_overlap_region:
        cmd += " --correction"
    if interleaved:
        cmd += " --interleaved_in"

    # adapters
    if trim_adapters:
        if adapter_file is None:
            # default Illumina TruSeq adapters
            cmd += " --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
            cmd += " --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
        else:
            cmd += f" --adapter_fasta {adapter_file}"
    else:
        cmd += " --disable_adapter_trimming"

    # quality
    if quality_trim:
        cmd += " --cut_tail"
        cmd += f" --cut_tail_window_size {int(window_size)}"
        cmd += f" --cut_tail_mean_quality {int(quality_cutoff)}"
    else:
        cmd += " --disable_quality_filtering"

    # log
    if log_directory is None:
        log_directory = tempfile.mkdtemp()
    make_dir(log_directory)
    if name is None:
        if merged.endswith(".gz"):
            name = ".".join(os.path.basename(merged).split(".")[:-2])
        else:
            name = ".".join(os.path.basename(merged).split(".")[:-1])
    cmd += f" --html '{os.path.join(log_directory, name)}_fastp-report.html'"
    cmd += f" --json '{os.path.join(log_directory, name)}_fastp-report.json'"

    # additional CLI args
    if additional_args is not None:
        cmd += f" {additional_args}"

    # merge reads
    p = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    if debug:
        print(stdout.decode())
        print(stderr.decode())
    if p.returncode != 0:
        err = f"Error merging reads with fastp: {stderr.decode()}"
        raise ValueError(err)
    return merged


def _print_files_to_merge(files: Iterable[str], name: str = "file") -> None:
    num_files = len(files)
    plural = f"{name}s" if num_files > 1 else name
    print("")
    print(f"found {num_files} {plural} to merge:")
    if num_files < 6:
        for f in files:
            print(f"  {f}")
    else:
        for f in files[:5]:
            print(f"  {f}")
        print(f"  ... and {num_files - 5} more")
    print("")
