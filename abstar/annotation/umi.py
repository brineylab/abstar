# Copyright (c) 2024 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

import os
import shutil
import tempfile
from typing import Iterable, Optional, Union

import abutils
from abutils import Sequence
from abutils.tl import PairwiseAlignment

__all__ = ["parse_umis"]


class UMI:
    """
    Class for parsing UMIs from sequence data.
    """

    def __init__(
        self,
        sequence: Sequence,
        pattern: str,
        length: Optional[int] = None,
        ignore_strand: bool = False,
        # extra_length_for_alignment: int = 25,
    ):
        """
        Parameters
        ----------
        sequence : str or Sequence
            Raw input sequence

        pattern : str
            Pattern for identifying the UMI

        length : int, default=None
            Length of the UMI

        ignore_strand : bool, default=False
            Ignore the strand of the sequence

        extra_length_for_alignment : int, default=25
            Extra length added to the query sequence during pattern alignment. This allows for some ambiguity
            in the location of the UMI without needing to align the complete sequence (saving compute time).

        """
        self.length = length
        self.pattern = pattern if pattern is None else pattern.strip().upper()
        self.ignore_strand = ignore_strand
        # self.extra_length_for_alignment = extra_length_for_alignment
        self.sequence = self.process_sequence(sequence)
        self._leading_aln = None
        self._trailing_aln = None
        self._umi = None
        self._num_mismatches = None

    @property
    def leading(self) -> Optional[str]:
        """
        Returns the leading conserved region of the pattern, if present.
        """
        if self.pattern is not None:
            leading = self.pattern.split("[UMI]")[0]
            return leading if leading else None

    @property
    def trailing(self) -> Optional[str]:
        """
        Returns the trailing conserved region of the pattern, if present.
        """
        if self.pattern is not None:
            trailing = self.pattern.split("[UMI]")[-1]
            return trailing if trailing else None

    @property
    def leading_aln(self) -> Optional[PairwiseAlignment]:
        """
        Returns the leading conserved region of the pattern, if present.
        """
        if self.pattern is not None:
            if self._leading_aln is None:
                if self.leading is not None:
                    self._leading_aln = self.align(self.leading, self.sequence)
            return self._leading_aln

    @property
    def trailing_aln(self) -> Optional[PairwiseAlignment]:
        """
        Returns the trailing conserved region of the pattern, if present.
        """
        if self.pattern is not None:
            if self._trailing_aln is None:
                if self.trailing is not None:
                    self._trailing_aln = self.align(self.trailing, self.sequence)
            return self._trailing_aln

    @property
    def umi(self) -> Optional[str]:
        """
        Returns the UMI, if present.
        """
        if self._umi is None:
            self._umi = self.get_umi()
        return self._umi

    @property
    def num_mismatches(self) -> Optional[int]:
        """
        Returns the number of mismatches between the pattern and sequence.
        """
        if self._num_mismatches is None:
            self._num_mismatches = self.get_mismatches()
        return self._num_mismatches

    def process_sequence(self, sequence: Sequence) -> str:
        """
        Processes the raw sequence:
          * compute the reverse complement if length is negative
          * truncate sequence prior to pattern alignment

        Parameters
        ----------
        sequence : str or Sequence
            Raw sequence

        Returns
        -------
        truncated : str
        """
        if self.pattern is not None:
            # length_for_aln = (
            #     len(self.pattern.replace("[UMI]", "")) + self.extra_length_for_alignment
            # )
            # length_for_aln += self.length if self.length is not None else 16
            # sequence = sequence[:length_for_aln]
            if self.length < 0 and not self.ignore_strand:
                sequence = abutils.tl.reverse_complement(sequence)
        return sequence

    def align(self, pattern: str, sequence: Sequence) -> PairwiseAlignment:
        """
        Performs local alignment between a pattern and sequence.
        """
        aln = abutils.tl.local_alignment(
            pattern, sequence, gap_open=-20, gap_extend=-20
        )
        return aln

    def get_mismatches(self) -> Optional[int]:
        """
        Returns the total number of mismatches between the pattern(s)
        and the input sequence.
        """
        if self.pattern is None:
            return 0
        total = len(self.pattern.replace("[UMI]", ""))
        matches = 0
        for aln in [self.leading_aln, self.trailing_aln]:
            if aln is not None:
                matches += aln.alignment_midline.count("|")
        return total - matches

    def get_umi(self) -> Optional[str]:
        """
        Parses UMI from the input sequence using the UMI pattern(s).
        """
        # length without pattern means we just slice
        if self.pattern is None:
            if self.length < 0:
                return self.sequence[self.length :]
            else:
                return self.sequence[: self.length]
        ## something like "ATGC[UMI]"
        ## we need to take `length` residues following the end of the leading alignment
        if self.trailing is None and self.leading is not None:
            start = self.leading_aln.target_end + 1
            end = start + self.length
            if end > len(self.sequence):
                return None
        ## something like "[UMI]ATGC"
        ## we need to take `length` residues preceeding the start of the trailing alignment
        elif self.leading is None and self.trailing is not None:
            end = self.trailing_aln.target_begin
            if self.length is not None:
                start = end - self.length
            else:
                start = 0
            ## if self.length is longer than the region preceeding the start of the
            ## trailing alignment, return None because we can't parse a complete UMI
            if start < 0:
                return None
        else:
            start = self.leading_aln.target_end + 1
            end = self.trailing_aln.target_begin
        return self.sequence[start:end]


def parse_umis(
    sequences: Union[str, Sequence, Iterable, dict],
    output_file: Optional[str] = None,
    pattern: Union[str, Iterable, None] = None,
    length: Union[int, Iterable, None] = None,
    allowed_mismatches: int = 1,
    # extra_length_for_alignment: int = 25,
    ignore_strand: bool = False,
    fmt: str = "fasta",
    sequence_key: str = "sequence_input",
    id_key: str = "sequence_id",
) -> str:
    """
    Parses unique molecular identifiers (UMIs) from bulk NGS sequencing data.

    Parameters
    ----------
    sequences : str or iterable
        Can be one of several things:
            1. path to a FASTA- or FASTQ-formatted file
            2. a single sequence, as a string or a ``Sequence`` object
            3. a list of ``Sequence`` objects or anything accepted by ``abutils.Sequence()``

    output_file : str, default=None
        Path to an output file. If `sequences` is not a file, and `output_file` is not
        provided `sequences` are converted to ``abutils.Sequence`` objects and the UMI
        is added to the ``"umi"`` field of ``sequence.annotations``. If `sequences` is
        a file and `output` is not provided, UMI parsing is done in-place and the
        `sequences` file is updated. UMIs are added to the end of each sequence ID with
        an underscore (``"_"``) separating the name and UMI. Sequences for which at least
        one UMI was not found are not included in the output.

    pattern : str or iterable, default=None
        Pattern (or iterable of patterns) for identifying the location of the UMI,
        or the name of a built-in pattern. Built-in options include ``"smartseq-human-bcr"``.
        Patterns may optionally contain leading and/or trailing conserved regions, with
        the UMI position within the pattern represented by ``"[UMI]"``. As an example,
        one of the built-in patterns for SmartSeq-BCR UMIs is::

            "[UMI]TCAGCGGGAAGACATT"

        which would be a UMI sequence followed immediately by ``"TCAGCGGGAAGACATT"``.
        By default, the pattern is matched to the 5' -> 3' strand. This allows
        users to more easily construct patterns from their amplification primers without
        needing to worry about reverse-complementing patterns for UMIs at the 3' end of the
        input sequence. To override this, set `ignore_strand` to ``True``. If `pattern`
        is not provided, UMIs will be parsed using only `length`, starting at the start
        (if `length` is positive) or end (if `length` is negative) of the sequence.

        .. note::

            The UMIs for all `patterns` that meet the `allowed_mismatches` criteria in
            the conserved leading/trailing portions will be concatenated into the final
            UMI. This allows the use of multiple `patterns` for either:
              - different patterns designed to match one of several different types of
                sequences in a heterogeneous sample. For example, if heavy, kappa and
                lambda primers each have different conserved regions flanking the UMI
                and the input file contains a mix of heavy, kappa and lambda chains,
                supplying all patterns (assuming they're sufficiently different from
                each other) will allow parsing of UMIs from all chains.
              - different patterns for sequences that contain multiple UMIs, either
                at opposite ends of the sequence or on the same end of the sequence
                but in different locations and with different conserved flanking regions.

        If multiple UMIs are parsed from a single sequence, they will be concatenated with
        ``"+"`` separating the UMIs.

    length : int or iterable, default=None
        Length of the UMI sequence, or iterable of UMI lengths. If `length` is provided
        and pattern is not, the UMI will be sliced directly from the start (if `length` is
        positive) or the end (if `length` is negative) of the sequencing without attempting
        to match UMI-flanking sequence patterns. If multiple lengths are
        provided, there must be an equal number of `patterns`, and they should be in the
        same order (the first `pattern` should correspond to the first UMI `length`). If
        `length` is positive, the UMI will be parsed from the start of the sequence. If
        `length` is negative, the UMI will be parsed from the end of the sequence. If multiple
        `patterns` are provided with a single `length`, that `length` will be used for all
        `patterns`. Required if `pattern` does not have a conserved trailing region.
        If `length` is not provided and `pattern` has a trailing conserved region but not a
        leading conserved region, the entire portion of `sequence` preceeding the trailing
        region will be parsed as the UMI. Ignored if both leading and trailing sequences
        are present in `pattern`, as the entire region between the conserved flanking regions
        will be parsed as the UMI regardless of what `length` is provided.

    allowed_mismatches : int, default=None
        Total number of mismatches allowed when aligning the conserved flanking regions
        of `pattern`. Mismatches are calculated by subtracting the number of identically
        matching alignment positions from the total length of the conserved flanking region.
        If not provided, the default is to allow ``1`` mismatch.

    extra_length_for_alignment : int, default=25
        To speed alignment and avoid coincidental matches in regions of the sequence that are
        not logical locations for UMIs, only a fragment of each input sequences (sliced from
        the start or end of the seqeunce, depending on `length`) is aligned against conserved
        `pattern` regions. The size of this fragment is computed by summing the total length of
        conserved flanking regions, the UMI `length`, and `extra_length_for_alignment`. The default
        is ``25``, which should be sufficient for most cases, however, `extra_length_for_alignment`
        may need to be increased if the UMI is not located close to the start or end of the
        input seqeunce.

    ignore_strand : bool, default=False
        If ``True``, patterns matching to the end of the input sequence (those with a negative
        `length`) will be aligned to the unmodified input sequence rather than the reverse
        complement. This means that the supplied `pattern` must be in the 3' -> 5' orientation.

    fmt : str, default="fasta"
        Format of the output file. Only used if `output_file` is ``True``. Options are ``"fasta"``
        and ``"fastq"``.

        .. warning::

            If selecting ``"fastq"`` output when providing a FASTA-formatted input file,
            the resulting FASTQ file may be malformed as there is no quality information
            in the input file.

    sequence_key : str, default="raw_input"
        Field in `sequences` from which the UMI should be parsed. Only used if `sequences`
        are ``abutils.Sequence`` objects. Default is ``"raw_input"``.  If `sequence_key` is
        not found, or is ``None`` ``sequence.sequence`` will be used.

    id_key : str, default="sequence_id"
        Field in `sequences` containing the sequence ID. Only used if `sequences` are
        ``abutils.Sequence`` objects. Default is ``"sequence_id"``. If `id_key` is not found
        or is ``None`` ``sequence.id`` will be used.


    Returns
    -------
    output_file : str
        Path to the output file, with a format (``"fasta"`` or ``"fastq"``) matching `fmt`.

    """
    # lengths and patterns
    if pattern is None:
        if isinstance(length, int):
            patterns = [None]
            lengths = [length]
        else:
            patterns = [None] * len(length)
            lengths = length
    elif pattern is not None:
        if isinstance(pattern, str):
            if pattern.lower() in BUILTIN_PATTERNS:
                name = pattern.lower()
                patterns = _get_patterns(name)
                length = _get_lengths(name)
                if allowed_mismatches is None:
                    allowed_mismatches = _get_allowed_mismatches(name)
            else:
                patterns = [pattern]
        else:
            patterns = pattern
        if isinstance(length, int) or length is None:
            lengths = [length] * len(patterns)
        else:
            lengths = length
        if len(patterns) != len(lengths):
            err = "\nERROR: if pattern and length are both iterables,"
            err += " they must be the same length.\n"
            err += f"  - pattern has {len(patterns)} elements: {', '.join(patterns)}\n"
            err += f"  - length has {len(lengths)} elements: {', '.join(lengths)}\n"
            raise ValueError(err)
    # validate format
    fmt = fmt.lower()
    if fmt not in ["fasta", "fastq"]:
        err = "\nERROR: format must be either 'fasta' or 'fastq'. "
        err += f"You provided {fmt}.\n"
        raise ValueError(err)
    # input processing
    if isinstance(sequences, str) and os.path.isfile(sequences):
        output = _parse_umis_from_file(
            input_file=sequences,
            patterns=patterns,
            lengths=lengths,
            output_file=output_file,
            allowed_mismatches=allowed_mismatches,
            # extra_length_for_alignment=extra_length_for_alignment,
            ignore_strand=ignore_strand,
            fmt=fmt,
        )
    elif isinstance(sequences, (str, Sequence)):
        output = _parse_umis_from_single_sequence(
            sequence=sequences,
            patterns=patterns,
            lengths=lengths,
            allowed_mismatches=allowed_mismatches,
            # extra_length_for_alignment=extra_length_for_alignment,
            ignore_strand=ignore_strand,
        )
    else:
        output = _parse_umis_from_sequences(
            sequences=sequences,
            patterns=patterns,
            lengths=lengths,
            output_file=output_file,
            allowed_mismatches=allowed_mismatches,
            # extra_length_for_alignment=extra_length_for_alignment,
            ignore_strand=ignore_strand,
            fmt=fmt,
            sequence_key=sequence_key,
            id_key=id_key,
        )
    return output


def _parse_umis_from_file(
    input_file: str,
    patterns: Iterable,
    lengths: Iterable,
    output_file: Optional[str] = None,
    allowed_mismatches: int = 1,
    # extra_length_for_alignment: int = 25,
    ignore_strand: bool = False,
    separator: str = "_",
    fmt: str = "fasta",
) -> str:
    """
    Parses UMIs from a FASTA- or FASTQ-formatted file.
    """
    # set up output file
    inplace = False
    if output_file is None:
        output_file = tempfile.NamedTemporaryFile(delete=False).name
        inplace = True
    else:
        abutils.io.make_dir(os.path.dirname(output_file))
    with open(output_file, "w") as ofile:
        for s in abutils.io.parse_fastx(input_file):
            # with open(input_file, "r") as ifile:
            #     for s in SeqIO.parse(ifile, fmt.lower()):
            #         s = abutils.Sequence(s)
            umi_list = []
            for pattern, length in zip(patterns, lengths):
                u = UMI(
                    sequence=s,
                    pattern=pattern,
                    length=length,
                    # extra_length_for_alignment=extra_length_for_alignment,
                    ignore_strand=ignore_strand,
                )
                if u.num_mismatches <= allowed_mismatches:
                    umi_list.append(u.umi)
            # if we couldn't find a UMI, skip the sequence
            umi_list = [u for u in umi_list if u is not None]
            if not umi_list:
                continue
            umi = "+".join(umi_list)
            s.id = f"{s.id}{separator}{umi}"
            if fmt.lower() == "fastq":
                fstring = s.fastq + "\n"
            else:
                fstring = s.fasta + "\n"
            ofile.write(fstring)
    if inplace:
        output_file = shutil.move(output_file, input_file)
    return output_file


def _parse_umis_from_sequences(
    sequences: Iterable,
    patterns: Iterable,
    lengths: Iterable,
    output_file: str,
    allowed_mismatches: int = 1,
    # extra_length_for_alignment: int = 25,
    ignore_strand: bool = False,
    fmt: str = "fasta",
    sequence_key: str = "sequence_input",
    id_key: str = "sequence_id",
) -> str:
    """
    Parses UMIs from a list of sequences.
    """
    inplace = True if output_file is None else False
    # output is a list of either abutils.Sequence objects (if inplace == True)
    # or a list of FAST[AQ]-formatted strings (if inplace == False)
    output = []
    for s in sequences:
        if isinstance(s, Sequence):
            if sequence_key in s:
                s.sequence = s[sequence_key]
            if id_key in s:
                s.id = s[id_key]
        else:
            s = Sequence(s)
        umi_list = []
        for pattern, length in zip(patterns, lengths):
            u = UMI(
                sequence=s,
                pattern=pattern,
                length=length,
                # extra_length_for_alignment=extra_length_for_alignment,
                ignore_strand=ignore_strand,
            )
            if u.num_mismatches <= allowed_mismatches:
                umi_list.append(u.umi)
        # if we couldn't find a UMI, skip the sequence
        umi_list = [u for u in umi_list if u is not None]
        if not umi_list:
            continue
        umi = "+".join(umi_list)
        if inplace:
            s["umi"] = umi
            output.append(s)
        else:
            s.id = f"{s.id}_{umi}"
            if fmt.lower() == "fastq":
                output.append(s.fastq)
            else:
                output.append(s.fasta)
    if inplace:
        return output
    else:
        with open(output_file, "w") as ofile:
            ofile.write("\n".join(output))
        return output_file


def _parse_umis_from_single_sequence(
    sequence: Union[str, Sequence, Iterable, dict],
    patterns: Iterable,
    lengths: Iterable,
    allowed_mismatches: int = 1,
    # extra_length_for_alignment: int = 25,
    ignore_strand: bool = False,
    sequence_key: str = "sequence_input",
    id_key: str = "sequence_id",
) -> str:
    """
    Parses UMIs from a single sequence.
    """
    # process input
    s = Sequence(sequence)
    if sequence_key in s:
        s.sequence = s[sequence_key]
    if id_key in s:
        s.id = s[id_key]
    # find UMIs
    umi_list = []
    for pattern, length in zip(patterns, lengths):
        u = UMI(
            sequence=s,
            pattern=pattern,
            length=length,
            # extra_length_for_alignment=extra_length_for_alignment,
            ignore_strand=ignore_strand,
        )
        if u.num_mismatches <= allowed_mismatches:
            umi_list.append(u.umi)
    # concatenate all of the UMI segments with "+"
    umi = "+".join([u for u in umi_list if u is not None])
    # if we couldn't find a UMI (UMI is an empty string), return None
    if not umi:
        return None
    return umi


def _get_patterns(name):
    p = BUILTIN_PATTERNS.get(name, None)
    if p is not None:
        return p.get("pattern", None)


def _get_lengths(name):
    p = BUILTIN_PATTERNS.get(name, None)
    if p is not None:
        return p.get("length", None)


def _get_allowed_mismatches(name):
    p = BUILTIN_PATTERNS.get(name, None)
    if p is not None:
        return p.get("allowed_mismatches", None)


# ---------------------------
#        PATTERNS
# ---------------------------

# Takara's Smart-Seq Human BCR kit
# https://www.takarabio.com/products/next-generation-sequencing/immune-profiling/human-repertoire/smart-seq-human-bcr-with-umis
SMARTSEQ_HUMAN_BCR_PATTERN = [
    "[UMI]TCAGCGGGAAGACATT",
    "[UMI]GGGCGGATGGACTACC",
    "[UMI]CCGATGGGCGCTTGGT",
    "[UMI]GGAACACATGCGGAGC",
]
SMARTSEQ_HUMAN_BCR_LENGTH = 12
SMARTSEQ_HUMAN_BCR_MISMATCH = 2


BUILTIN_PATTERNS = {
    "smartseq-human-bcr": {
        "pattern": SMARTSEQ_HUMAN_BCR_PATTERN,
        "length": SMARTSEQ_HUMAN_BCR_LENGTH,
        "allowed_mismatches": SMARTSEQ_HUMAN_BCR_MISMATCH,
    },
}


# import os
# import shutil
# import sys
# import tempfile
# from typing import Union, Iterable, Optional

# from Bio import SeqIO
# from Bio.SeqRecord import SeqRecord

# import abutils


# class UMI:
#     """ """

#     def __init__(
#         self,
#         sequence: Union[str, abutils.Sequence, SeqRecord],
#         pattern: str,
#         length: Optional[int] = None,
#         ignore_strand: bool = False,
#         extra_length_for_alignment: int = 25,
#     ):
#         """ """
#         self.length = length
#         self.pattern = pattern.strip().upper()
#         self.ignore_strand = ignore_strand
#         self.extra_length_for_alignment = extra_length_for_alignment
#         self.sequence = self.process_sequence(sequence)
#         self._leading_aln = None
#         self._trailing_aln = None
#         self._umi = None
#         self._num_mismatches = None

#     @property
#     def leading(self) -> Optional[str]:
#         leading = self.pattern.split("[UMI]")[0]
#         return leading if leading else None

#     @property
#     def trailing(self) -> Optional[str]:
#         trailing = self.pattern.split("[UMI]")[-1]
#         return trailing if trailing else None

#     @property
#     def leading_aln(self) -> abutils.SSWAlignment:
#         if self._leading_aln is None:
#             if self.leading is not None:
#                 self._leading_aln = self.align(self.leading, self.sequence)
#         return self._leading_aln

#     @property
#     def trailing_aln(self) -> abutils.SSWAlignment:
#         if self._trailing_aln is None:
#             if self.trailing is not None:
#                 self._trailing_aln = self.align(self.trailing, self.sequence)
#         return self._trailing_aln

#     @property
#     def umi(self) -> Optional[str]:
#         if self._umi is None:
#             self._umi = self.get_umi()
#         return self._umi

#     @property
#     def num_mismatches(self) -> int:
#         """
#         Returns the number of mismatches between the pattern and sequence.
#         """
#         if self._num_mismatches is None:
#             self._num_mismatches = self.get_mismatches()
#         return self._num_mismatches

#     def process_sequence(self, sequence) -> str:
#         """
#         Processes the raw sequence:
#           * reverse complement if length is negative
#           * truncate sequence prior to pattern alignment

#         Parameters
#         ----------
#         sequence : str, Sequence, or SeqRecord
#             Raw sequence

#         Returns
#         -------
#         truncated : str
#         """
#         s = abutils.Sequence(sequence)
#         length_for_aln = (
#             len(self.pattern.replace("[UMI]", "")) + self.extra_length_for_alignment
#         )
#         if self.length is not None:
#             length_for_aln += self.length
#         if self.length < 0 and not self.ignore_strand:
#             return s.reverse_complement[:length_for_aln]
#         else:
#             return s[:length_for_aln]

#     def align(
#         self, pattern: str, sequence: Union[abutils.Sequence, str]
#     ) -> abutils.SSWAlignment:
#         """
#         Performs local alignment between a pattern and sequence.
#         """
#         aln = abutils.tl.local_alignment(
#             pattern, sequence, gap_open=-20, gap_extend=-20
#         )
#         return aln

#     def get_mismatches(self) -> Optional[int]:
#         """
#         Returns the total number of mismatches between the pattern(s)
#         and the input sequence.
#         """
#         total = len(self.pattern.replace("[UMI]", ""))
#         matches = 0
#         for aln in [self.leading_aln, self.trailing_aln]:
#             if aln is not None:
#                 matches += aln.alignment_midline.count("|")
#         return total - matches

#     def get_umi(self) -> Optional[str]:
#         """
#         Parses UMI from the input sequence using the UMI pattern(s).
#         """
#         ## something like "ATGC[UMI]"
#         ## we need to take `length` residues following the end of the leading alignment
#         if self.trailing is None and self.leading is not None:
#             start = self.leading_aln.target_end + 1
#             end = start + self.length
#             if end > len(self.sequence):
#                 return None
#         ## something like "[UMI]ATGC"
#         ## we need to take `length` residues preceeding the start of the trailing alignment
#         elif self.leading is None and self.trailing is not None:
#             end = self.trailing_aln.target_begin
#             if self.length is not None:
#                 start = end - self.length
#             else:
#                 start = 0
#             ## if self.length is longer than the region preceeding the start of the
#             ## trailing alignment, return None because we can't parse a complete UMI
#             if start < 0:
#                 return None
#         else:
#             start = self.leading_aln.target_end + 1
#             end = self.trailing_aln.target_begin
#         return self.sequence[start:end]


# def parse_umis(
#     sequences: Union[str, Iterable],
#     output_file: Optional[str] = None,
#     pattern: Union[str, Iterable, None] = None,
#     length: Union[int, Iterable, None] = None,
#     allowed_mismatches: Optional[int] = None,
#     extra_length_for_alignment: int = 25,
#     ignore_strand: bool = False,
#     fmt: str = "fasta",
#     sequence_key: str = "raw_input",
#     id_key: str = "sequence_id",
# ):
#     """
#     Parses unique molecular identifiers (UMIs) from bulk NGS sequencing data.

#     Parameters
#     ----------
#     sequences : str or iterable
#         Can be one of several things:
#             1. path to a FASTA- or FASTQ-formatted file
#             2. a list of BioPython ``SeqRecord`` objects
#             3. a list of abutils ``Sequence`` objects
#             4. a list of lists/tuples, of the format ``[sequence_id, sequence]``

#     output_file : str, default=None
#         Path to an output file. Required if `sequences` is not a
#         file. If `sequences` is a file and `output` is not provided, UMI parsing is
#         done in-place and the `sequences` file is updated.

#     pattern : str or iterable, default=None
#         Pattern (or iterable of patterns) for identifying the location of the UMI,
#         or the name of a built-in pattern. Built-in options include ``"smartseq-human-bcr"``.
#         Patterns may optionally contain leading and/or trailing conserved regions, with
#         the UMI position within the pattern represented by ``"[UMI]"``. As an example,
#         the built-in pattern for SmartSeq-BCR UMIs is::

#             "[UMI]TCAGCGGGAAGACATT"

#         which would be a UMI sequence followed immediately by ``"TCAGCGGGAAGACATT"``.
#         By default, the pattern is matched on the 5' -> 3' strand. This allows
#         users to more easily construct patterns from their amplification primers without
#         needing to worry about reverse-complementing patterns for UMIs at the 3' end of the
#         input sequence. To override this, set `ignore_strand` to ``True``. If `pattern`
#         is not provided, UMIs will be parsed using only `length`, starting at the start
#         or end of the sequence.

#         .. note::

#             The UMIs for all `patterns` that meet the `allowed_mismatches` criteria in
#             the conserved leading/trailing portions will be concatenated into the final
#             UMI. This allows the use of multiple `patterns` for either:
#               * different patterns designed to match one of several different types of
#                 sequences in a heterogeneous sample. For example, if heavy, kappa and
#                 lambda primers each have different conserved regions flanking the UMI
#                 and the input file contains a mix of heavy, kappa and lambda chains,
#                 supplying all patterns (assuming they're sufficiently different from
#                 each other) will allow parsing of UMIs from all chains.
#               * different patterns for sequences that contain multiple UMIs, either
#                 at opposite ends of the sequence or on the same end of the sequence
#                 but in different locations and with different conserved flanking regions.

#     length : int or iterable, default=None
#         Length of the UMI sequence, or iterable of UMI lengths. If multiple lengths are
#         provided, there must be an equal number of `patterns`, and they must be in the
#         same order (the first `pattern` should correspond to the first UMI `length`). If
#         `length` is positive, the UMI will be parsed from the start of the sequence. If
#         `length` is negative, the UMI will be parsed from the end of the sequence. If multiple
#         `patterns` are provided with a single `length`, that `length` will be used for all
#         `patterns`. Required if `pattern` does not have a conserved trailing region.
#         If `length` is not provided and a trailing conserved region is present in `pattern`,
#         the entire portion of `sequence` preceeding the trailing region will be parsed as the
#         UMI. Ignored if both leading and trailing sequences are present in `pattern`,
#         as the entire region between the conserved flanking regions will be parsed as
#         the UMI regardless of what `length` is provided.


#     """
#     # lengths and patterns
#     if isinstance(pattern, str):
#         if pattern.lower() in BUILTIN_PATTERNS:
#             name = pattern.lower()
#             patterns = _get_patterns(name)
#             length = _get_lengths(name)
#             if allowed_mismatches is None:
#                 allowed_mismatches = _get_allowed_mismatches(name)
#         else:
#             patterns = [pattern]
#     else:
#         patterns = pattern
#     if isinstance(length, int) or length is None:
#         lengths = [length] * len(patterns)
#     else:
#         lengths = length
#     if len(patterns) != len(lengths):
#         err = "\nERROR: if pattern and length are both iterables,"
#         err += " they must be the same length.\n"
#         err += f"  - pattern has {len(patterns)} elements: {', '.join(patterns)}\n"
#         err += f"  - length has {len(lengths)} elements: {', '.join(lengths)}\n"
#         print(err)
#         sys.exit()
#     # mismatches
#     if allowed_mismatches is None:
#         allowed_mismatches = 1
#     # validate format
#     fmt = fmt.lower()
#     if fmt not in ["fasta", "fastq"]:
#         err = "\nERROR: format must be either 'fasta' or 'fastq'. "
#         err += f"You provided {fmt}.\n"
#         print(err)
#         sys.exit()
#     # input processing
#     if isinstance(sequences, str) and os.path.isfile(sequences):
#         output_file = _parse_umis_from_file(
#             input_file=sequences,
#             patterns=patterns,
#             lengths=lengths,
#             output_file=output_file,
#             allowed_mismatches=allowed_mismatches,
#             extra_length_for_alignment=extra_length_for_alignment,
#             ignore_strand=ignore_strand,
#             fmt=fmt,
#         )
#     else:
#         if output_file is None:
#             err = "\nERROR: if input is not a FASTA- or FASTQ-formatted file, "
#             err += f"output_file must be provided.\n"
#             print(err)
#             sys.exit()
#         output_file = _parse_umis_from_sequences(
#             sequences=sequences,
#             patterns=patterns,
#             lengths=lengths,
#             output_file=output_file,
#             allowed_mismatches=allowed_mismatches,
#             extra_length_for_alignment=extra_length_for_alignment,
#             ignore_strand=ignore_strand,
#             fmt=fmt,
#             sequence_key=sequence_key,
#             id_key=id_key,
#         )


# def _parse_umis_from_file(
#     input_file: str,
#     patterns: Iterable,
#     lengths: Iterable,
#     output_file: Optional[str] = None,
#     allowed_mismatches: int = 1,
#     extra_length_for_alignment: int = 25,
#     ignore_strand: bool = False,
#     fmt: str = "fasta",
# ) -> str:
#     """
#     Parses UMIs from a FASTA- or FASTQ-formatted file.
#     """
#     # set up output file
#     inplace = False
#     if output_file is None:
#         output_file = tempfile.NamedTemporaryFile(delete=False).name
#         inplace = True
#     with open(output_file, "w") as ofile:
#         with open(input_file, "r") as ifile:
#             for s in SeqIO.parse(ifile, fmt.lower()):
#                 s = abutils.Sequence(s)
#                 umi_list = []
#                 for p, l in zip(patterns, lengths):
#                     u = UMI(
#                         sequence=s,
#                         pattern=p,
#                         length=l,
#                         extra_length_for_alignment=extra_length_for_alignment,
#                         ignore_strand=ignore_strand,
#                     )
#                     if u.num_mismatches <= allowed_mismatches:
#                         umi_list.append(u.umi)
#                 # if we couldn't find a UMI, skip the sequence
#                 umi_list = [u for u in umi_list if u is not None]
#                 if not umi_list:
#                     continue
#                 umi = "+".join(umi_list)
#                 s.id = f"{s.id}_{umi}"
#                 if fmt.lower() == "fastq":
#                     fstring = s.fastq + "\n"
#                 else:
#                     fstring = s.fasta + "\n"
#                 ofile.write(fstring)
#     if inplace:
#         output_file = shutil.move(output_file, input_file)
#     return output_file


# def _parse_umis_from_sequences(
#     sequences: Iterable,
#     patterns: Iterable,
#     lengths: Iterable,
#     output_file: str,
#     allowed_mismatches: int = 1,
#     extra_length_for_alignment: int = 25,
#     ignore_strand: bool = False,
#     fmt: str = "fasta",
#     sequence_key: str = "raw_input",
#     id_key: str = "sequence_id",
# ) -> str:
#     """
#     Parses UMIs from a list of sequences.
#     """
#     with open(output_file, "w") as ofile:
#         for s in sequences:
#             if isinstance(s, abutils.Sequence):
#                 if sequence_key in s:
#                     s.sequence = s[sequence_key]
#                 if id_key in s:
#                     s.id = s[id_key]
#             s = abutils.Sequence(s)
#             umi_list = []
#             for p, l in zip(patterns, lengths):
#                 u = UMI(
#                     sequence=s,
#                     pattern=p,
#                     length=l,
#                     extra_length_for_alignment=extra_length_for_alignment,
#                     ignore_strand=ignore_strand,
#                 )
#                 if u.num_mismatches <= allowed_mismatches:
#                     umi_list.append(u.umi)
#             # if we couldn't find a UMI, skip the sequence
#             umi_list = [u for u in umi_list if u is not None]
#             if not umi_list:
#                 continue
#             umi = "+".join(umi_list)
#             s.id = f"{s.id}_{umi}"
#             if fmt.lower() == "fastq":
#                 fstring = s.fastq + "\n"
#             else:
#                 fstring = s.fasta + "\n"
#             ofile.write(fstring)
#     return output_file


# def _get_patterns(name):
#     p = BUILTIN_PATTERNS.get(name, None)
#     if p is not None:
#         return p.get("pattern", None)


# def _get_lengths(name):
#     p = BUILTIN_PATTERNS.get(name, None)
#     if p is not None:
#         return p.get("length", None)


# def _get_allowed_mismatches(name):
#     p = BUILTIN_PATTERNS.get(name, None)
#     if p is not None:
#         return p.get("allowed_mismatches", None)


# # ---------------------------
# #        PATTERNS
# # ---------------------------

# # Takara's Smart-Seq Human BCR kit
# # https://www.takarabio.com/products/next-generation-sequencing/immune-profiling/human-repertoire/smart-seq-human-bcr-with-umis
# SMARTSEQ_HUMAN_BCR_PATTERN = [
#     "[UMI]TCAGCGGGAAGACATT",
#     "[UMI]GGGCGGATGGACTACC",
#     "[UMI]CCGATGGGCGCTTGGT",
#     "[UMI]GGAACACATGCGGAGC",
# ]
# SMARTSEQ_HUMAN_BCR_LENGTH = 12
# SMARTSEQ_HUMAN_BCR_MISMATCH = 2


# BUILTIN_PATTERNS = {
#     "smartseq-human-bcr": {
#         "pattern": SMARTSEQ_HUMAN_BCR_PATTERN,
#         "length": SMARTSEQ_HUMAN_BCR_LENGTH,
#         "allowed_mismatches": SMARTSEQ_HUMAN_BCR_MISMATCH,
#     },
# }
