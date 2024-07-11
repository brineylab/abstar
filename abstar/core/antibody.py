#!/usr/bin/env python
# filename: antibody.py

#
# Copyright (c) 2024 Bryan Briney
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


from dataclasses import dataclass, field
from typing import Iterable, Optional, Union

from ..utils.mixins import LoggingMixin


@dataclass
class Antibody(LoggingMixin):
    """
    Class for storing antibody annotation data.

    Includes the LoggingMixin, which provides methods for logging and exception handling.

    """

    # most useful info up front
    sequence_id: str = None
    v_gene: str = None
    d_gene: str = None
    j_gene: str = None
    c_gene: str = None
    cdr3_length: int = None
    junction_aa: str = None
    v_identity: float = None
    v_identity_aa: float = None
    d_identity: float = None
    d_identity_aa: float = None
    j_identity: float = None
    j_identity_aa: float = None
    productive: bool = True
    complete_vdj: bool = False

    # everything else
    sequence: str = None
    germline: str = None
    sequence_aa: str = None
    germline_aa: str = None
    sequence_alignment: str = None
    germline_alignment: str = None
    sequence_alignment_aa: str = None
    germline_alignment_aa: str = None
    umi: str = None
    quality: str = None
    locus: str = None
    species: str = None
    germline_database: str = None
    sequence_input: str = None
    sequence_oriented: str = None
    rev_comp: bool = False
    productivity_issues: list = field(default_factory=list)
    stop_codon: bool = False
    v_call: str = None
    v_score: float = None
    v_support: float = None
    v_cigar: str = None
    v_sequence: str = None
    v_germline: str = None
    v_mutations: str = None
    v_mutations_aa: str = None
    v_mutation_count: int = None
    v_mutation_count_aa: int = None
    v_insertions: str = None
    v_deletions: str = None
    v_frameshift: bool = False
    frame: str = None
    d_call: str = None
    d_score: float = None
    d_support: float = None
    d_cigar: str = None
    d_sequence: str = None
    d_germline: str = None
    d_frame: str = None
    j_call: str = None
    j_score: float = None
    j_support: float = None
    j_cigar: str = None
    j_sequence: str = None
    j_germline: str = None
    c_call: str = None
    c_score: float = None
    c_identity: float = None
    c_identity_aa: float = None
    c_cigar: str = None
    c_sequence: str = None
    c_germline: str = None
    np1: str = None
    np2: str = None
    np1_length: int = None
    np2_length: int = None
    fwr1: str = None
    fwr1_aa: str = None
    cdr1: str = None
    cdr1_aa: str = None
    fwr2: str = None
    fwr2_aa: str = None
    cdr2: str = None
    cdr2_aa: str = None
    fwr3: str = None
    fwr3_aa: str = None
    cdr3: str = None
    cdr3_aa: str = None
    fwr4: str = None
    fwr4_aa: str = None
    junction: str = None
    v_sequence_start: int = None
    v_sequence_end: int = None
    v_germline_start: int = None
    v_germline_end: int = None
    j_sequence_start: int = None
    j_sequence_end: int = None
    j_germline_start: int = None
    j_germline_end: int = None

    def __post_init__(self):
        # establish the list of AIRR output fields
        self.airr_fields = list(self.__dict__.keys())

        # initialize the LoggingMixin
        super().__init__()

    def to_dict(
        self,
        include: Optional[Union[Iterable, str]] = None,
        exclude: Optional[Union[Iterable, str]] = None,
    ) -> dict:
        """
        Convert the Antibody object to a dictionary of annotations.

        Parameters:
        ----------
        include : Iterable or str, default: None
            Fields to include in the dictionary, in addition to the default fields.

        exclude : Iterable or str, default: None
            Fields to exclude from the dictionary.

        Returns:
        --------
        dict: The dictionary representation of the antibody.

        """
        airr_fields = self.airr_fields

        # excluded fields
        if exclude is not None:
            if isinstance(exclude, str):
                exclude = [exclude]
            airr_fields = [a for a in airr_fields if a not in exclude]

        # included fields
        if include is not None:
            if isinstance(include, str):
                include = [include]
            airr_fields.extend(include)

        return {k: self.__dict__.get(k, None) for k in airr_fields}
