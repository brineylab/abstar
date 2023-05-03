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


from __future__ import absolute_import, division, print_function, unicode_literals

import math
import os
import traceback
from typing import Optional, Union, Iterable

import parasail

from Bio import SeqIO
from Bio.Seq import Seq

from abutils.core.sequence import Sequence
from abutils.utils.alignment import global_alignment, local_alignment
from abutils.utils.codons import codon_lookup
from abutils.utils.decorators import lazy_property

# from .antibody import Antibody
from ..utils import indels
from ..utils.mixins import LoggingMixin


class GermlineSegment(LoggingMixin):
    """
    docstring for Germline

    Parameters
    ----------
    raw : str
        Raw germline assignment, which includes the species name.

    species : str
        Species from which the germline gene was derived. Choices include
        'human', 'macaque', mouse', rabbit'.

    db_name : str
        Name of the germline database, as a string.

    receptor : str
        Receptor name. Options are ``"bcr"`` and ``"tcr"``.

    score : int or float, default=None
        Score of the top germline alignment. Optional. If ``score`` is not
        provided, the realignment score will be used instead.

    strand : str, default=None
        Strand of the alignment. Options are ``'+'`` for the positive strand
        (meaning the germline gene and the query sequence are in the same orientation)
        or ``'-'`` for the negative strand (query sequence is the reverse complement of
        the germline gene).

    others : list(GermlineSegment), default=None
        A list of additional high scoring germline genes. Each element of the list
        should be a ``GermlineSegment``.

    assigner_name : str, default=None
        The assigner name. Will be converted to lowercase. If not provided,
        `assigner_name` will be set to ``'unknown'``.

    initialize_log : bool, default=True
        Whether or not to initialize logging.

    """

    def __init__(
        self,
        raw: str,
        species: str,
        db_name: str,
        receptor: str,
        score: Optional[Union[float, int]] = None,
        strand: Optional[str] = None,
        others: Optional[Iterable] = None,
        assigner_name: Optional[str] = None,
        initialize_log: bool = True,
    ):
        super(GermlineSegment, self).__init__()
        LoggingMixin.__init__(self)
        self.raw_assignment = raw
        self.full = raw.split("__")[0]
        self.species = species.lower()
        self.db_name = db_name
        self.receptor = receptor.lower()
        self.assigner_score = score
        self.strand = strand
        self.others = others
        self.gene_type = self.full[3].upper()
        self.assigner = (
            assigner_name.lower() if assigner_name is not None else "unknown"
        )
        self._family = None
        self._gene = None
        self._chain = None

        # Optional properties for assigners to populate.
        # If populated by an assigner, they will be used to
        # force realignment to these parameters.
        #
        # Note that the query_start/query_end positions should
        # be 0-indexed and apply to the full query sequence. So if the sequence
        # is truncated following V-gene assignment and the truncated
        # sequence is submitted for J-gene assignment, make sure to
        # adjust the positions for the J-gene so that the numberings apply to
        # the complete input sequence, not the truncated sequence that was
        # used to identify the J-gene.
        #
        # Also note that because the end position of Python's slicing operator
        # is exclusive, to get the aligned sequence from the oriented_input, you
        # need to add one to the query_end (so that the slice includes query_end)
        self.query_start = None
        self.query_end = None
        self.germline_start = None
        self.germline_end = None

        # initialize log
        if initialize_log:
            self.initialize_log()

        # These properties are populated by AbStar.
        # Assigners don't need to populate these (and they'll be overwritten
        # by AbStar if an assigner does populate them).
        self.score = None
        self.realignment = None
        self.raw_query = None
        self.raw_germline = None
        self.query_alignment = None
        self.germline_alignment = None
        self.alignment_midline = None
        self.alignment_length = None
        self.alignment_reading_frame = None  # 0-based, so if alignment start is the start of a codon, reading frame will be 0
        self.imgt_germline = None
        self.imgt_gapped_alignment = None
        self.imgt_nt_positions = []
        self.imgt_aa_positions = []
        self._imgt_position_from_raw = {}
        self._raw_position_from_imgt = {}
        self._initial_correct_imgt_nt_position_from_imgt = None
        self._initial_correct_imgt_aa_position_from_imgt = None
        self._correct_imgt_nt_position_from_imgt = None
        self._correct_imgt_aa_position_from_imgt = None
        self.fs_indel_adjustment = 0
        self.nfs_indel_adjustment = 0
        self.has_insertion = "no"
        self.has_deletion = "no"
        self._insertions = None
        self._deletions = None
        self.coding_region = None
        self.aa_sequence = None
        self.regions = None

    @property
    def family(self):
        if self._family is None:
            if "-" in self.full:
                self._family = self.full.split("-")[0]
        return self._family

    @family.setter
    def family(self, family):
        self._family = family

    @property
    def gene(self):
        if self._gene is None:
            if "*" in self.full:
                self._gene = self.full.split("*")[0]
        return self._gene

    @gene.setter
    def gene(self, gene):
        self._gene = gene

    @property
    def chain(self):
        if self._chain is None:
            c = {
                "H": "heavy",
                "K": "kappa",
                "L": "lambda",
                "A": "alpha",
                "B": "beta",
                "G": "gamma",
                "D": "delta",
            }
            self._chain = c.get(self.full[2], None)
        return self._chain

    @property
    def insertions(self):
        if self._insertions is None:
            return []
        return self._insertions

    @insertions.setter
    def insertions(self, insertions):
        self._insertions = insertions

    @property
    def deletions(self):
        if self._deletions is None:
            return []
        return self._deletions

    @deletions.setter
    def deletions(self, deletions):
        self._deletions = deletions

    def correct_imgt_nt_position_from_imgt(self, position):
        if self._correct_imgt_nt_position_from_imgt is not None:
            p = self._correct_imgt_nt_position_from_imgt.get(position, None)
            if p is not None:
                return p
        if self._initial_correct_imgt_nt_position_from_imgt is not None:
            p = self._initial_correct_imgt_nt_position_from_imgt.get(position, None)
            if p is not None:
                return p
        return None

    def correct_imgt_aa_position_from_imgt(self, position):
        if self._correct_imgt_aa_position_from_imgt is not None:
            p = self._correct_imgt_aa_position_from_imgt.get(position, None)
            if p is not None:
                return p
        if self._initial_correct_imgt_aa_position_from_imgt is not None:
            p = self._initial_correct_imgt_aa_position_from_imgt.get(position, None)
            if p is not None:
                return p
        return None

    def initialize_log(self):
        log = []
        log.append("GERMLINE: {}".format(self.full))
        self._log = log

    def realign_germline(
        self,
        antibody,
        query_start: Optional[int] = None,
        query_end: Optional[int] = None,
    ) -> None:
        """
        Due to restrictions on the available scoring parameters in BLASTn, incorrect truncation
        of the v-gene alignment can occur. This function re-aligns the query sequence with
        the identified germline variable gene using more appropriate alignment parameters.

        Parameters
        ----------
        antibody : Antibody
            The ``Antibody`` object to which this ``GermlineSegment`` object is attached.

        query_start: int, default=None
            Position in the input sequence at which the re-alignment should start. If not provided,
            the start of the input sequence is used.

        query_end: int, default=None
            Position in the input sequence at which the re-alignment should end. If not provided,
            the end of the input sequence is used.

        """
        oriented_input = antibody.oriented_input
        germline_seq = self._get_germline_sequence_for_realignment(antibody)
        if germline_seq is None:
            antibody.log("GET GERMLINE SEQUENCE ERROR")
            antibody.log("RAW ASSIGNMENT:", self.raw_assignment)
            antibody.log("GERMLINE GENE:", self.full)
        aln_params = self._realignment_scoring_params(self.gene_type)
        # if the alignment start/end positions have been annotated by the assigner,
        # force re-alignment using those parameters
        if all(
            [
                x is not None
                for x in [
                    self.query_start,
                    self.query_end,
                    self.germline_start,
                    self.germline_end,
                ]
            ]
        ):
            query = oriented_input.sequence[self.query_start : self.query_end]
            germline = germline_seq[self.germline_start : self.germline_end]
            alignment = global_alignment(query, germline, **aln_params)
        # use local alignment to determine alignment start/end positions if
        # they haven't already been determined by the assigner
        else:
            query = oriented_input.sequence[query_start:query_end]
            if len(query) == 0:
                antibody.log("GERMLINE REALIGNMENT ERROR: query sequence is empty")
                antibody.log("REALIGNMENT QUERY SEQUENCE:", query)
                antibody.log("QUERY START:", query_start)
                antibody.log("QUERY END:", query_end)
                return
            alignment = local_alignment(query, germline_seq, **aln_params)
            # fix for a fairly rare edge case where coincidental matching to 2-3 residues at the extreme
            # 3' end of K/L germline V genes can result in incorrect identification of the
            # end of the V gene region (making the called V region far too long and, in some cases,
            # extending beyond the junction and into FR4). What we do here is drop the last 2 nucleotides
            # of the aligned germline and re-align to see whether that substantialy truncates the resulting
            # alignment (by at least 2 additional nucleotides). If so, we use the new alignment instead.
            if self.gene_type == "V" and self.chain in ["kappa", "lambda"]:
                germline_trunc = germline_seq[: alignment.target_end - 2]
                alignment2 = local_alignment(query, germline_trunc, **aln_params)
                if alignment.query_end - alignment2.query_end >= 4:
                    alignment2.raw_target = (
                        alignment.raw_target
                    )  # swap the truncated germline with the full one
                    alignment = alignment2
            # occasionally, a small number of mismatches can cause a sequence which encodes a
            # complete germline gene region to be incorrectly truncated at the 5' end (V genes) or 3' end (J genes)
            # this can cause a problem when synthesizing mAbs for recombinant expression, as the annotated VDJ
            # sequence is incomplete
            #
            # here we check to see if the alignment was truncated and whether there are additional bases in the
            # input sequence that should be included in the alignment
            # if so, we force the alignment to extend to the entire germline gene region
            if self.gene_type == "V":
                if alignment.target_begin > 0:
                    # the input sequence should contain at least enough truncated 5' residues to reach the
                    # start of the germline V gene segment
                    if alignment.query_begin >= alignment.target_begin:
                        antibody.log("\nFORCING FULL-LENGTH V-GENE REALIGNMENT")
                        antibody.log("ORIGINAL RAW QUERY:", alignment.raw_query)
                        antibody.log("ORIGINAL RAW GERMLINE:", alignment.raw_target)
                        antibody.log(
                            "ORIGINAL ALIGNED QUERY: ", alignment.aligned_query
                        )
                        antibody.log(
                            "ORIGINAL ALIGNED GERMLINE:", alignment.aligned_target
                        )
                        antibody.log("ORIGINAL QUERY START: ", alignment.query_begin)
                        antibody.log("ORIGINAL QUERY END: ", alignment.query_end)
                        antibody.log("ORIGINAL GERMLINE START:", alignment.target_begin)
                        antibody.log("ORIGINAL GERMLINE END: ", alignment.target_end)
                        antibody.log(
                            "NUMBER OF TRUNCATED RESIDUES:", alignment.target_begin
                        )

                        trunc_length = alignment.target_begin
                        # get the truncated portion of the query sequence
                        # and prepend to the aligned query sequence
                        query_truncation = alignment.raw_query[
                            alignment.query_begin - trunc_length : alignment.query_begin
                        ]
                        new_aligned_query = query_truncation + alignment.aligned_query
                        alignment.aligned_query = new_aligned_query
                        # update the query start position
                        q_begin = alignment.query_begin - trunc_length
                        alignment.query_begin = q_begin
                        # get the truncated portion of the target (germline) sequence
                        # and prepend to the aligned target sequence
                        target_truncation = alignment.raw_target[
                            alignment.target_begin
                            - trunc_length : alignment.target_begin
                        ]
                        new_aligned_target = (
                            target_truncation + alignment.aligned_target
                        )
                        alignment.aligned_target = new_aligned_target
                        # update the target (germline) start position
                        t_begin = 0
                        alignment.target_begin = t_begin

                        antibody.log("NEW ALIGNED QUERY: ", alignment.aligned_query)
                        antibody.log("NEW ALIGNED GERMLINE:", alignment.aligned_target)

            if self.gene_type == "J":
                if alignment.target_end + 1 < len(alignment.raw_target):
                    # the input sequence should contain at least enough truncated 3' residues to reach the
                    # end of the germline J gene segment
                    query_truncation_length = len(alignment.raw_query) - (
                        alignment.query_end + 1
                    )
                    target_truncation_length = len(alignment.raw_target) - (
                        alignment.target_end + 1
                    )
                    if query_truncation_length >= target_truncation_length:
                        antibody.log("\nFORCING FULL-LENGTH V-GENE REALIGNMENT")
                        antibody.log("ORIGINAL RAW QUERY:", alignment.raw_query)
                        antibody.log("ORIGINAL RAW GERMLINE:", alignment.raw_target)
                        antibody.log(
                            "ORIGINAL ALIGNED QUERY: ", alignment.aligned_query
                        )
                        antibody.log(
                            "ORIGINAL ALIGNED GERMLINE:", alignment.aligned_target
                        )
                        antibody.log("ORIGINAL QUERY START: ", alignment.query_begin)
                        antibody.log("ORIGINAL QUERY END: ", alignment.query_end)
                        antibody.log("ORIGINAL GERMLINE START:", alignment.target_begin)
                        antibody.log("ORIGINAL GERMLINE END: ", alignment.target_end)
                        antibody.log(
                            "NUMBER OF TRUNCATED RESIDUES:", target_truncation_length
                        )

                        # get the truncated portion of the query sequence
                        # and append to the aligned query sequence
                        query_truncation = alignment.raw_query[
                            alignment.query_end
                            + 1 : alignment.query_end
                            + target_truncation_length
                            + 1
                        ]
                        new_aligned_query = alignment.aligned_query + query_truncation
                        alignment.aligned_query = new_aligned_query
                        # update the query end position
                        q_end = alignment.query_end + target_truncation_length
                        alignment.query_end = q_end
                        # get the truncated portion of the target (germline) sequence
                        # and append to the aligned target sequence
                        target_truncation = alignment.raw_target[
                            alignment.target_end
                            + 1 : alignment.target_end
                            + target_truncation_length
                            + 1
                        ]
                        new_aligned_target = (
                            alignment.aligned_target + target_truncation
                        )
                        alignment.aligned_target = new_aligned_target
                        # update the target (germline) end position
                        t_end = alignment.target_end + target_truncation_length
                        alignment.target_end = t_end

                        antibody.log("NEW ALIGNED QUERY: ", alignment.aligned_query)
                        antibody.log("NEW ALIGNED GERMLINE:", alignment.aligned_target)
        if alignment:
            self._process_realignment(antibody, alignment, query_start)
        else:
            antibody.log("GERMLINE REALIGNMENT ERROR")
            antibody.log("REALIGNMENT QUERY SEQUENCE:", query)
            antibody.log("QUERY START:", query_start)
            antibody.log("QUERY END:", query_end)

    def gapped_imgt_realignment(self):
        """
        Aligns to gapped IMGT germline sequence. Used to determine
        IMGT-formatted position numberings so that identifying
        antibody regions is simplified.
        """
        self.imgt_germline = get_imgt_germlines(
            self.db_name,
            gene_type=self.gene_type,
            receptor=self.receptor,
            gene=self.full,
        )
        query = self.germline_alignment.replace("-", "")
        aln_params = self._realignment_scoring_params(self.gene_type)
        aln_params["gap_open"] = -11
        # aln_matrix = self._get_gapped_imgt_substitution_matrix()
        self.imgt_gapped_alignment = local_alignment(
            query,
            self.imgt_germline.gapped_nt_sequence,
            # matrix=aln_matrix,
            **aln_params,
        )
        self.alignment_reading_frame = (
            (2 * (self.imgt_gapped_alignment.target_begin % 3)) % 3
        ) + (
            self.imgt_germline.coding_start - 1
        )  # IMGT coding start is 1-based
        self.coding_region = self._get_coding_region()
        self.aa_sequence = self._get_aa_sequence()
        try:
            self._imgt_numbering()
        except:
            self.exception("IMGT NUMBERING", traceback.format_exc(), sep="\n")

    def _get_gapped_imgt_substitution_matrix(self):
        residues = "ACGTN."
        m = parasail.matrix_create(residues, 3, -2)
        for i in range(len(residues)):
            m[5, i] = -3
            m[i, 5] = -3
        return m

        # matrix = {}
        # residues = ["A", "C", "G", "T", "N", "."]
        # for r1 in residues:
        #     matrix[r1] = {}
        #     for r2 in residues:
        #         if r1 == r2:
        #             score = 3
        #         elif any([r1 == ".", r2 == "."]):
        #             score = -3
        #         else:
        #             score = -2
        #         matrix[r1][r2] = score
        # return matrix

    def get_imgt_position_from_raw(self, raw):
        return self._imgt_position_from_raw.get(raw, None)

    def get_raw_position_from_imgt(self, imgt):
        return self._raw_position_from_imgt.get(imgt, None)

    def _process_realignment(self, antibody, aln, query_start):
        self.realignment = aln
        self.score = aln.score
        self.raw_query = aln.raw_query
        self.raw_germline = aln.raw_target
        self.query_alignment = aln.aligned_query
        self.germline_alignment = aln.aligned_target
        self.alignment_midline = "".join(
            [
                "|" if q == g else " "
                for q, g in zip(aln.aligned_query, aln.aligned_target)
            ]
        )
        # only update alignment start/end positions if not already annotated by aligner
        if any(
            [
                x is None
                for x in [
                    self.query_start,
                    self.query_end,
                    self.germline_start,
                    self.germline_end,
                ]
            ]
        ):
            offset = query_start if query_start is not None else 0
            self.query_start = aln.query_begin + offset
            self.query_end = aln.query_end + offset
            self.germline_start = aln.target_begin
            self.germline_end = aln.target_end
        self._fix_ambigs(antibody)
        self._find_indels(antibody)

    def _get_germline_sequence_for_realignment(self, antibody):
        """
        Identifies the appropriate germline variable gene from a database of all
        germline variable genes.

        Returns:
        --------

            str: Germline sequence, or ``None`` if the requested germline gene could not be found.
        """
        # mod_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        # db_file = os.path.join(mod_dir, 'assigners/germline_dbs/{}_{}.fasta'.format(self.species.lower(), self.gene_type))
        germ_dir = get_germline_database_directory(self.db_name, self.receptor)
        antibody.log("GERMLINE DIRECTORY: {}".format(germ_dir))
        db_file = os.path.join(
            germ_dir, "ungapped/{}.fasta".format(self.gene_type.lower())
        )
        antibody.log("GERMLINE DB FILE: {}".format(db_file))
        try:
            for s in SeqIO.parse(open(db_file), "fasta"):
                if s.id == self.raw_assignment:
                    return str(s.seq)
            # TODO: log that the germline gene wasn't found in the database file
            return None
        except:
            # TODO: log that the germline database file couldn't be found
            return None

    def _imgt_numbering(self):
        aln_start = self.query_start
        aln_pos = 0
        imgt_start = self.imgt_gapped_alignment.target_begin + 1
        imgt_pos = imgt_start
        raw_position_from_imgt = {}
        imgt_position_from_raw = {}

        # imgt_start_offset is for J-genes only. Since the first position of the gapped IMGT
        # V-gene is the first position of the antibody seqeunce, IMGT numbering of the
        # V-gene should start at position 1. Not true with J-genes.
        # When processing V-genes, imgt_start_offset will always be 0.
        imgt_start_offset = self._get_imgt_start_offset()

        # Because we're iterating over the alignment but also want to track the raw (oriented_input)
        # query position, we need to adjust the alignment numbering in case there's a deletion
        # in the query sequence. With a deletion, the alignment position increases, but the raw
        # position shouldn't. query_del_adjustment allows us to adjust the alignment position so that
        # the raw query position is accurately tracked.
        query_del_adjustment = 0

        for gl in self.imgt_germline.gapped_nt_sequence[imgt_start:]:
            # only compute IMGT AA position numbers for V-genes. J-gene numbering is dependent on the CDR3 length
            # (due to IMGT's CDR3 naming conventions), so the junction must be annotated before J-gene AA position
            # numbers can be assigned.
            # if self.gene_type == 'V':

            # check to see if the position we're looking at is the start of a codon
            is_codon_start = (imgt_pos + imgt_start_offset) % 3 == 1

            # start by adding the aa position number (only if we're at the start of a codon)
            if is_codon_start:
                codon = self.germline_alignment[aln_pos : aln_pos + 3]
                if codon.count("-") >= 2:
                    self.imgt_aa_positions.append(None)
                else:
                    self.imgt_aa_positions.append(
                        (imgt_pos + imgt_start_offset + 2) / 3
                    )

            # If the gapped IMGT germline is '.' (indicating a gap introduced by IMGT for numbering purposes),
            # we only need to increment the IMGT position and indicate the lack of sequence at the IMGT position.
            if gl == ".":
                raw_position_from_imgt[imgt_pos + imgt_start_offset] = None
                imgt_pos += 1
                continue

            # If there's a gap in the query alignment (deletion in the query sequence)
            # there's no equivalent IMGT position in the query.
            if self.query_alignment[aln_pos] == "-":
                raw_position_from_imgt[imgt_pos + imgt_start_offset] = None
                aln_pos += 1
                imgt_pos += 1
                query_del_adjustment += 1
                continue

            # If there's a gap in the germline alignment (insertion in the query sequence)
            # we need to increment the alignment position until the insertion is finished.
            # Since there isn't a germline IMGT position that's equivalent to the query position,
            # we don't want to increment the IMGT position or continue iterating through the IMGT
            # germline sequence.
            if self.germline_alignment[aln_pos] == "-":
                self.log("INFO: Found an insertion in the query sequence!")
                while self.germline_alignment[aln_pos] == "-":
                    imgt_position_from_raw[
                        aln_pos + self.query_start - query_del_adjustment
                    ] = None
                    self.imgt_nt_positions.append(None)
                    aln_pos += 1

            # if the gapped IMGT germline isn't '.' and there's not an insertion in the query
            # sequence (or we've already iterated past it), we record both the IMGT position
            # and the raw (oriented_input) position and increment both the aligned and IMGt positions.
            raw_position_from_imgt[imgt_pos + imgt_start_offset] = (
                aln_pos + aln_start - query_del_adjustment
            )
            imgt_position_from_raw[aln_pos + aln_start - query_del_adjustment] = (
                imgt_pos + imgt_start_offset
            )
            self.imgt_nt_positions.append(imgt_pos + imgt_start_offset)
            aln_pos += 1
            imgt_pos += 1
            if aln_pos >= len(self.germline_alignment):
                break
        self._raw_position_from_imgt = raw_position_from_imgt
        self._imgt_position_from_raw = imgt_position_from_raw
        if self.insertions or self.deletions:
            self._calculate_imgt_indel_positions()

    def _get_imgt_start_offset(self):
        if self.gene_type == "V":
            return 0
        # find the start of FR4 in the IMGT gapped germline gene
        end_res = "W" if self.chain == "heavy" else "F"
        for i, res in enumerate(self.imgt_germline.ungapped_aa_sequence):
            if (
                res == end_res
                and end_res not in self.imgt_germline.ungapped_aa_sequence[i + 1 :]
            ):
                nts_from_start_to_fr4 = (self.imgt_germline.coding_start) + (i * 3)
                break
        # the IMGT start offset is the conserved FR4 start position (352)
        # minus the number of nts from the start of the germline gene to FR4
        return 352 - nts_from_start_to_fr4

    def _get_coding_region(self):
        coding_region = self.query_alignment[self.alignment_reading_frame :].replace(
            "-", ""
        )
        truncation = len(coding_region) % 3
        if truncation > 0:
            coding_region = coding_region[:-truncation]
        return coding_region

    def _get_aa_sequence(self):
        return Seq(self.coding_region).translate()

    def _fix_ambigs(self, antibody):
        """
        Fixes ambiguous nucleotides by replacing them with the germline nucleotide.
        """
        self.query_alignment = "".join(
            [
                q if q.upper() != "N" else g
                for q, g in zip(self.query_alignment, self.germline_alignment)
            ]
        )
        # don't forget to also correct ambigs in the oriented_input sequence
        antibody.oriented_input.sequence = (
            antibody.oriented_input.sequence[: self.query_start]
            + self.query_alignment.replace("-", "")
            + antibody.oriented_input.sequence[self.query_end + 1 :]
        )

    def _indel_check(self):
        if any(["-" in self.query_alignment, "-" in self.germline_alignment]):
            return True
        return False

    def _find_indels(self, antibody):
        """
        Identifies and annotates indels in the query sequence.
        """
        if self._indel_check():
            self.insertions = indels.find_insertions(antibody, self)
            if self.insertions:
                # only set self.has_insertion to 'yes' if the sequence has an in-frame insertion
                if "yes" in [ins["in frame"] for ins in self.insertions]:
                    self.has_insertion = "yes"
            self.deletions = indels.find_deletions(antibody, self)
            if self.deletions:
                # only set self.has_deletion to 'yes' if the sequence has an in-frame deletion
                if "yes" in [deletion["in frame"] for deletion in self.deletions]:
                    self.has_deletion = "yes"

    def _calculate_imgt_indel_positions(self):
        if self.insertions:
            for i in self.insertions:
                # since there's possibly not a direct IMGT correlate to the position
                # at the start of the insertion, need to find the closest IMGT correlate
                raw_pos = i.raw_position
                while self.get_imgt_position_from_raw(raw_pos) is None:
                    raw_pos -= 1
                i.imgt_position = self.get_imgt_position_from_raw(raw_pos)
                i.imgt_codon = int(math.ceil(i.imgt_position / 3.0))
        if self.deletions:
            for d in self.deletions:
                # since there's possibly not a direct IMGT correlate to the position
                # at the start of the deletion, need to find the closest IMGT correlate
                raw_pos = d.raw_position
                while self.get_imgt_position_from_raw(raw_pos) is None:
                    raw_pos -= 1
                d.imgt_position = self.get_imgt_position_from_raw(raw_pos)
                d.imgt_codon = int(math.ceil(d.imgt_position / 3.0))

    @staticmethod
    def _realignment_scoring_params(gene):
        """
        Returns realignment scoring paramaters for a given gene type.

        Args:

            gene (str): the gene type ('V', 'D', or 'J')


        Returns:

            dict: realignment scoring parameters
        """
        scores = {
            "V": {"match": 3, "mismatch": -2, "gap_open": -22, "gap_extend": -1},
            "D": {"match": 3, "mismatch": -2, "gap_open": -22, "gap_extend": -1},
            "J": {"match": 3, "mismatch": -2, "gap_open": -22, "gap_extend": -1},
        }
        return scores[gene]


def get_germline_database_directory(species, receptor="bcr"):
    species = species.lower()
    receptor = receptor.lower()
    addon_dir = os.path.expanduser(f"~/.abstar/germline_dbs/{receptor}")
    if os.path.isdir(addon_dir):
        if species.lower() in [os.path.basename(d[0]) for d in os.walk(addon_dir)]:
            return os.path.join(addon_dir, species.lower())
    mod_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(mod_dir, f"assigners/germline_dbs/{receptor}/{species}")


def get_imgt_germlines(db_name, gene_type, receptor="bcr", gene=None):
    """
    Returns one or more IMGTGermlineGene objects that each contain a single IMGT-gapped germline gene.

    Args:
    -----

        species (str): Species for which the germline genes should be obtained.

        gene_type (str): Options are 'V', 'D', and 'J'.

        receptor (str): Options are ``'bcr'`` and ``'tcr'``.

        gene (str): Full name of a germline gene (using IMGT-style names, like IGHV1-2*02).
                    If provided, a single ``IMGTGermlineGene`` object will be returned, or None if the
                    specified gene could not be found. If not provided, a list of ``IMGTGermlineGene``
                    objects for all germline genes matching the ``species`` and ``gene_type`` will be returned.

    Returns:
    --------

        IMGTGermlineGene: a single ``IMGTGermlineGene`` object (if ``gene`` is provided) or a list of
                          ``IMGTGermlineGene`` objects. If no sequences match the criteria or if a germline
                          database for the requested species is not found, ``None`` is returned.
    """
    # mod_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    # db_file = os.path.join(mod_dir, 'assigners/germline_dbs/imgt_gapped/{}_{}_imgt-gapped.fasta'.format(species, gene_type))
    germ_dir = get_germline_database_directory(db_name, receptor.lower())
    db_file = os.path.join(germ_dir, "imgt_gapped/{}.fasta".format(gene_type.lower()))
    try:
        germs = [IMGTGermlineGene(g) for g in SeqIO.parse(open(db_file, "r"), "fasta")]
    except:
        # TODO: log that the germline database file couldn't be found

        # print('Could not locate the IMGT germline database file ({}).'.format(db_file))
        # print(traceback.format_exc())

        return None
    if gene is None:
        return germs
    try:
        return [g for g in germs if g.name == gene][0]
    except IndexError:
        # print('Could not locate the IMGT germline gene ({}).'.format(gene))
        # print(traceback.format_exc())

        return None


def get_germlines(db_name, gene_type, receptor="bcr", chain=None, gene=None):
    """
    Returns one or more IMGTGermlineGene objects that each contain a single IMGT-gapped germline gene.

    Args:
    -----

        species (str): Species for which the germline genes should be obtained.

        gene_type (str): Options are ``'V'``, ``'D'``, and ``'J'``.

        receptor (str): Options are ``'bcr'`` and ``'tcr'``.

        chain (str): Options are ``'heavy'``, ``'kappa'``, and ``'lambda'`` for BCRs and ``'alpha'``,
        ``'beta'``, ``'gamma'``, and ``'delta'`` for TCRs. If not provided, germline sequences from
        all chains are returned.

        gene (str): Full name of a germline gene (using IMGT-style names, like ``'IGHV1-2*02'``).
        If provided, a single ``IMGTGermlineGene`` object will be returned, or ``None`` if the
        specified gene could not be found. If not provided, a list of ``IMGTGermlineGene`` objects
        for all germline genes matching the ``species``, ``receptor`` and ``gene_type`` will be returned.

    Returns:
    --------

        IMGTGermlineGene: a single ``IMGTGermlineGene`` object (if ``gene`` is provided) or a list of
                          ``IMGTGermlineGene`` objects. If no sequences match the criteria or if a germline
                          database for the requested species is not found, ``None`` is returned.
    """
    germs = get_imgt_germlines(db_name, gene_type, receptor=receptor.lower(), gene=gene)
    if germs is None:
        return germs
    if chain is not None:
        c = chain[0].upper()
        germs = [g for g in germs if g.name[2] == c]
    if len(germs) == 1:
        return germs[0]
    return germs


class IMGTGermlineGene(object):
    """docstring for IMGTGermlineGene"""

    # species_lookup exists to translate IMGT names to AbStar germline DB names.
    # if the IMGT species isn't in species_lookup, AbStar will try to use the IMGT
    # species name directly.
    species_lookup = {"homo sapiens": "human"}

    def __init__(self, sequence, species=None):
        self.raw_sequence = Sequence(str(sequence.seq), id=sequence.description)
        self._species = species
        self.gapped_nt_sequence = self.raw_sequence.sequence
        self.ungapped_nt_sequence = self.gapped_nt_sequence.replace(".", "")

    @lazy_property
    def accession(self):
        return self.raw_sequence.id.split("|")[0].strip()

    @lazy_property
    def name(self):
        return self.raw_sequence.id.split("|")[1].strip()

    @property
    def species(self):
        if self._species is None:
            species = self.raw_sequence.id.split("|")[2].strip().lower()
            self._species = self.species_lookup.get(species, species)
        return self._species

    @lazy_property
    def functionality(self):
        return self.raw_sequence.id.split("|")[3].strip()

    @lazy_property
    def gene_type(self):
        return self.raw_sequence.id.split("|")[4].strip()[0].upper()

    @lazy_property
    def coding_start(self):
        # NOTE: uses 1-based indexing, so need to adjust if using for slicing
        return int(self.raw_sequence.id.split("|")[7].strip())

    @lazy_property
    def nt_length(self):
        return int(self.raw_sequence.id.split("|")[12].strip().split("+")[0])

    @lazy_property
    def gap_length(self):
        return int(
            self.raw_sequence.id.split("|")[12].strip().split("+")[1].split("=")[0]
        )

    @lazy_property
    def total_length(self):
        return int(self.raw_sequence.id.split("|")[12].strip().split("=")[1])

    @lazy_property
    def partial(self):
        p = []
        partial = self.raw_sequence.id.split("|")[13]
        if "3'" in partial:
            p.append("3'")
        if "5'" in partial:
            p.append("5'")
        return p

    @lazy_property
    def is_rev_comp(self):
        if self.raw_sequence.id.split("|")[14].strip() == "rev-compl":
            is_rev_comp = True
        else:
            is_rev_comp = False
        return is_rev_comp

    @lazy_property
    def gapped_aa_sequence(self):
        res = []
        coding = self.gapped_nt_sequence[self.coding_start - 1 :]
        for codon in (coding[pos : pos + 3] for pos in range(0, len(coding), 3)):
            if len(codon) != 3:
                continue
            if codon == "...":
                res.append(".")
            else:
                res.append(codon_lookup.get(codon, "X"))
        return "".join(res)

    @lazy_property
    def ungapped_aa_sequence(self):
        return self.gapped_aa_sequence.replace(".", "")
