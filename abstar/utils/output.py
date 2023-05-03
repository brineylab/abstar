#!/usr/bin/python
# filename: output.py

#
# Copyright (c) 2015 Bryan Briney
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

import collections
import json
import os
import traceback
import uuid
import pandas as pd

from abutils.utils import log
from abstar.utils.parquet_schema import schema

from .cigar import make_cigar


OUTPUT_SEPARATORS = {"airr": "\t", "imgt": ",", "tabular": ","}
OUTPUT_EXTENSIONS = {"airr": ".tsv", "json": ".json", "imgt": ".csv", "tabular": ".csv"}


def get_abstar_results(antibodies, pretty=False, padding=True, raw=False, keys=None):
    return [AbstarResult(ab, pretty, padding, raw, keys) for ab in antibodies]


def get_abstar_result(ab, pretty=False, padding=True, raw=False, keys=None):
    return AbstarResult(ab, pretty, padding, raw, keys)


class AbstarResult(object):
    """docstring for AbstarOutput"""

    def __init__(self, antibody, pretty, padding, raw, keys):
        super(AbstarResult, self).__init__()
        self.antibody = antibody
        self.pretty = pretty
        self.padding = padding
        self.raw = raw
        self.keys = keys
        # property vars
        self._json_output = None
        self._imgt_output = None
        self._tabular_output = None
        self._airr_output = None
        self._imgt_header = None
        self._tabular_header = None
        self._airr_header = None

    @property
    def json_output(self):
        if self._json_output is None:
            try:
                self._json_output = self._build_json_output()
            except:
                self.antibody.exception(
                    "JSON CREATION EXCEPTION", traceback.format_exc()
                )
                self._json_output = None
        return self._json_output

    @json_output.setter
    def json_output(self, json):
        self._json_output = json

    @property
    def imgt_output(self):
        if self._imgt_output is None:
            # TODO: build IMGT-formatted output
            pass
        return self._imgt_output

    @imgt_output.setter
    def imgt_output(self, imgt):
        self._imgt_output = imgt

    @property
    def tabular_output(self):
        if self._tabular_output is None:
            self._tabular_output = self._build_tabular_output()
        return self._tabular_output

    @tabular_output.setter
    def tabular_output(self, tabular):
        self._tabular_output = tabular

    @property
    def airr_output(self):
        if self._airr_output is None:
            self._airr_output = self._build_airr_output()
        return self._airr_output

    @airr_output.setter
    def airr_output(self, airr):
        self._airr_output = airr

    # @property
    # def parquet_output(self):
    #     if self._parquet_output is None:
    #         self._parquet_output = self._build_parquet_output()
    #     return self._parquet_output

    # @parquet_output.setter
    # def parquet_output(self, parquet):
    #     self._parquet_output = parquet

    def _build_json_output(self, raw=False):
        """
        Assembles AbAnalyze output in JSON format.

        Input is a VDJ object.

        Output is a JSON-formatted output string that should be suitable
        for writing to an output file.
        """
        d_info = {}
        mut_count_nt = (
            self.antibody.v.nt_mutations.count + self.antibody.j.nt_mutations.count
        )
        mut_count_aa = (
            self.antibody.v.aa_mutations.count + self.antibody.j.aa_mutations.count
        )
        div_muts_nt = {}
        nt_identity = {
            "v": self.antibody.v.nt_identity,
            "j": self.antibody.j.nt_identity,
        }
        assigner_scores = {
            "v": self.antibody.v.assigner_score,
            "j": self.antibody.j.assigner_score,
        }
        align_info = {
            "v_start": self.antibody.v.germline_start,
            "v_end": self.antibody.v.germline_end,
            "j_start": self.antibody.j.germline_start,
            "j_end": self.antibody.j.germline_end,
        }
        germ_alignments_nt = {
            "var": {
                "query": self.antibody.v.query_alignment,
                "germ": self.antibody.v.germline_alignment,
                "midline": self.antibody.v.alignment_midline,
            },
            "join": {
                "query": self.antibody.j.query_alignment,
                "germ": self.antibody.j.germline_alignment,
                "midline": self.antibody.j.alignment_midline,
            },
        }
        exo_trim = {
            "var_3": len(self.antibody.v.raw_germline)
            - (self.antibody.v.germline_end + 1),
            "join_5": self.antibody.j.germline_start,
        }

        try:
            isotype = self.antibody.isotype.isotype
            isotype_score = self.antibody.isotype.score
            isotype_alignment = {
                "query": self.antibody.isotype.alignment.aligned_query,
                "midline": self.antibody.isotype.alignment.alignment_midline,
                "isotype": self.antibody.isotype.alignment.aligned_target,
            }
        except AttributeError:
            isotype = "unknown"
            isotype_score = ""
            isotype_alignment = {}

        if self.antibody.d:
            d_info = {
                "full": self.antibody.d.full,
                "fam": self.antibody.d.family,
                "gene": self.antibody.d.gene,
                "score": self.antibody.d.score,
                "assigner_score": self.antibody.d.assigner_score,
                # 'frame': self.antibody.d.reading_frame,
                "others": [
                    {"full": o.full, "assigner_score": o.assigner_score}
                    for o in self.antibody.d.others
                ],
            }
            assigner_scores["d"] = self.antibody.d.assigner_score
            # mut_count_nt += self.antibody.d.nt_mutations.count
            # mut_count_aa += self.antibody.d.aa_mutations.count
            # div_muts_nt = {'num': self.antibody.d.nt_mutations.count,
            #                'muts': [{'loc': m['pos'],
            #                          'mut': m['mut']} for m in self.antibody.d.nt_mutations.all_mutations]}
            # nt_identity['d'] = self.antibody.d.nt_identity
            align_info["d_start"] = self.antibody.d.germline_start
            align_info["d_end"] = self.antibody.d.germline_end
            germ_alignments_nt["div"] = {
                "query": self.antibody.d.query_alignment,
                "germ": self.antibody.d.germline_alignment,
                "midline": self.antibody.d.alignment_midline,
            }
            junc = {
                "v_nt": self.antibody.junction.v_nt,
                "n1_nt": self.antibody.junction.n1_nt,
                "d_nt": self.antibody.junction.d_nt,
                "n2_nt": self.antibody.junction.n2_nt,
                "j_nt": self.antibody.junction.j_nt,
                #                     'd_cdr3_pos': {'start': self.antibody.junction.d_start_position_nt,
                #                                    'end': self.antibody.junction.d_end_position_nt},
                "d_dist_from_cdr3_start": self.antibody.junction.d_dist_from_cdr3_start_nt,
                "d_dist_from_cdr3_end": self.antibody.junction.d_dist_from_cdr3_end_nt,
            }
            # germ_junc = {'n1_nt': self.antibody.germ_junction.n1_nt,
            #              'n2_nt': self.antibody.germ_junction.n2_nt,
            #              'd_nt': self.antibody.germ_junction.d_nt,
            #              'd_cdr3_pos': {'start': self.antibody.germ_junction.d_start_position_nt,
            #                             'end': self.antibody.germ_junction.d_end_position_nt},
            #              'd_dist_from_cdr3_start': self.antibody.germ_junction.d_dist_from_cdr3_start_nt,
            #              'd_dist_from_cdr3_end': self.antibody.germ_junction.d_dist_from_cdr3_end_nt}
            exo_trim["div_5"] = self.antibody.d.germline_start
            exo_trim["div_3"] = len(self.antibody.d.raw_germline) - (
                self.antibody.d.germline_end + 1
            )

        else:
            junc = {
                "v_nt": self.antibody.junction.v_nt,
                "n_nt": self.antibody.junction.n_nt,
                "j_nt": self.antibody.junction.j_nt,
            }

        output = collections.OrderedDict(
            [
                ("seq_id", self.antibody.id),
                ("uid", self.antibody.uid),
                ("uaid", self.antibody.uid),
                ("chain", self.antibody.chain),
                (
                    "v_gene",
                    {
                        "full": self.antibody.v.full,
                        "fam": self.antibody.v.family,
                        "gene": self.antibody.v.gene,
                        "score": self.antibody.v.score,
                        "assigner_score": self.antibody.v.assigner_score,
                        "others": [
                            {"full": o.full, "assigner_score": o.assigner_score}
                            for o in self.antibody.v.others
                        ],
                    },
                ),
                ("d_gene", d_info),
                (
                    "j_gene",
                    {
                        "full": self.antibody.j.full,
                        "gene": self.antibody.j.gene,
                        "score": self.antibody.j.score,
                        "assigner_score": self.antibody.j.assigner_score,
                        "others": [
                            {"full": o.full, "assigner_score": o.assigner_score}
                            for o in self.antibody.j.others
                        ]
                        #                                    for germ, score in zip(self.antibody.j.all_germlines[1:], self.antibody.j.all_scores[1:])]
                    },
                ),
                ("assigner_scores", assigner_scores),
                ("vdj_assigner", self.antibody.v.assigner),
                ("isotype", isotype),
                ("isotype_score", isotype_score),
                ("isotype_alignment", isotype_alignment),
                ("nt_identity", nt_identity),
                (
                    "aa_identity",
                    {
                        "v": self.antibody.v.aa_identity,
                        "j": self.antibody.j.aa_identity,
                    },
                ),
                ("junc_len", len(self.antibody.junction.junction_aa)),
                ("cdr3_len", len(self.antibody.junction.cdr3_aa)),
                ("vdj_nt", self.antibody.vdj_nt),
                ("gapped_vdj_nt", self.antibody.gapped_vdj_nt),
                ("fr1_nt", self.antibody.v.regions.nt_seqs["FR1"]),
                ("cdr1_nt", self.antibody.v.regions.nt_seqs["CDR1"]),
                ("fr2_nt", self.antibody.v.regions.nt_seqs["FR2"]),
                ("cdr2_nt", self.antibody.v.regions.nt_seqs["CDR2"]),
                ("fr3_nt", self.antibody.v.regions.nt_seqs["FR3"]),
                ("cdr3_nt", self.antibody.junction.cdr3_nt),
                ("fr4_nt", self.antibody.j.regions.nt_seqs["FR4"]),
                ("vdj_germ_nt", self.antibody.vdj_germ_nt),
                ("gapped_vdj_germ_nt", self.antibody.gapped_vdj_germ_nt),
                # ('fr1_germ_nt', self.antibody.v.regions.germline_nt_seqs['FR1']),
                # ('cdr1_germ_nt', self.antibody.v.regions.germline_nt_seqs['CDR1']),
                # ('fr2_germ_nt', self.antibody.v.regions.germline_nt_seqs['FR2']),
                # ('cdr2_germ_nt', self.antibody.v.regions.germline_nt_seqs['CDR2']),
                # ('fr3_germ_nt', self.antibody.v.regions.germline_nt_seqs['FR3']),
                # ('fr4_germ_nt', self.antibody.j.regions.germline_nt_seqs['FR4']),
                ("junc_nt", self.antibody.junction.junction_nt),
                (
                    "region_len_nt",
                    {
                        "fr1": len(self.antibody.v.regions.nt_seqs["FR1"]),
                        "cdr1": len(self.antibody.v.regions.nt_seqs["CDR1"]),
                        "fr2": len(self.antibody.v.regions.nt_seqs["FR2"]),
                        "cdr2": len(self.antibody.v.regions.nt_seqs["CDR2"]),
                        "fr3": len(self.antibody.v.regions.nt_seqs["FR3"]),
                        "cdr3": len(self.antibody.junction.cdr3_nt),
                        "fr4": len(self.antibody.j.regions.nt_seqs["FR4"]),
                    },
                ),
                (
                    "var_muts_nt",
                    {
                        "num": self.antibody.v.nt_mutations.count,
                        "muts": [
                            m.json_formatted for m in self.antibody.v.nt_mutations
                        ],
                    },
                ),
                # ('div_muts_nt', div_muts_nt),
                (
                    "join_muts_nt",
                    {
                        "num": self.antibody.j.nt_mutations.count,
                        "muts": [
                            m.json_formatted for m in self.antibody.j.nt_mutations
                        ],
                    },
                ),
                ("mut_count_nt", mut_count_nt),
                ("vdj_aa", self.antibody.vdj_aa),
                ("fr1_aa", self.antibody.v.regions.aa_seqs["FR1"]),
                ("cdr1_aa", self.antibody.v.regions.aa_seqs["CDR1"]),
                ("fr2_aa", self.antibody.v.regions.aa_seqs["FR2"]),
                ("cdr2_aa", self.antibody.v.regions.aa_seqs["CDR2"]),
                ("fr3_aa", self.antibody.v.regions.aa_seqs["FR3"]),
                ("cdr3_aa", self.antibody.junction.cdr3_aa),
                ("fr4_aa", self.antibody.j.regions.aa_seqs["FR4"]),
                ("vdj_germ_aa", self.antibody.vdj_germ_aa),
                # ('fr1_germ_aa', self.antibody.v.regions.germline_aa_seqs['FR1']),
                # ('cdr1_germ_aa', self.antibody.v.regions.germline_aa_seqs['CDR1']),
                # ('fr2_germ_aa', self.antibody.v.regions.germline_aa_seqs['FR2']),
                # ('cdr2_germ_aa', self.antibody.v.regions.germline_aa_seqs['CDR2']),
                # ('fr3_germ_aa', self.antibody.v.regions.germline_aa_seqs['FR3']),
                # ('fr4_germ_aa', self.antibody.j.regions.germline_aa_seqs['FR4']),
                ("junc_aa", self.antibody.junction.junction_aa),
                (
                    "region_len_aa",
                    {
                        "fr1": len(self.antibody.v.regions.aa_seqs["FR1"]),
                        "cdr1": len(self.antibody.v.regions.aa_seqs["CDR1"]),
                        "fr2": len(self.antibody.v.regions.aa_seqs["FR2"]),
                        "cdr2": len(self.antibody.v.regions.aa_seqs["CDR2"]),
                        "fr3": len(self.antibody.v.regions.aa_seqs["FR3"]),
                        "cdr3": len(self.antibody.junction.cdr3_aa),
                        "fr4": len(self.antibody.j.regions.aa_seqs["FR4"]),
                    },
                ),
                (
                    "var_muts_aa",
                    {
                        "num": self.antibody.v.aa_mutations.count,
                        "muts": [
                            m.json_formatted for m in self.antibody.v.aa_mutations
                        ],
                    },
                ),
                (
                    "join_muts_aa",
                    {
                        "num": self.antibody.j.aa_mutations.count,
                        "muts": [
                            m.json_formatted for m in self.antibody.j.aa_mutations
                        ],
                    },
                ),
                ("v_ins", [i.json_formatted for i in self.antibody.v.insertions]),
                ("v_del", [i.json_formatted for i in self.antibody.v.deletions]),
                ("j_ins", [i.json_formatted for i in self.antibody.j.insertions]),
                ("j_del", [i.json_formatted for i in self.antibody.j.deletions]),
                (
                    "region_muts_nt",
                    {
                        "fr1": {
                            "num": len(self.antibody.nt_mutations.in_region("FR1")),
                            "muts": [
                                m.json_formatted
                                for m in self.antibody.nt_mutations.in_region("FR1")
                            ],
                        },
                        "cdr1": {
                            "num": len(self.antibody.nt_mutations.in_region("CDR1")),
                            "muts": [
                                m.json_formatted
                                for m in self.antibody.nt_mutations.in_region("CDR1")
                            ],
                        },
                        "fr2": {
                            "num": len(self.antibody.nt_mutations.in_region("FR2")),
                            "muts": [
                                m.json_formatted
                                for m in self.antibody.nt_mutations.in_region("FR2")
                            ],
                        },
                        "cdr2": {
                            "num": len(self.antibody.nt_mutations.in_region("CDR2")),
                            "muts": [
                                m.json_formatted
                                for m in self.antibody.nt_mutations.in_region("CDR2")
                            ],
                        },
                        "fr3": {
                            "num": len(self.antibody.nt_mutations.in_region("FR3")),
                            "muts": [
                                m.json_formatted
                                for m in self.antibody.nt_mutations.in_region("FR3")
                            ],
                        },
                        "fr4": {
                            "num": len(self.antibody.nt_mutations.in_region("FR4")),
                            "muts": [
                                m.json_formatted
                                for m in self.antibody.nt_mutations.in_region("FR4")
                            ],
                        },
                    },
                ),
                (
                    "region_muts_aa",
                    {
                        "fr1": {
                            "num": len(self.antibody.aa_mutations.in_region("FR1")),
                            "muts": [
                                m.json_formatted
                                for m in self.antibody.aa_mutations.in_region("FR1")
                            ],
                        },
                        "cdr1": {
                            "num": len(self.antibody.aa_mutations.in_region("CDR1")),
                            "muts": [
                                m.json_formatted
                                for m in self.antibody.aa_mutations.in_region("CDR1")
                            ],
                        },
                        "fr2": {
                            "num": len(self.antibody.aa_mutations.in_region("FR2")),
                            "muts": [
                                m.json_formatted
                                for m in self.antibody.aa_mutations.in_region("FR2")
                            ],
                        },
                        "cdr2": {
                            "num": len(self.antibody.aa_mutations.in_region("CDR2")),
                            "muts": [
                                m.json_formatted
                                for m in self.antibody.aa_mutations.in_region("CDR2")
                            ],
                        },
                        "fr3": {
                            "num": len(self.antibody.aa_mutations.in_region("FR3")),
                            "muts": [
                                m.json_formatted
                                for m in self.antibody.aa_mutations.in_region("FR3")
                            ],
                        },
                        "fr4": {
                            "num": len(self.antibody.aa_mutations.in_region("FR4")),
                            "muts": [
                                m.json_formatted
                                for m in self.antibody.aa_mutations.in_region("FR4")
                            ],
                        },
                    },
                ),
                ("prod", "yes" if self.antibody.productivity.is_productive else "no"),
                (
                    "productivity_issues",
                    ", ".join(self.antibody.productivity.productivity_issues),
                ),
                (
                    "junction_in_frame",
                    "yes" if self.antibody.junction.in_frame else "no",
                ),
                ("raw_input", self.antibody.raw_input.sequence),
                ("oriented_input", self.antibody.oriented_input.sequence),
                ("strand", self.antibody.strand),
                ("germ_alignments_nt", germ_alignments_nt),
                ("exo_trimming", exo_trim),
                ("junc_nt_breakdown", junc),
                ("germline_database", self.antibody.germ_db),
                ("vdj_assigner", self.antibody.v.assigner),
                ("species", self.antibody.species),
                (
                    "align_info",
                    align_info,
                ),  # TODO!!  Add things like V/D/J start and end positions, etc.
            ]
        )

        if self.keys is not None:
            output = collections.OrderedDict(
                [(k, v) for k, v in output.items() if k in self.keys]
            )

        if self.padding:
            output["padding"] = ["n" * 100] * 10

        # remove empty entries to save MongoDB space
        for i in list(output.keys()):
            if output[i] == "":
                del output[i]
            elif output[i] == []:
                del output[i]
            elif output[i] == {}:
                del output[i]
            elif output[i] == None:
                del output[i]
        if self.raw:
            return output
        if self.pretty:
            return json.dumps(output, indent=4)
        else:
            return json.dumps(output)

    def _build_tabular_output(self):
        # pre-fill the isotype info
        try:
            isotype = self.antibody.isotype.isotype
        except AttributeError:
            isotype = "unknown"
        # pre-fill the D-gene info
        if self.antibody.d:
            d_full = self.antibody.d.full
            d_gene = self.antibody.d.gene
        else:
            d_full = "-"
            d_gene = "-"
        output = collections.OrderedDict(
            [
                ("seq_id", self.antibody.id),
                ("uid", self.antibody.uid),
                ("chain", self.antibody.chain),
                (
                    "productive",
                    "yes" if self.antibody.productivity.is_productive else "no",
                ),
                ("v_full", self.antibody.v.full),
                ("v_gene", self.antibody.v.gene),
                # ('v_allele', vdj.v.top_germline.split('*')[-1]),
                ("d_full", d_full),
                ("d_gene", d_gene),
                ("j_full", self.antibody.j.full),
                ("j_gene", self.antibody.j.gene),
                # ('j_allele', vdj.j.top_germline.split('*')[-1]),
                ("junction_aa", self.antibody.junction.junction_aa),
                ("junction_nt", self.antibody.junction.junction_nt),
                ("cdr3_length", str(len(self.antibody.junction.cdr3_aa))),
                ("fr1_aa", self.antibody.v.regions.aa_seqs["FR1"]),
                ("fr2_aa", self.antibody.v.regions.aa_seqs["FR2"]),
                ("fr3_aa", self.antibody.v.regions.aa_seqs["FR3"]),
                ("fr4_aa", self.antibody.j.regions.aa_seqs["FR4"]),
                ("cdr1_aa", self.antibody.v.regions.aa_seqs["CDR1"]),
                ("cdr2_aa", self.antibody.v.regions.aa_seqs["CDR2"]),
                ("cdr3_nt", self.antibody.junction.cdr3_nt),
                ("cdr3_aa", self.antibody.junction.cdr3_aa),
                ("v_start", str(self.antibody.v.germline_start)),
                ("vdj_nt", self.antibody.vdj_nt),
                ("vdj_aa", self.antibody.vdj_aa),
                (
                    "var_muts_nt",
                    "|".join(
                        [m.abstar_formatted for m in self.antibody.v.nt_mutations]
                    ),
                ),
                (
                    "var_muts_aa",
                    "|".join(
                        [m.abstar_formatted for m in self.antibody.v.aa_mutations]
                    ),
                ),
                ("var_identity_nt", str(self.antibody.v.nt_identity)),
                ("var_identity_aa", str(self.antibody.v.aa_identity)),
                ("var_mut_count_nt", str(self.antibody.v.nt_mutations.count)),
                ("var_mut_count_aa", str(self.antibody.v.aa_mutations.count)),
                (
                    "var_ins",
                    "|".join([i.abstar_formatted for i in self.antibody.v.insertions]),
                ),
                (
                    "var_del",
                    "|".join([d.abstar_formatted for d in self.antibody.v.deletions]),
                ),
                ("isotype", isotype),
                ("species", self.antibody.species),
                ("raw_input", self.antibody.raw_input.sequence),
            ]
        )
        return ",".join(output.values())

    def _build_airr_output(self):
        # format specification: https://docs.airr-community.org/en/latest/datarep/rearrangements.html
        try:
            c_call = self.antibody.isotype.isotype
        except AttributeError:
            c_call = ""
        if self.antibody.d:
            d_call = self.antibody.d.full
            d_gene = self.antibody.d.gene
            d_cigar = make_cigar(self.antibody.d)
            d_score = str(self.antibody.d.score)
            # d_identity = str(self.antibody.d.nt_identity / 100.)
        else:
            d_call = ""
            d_gene = ""
            d_cigar = ""
            d_score = ""
            # d_identity = ''
        output = collections.OrderedDict(
            [
                ("sequence_id", self.antibody.id),
                ("umi", self.antibody.uid),
                ("sequence", self.antibody.vdj_nt),
                ("sequence_aa", self.antibody.vdj_aa),
                ("sequence_input", self.antibody.raw_input.sequence),
                ("rev_comp", "True" if self.antibody.v.strand == "-" else "False"),
                (
                    "productive",
                    "True" if self.antibody.productivity.is_productive else "False",
                ),
                (
                    "productivity_issues",
                    "|".join(self.antibody.productivity.productivity_issues),
                ),
                ("stop_codon", "True" if "*" in self.antibody.vdj_aa else "False"),
                ("v_call", self.antibody.v.full),
                ("v_gene", self.antibody.v.gene),
                ("v_score", str(self.antibody.v.assigner_score)),
                ("v_identity", str(self.antibody.v.nt_identity / 100.0)),
                ("v_identity_aa", str(self.antibody.v.aa_identity / 100.0)),
                ("d_call", d_call),
                ("d_gene", d_gene),
                ("d_score", d_score),
                # ('d_identity', d_identity),
                ("j_call", self.antibody.j.full),
                ("j_gene", self.antibody.j.gene),
                ("j_score", str(self.antibody.j.assigner_score)),
                ("j_identity", str(self.antibody.j.nt_identity / 100.0)),
                ("j_identity_aa", str(self.antibody.j.aa_identity / 100.0)),
                ("c_call", c_call),
                ("cdr3_length", str(len(self.antibody.junction.cdr3_aa))),
                ("fwr1_aa", self.antibody.v.regions.aa_seqs["FR1"]),
                ("fwr2_aa", self.antibody.v.regions.aa_seqs["FR2"]),
                ("fwr3_aa", self.antibody.v.regions.aa_seqs["FR3"]),
                ("fwr4_aa", self.antibody.j.regions.aa_seqs["FR4"]),
                ("cdr1_aa", self.antibody.v.regions.aa_seqs["CDR1"]),
                ("cdr2_aa", self.antibody.v.regions.aa_seqs["CDR2"]),
                ("cdr3_aa", self.antibody.junction.cdr3_aa),
                ("sequence_alignment", self.antibody.gapped_vdj_nt),
                ("germline_alignment", self.antibody.gapped_vdj_germ_nt),
                (
                    "v_mutations",
                    "|".join(
                        [m.abstar_formatted for m in self.antibody.v.nt_mutations]
                    ),
                ),
                (
                    "v_mutations_aa",
                    "|".join(
                        [m.abstar_formatted for m in self.antibody.v.aa_mutations]
                    ),
                ),
                ("v_mutation_count", str(self.antibody.v.nt_mutations.count)),
                ("v_mutation_count_aa", str(self.antibody.v.aa_mutations.count)),
                (
                    "v_insertions",
                    "|".join([i.abstar_formatted for i in self.antibody.v.insertions])
                    if self.antibody.v.insertions
                    else "",
                ),
                (
                    "v_deletions",
                    "|".join([i.abstar_formatted for i in self.antibody.v.deletions])
                    if self.antibody.v.deletions
                    else "",
                ),
                ("junction", self.antibody.junction.junction_nt),
                ("junction_aa", self.antibody.junction.junction_aa),
                (
                    "junction_in_frame",
                    "True" if self.antibody.junction.in_frame else "False",
                ),
                ("junction_v", self.antibody.junction.v_nt),
                ("junction_d", self.antibody.junction.d_nt),
                ("junction_j", self.antibody.junction.j_nt),
                ("junction_n", self.antibody.junction.n_nt),
                ("junction_n1", self.antibody.junction.n1_nt),
                ("junction_n2", self.antibody.junction.n2_nt),
                ("junction_v_aa", self.antibody.junction.v_aa),
                ("junction_d_aa", self.antibody.junction.d_aa),
                ("junction_j_aa", self.antibody.junction.j_aa),
                ("junction_n_aa", self.antibody.junction.n_aa),
                ("junction_n1_aa", self.antibody.junction.n1_aa),
                ("junction_n2_aa", self.antibody.junction.n2_aa),
                ("junction_germ_v", self.antibody.junction.v_germ_nt),
                ("junction_germ_d", self.antibody.junction.d_germ_nt),
                ("junction_germ_j", self.antibody.junction.j_germ_nt),
                ("junction_germ_v_aa", self.antibody.junction.v_germ_aa),
                ("junction_germ_d_aa", self.antibody.junction.d_germ_aa),
                ("junction_germ_j_aa", self.antibody.junction.j_germ_aa),
                ("isotype", c_call),
                ("locus", self.antibody.v.gene[:3].upper()),
                ("v_cigar", make_cigar(self.antibody.v)),
                ("d_cigar", d_cigar),
                ("j_cigar", make_cigar(self.antibody.j)),
                ("species", self.antibody.species),
                ("germline_database", self.antibody.germ_db),
            ]
        )
        return "\t".join([v if v is not None else "" for v in output.values()])


def get_output_suffix(output_format):
    return OUTPUT_EXTENSIONS[output_format.lower()]


def get_output_separator(output_format):
    return OUTPUT_SEPARATORS[output_format.lower()]


def get_parquet_dtypes(output_format):
    if output_format.lower() == "tabular":
        dtypes = {
            "var_ins": "object",
            "var_del": "object",
            "var_muts_nt": "object",
            "var_muts_aa": "object",
        }
    elif output_format == "json":
        dtypes = {
            "seq_id": object,
            "chain": object,
            "v_gene": object,
            "d_gene": object,
            "j_gene": object,
            "assigner_scores": object,
            "vdj_assigner": object,
            "isotype": object,
            "isotype_score": int,
            "isotype_alignment": object,
            "nt_identity": object,
            "aa_identity": object,
            "junc_len": int,
            "cdr3_len": int,
            "vdj_nt": object,
            "gapped_vdj_nt": object,
            "fr1_nt": object,
            "cdr1_nt": object,
            "fr2_nt": object,
            "cdr2_nt": object,
            "fr3_nt": object,
            "cdr3_nt": object,
            "fr4_nt": object,
            "vdj_germ_nt": object,
            "gapped_vdj_germ_nt": object,
            "junc_nt": object,
            "region_len_nt": object,
            "var_muts_nt": object,
            "join_muts_nt": object,
            "mut_count_nt": int,
            "vdj_aa": object,
            "fr1_aa": object,
            "cdr1_aa": object,
            "fr2_aa": object,
            "cdr2_aa": object,
            "fr3_aa": object,
            "cdr3_aa": object,
            "fr4_aa": object,
            "vdj_germ_aa": object,
            "junc_aa": object,
            "region_len_aa": object,
            "var_muts_aa": object,
            "join_muts_aa": object,
            "region_muts_nt": object,
            "region_muts_aa": object,
            "prod": object,
            "productivity_issues": object,
            "junction_in_frame": object,
            "raw_input": object,
            "oriented_input": object,
            "strand": object,
            "germ_alignments_nt": object,
            "exo_trimming": object,
            "junc_nt_breakdown": object,
            "germline_database": object,
            "species": object,
            "align_info": object,
            "j_del": object,
            "v_del": object,
            "j_ins": object,
            "v_ins": object,
        }
    else:
        dtypes = {}
    return dtypes


def get_header(output_type):
    if output_type == "tabular":
        return ",".join(TABULAR_HEADER)
    if output_type == "imgt":
        return ",".join(IMGT_HEADER)
    if output_type == "airr":
        return "\t".join(AIRR_HEADER)
    return None


def get_output(result, output_type):
    if output_type.lower() == "json":
        return result.json_output
    elif output_type.lower() == "imgt":
        return result.imgt_output
    elif output_type.lower() == "tabular":
        return result.tabular_output
    elif output_type.lower() == "airr":
        return result.airr_output
    else:
        return result.json_output


# def write_output(outputs, outfiles):
#     for _outputs, outfile in zip(outputs, outfiles):
#         with open(outfile, 'w') as f:
#             f.write('\n'.join(_outputs))


def write_output(output_dict, output_dir, output_prefix, write_parquet: bool):
    output_file_dict = {}
    for fmt in output_dict.keys():
        subdir = os.path.join(output_dir, fmt)

        if fmt == "json" and write_parquet:
            output_name = output_prefix + ".snappy.parquet"
            output_file = os.path.join(subdir, output_name)
            dtypes = get_parquet_dtypes(fmt)
            df = pd.DataFrame.from_records(
                [json.loads(line) for line in output_dict[fmt]]
            )
            df = df.reindex(columns=dtypes).astype(dtypes)
            df.to_parquet(
                output_file,
                index=False,
                engine="pyarrow",
                compression="snappy",
                schema=schema,
            )
        else:
            output_name = output_prefix + get_output_suffix(fmt)
            output_file = os.path.join(subdir, output_name)
            with open(output_file, "w") as f:
                f.write("\n".join(output_dict[fmt]))
                f.write("\n")

        output_file_dict[fmt] = output_file
    return output_file_dict


def build_output(vdjs, output_type, pretty, padding):
    logger = log.get_logger()
    try:
        vdjs = [vdj for vdj in vdjs if vdj.rearrangement]
        if output_type.lower() == "json":
            output = []
            for vdj in vdjs:
                try:
                    output.append(_json_output(vdj, pretty, padding))
                except AttributeError:
                    logger.debug("OUTPUT ERROR: {}".format(vdj.id))
        elif output_type.lower() == "imgt":
            header, firstvals = _imgt_summary_output(vdjs[0], header=True)
            output = [
                header,
                firstvals,
            ]
            for vdj in vdjs[1:]:
                try:
                    output.append(_imgt_summary_output(vdj))
                except AttributeError:
                    logger.debug("OUTPUT ERROR: {}".format(vdj.id))
        elif output_type.lower() == "hadoop":
            output = []
            for vdj in vdjs:
                try:
                    output.append(_hadoop_tabular_output(vdj))
                except AttributeError:
                    logger.debug("OUTPUT ERROR: {}".format(vdj.id))
        return output
    except:
        logger.debug(
            "FILE-LEVEL OUTPUT ERROR: sequences {} - {}, output_type = {}".format(
                vdjs[0].id, vdjs[-1].id, output_type
            )
        )
        logger.debug(traceback.format_exc())


def as_dict(vdjs):
    output = []
    for vdj in vdjs:
        output.append(_json_output(vdj, False, False, raw=True))
    return output


def output_func(output_type):
    outputs = {
        "json": _json_output,
        "imgt": _imgt_summary_output,
        "hadoop": _hadoop_tabular_output,
    }
    return outputs[output_type]


def _json_output(vdj, pretty, padding, raw=False):
    """
    Assembles AbAnalyze output in JSON format.

    Input is a VDJ object.

    Output is a JSON-formatted output string that should be suitable
    for writing to an output file.
    """
    d_info = {}
    mut_count_nt = vdj.v.nt_mutations.mutation_count + vdj.j.nt_mutations.mutation_count
    div_muts_nt = {}
    nt_identity = {
        "v": vdj.v.nt_mutations.germline_identity,
        "j": vdj.j.nt_mutations.germline_identity,
    }
    align_info = {
        "v_start": vdj.v.germline_start,
        "v_end": vdj.v.germline_end,
        "j_start": vdj.j.germline_start,
        "j_end": vdj.j.germline_end,
    }
    germ_alignments_nt = {
        "var": {
            "query": vdj.v.query_alignment,
            "germ": vdj.v.germline_alignment,
            "midline": vdj.v.alignment_midline,
        },
        "join": {
            "query": vdj.j.query_alignment,
            "germ": vdj.j.germline_alignment,
            "midline": vdj.j.alignment_midline,
        },
    }
    exo_trim = {
        "var_3": len(vdj.v.germline_seq) - vdj.v.germline_end,
        "join_5": vdj.j.germline_start,
    }

    try:
        isotype = vdj.isotype.isotype
        isotype_score = vdj.isotype.score
        isotype_alignment = {
            "query": vdj.isotype.alignment.aligned_query,
            "midline": vdj.isotype.alignment.alignment_midline,
            "isotype": vdj.isotype.alignment.aligned_target,
        }
    except AttributeError:
        isotype = "unknown"
        isotype_score = ""
        isotype_alignment = {}

    if vdj.d:
        d_info = {
            "full": vdj.d.top_germline,
            "fam": vdj.d.top_germline.split("-")[0] if vdj.d.top_germline else None,
            "gene": vdj.d.top_germline.split("*")[0] if vdj.d.top_germline else None,
            "score": vdj.d.top_score,
            "frame": vdj.d.reading_frame,
            "others": [
                {"full": germ, "score": score}
                for germ, score in zip(vdj.d.all_germlines[1:], vdj.d.all_scores[1:])
            ],
        }
        mut_count_nt += vdj.d.nt_mutations.mutation_count
        div_muts_nt = {
            "num": vdj.d.nt_mutations.mutation_count,
            "muts": [
                {"loc": m["pos"], "mut": m["mut"]}
                for m in vdj.d.nt_mutations.all_mutations
            ],
        }
        nt_identity["d"] = vdj.d.nt_mutations.germline_identity
        align_info["d_start"] = vdj.d.germline_start
        align_info["d_end"] = vdj.d.germline_end
        germ_alignments_nt["div"] = {
            "query": vdj.d.query_alignment,
            "germ": vdj.d.germline_alignment,
            "midline": vdj.d.alignment_midline,
        }
        junc = {
            "n1_nt": vdj.junction.n1_nt,
            "n2_nt": vdj.junction.n2_nt,
            "d_nt": vdj.junction.d_nt,
            "d_cdr3_pos": {
                "start": vdj.junction.d_start_position_nt,
                "end": vdj.junction.d_end_position_nt,
            },
            "d_dist_from_cdr3_start": vdj.junction.d_dist_from_cdr3_start_nt,
            "d_dist_from_cdr3_end": vdj.junction.d_dist_from_cdr3_end_nt,
        }
        germ_junc = {
            "n1_nt": vdj.germ_junction.n1_nt,
            "n2_nt": vdj.germ_junction.n2_nt,
            "d_nt": vdj.germ_junction.d_nt,
            "d_cdr3_pos": {
                "start": vdj.germ_junction.d_start_position_nt,
                "end": vdj.germ_junction.d_end_position_nt,
            },
            "d_dist_from_cdr3_start": vdj.germ_junction.d_dist_from_cdr3_start_nt,
            "d_dist_from_cdr3_end": vdj.germ_junction.d_dist_from_cdr3_end_nt,
        }
        exo_trim["div_5"] = vdj.d.germline_start
        exo_trim["div_3"] = vdj.d.germline_end

    else:
        junc = {"n_nt": vdj.junction.n1_nt}
        germ_junc = {"n_nt": vdj.germ_junction.n1_nt}

    output = collections.OrderedDict(
        [
            ("seq_id", vdj.id),
            ("uaid", vdj.uaid),
            ("chain", vdj.chain),
            (
                "v_gene",
                {
                    "full": vdj.v.top_germline,
                    "fam": vdj.v.top_germline.split("-")[0],
                    "gene": vdj.v.top_germline.split("*")[0],
                    "score": vdj.v.top_score,
                    "others": [
                        {
                            "full": germ,
                            # 'fam': germ.split('-')[0],
                            # 'gene': germ.split('*')[0],
                            "score": score,
                        }
                        for germ, score in zip(
                            vdj.v.all_germlines[1:], vdj.v.all_scores[1:]
                        )
                    ],
                },
            ),
            ("d_gene", d_info),
            (
                "j_gene",
                {
                    "full": vdj.j.top_germline,
                    "gene": vdj.j.top_germline.split("*")[0],
                    "score": vdj.j.top_score,
                    "others": [
                        {
                            "full": germ,
                            # 'fam': germ.split('*')[0],
                            "score": score,
                        }
                        for germ, score in zip(
                            vdj.j.all_germlines[1:], vdj.j.all_scores[1:]
                        )
                    ],
                },
            ),
            ("bitscores", {"v": vdj.v.top_bitscore, "j": vdj.j.top_bitscore}),
            ("e_values", {"v": vdj.v.top_evalue, "j": vdj.j.top_evalue}),
            ("isotype", isotype),
            ("isotype_score", isotype_score),
            ("isotype_alignment", isotype_alignment),
            ("nt_identity", nt_identity),
            (
                "aa_identity",
                {
                    "v": vdj.v.aa_mutations.germline_identity,
                    "j": vdj.j.aa_mutations.germline_identity,
                },
            ),
            ("junc_len", len(vdj.junction.junction_aa)),
            ("cdr3_len", len(vdj.junction.cdr3_aa)),
            ("vdj_nt", vdj.vdj_nt),
            ("gapped_vdj_nt", vdj.gapped_vdj_nt),
            ("fr1_nt", vdj.v.regions.nt_seqs["FR1"]),
            ("cdr1_nt", vdj.v.regions.nt_seqs["CDR1"]),
            ("fr2_nt", vdj.v.regions.nt_seqs["FR2"]),
            ("cdr2_nt", vdj.v.regions.nt_seqs["CDR2"]),
            ("fr3_nt", vdj.v.regions.nt_seqs["FR3"]),
            ("cdr3_nt", vdj.junction.cdr3_nt),
            ("fr4_nt", vdj.j.regions.nt_seqs["FR4"]),
            ("vdj_germ_nt", vdj.vdj_germ_nt),
            ("gapped_vdj_germ_nt", vdj.gapped_vdj_germ_nt),
            ("fr1_germ_nt", vdj.v.regions.germline_nt_seqs["FR1"]),
            ("cdr1_germ_nt", vdj.v.regions.germline_nt_seqs["CDR1"]),
            ("fr2_germ_nt", vdj.v.regions.germline_nt_seqs["FR2"]),
            ("cdr2_germ_nt", vdj.v.regions.germline_nt_seqs["CDR2"]),
            ("fr3_germ_nt", vdj.v.regions.germline_nt_seqs["FR3"]),
            ("fr4_germ_nt", vdj.j.regions.germline_nt_seqs["FR4"]),
            ("junc_nt", vdj.junction.junction_nt),
            (
                "region_len_nt",
                {
                    "fr1": int(vdj.v.regions.nt_lengths["FR1"]),
                    "cdr1": int(vdj.v.regions.nt_lengths["CDR1"]),
                    "fr2": int(vdj.v.regions.nt_lengths["FR2"]),
                    "cdr2": int(vdj.v.regions.nt_lengths["CDR2"]),
                    "fr3": int(vdj.v.regions.nt_lengths["FR3"]),
                    "cdr3": len(vdj.junction.cdr3_nt),
                    "fr4": int(vdj.j.regions.nt_lengths["FR4"]),
                },
            ),
            (
                "var_muts_nt",
                {
                    "num": vdj.v.nt_mutations.mutation_count,
                    "muts": [
                        {"loc": m["pos"], "mut": m["mut"]}
                        for m in vdj.v.nt_mutations.all_mutations
                    ],
                },
            ),
            ("div_muts_nt", div_muts_nt),
            (
                "join_muts_nt",
                {
                    "num": vdj.j.nt_mutations.mutation_count,
                    "muts": [
                        {"loc": m["pos"], "mut": m["mut"]}
                        for m in vdj.j.nt_mutations.all_mutations
                    ],
                },
            ),
            ("mut_count_nt", mut_count_nt),
            ("vdj_aa", vdj.vdj_aa),
            ("fr1_aa", vdj.v.regions.aa_seqs["FR1"]),
            ("cdr1_aa", vdj.v.regions.aa_seqs["CDR1"]),
            ("fr2_aa", vdj.v.regions.aa_seqs["FR2"]),
            ("cdr2_aa", vdj.v.regions.aa_seqs["CDR2"]),
            ("fr3_aa", vdj.v.regions.aa_seqs["FR3"]),
            ("cdr3_aa", vdj.junction.cdr3_aa),
            ("fr4_aa", vdj.j.regions.aa_seqs["FR4"]),
            ("vdj_germ_aa", vdj.vdj_germ_aa),
            ("fr1_germ_aa", vdj.v.regions.germline_aa_seqs["FR1"]),
            ("cdr1_germ_aa", vdj.v.regions.germline_aa_seqs["CDR1"]),
            ("fr2_germ_aa", vdj.v.regions.germline_aa_seqs["FR2"]),
            ("cdr2_germ_aa", vdj.v.regions.germline_aa_seqs["CDR2"]),
            ("fr3_germ_aa", vdj.v.regions.germline_aa_seqs["FR3"]),
            ("fr4_germ_aa", vdj.j.regions.germline_aa_seqs["FR4"]),
            ("junc_aa", vdj.junction.junction_aa),
            (
                "region_len_aa",
                {
                    "fr1": int(vdj.v.regions.aa_lengths["FR1"]),
                    "cdr1": int(vdj.v.regions.aa_lengths["CDR1"]),
                    "fr2": int(vdj.v.regions.aa_lengths["FR2"]),
                    "cdr2": int(vdj.v.regions.aa_lengths["CDR2"]),
                    "fr3": int(vdj.v.regions.aa_lengths["FR3"]),
                    "cdr3": len(vdj.junction.cdr3_aa),
                    "fr4": int(vdj.j.regions.aa_lengths["FR4"]),
                },
            ),
            (
                "var_muts_aa",
                {
                    "num": vdj.v.aa_mutations.mutation_count,
                    "muts": [
                        {"loc": m["pos"], "mut": m["mut"]}
                        for m in vdj.v.aa_mutations.all_mutations
                    ],
                },
            ),
            (
                "v_ins",
                [
                    {
                        "loc": i["pos"],
                        "len": i["len"],
                        "seq": i["seq"],
                        "in_frame": i["in frame"],
                    }
                    for i in vdj.v.insertions
                ],
            ),
            (
                "v_del",
                [
                    {
                        "loc": d["pos"],
                        "len": d["len"],
                        "seq": d["seq"],
                        "in_frame": d["in frame"],
                    }
                    for d in vdj.v.deletions
                ],
            ),
            (
                "j_ins",
                [
                    {
                        "loc": i["pos"],
                        "len": i["len"],
                        "seq": i["seq"],
                        "in_frame": i["in frame"],
                    }
                    for i in vdj.j.insertions
                ],
            ),
            (
                "j_del",
                [
                    {
                        "loc": d["pos"],
                        "len": d["len"],
                        "seq": d["seq"],
                        "in_frame": d["in frame"],
                    }
                    for d in vdj.j.deletions
                ],
            ),
            (
                "region_muts_nt",
                {
                    "fr1": {
                        "num": vdj.v.nt_mutations.region_mutation_count["FR1"],
                        "muts": [
                            {"loc": m["pos"], "mut": m["mut"]}
                            for m in vdj.v.nt_mutations.region_mutations["FR1"]
                        ],
                    },
                    "cdr1": {
                        "num": vdj.v.nt_mutations.region_mutation_count["CDR1"],
                        "muts": [
                            {"loc": m["pos"], "mut": m["mut"]}
                            for m in vdj.v.nt_mutations.region_mutations["CDR1"]
                        ],
                    },
                    "fr2": {
                        "num": vdj.v.nt_mutations.region_mutation_count["FR2"],
                        "muts": [
                            {"loc": m["pos"], "mut": m["mut"]}
                            for m in vdj.v.nt_mutations.region_mutations["FR2"]
                        ],
                    },
                    "cdr2": {
                        "num": vdj.v.nt_mutations.region_mutation_count["CDR2"],
                        "muts": [
                            {"loc": m["pos"], "mut": m["mut"]}
                            for m in vdj.v.nt_mutations.region_mutations["CDR2"]
                        ],
                    },
                    "fr3": {
                        "num": vdj.v.nt_mutations.region_mutation_count["FR3"],
                        "muts": [
                            {"loc": m["pos"], "mut": m["mut"]}
                            for m in vdj.v.nt_mutations.region_mutations["FR3"]
                        ],
                    },
                    "fr4": {
                        "num": vdj.j.nt_mutations.region_mutation_count["FR4"],
                        "muts": [
                            {"loc": m["pos"], "mut": m["mut"]}
                            for m in vdj.j.nt_mutations.region_mutations["FR4"]
                        ],
                    },
                },
            ),
            (
                "region_muts_aa",
                {
                    "fr1": {
                        "num": vdj.v.aa_mutations.region_mutation_count["FR1"],
                        "muts": [
                            {"loc": m["pos"], "mut": m["mut"]}
                            for m in vdj.v.aa_mutations.region_mutations["FR1"]
                        ],
                    },
                    "cdr1": {
                        "num": vdj.v.aa_mutations.region_mutation_count["CDR1"],
                        "muts": [
                            {"loc": m["pos"], "mut": m["mut"]}
                            for m in vdj.v.aa_mutations.region_mutations["CDR1"]
                        ],
                    },
                    "fr2": {
                        "num": vdj.v.aa_mutations.region_mutation_count["FR2"],
                        "muts": [
                            {"loc": m["pos"], "mut": m["mut"]}
                            for m in vdj.v.aa_mutations.region_mutations["FR2"]
                        ],
                    },
                    "cdr2": {
                        "num": vdj.v.aa_mutations.region_mutation_count["CDR2"],
                        "muts": [
                            {"loc": m["pos"], "mut": m["mut"]}
                            for m in vdj.v.aa_mutations.region_mutations["CDR2"]
                        ],
                    },
                    "fr3": {
                        "num": vdj.v.aa_mutations.region_mutation_count["FR3"],
                        "muts": [
                            {"loc": m["pos"], "mut": m["mut"]}
                            for m in vdj.v.aa_mutations.region_mutations["FR3"]
                        ],
                    },
                    "fr4": {
                        "num": vdj.j.aa_mutations.region_mutation_count["FR4"],
                        "muts": [
                            {"loc": m["pos"], "mut": m["mut"]}
                            for m in vdj.j.aa_mutations.region_mutations["FR4"]
                        ],
                    },
                },
            ),
            ("prod", vdj.productive),
            ("junction_in_frame", "yes" if vdj.junction.in_frame else "no"),
            ("raw_input", vdj.raw_input),
            ("raw_query", vdj.raw_query),
            ("strand", vdj.strand),
            (
                "codons",
                {
                    "vdj": vdj.codons.vdj_codons,
                    "vdj_regions": vdj.codons.vdj_codon_regions,
                    "v": vdj.codons.v_codons,
                    "v_germ": vdj.codons.v_germ_codons,
                    "j": vdj.codons.j_codons,
                    "j_germ": vdj.codons.j_germ_codons,
                },
            ),
            (
                "gapped_codons",
                {
                    "vdj": vdj.gapped_codons.vdj_codons,
                    "vdj_regions": vdj.gapped_codons.vdj_codon_regions,
                    "v": vdj.gapped_codons.v_codons,
                    "v_germ": vdj.gapped_codons.v_germ_codons,
                    "j": vdj.gapped_codons.j_codons,
                    "j_germ": vdj.gapped_codons.j_germ_codons,
                },
            ),
            ("germ_alignments_nt", germ_alignments_nt),
            ("exo_trim", exo_trim),
            ("junc", junc),
            ("germ_junc", germ_junc),
            (
                "align_info",
                align_info,
            ),  # TODO!!  Add things like V/D/J start and end positions, etc.
            ("vdj_region_string", vdj.vdj_region_string),
            ("gapped_vdj_region_string", vdj.gapped_vdj_region_string),
            ("germline_database", vdj.germ_db),
            ("species", vdj.species),
        ]
    )

    if padding:
        output["padding"] = ["n" * 100] * 10

    # remove empty entries to save MongoDB space
    for i in output.keys():
        if output[i] == "":
            del output[i]
        elif output[i] == []:
            del output[i]
        elif output[i] == {}:
            del output[i]
        elif output[i] == None:
            del output[i]
    if raw:
        return output
    if pretty:
        return json.dumps(output, indent=4)
    else:
        return json.dumps(output)


def _imgt_summary_output(vdj, header=False):
    output = collections.OrderedDict(
        [
            ("Sequence number", uuid.uuid4()),
            ("Sequence ID", vdj.id),
            (
                "Functionality",
                "productive" if vdj.productive == "yes" else "non-productive",
            ),
            ("V-GENE and allele", vdj.v.top_germline),
            ("V-REGION score", vdj.v.top_score),
            ("V-REGION identity %", round(vdj.v.nt_mutations.germline_identity, 2)),
            ("V-REGION identity nt", _get_iden_nt(vdj.v)),
            (
                "V-REGION identity % (with ins/del events)",
                round(vdj.v.nt_mutations.germline_identity, 2),
            ),
            ("V-REGION identity nt (with ins/del events)", _get_iden_nt(vdj.v)),
            ("J-GENE and allele", vdj.j.top_germline),
            ("J-REGION score", vdj.j.top_score),
            ("J-REGION identity %", round(vdj.j.nt_mutations.germline_identity, 2)),
            ("J-REGION identity nt", _get_iden_nt(vdj.j)),
            ("D-GENE and allele", vdj.d.top_germline if vdj.d is not None else ""),
            ("D-REGION reading frame", vdj.d.reading_frame if vdj.d else ""),
            ("CDR1-IMGT length", vdj.v.regions.aa_lengths["CDR1"]),
            ("CDR2-IMGT length", vdj.v.regions.aa_lengths["CDR2"]),
            ("CDR3-IMGT length", len(vdj.junction.cdr3_aa)),
            (
                "CDR-IMGT lengths",
                "{}.{}.{}".format(
                    vdj.v.regions.aa_lengths["CDR1"],
                    vdj.v.regions.aa_lengths["CDR2"],
                    len(vdj.junction.cdr3_aa),
                ),
            ),
            (
                "FR-IMGT lengths",
                "[{}.{}.{}.{}]".format(
                    vdj.v.regions.aa_lengths["FR1"],
                    vdj.v.regions.aa_lengths["FR2"],
                    vdj.v.regions.aa_lengths["FR3"],
                    vdj.j.regions.aa_lengths["FR4"],
                ),
            ),
            ("AA JUNCTION", vdj.junction.junction_aa),
            (
                "JUNCTION frame",
                "in-frame"
                if vdj.junction.junction_aa[0] == "C"
                and vdj.junction.junction_aa[-1] in ["W", "F"]
                else "out-of-frame",
            ),
            ("Orientation", "+" if vdj.strand == "plus" else "-"),
            ("Functionality comment", ""),
            ("V-REGION potential ins/del", ""),
            ("J-GENE and allele comment", ""),
            ("V-REGION insertions", _get_imgt_indel_string(vdj.v, "ins")),
            ("V-REGION deletions", _get_imgt_indel_string(vdj.v, "del")),
            ("Sequence", vdj.vdj_nt),
        ]
    )
    data = "\t".join([str(v) for v in output.values()])
    if not header:
        return data
    headers = "\t".join(output.keys())
    return headers, data


def _get_iden_nt(v):
    v_len = len(v.query_alignment.replace("-", ""))
    muts = v.nt_mutations.mutation_count
    return "{}/{} nt".format(v_len - muts, v_len)


def _get_imgt_indel_string(v, indel_type, hadoop=False):
    indel_list = []
    if indel_type == "ins":
        try:
            indels = v.insertions
        except KeyError:
            return ""
    if indel_type == "del":
        try:
            indels = v.deletions
        except KeyError:
            return ""
    for i in indels:
        indel_list.append("{}|{}|{}".format(i["pos"], i["len"], i["seq"]))
    if hadoop:
        return " ".join(indel_list)
    return ", ".join(indel_list)


def _hadoop_tabular_output(vdj):
    output = collections.OrderedDict(
        [
            ("seq_uuid", uuid.uuid4()),
            ("seq_id", vdj.id),
            ("chain", vdj.chain),
            ("productive", vdj.productive),
            ("vgene_fam", vdj.v.top_germline.split("-")[0]),
            ("vgene_gene", vdj.v.top_germline.split("*")[0]),
            ("vgene_allele", vdj.v.top_germline.split("*")[-1]),
            # ('dgene_fam', '' if not vdj.d or vdj.d.top_germline == None else vdj.d.top_germline.split('-')[0]),
            # ('dgene_gene', '' if not vdj.d or vdj.d.top_germline == None else vdj.d.top_germline.split('*')[0]),
            # ('dgene_allele', '' if not vdj.d or vdj.d.top_germline == None else vdj.d.top_germline.split('*')[-1]),
            ("jgene_gene", vdj.j.top_germline.split("*")[0]),
            ("jgene_allele", vdj.j.top_germline.split("*")[-1]),
            ("junction_aa", vdj.junction.junction_aa),
            ("junction_nt", vdj.junction.junction_nt),
            ("cdr3_length", len(vdj.junction.cdr3_aa)),
            ("fr1_aa", vdj.v.regions.aa_seqs["FR1"]),
            ("fr2_aa", vdj.v.regions.aa_seqs["FR2"]),
            ("fr3_aa", vdj.v.regions.aa_seqs["FR3"]),
            ("fr4_aa", vdj.j.regions.aa_seqs["FR4"]),
            ("cdr1_aa", vdj.v.regions.aa_seqs["CDR1"]),
            ("cdr2_aa", vdj.v.regions.aa_seqs["CDR2"]),
            ("cdr3_aa", vdj.junction.cdr3_aa),
            ("vdj_nt", vdj.vdj_nt),
            ("vj_aa", vdj.vdj_aa),
            ("var_ins", _get_imgt_indel_string(vdj.v, "ins", hadoop=True)),
            ("var_del", _get_imgt_indel_string(vdj.v, "del", hadoop=True)),
        ]
    )
    return ",".join([str(v) for v in output.values()])


IMGT_HEADER = [
    "Sequence number",
    "Sequence ID",
    "Functionality",
    "V-GENE and allele",
    "V-REGION score",
    "V-REGION identity %",
    "V-REGION identity nt",
    "V-REGION identity %with ins/del events)",
    "V-REGION identity ntwith ins/del events)",
    "J-GENE and allele",
    "J-REGION score",
    "J-REGION identity %",
    "J-REGION identity nt",
    "D-GENE and allele",
    "D-REGION reading frame",
    "CDR1-IMGT length",
    "CDR2-IMGT length",
    "CDR3-IMGT length",
    "CDR-IMGT lengths",
    "FR-IMGT lengths",
    "AA JUNCTION",
    "JUNCTION frame",
    "Orientation",
    "Functionality comment",
    "V-REGION potential ins/del",
    "J-GENE and allele comment",
    "V-REGION insertions",
    "V-REGION deletions",
    "Sequence",
]

TABULAR_HEADER = [
    "seq_id",
    "uid",
    "chain",
    "productive",
    "v_full",
    "v_gene",
    "d_full",
    "d_gene",
    "j_full",
    "j_gene",
    "junction_aa",
    "junction_nt",
    "cdr3_length",
    "fr1_aa",
    "fr2_aa",
    "fr3_aa",
    "fr4_aa",
    "cdr1_aa",
    "cdr2_aa",
    "cdr3_nt",
    "cdr3_aa",
    "v_start",
    "vdj_nt",
    "vdj_aa",
    "var_muts_nt",
    "var_muts_aa",
    "var_identity_nt",
    "var_identity_aa",
    "var_mut_count_nt",
    "var_mut_count_aa",
    "var_ins",
    "var_del",
    "isotype",
    "species",
    "raw_input",
]

AIRR_HEADER = [
    "sequence_id",
    "umi",
    "sequence",
    "sequence_aa",
    "sequence_input",
    "rev_comp",
    "productive",
    "productivity_issues",
    "stop_codon",
    "v_call",
    "v_gene",
    "v_score",
    "v_identity",
    "v_identity_aa",
    "d_call",
    "d_gene",
    "d_score",
    "j_call",
    "j_gene",
    "j_score",
    "j_identity",
    "j_identity_aa",
    "c_call",
    "cdr3_length",
    "fwr1_aa",
    "fwr2_aa",
    "fwr3_aa",
    "fwr4_aa",
    "cdr1_aa",
    "cdr2_aa",
    "cdr3_aa",
    "sequence_alignment",
    "germline_alignment",
    "v_mutations",
    "v_mutations_aa",
    "v_mutation_count",
    "v_mutation_count_aa",
    "v_insertions",
    "v_deletions",
    "junction",
    "junction_aa",
    "junction_in_frame",
    "junction_v",
    "junction_d",
    "junction_j",
    "junction_n",
    "junction_n1",
    "junction_n2",
    "junction_v_aa",
    "junction_d_aa",
    "junction_j_aa",
    "junction_n_aa",
    "junction_n1_aa",
    "junction_n2_aa",
    "junction_germ_v",
    "junction_germ_d",
    "junction_germ_j",
    "junction_germ_v_aa",
    "junction_germ_d_aa",
    "junction_germ_j_aa",
    "isotype",
    "locus",
    "v_cigar",
    "d_cigar",
    "j_cigar",
    "species",
    "germline_database",
]

