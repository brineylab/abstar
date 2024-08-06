import pyarrow as pa


# Add pyarrow schema to use when writing JSON-formatted abstar output to parquet files.
# Adding this explicitly defined schema attributes allows for nested data structures to be queried more efficiently.
schema = pa.schema(
    [
        pa.field("seq_id", pa.string()),
        pa.field("chain", pa.string()),
        pa.field(
            "v_gene",
            pa.struct(
                [
                    pa.field("full", pa.string()),
                    pa.field("fam", pa.string()),
                    pa.field("gene", pa.string()),
                    pa.field("score", pa.int64()),
                    pa.field("assigner_score", pa.float64()),
                    pa.field(
                        "others",
                        pa.list_(
                            pa.struct(
                                [
                                    pa.field("full", pa.string()),
                                    pa.field("assigner_score", pa.float64()),
                                ]
                            )
                        ),
                    ),
                ]
            ),
        ),
        pa.field(
            "d_gene",
            pa.struct(
                [
                    pa.field("full", pa.string()),
                    pa.field("fam", pa.string()),
                    pa.field("gene", pa.string()),
                    pa.field("score", pa.int64()),
                    pa.field("assigner_score", pa.int64()),
                    pa.field(
                        "others",
                        pa.list_(
                            pa.struct(
                                [
                                    pa.field("full", pa.string()),
                                    pa.field("assigner_score", pa.int64()),
                                ]
                            )
                        ),
                    ),
                ]
            ),
        ),
        pa.field(
            "j_gene",
            pa.struct(
                [
                    pa.field("full", pa.string()),
                    pa.field("gene", pa.string()),
                    pa.field("score", pa.int64()),
                    pa.field("assigner_score", pa.float64()),
                    pa.field(
                        "others",
                        pa.list_(
                            pa.struct(
                                [
                                    pa.field("full", pa.string()),
                                    pa.field("assigner_score", pa.float64()),
                                ]
                            )
                        ),
                    ),
                ]
            ),
        ),
        pa.field(
            "assigner_scores",
            pa.struct(
                [
                    pa.field("v", pa.float64()),
                    pa.field("d", pa.int64()),
                    pa.field("j", pa.float64()),
                ]
            ),
        ),
        pa.field("vdj_assigner", pa.string()),
        pa.field("isotype", pa.string()),
        pa.field("isotype_score", pa.int64()),
        pa.field(
            "isotype_alignment",
            pa.struct(
                [
                    pa.field("query", pa.string()),
                    pa.field("midline", pa.string()),
                    pa.field("isotype", pa.string()),
                ]
            ),
        ),
        pa.field(
            "nt_identity",
            pa.struct(
                [
                    pa.field("v", pa.float64()),
                    pa.field("j", pa.float64()),
                ]
            ),
        ),
        pa.field(
            "aa_identity",
            pa.struct(
                [
                    pa.field("v", pa.float64()),
                    pa.field("j", pa.float64()),
                ]
            ),
        ),
        pa.field("junc_len", pa.int64()),
        pa.field("cdr3_len", pa.int64()),
        pa.field("vdj_nt", pa.string()),
        pa.field("gapped_vdj_nt", pa.string()),
        pa.field("fr1_nt", pa.string()),
        pa.field("cdr1_nt", pa.string()),
        pa.field("fr2_nt", pa.string()),
        pa.field("cdr2_nt", pa.string()),
        pa.field("fr3_nt", pa.string()),
        pa.field("cdr3_nt", pa.string()),
        pa.field("fr4_nt", pa.string()),
        pa.field("vdj_germ_nt", pa.string()),
        pa.field("gapped_vdj_germ_nt", pa.string()),
        pa.field("junc_nt", pa.string()),
        pa.field(
            "region_len_nt",
            pa.struct(
                [
                    pa.field("fr1", pa.int64()),
                    pa.field("cdr1", pa.int64()),
                    pa.field("fr2", pa.int64()),
                    pa.field("cdr2", pa.int64()),
                    pa.field("fr3", pa.int64()),
                    pa.field("cdr3", pa.int64()),
                    pa.field("fr4", pa.int64()),
                ]
            ),
        ),
        pa.field(
            "var_muts_nt",
            pa.struct(
                [
                    pa.field("num", pa.int64()),
                    pa.field(
                        "muts",
                        pa.list_(
                            pa.struct(
                                [
                                    pa.field("was", pa.string()),
                                    pa.field("is", pa.string()),
                                    pa.field("raw_position", pa.string()),
                                    pa.field("position", pa.string()),
                                    pa.field("codon", pa.string()),
                                ]
                            )
                        ),
                    ),
                ]
            ),
        ),
        pa.field(
            "join_muts_nt",
            pa.struct(
                [
                    pa.field("num", pa.int64()),
                    pa.field(
                        "muts",
                        pa.list_(
                            pa.struct(
                                [
                                    pa.field("was", pa.string()),
                                    pa.field("is", pa.string()),
                                    pa.field("raw_position", pa.string()),
                                    pa.field("position", pa.string()),
                                    pa.field("codon", pa.string()),
                                ]
                            )
                        ),
                    ),
                ]
            ),
        ),
        pa.field("mut_count_nt", pa.int64()),
        pa.field("vdj_aa", pa.string()),
        pa.field("fr1_aa", pa.string()),
        pa.field("cdr1_aa", pa.string()),
        pa.field("fr2_aa", pa.string()),
        pa.field("cdr2_aa", pa.string()),
        pa.field("fr3_aa", pa.string()),
        pa.field("cdr3_aa", pa.string()),
        pa.field("fr4_aa", pa.string()),
        pa.field("vdj_germ_aa", pa.string()),
        pa.field("junc_aa", pa.string()),
        pa.field(
            "region_len_aa",
            pa.struct(
                [
                    pa.field("fr1", pa.int64()),
                    pa.field("cdr1", pa.int64()),
                    pa.field("fr2", pa.int64()),
                    pa.field("cdr2", pa.int64()),
                    pa.field("fr3", pa.int64()),
                    pa.field("cdr3", pa.int64()),
                    pa.field("fr4", pa.int64()),
                ]
            ),
        ),
        pa.field(
            "var_muts_aa",
            pa.struct(
                [
                    pa.field("num", pa.int64()),
                    pa.field(
                        "muts",
                        pa.list_(
                            pa.struct(
                                [
                                    pa.field("was", pa.string()),
                                    pa.field("is", pa.string()),
                                    pa.field("raw_position", pa.string()),
                                    pa.field("position", pa.string()),
                                    pa.field("codon", pa.string()),
                                ]
                            )
                        ),
                    ),
                ]
            ),
        ),
        pa.field(
            "join_muts_aa",
            pa.struct(
                [
                    pa.field("num", pa.int64()),
                    pa.field(
                        "muts",
                        pa.list_(
                            pa.struct(
                                [
                                    pa.field("was", pa.string()),
                                    pa.field("is", pa.string()),
                                    pa.field("raw_position", pa.string()),
                                    pa.field("position", pa.string()),
                                    pa.field("codon", pa.string()),
                                ]
                            )
                        ),
                    ),
                ]
            ),
        ),
        pa.field(
            "region_muts_nt",
            pa.struct(
                [
                    pa.field(
                        "fr1",
                        pa.struct(
                            [
                                pa.field("num", pa.int64()),
                                pa.field(
                                    "muts",
                                    pa.list_(
                                        pa.struct(
                                            [
                                                pa.field("was", pa.string()),
                                                pa.field("is", pa.string()),
                                                pa.field("raw_position", pa.string()),
                                                pa.field("position", pa.string()),
                                                pa.field("codon", pa.string()),
                                            ]
                                        )
                                    ),
                                ),
                            ]
                        ),
                    ),
                    pa.field(
                        "cdr1",
                        pa.struct(
                            [
                                pa.field("num", pa.int64()),
                                pa.field(
                                    "muts",
                                    pa.list_(
                                        pa.struct(
                                            [
                                                pa.field("was", pa.string()),
                                                pa.field("is", pa.string()),
                                                pa.field("raw_position", pa.string()),
                                                pa.field("position", pa.string()),
                                                pa.field("codon", pa.string()),
                                            ]
                                        )
                                    ),
                                ),
                            ]
                        ),
                    ),
                    pa.field(
                        "fr2",
                        pa.struct(
                            [
                                pa.field("num", pa.int64()),
                                pa.field(
                                    "muts",
                                    pa.list_(
                                        pa.struct(
                                            [
                                                pa.field("was", pa.string()),
                                                pa.field("is", pa.string()),
                                                pa.field("raw_position", pa.string()),
                                                pa.field("position", pa.string()),
                                                pa.field("codon", pa.string()),
                                            ]
                                        )
                                    ),
                                ),
                            ]
                        ),
                    ),
                    pa.field(
                        "cdr2",
                        pa.struct(
                            [
                                pa.field("num", pa.int64()),
                                pa.field(
                                    "muts",
                                    pa.list_(
                                        pa.struct(
                                            [
                                                pa.field("was", pa.string()),
                                                pa.field("is", pa.string()),
                                                pa.field("raw_position", pa.string()),
                                                pa.field("position", pa.string()),
                                                pa.field("codon", pa.string()),
                                            ]
                                        )
                                    ),
                                ),
                            ]
                        ),
                    ),
                    pa.field(
                        "fr3",
                        pa.struct(
                            [
                                pa.field("num", pa.int64()),
                                pa.field(
                                    "muts",
                                    pa.list_(
                                        pa.struct(
                                            [
                                                pa.field("was", pa.string()),
                                                pa.field("is", pa.string()),
                                                pa.field("raw_position", pa.string()),
                                                pa.field("position", pa.string()),
                                                pa.field("codon", pa.string()),
                                            ]
                                        )
                                    ),
                                ),
                            ]
                        ),
                    ),
                    pa.field(
                        "cdr3",
                        pa.struct(
                            [
                                pa.field("num", pa.int64()),
                                pa.field(
                                    "muts",
                                    pa.list_(
                                        pa.struct(
                                            [
                                                pa.field("was", pa.string()),
                                                pa.field("is", pa.string()),
                                                pa.field("raw_position", pa.string()),
                                                pa.field("position", pa.string()),
                                                pa.field("codon", pa.string()),
                                            ]
                                        )
                                    ),
                                ),
                            ]
                        ),
                    ),
                    pa.field(
                        "fr4",
                        pa.struct(
                            [
                                pa.field("num", pa.int64()),
                                pa.field(
                                    "muts",
                                    pa.list_(
                                        pa.struct(
                                            [
                                                pa.field("was", pa.string()),
                                                pa.field("is", pa.string()),
                                                pa.field("raw_position", pa.string()),
                                                pa.field("position", pa.string()),
                                                pa.field("codon", pa.string()),
                                            ]
                                        )
                                    ),
                                ),
                            ]
                        ),
                    ),
                ]
            ),
        ),
        pa.field(
            "region_muts_aa",
            pa.struct(
                [
                    pa.field(
                        "fr1",
                        pa.struct(
                            [
                                pa.field("num", pa.int64()),
                                pa.field(
                                    "muts",
                                    pa.list_(
                                        pa.struct(
                                            [
                                                pa.field("was", pa.string()),
                                                pa.field("is", pa.string()),
                                                pa.field("raw_position", pa.string()),
                                                pa.field("position", pa.string()),
                                                pa.field("codon", pa.string()),
                                            ]
                                        )
                                    ),
                                ),
                            ]
                        ),
                    ),
                    pa.field(
                        "cdr1",
                        pa.struct(
                            [
                                pa.field("num", pa.int64()),
                                pa.field(
                                    "muts",
                                    pa.list_(
                                        pa.struct(
                                            [
                                                pa.field("was", pa.string()),
                                                pa.field("is", pa.string()),
                                                pa.field("raw_position", pa.string()),
                                                pa.field("position", pa.string()),
                                                pa.field("codon", pa.string()),
                                            ]
                                        )
                                    ),
                                ),
                            ]
                        ),
                    ),
                    pa.field(
                        "fr2",
                        pa.struct(
                            [
                                pa.field("num", pa.int64()),
                                pa.field(
                                    "muts",
                                    pa.list_(
                                        pa.struct(
                                            [
                                                pa.field("was", pa.string()),
                                                pa.field("is", pa.string()),
                                                pa.field("raw_position", pa.string()),
                                                pa.field("position", pa.string()),
                                                pa.field("codon", pa.string()),
                                            ]
                                        )
                                    ),
                                ),
                            ]
                        ),
                    ),
                    pa.field(
                        "cdr2",
                        pa.struct(
                            [
                                pa.field("num", pa.int64()),
                                pa.field(
                                    "muts",
                                    pa.list_(
                                        pa.struct(
                                            [
                                                pa.field("was", pa.string()),
                                                pa.field("is", pa.string()),
                                                pa.field("raw_position", pa.string()),
                                                pa.field("position", pa.string()),
                                                pa.field("codon", pa.string()),
                                            ]
                                        )
                                    ),
                                ),
                            ]
                        ),
                    ),
                    pa.field(
                        "fr3",
                        pa.struct(
                            [
                                pa.field("num", pa.int64()),
                                pa.field(
                                    "muts",
                                    pa.list_(
                                        pa.struct(
                                            [
                                                pa.field("was", pa.string()),
                                                pa.field("is", pa.string()),
                                                pa.field("raw_position", pa.string()),
                                                pa.field("position", pa.string()),
                                                pa.field("codon", pa.string()),
                                            ]
                                        )
                                    ),
                                ),
                            ]
                        ),
                    ),
                    pa.field(
                        "cdr3",
                        pa.struct(
                            [
                                pa.field("num", pa.int64()),
                                pa.field(
                                    "muts",
                                    pa.list_(
                                        pa.struct(
                                            [
                                                pa.field("was", pa.string()),
                                                pa.field("is", pa.string()),
                                                pa.field("raw_position", pa.string()),
                                                pa.field("position", pa.string()),
                                                pa.field("codon", pa.string()),
                                            ]
                                        )
                                    ),
                                ),
                            ]
                        ),
                    ),
                    pa.field(
                        "fr4",
                        pa.struct(
                            [
                                pa.field("num", pa.int64()),
                                pa.field(
                                    "muts",
                                    pa.list_(
                                        pa.struct(
                                            [
                                                pa.field("was", pa.string()),
                                                pa.field("is", pa.string()),
                                                pa.field("raw_position", pa.string()),
                                                pa.field("position", pa.string()),
                                                pa.field("codon", pa.string()),
                                            ]
                                        )
                                    ),
                                ),
                            ]
                        ),
                    ),
                ]
            ),
        ),
        pa.field("prod", pa.string()),
        pa.field("productivity_issues", pa.string()),
        pa.field("junction_in_frame", pa.string()),
        pa.field("raw_input", pa.string()),
        pa.field("oriented_input", pa.string()),
        pa.field("strand", pa.string()),
        pa.field(
            "germ_alignments_nt",
            pa.struct(
                [
                    pa.field(
                        "var",
                        pa.struct(
                            [
                                pa.field("query", pa.string()),
                                pa.field("germ", pa.string()),
                                pa.field("midline", pa.string()),
                            ]
                        ),
                    ),
                    pa.field(
                        "join",
                        pa.struct(
                            [
                                pa.field("query", pa.string()),
                                pa.field("germ", pa.string()),
                                pa.field("midline", pa.string()),
                            ]
                        ),
                    ),
                    pa.field(
                        "div",
                        pa.struct(
                            [
                                pa.field("query", pa.string()),
                                pa.field("germ", pa.string()),
                                pa.field("midline", pa.string()),
                            ]
                        ),
                    ),
                ]
            ),
        ),
        pa.field(
            "exo_trimming",
            pa.struct(
                [
                    pa.field("var_3", pa.int64()),
                    pa.field("join_5", pa.int64()),
                    pa.field("div_5", pa.int64()),
                    pa.field("div_3", pa.int64()),
                ]
            ),
        ),
        pa.field(
            "junc_nt_breakdown",
            pa.struct(
                [
                    pa.field("v_nt", pa.string()),
                    pa.field("n_nt", pa.string()),
                    pa.field("n1_nt", pa.string()),
                    pa.field("d_nt", pa.string()),
                    pa.field("n2_nt", pa.string()),
                    pa.field("j_nt", pa.string()),
                    pa.field("d_dist_from_cdr3_start", pa.int64()),
                    pa.field("d_dist_from_cdr3_end", pa.int64()),
                ]
            ),
        ),
        pa.field("germline_database", pa.string()),
        pa.field("species", pa.string()),
        pa.field(
            "align_info",
            pa.struct(
                [
                    pa.field("v_start", pa.int64()),
                    pa.field("v_end", pa.int64()),
                    pa.field("d_start", pa.int64()),
                    pa.field("d_end", pa.int64()),
                    pa.field("j_start", pa.int64()),
                    pa.field("j_end", pa.int64()),
                ]
            ),
        ),
        pa.field(
            "j_del",
            pa.list_(
                pa.struct(
                    [
                        pa.field("in_frame", pa.string()),
                        pa.field("length", pa.int64()),
                        pa.field("sequence", pa.string()),
                        pa.field("position", pa.string()),
                        pa.field("codon", pa.string()),
                    ]
                )
            ),
        ),
        pa.field(
            "v_del",
            pa.list_(
                pa.struct(
                    [
                        pa.field("in_frame", pa.string()),
                        pa.field("length", pa.int64()),
                        pa.field("sequence", pa.string()),
                        pa.field("position", pa.string()),
                        pa.field("codon", pa.string()),
                    ]
                )
            ),
        ),
        pa.field(
            "j_ins",
            pa.list_(
                pa.struct(
                    [
                        pa.field("in_frame", pa.string()),
                        pa.field("length", pa.int64()),
                        pa.field("sequence", pa.string()),
                        pa.field("position", pa.string()),
                        pa.field("codon", pa.string()),
                    ]
                )
            ),
        ),
        pa.field(
            "v_ins",
            pa.list_(
                pa.struct(
                    [
                        pa.field("in_frame", pa.string()),
                        pa.field("length", pa.int64()),
                        pa.field("sequence", pa.string()),
                        pa.field("position", pa.string()),
                        pa.field("codon", pa.string()),
                    ]
                )
            ),
        ),
    ]
)
