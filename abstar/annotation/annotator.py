# Copyright (c) 2024 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

import os
import traceback
from typing import Iterable, Optional, Union

import abutils

# import pandas as pd
import polars as pl

from .antibody import Antibody
from .germline import (
    get_germline,
    process_cgene_alignment,
    process_dgene_alignment,
    process_jgene_alignment,
    process_vgene_alignment,
    realign_germline,
    reassign_dgene,
)
from .indels import annotate_deletions, annotate_insertions
from .mask import (
    generate_cdr_mask,
    generate_gene_segment_mask,
    generate_nongermline_mask,
)
from .mutations import annotate_c_mutations, annotate_v_mutations
from .positions import get_gapped_sequence, get_ungapped_position_from_aligned
from .productivity import assess_productivity
from .regions import get_region_sequence, identify_cdr3_regions
from .schema import OUTPUT_SCHEMA
from .umi import parse_umis


def annotate(
    input_file: str,
    output_directory: str,
    germline_database: str,
    log_directory: Optional[str] = None,
    umi_pattern: Optional[str] = None,
    umi_length: Optional[int] = None,
    debug: bool = False,
) -> Iterable[Union[str, Optional[str], Optional[str]]]:
    """
    Annotates a Parquet file of V(D)J-assigned antibody sequences.

    Parameters
    ----------
    input_file : str
        Path to the input file, in Parquet format. The following fields are required:

            - sequence_id
            - sequence_input
            - quality
            - rev_comp
            - v_call
            - v_support
            - d_call
            - d_support
            - j_call
            - j_support
            - c_call
            - c_support

        although any of the values (aside from sequence_id, sequence_input, and rev_comp)
        can be None.

    output_directory : str
        Path to the directory where the annotated Parquet file will be written.

    germline_database : str
        Name of the germline database to use for annotation.

    log_directory : Optional[str], default=None
        Path to the directory where the log files will be written.

    umi_pattern : Optional[str], default=None
        Pattern to match for parsing the UMI sequence, or name of a built-in pattern.

    umi_length : Optional[int], default=None
        Length of the UMI sequence.

    debug : bool, default=False
        Whether to write log files for failed and succeeded sequences.

    Returns
    -------
    output_file : str
        Path to the output Parquet file containing annotated sequences.

    failed_logfile : Optional[str]
        Path to the log file for failed sequences. Only returned if
        `log_directory` is provided

    succeeded_logfile : Optional[str]
        Path to the log file for succeeded sequences. Only returned if
        `debug=True` and `log_directory` is provided.

    """
    # load the input parquet file
    df = pl.read_parquet(input_file)

    # do annotations
    annotated = []
    for r in df.iter_rows(named=True):
        ab = Antibody(**r)
        try:
            ab = annotate_single_sequence(
                ab=ab,
                germline_database=germline_database,
                umi_pattern=umi_pattern,
                umi_length=umi_length,
                debug=debug,
            )
        except Exception as e:
            ab.exception("ANNOTATION EXCEPTION", traceback.format_exc())
        annotated.append(ab)

    # gather the results
    failed = [a for a in annotated if a.exceptions]
    succeeded = [a for a in annotated if not a.exceptions]
    succeeded_df = pl.DataFrame([s.to_dict() for s in succeeded], schema=OUTPUT_SCHEMA)

    # write logs
    basename = os.path.basename(input_file)
    if log_directory:
        failed_logfile = os.path.join(log_directory, f"{basename}.failed")
        with open(failed_logfile, "w") as f:
            for fail in failed:
                f.write(fail.format_log())
        if debug:
            succeeded_logfile = os.path.join(log_directory, f"{basename}.succeeded")
            with open(succeeded_logfile, "w") as f:
                for succ in succeeded:
                    f.write(succ.format_log())

    # write outputs
    output_file = os.path.join(output_directory, f"{basename}_annotated.parquet")
    succeeded_df.write_parquet(output_file)

    # returns
    if log_directory:
        if debug:
            return output_file, failed_logfile, succeeded_logfile
        else:
            return output_file, failed_logfile, None
    else:
        return output_file, None, None


def annotate_single_sequence(
    ab: Antibody,
    germline_database: str,
    umi_pattern: Optional[str] = None,
    umi_length: Optional[int] = None,
    debug: bool = False,
):
    # log information from the assigner
    ab.log("=" * (len(str(ab.sequence_id)) + 15))
    ab.log(" SEQUENCE ID:", ab.sequence_id)
    ab.log("=" * (len(str(ab.sequence_id)) + 15) + "\n")
    ab.log(f">{ab.sequence_id}\n{ab.sequence_input}\n")
    ab.log(f"{ab.v_call} | {ab.d_call} | {ab.j_call}")
    ab.log("REV COMP:", ab.rev_comp)

    # germline database, locus, and receptor type
    ab.germline_database = germline_database
    ab.locus = ab.v_call[:3].upper()
    ab.receptor_type = "bcr" if ab.locus[:2] == "IG" else "tcr"
    ab.log("GERMLINE DATABASE:", ab.germline_database)
    ab.log("LOCUS:", ab.locus)

    # umi
    if any([umi_pattern is not None, umi_length is not None]):
        ab.umi = parse_umis(
            ab.sequence_input,
            pattern=umi_pattern,
            length=umi_length,
        )

    # orient the input sequence
    if ab.rev_comp:
        ab.sequence_oriented = abutils.tl.reverse_complement(ab.sequence_input)
    else:
        ab.sequence_oriented = ab.sequence_input
    ab.log("SEQUENCE ORIENTED:", ab.sequence_oriented)

    # check germline calls for species name
    # (mainly for mixed species databases like humouse)
    # need to keep the raw call for calls to get_germline(), which must include the species name if it's present
    raw_v_call = ab.v_call
    raw_j_call = ab.j_call
    raw_d_call = ab.d_call
    raw_c_call = ab.c_call
    if len(vsplit := ab.v_call.split("__")) > 1:
        ab.v_call, ab.species = vsplit
    else:
        ab.species = ab.germline_database
    if len(jsplit := ab.j_call.split("__")) > 1:
        ab.j_call, _ = jsplit
    if ab.d_call is not None:
        if len(dsplit := ab.d_call.split("__")) > 1:
            ab.d_call, _ = dsplit
    if ab.c_call is not None:
        if len(c_split := ab.c_call.split("__")) > 1:
            ab.c_call, _ = c_split
    ab.log("SPECIES:", ab.species)

    # get genes from calls
    ab.v_gene = ab.v_call.split("*")[0]
    ab.j_gene = ab.j_call.split("*")[0]
    if ab.d_call is not None:
        ab.d_gene = ab.d_call.split("*")[0]
    if ab.c_call is not None:
        ab.c_gene = ab.c_call.split("*")[0]

    ab.log("\n--------")
    ab.log(" V GENE")
    ab.log("--------\n")

    # v-gene realignment
    v_sg, v_loc = realign_germline(
        sequence=ab.sequence_oriented,
        germline_name=raw_v_call,  # must include species annotation if it's present in the germline database
        germdb_name=ab.germline_database,
        imgt_gapped=False,
        semiglobal_aln_params=ALIGNMENT_PARAMS,
        local_aln_params=ALIGNMENT_PARAMS,
    )
    ab.log("V GERMLINE:", v_sg.target.sequence)

    ab.log("SEMIGLOBAL ALIGNMENT:")
    ab.log(f"     QUERY: {v_sg.aligned_query}")
    ab.log(f"            {v_sg.alignment_midline}")
    ab.log(f"  GERMLINE: {v_sg.aligned_target}")
    ab.log("LOCAL ALIGNMENT:")
    ab.log(f"     QUERY: {v_loc.aligned_query}")
    ab.log(f"            {v_loc.alignment_midline}")
    ab.log(f"  GERMLINE: {v_loc.aligned_target}")

    ab = process_vgene_alignment(semiglobal_aln=v_sg, local_aln=v_loc, ab=ab)
    ab.log("SEMIGLOBAL QUERY START:", v_sg.query_begin)
    ab.log("SEMIGLOBAL GERMLINE START:", v_sg.target_begin)
    ab.log("LOCAL QUERY START:", v_loc.query_begin)
    ab.log("LOCAL GERMLINE START:", v_loc.target_begin)
    ab.log("V SEQUENCE START:", ab.v_sequence_start)
    ab.log("V SEQUENCE END:", ab.v_sequence_end)
    ab.log("V GERMLINE START:", ab.v_germline_start)
    ab.log("V GERMLINE END:", ab.v_germline_end)
    ab.log("V FRAME:", ab.frame)
    ab.log("V SEQUENCE:", ab.v_sequence)
    ab.log("V GERMLINE:", ab.v_germline)
    ab.log("V SEQUENCE AA:", ab.v_sequence_aa)
    ab.log("V GERMLINE AA:", ab.v_germline_aa)

    # semiglobal alignment of AA sequences
    # first need to translate the query and germline sequences
    # query gets truncated at the start of the germline alignment
    # and gets translated in its frame. The full germline is
    # translated and used for alignment.
    query_aa_sg = abutils.tl.translate(v_sg.query[v_sg.query_begin :], frame=ab.frame)
    germline_aa_sg = abutils.tl.translate(v_sg.target)
    v_sg_aa = abutils.tl.semiglobal_alignment(
        query=query_aa_sg, target=germline_aa_sg, **ALIGNMENT_PARAMS
    )

    ab.log("SEMIGLOBAL ALIGNMENT AA:")
    ab.log(f"     QUERY: {v_sg_aa.aligned_query}")
    ab.log(f"            {v_sg_aa.alignment_midline}")
    ab.log(f"  GERMLINE: {v_sg_aa.aligned_target}")

    # global alignment of sequence and germline
    # using the start/end positions determined parsed from the
    # semiglobal and local alignments
    v_global = abutils.tl.global_alignment(
        ab.v_sequence,
        ab.v_germline,
        **ALIGNMENT_PARAMS,
    )
    ab.log("GLOBAL ALIGNMENT:")
    ab.log(f"     QUERY: {v_global.aligned_query}")
    ab.log(f"            {v_global.alignment_midline}")
    ab.log(f"  GERMLINE: {v_global.aligned_target}")

    # global alignment of the AA sequences
    v_global_aa = abutils.tl.global_alignment(
        ab.v_sequence_aa,
        ab.v_germline_aa,
        **ALIGNMENT_PARAMS,
    )
    ab.log("GLOBAL ALIGNMENT AA:")
    ab.log(f"     QUERY: {v_global_aa.aligned_query}")
    ab.log(f"            {v_global_aa.alignment_midline}")
    ab.log(f"  GERMLINE: {v_global_aa.aligned_target}")

    # gapped V-gene germline
    ab.v_germline_gapped = get_germline(
        raw_v_call,
        ab.germline_database,
        imgt_gapped=True,
        exact_match=True,
        truncate_species=False,
    ).sequence
    ab.v_germline_gapped_aa = abutils.tl.translate(
        ab.v_germline_gapped, allow_dots=True
    )
    ab.log("GAPPED GERMLINE:", ab.v_germline_gapped)
    ab.log("GAPPED GERMLINE AA:", ab.v_germline_gapped_aa)
    # ab.v_germline_gapped = ab.v_germline_gapped
    # ab.v_germline_gapped_aa = ab.v_germline_gapped_aa

    # gapped V-gene sequence
    ab.v_sequence_gapped = get_gapped_sequence(
        aligned_sequence=v_global.aligned_query,
        aligned_germline=v_global.aligned_target,
        gapped_germline=ab.v_germline_gapped,
        germline_start=ab.v_germline_start,
    )
    ab.v_sequence_gapped_aa = abutils.tl.translate(
        ab.v_sequence_gapped, frame=ab.frame, allow_dots=True
    )
    ab.log("GAPPED SEQUENCE:", ab.v_sequence_gapped)
    ab.log("GAPPED SEQUENCE AA:", ab.v_sequence_gapped_aa)
    # ab.v_sequence_gapped = ab.v_sequence_gapped
    # ab.v_sequence_gapped_aa = ab.v_sequence_gapped_aa

    # insertions
    if "-" in v_loc.aligned_target:
        ab.v_insertions = annotate_insertions(
            aligned_sequence=v_global.aligned_query,
            aligned_germline=v_global.aligned_target,
            gapped_germline=ab.v_germline_gapped,
            germline_start=ab.v_germline_start,
        )
        ab.log("V INSERTIONS:", ab.v_insertions)
        insertions = ab.v_insertions.split("|")
        cumulative_ins_length = sum(
            [int(i.split(">")[0].split(":")[1]) for i in insertions]
        )
    else:
        cumulative_ins_length = 0

    # deletions
    if "-" in v_loc.aligned_query:
        ab.v_deletions = annotate_deletions(
            aligned_sequence=v_loc.aligned_query,
            aligned_germline=v_loc.aligned_target,
            gapped_germline=ab.v_germline_gapped,
            germline_start=ab.v_germline_start,
        )
        ab.log("V DELETIONS:", ab.v_deletions)
        deletions = ab.v_deletions.split("|")
        cumulative_del_length = sum(
            [int(d.split(">")[0].split(":")[1]) for d in deletions]
        )
    else:
        cumulative_del_length = 0

    # only raise a productivity issue if the combination of all indels causes a frameshift
    # this avoids incorrectly flagging as non-productive a sequence with compensatory out-of-frame indels
    if (cumulative_ins_length - cumulative_del_length) % 3:
        ab.productivity_issues.append("out-of-frame indel(s)")
        ab.v_frameshift = True

    ab.log("\n--------")
    ab.log(" J GENE")
    ab.log("--------\n")

    # j-gene realignment
    jquery = ab.sequence_oriented[ab.v_sequence_end :]
    if len(jquery) == 0:
        raise ValueError(
            f"Query sequence for J-gene realignment is empty for {ab.sequence_id}"
        )
    # use a smaller gap open penalty for the J-gene alignment than we use for others (V-gene and constant region)
    j_sg, j_loc = realign_germline(
        sequence=jquery,
        germline_name=raw_j_call,  # must include species annotation if it's present in the germline database
        germdb_name=ab.germline_database,
        imgt_gapped=False,
        semiglobal_aln_params=ALIGNMENT_PARAMS | {"gap_open": -20},
        local_aln_params=ALIGNMENT_PARAMS | {"gap_open": -15},
    )
    ab.log("SEMIGLOBAL ALIGNMENT:")
    ab.log(f"     QUERY: {j_sg.aligned_query}")
    ab.log(f"            {j_sg.alignment_midline}")
    ab.log(f"  GERMLINE: {j_sg.aligned_target}")
    ab.log("")
    ab.log("LOCAL QUERY:", j_loc.query.sequence)
    ab.log("LOCAL GERMLINE:", j_loc.target.sequence)
    ab.log("LOCAL ALIGNMENT:")
    ab.log(f"     QUERY: {j_loc.aligned_query}")
    ab.log(f"            {j_loc.alignment_midline}")
    ab.log(f"  GERMLINE: {j_loc.aligned_target}")

    ab = process_jgene_alignment(
        oriented_input=ab.sequence_oriented,
        v_sequence_end=ab.v_sequence_end,
        semiglobal_aln=j_sg,
        local_aln=j_loc,
        ab=ab,
    )
    ab.log("SEMIGLOBAL QUERY START:", j_sg.query_begin)
    ab.log("LOCAL QUERY START:", j_loc.query_begin)
    ab.log("J SEQUENCE START:", ab.j_sequence_start)
    ab.log("J SEQUENCE END:", ab.j_sequence_end)
    ab.log("J GERMLINE START:", ab.j_germline_start)
    ab.log("J GERMLINE END:", ab.j_germline_end)
    ab.log("J SEQUENCE:", ab.j_sequence)
    ab.log("J GERMLINE:", ab.j_germline)

    ab.log("\n--------")
    ab.log(" D GENE")
    ab.log("--------\n")

    # d-gene realignment
    dquery = ab.sequence_oriented[ab.v_sequence_end : ab.j_sequence_start]
    # if the assigner made a d-gene call
    if dquery and ab.d_call is not None:
        _, d_loc = realign_germline(
            sequence=dquery,
            germline_name=raw_d_call,  # must include species annotation if it's present in the germline database
            germdb_name=ab.germline_database,
            imgt_gapped=False,
            skip_semiglobal=True,
            local_aln_params=ALIGNMENT_PARAMS,
        )
    # if not, we can try again using local pairwise alignment (for IGH/TRA/TRD chains only)
    elif len(dquery) >= 5 and ab.locus in ["IGH", "TRA", "TRD"]:
        d_loc = reassign_dgene(
            sequence=dquery,
            germdb_name=ab.germline_database,
            aln_params=ALIGNMENT_PARAMS,
        )
    # maybe there's not a d-gene (or it's a light/TRB/TRG chain)
    else:
        ab.d_call = None
        ab.d_gene = None
        ab.d_score = None
        d_loc = None
    # process the d-gene alignment
    if d_loc is not None:
        ab = process_dgene_alignment(
            oriented_input=ab.sequence_oriented,
            v_sequence_end=ab.v_sequence_end,
            local_aln=d_loc,
            ab=ab,
        )
        ab.np1 = ab.sequence_oriented[ab.v_sequence_end : ab.d_sequence_start]
        ab.np2 = ab.sequence_oriented[ab.d_sequence_end : ab.j_sequence_start]
        ab.np1_length = len(ab.np1)
        ab.np2_length = len(ab.np2)
    else:
        ab.np1 = ab.sequence_oriented[ab.v_sequence_end : ab.j_sequence_start]
        ab.np1_length = len(ab.np1)

    if d_loc is not None:
        ab.log("LOCAL ALIGNMENT:")
        ab.log(f"     QUERY: {d_loc.aligned_query}")
        ab.log(f"            {d_loc.alignment_midline}")
        ab.log(f"  GERMLINE: {d_loc.aligned_target}")

        ab.log("LOCAL QUERY START:", d_loc.query_begin)
        ab.log("D SEQUENCE START:", ab.d_sequence_start)
        ab.log("D SEQUENCE END:", ab.d_sequence_end)
        ab.log("D GERMLINE START:", ab.d_germline_start)
        ab.log("D GERMLINE END:", ab.d_germline_end)
        ab.log("D FRAME:", ab.d_frame)
        ab.log("D SEQUENCE:", ab.d_sequence)
        ab.log("D GERMLINE:", ab.d_germline)
        ab.log("NP1 SEQUENCE:", ab.np1)
        ab.log("NP1 LENGTH:", ab.np1_length)
        ab.log("NP2 SEQUENCE:", ab.np2)
        ab.log("NP2 LENGTH:", ab.np2_length)

    ab.log("\n----------------")
    ab.log(" CONSTANT REGION")
    ab.log("----------------\n")

    # constant region realignment
    cquery = ab.sequence_oriented[ab.j_sequence_end :]
    # if the assigner made a c-gene call
    if cquery and ab.c_call is not None:
        c_sg, c_loc = realign_germline(
            sequence=cquery,
            germline_name=raw_c_call,  # must include species annotation if it's present in the germline database
            germdb_name=ab.germline_database,
            imgt_gapped=False,
            semiglobal_aln_params=ALIGNMENT_PARAMS,
            local_aln_params=ALIGNMENT_PARAMS,
            force_constant=True,
        )

        ab.log("SEMIGLOBAL ALIGNMENT:")
        ab.log(f"     QUERY: {c_sg.aligned_query}")
        ab.log(f"            {c_sg.alignment_midline}")
        ab.log(f"  GERMLINE: {c_sg.aligned_target}")
        ab.log("LOCAL ALIGNMENT:")
        ab.log(f"     QUERY: {c_loc.aligned_query}")
        ab.log(f"            {c_loc.alignment_midline}")
        ab.log(f"  GERMLINE: {c_loc.aligned_target}")

        ab = process_cgene_alignment(
            oriented_input=ab.sequence_oriented,
            j_sequence_end=ab.j_sequence_end,
            semiglobal_aln=c_sg,
            local_aln=c_loc,
            ab=ab,
        )
        ab.log("SEMIGLOBAL QUERY START:", c_sg.query_begin)
        ab.log("LOCAL QUERY START:", c_loc.query_begin)
        ab.log("C SEQUENCE START:", ab.c_sequence_start)
        ab.log("C SEQUENCE END:", ab.c_sequence_end)
        ab.log("C GERMLINE START:", ab.c_germline_start)
        ab.log("C GERMLINE END:", ab.c_germline_end)
        ab.log("C SEQUENCE:", ab.c_sequence)
        ab.log("C GERMLINE:", ab.c_germline)
        ab.log("C SEQUENCE AA:", ab.c_sequence_aa)
        ab.log("C GERMLINE AA:", ab.c_germline_aa)

        if all(
            [
                ab.c_sequence,
                ab.c_germline,
                ab.c_sequence_aa,
                ab.c_germline_aa,
            ]
        ):
            # global nucleotide alignment
            c_global = abutils.tl.global_alignment(
                ab.c_sequence,
                ab.c_germline,
                **ALIGNMENT_PARAMS,
            )
            ab.log("GLOBAL ALIGNMENT:")
            ab.log(f"     QUERY: {c_global.aligned_query}")
            ab.log(f"            {c_global.alignment_midline}")
            ab.log(f"  GERMLINE: {c_global.aligned_target}")

            # # verify that the AA sequences are not empty
            # if len(ab.v_sequence_aa) == 0:
            #     raise ValueError(f"V-gene AA sequence is empty for {ab.sequence_id}")
            # if len(ab.v_germline_aa) == 0:
            #     raise ValueError(f"V-gene germline AA sequence is empty for {ab.sequence_id}")

            # global alignment of the AA sequences
            c_global_aa = abutils.tl.global_alignment(
                ab.c_sequence_aa,
                ab.c_germline_aa,
                **ALIGNMENT_PARAMS,
            )
            ab.log("GLOBAL ALIGNMENT AA:")
            ab.log(f"     QUERY: {c_global_aa.aligned_query}")
            ab.log(f"            {c_global_aa.alignment_midline}")
            ab.log(f"  GERMLINE: {c_global_aa.aligned_target}")

            # gapped Constant region germline
            ab.c_germline_gapped = get_germline(
                raw_c_call,
                ab.germline_database,
                imgt_gapped=True,
                exact_match=True,
                force_constant=True,
                truncate_species=False,
            ).sequence
            ab.c_germline_gapped_aa = abutils.tl.translate(
                ab.c_germline_gapped, allow_dots=True
            )
            ab.log("GAPPED GERMLINE:", ab.c_germline_gapped)
            ab.log("GAPPED GERMLINE AA:", ab.c_germline_gapped_aa)
            # ab.v_germline_gapped = ab.v_germline_gapped
            # ab.v_germline_gapped_aa = ab.v_germline_gapped_aa

            # gapped Constant region sequence
            ab.c_sequence_gapped = get_gapped_sequence(
                aligned_sequence=c_global.aligned_query,
                aligned_germline=c_global.aligned_target,
                gapped_germline=ab.c_germline_gapped,
                germline_start=ab.c_germline_start,
            )
            ab.c_sequence_gapped_aa = abutils.tl.translate(
                ab.c_sequence_gapped, frame=ab.frame, allow_dots=True
            )

            # nucleotide mutations
            complete_c_germline = c_sg.target.sequence
            ab = annotate_c_mutations(
                aligned_sequence=c_global.aligned_query,
                aligned_germline=c_global.aligned_target,
                gapped_germline=complete_c_germline,
                germline_start=ab.c_germline_start,
                is_aa=False,
                ab=ab,
            )
            ab.log("CONSTANT REGION MUTATIONS:", ab.c_mutations)
            ab.log("CONSTANT REGION MUTATION COUNT:", ab.c_mutation_count)

            # amino acid mutations
            complete_c_germline_aa = abutils.tl.translate(complete_c_germline)
            ab = annotate_c_mutations(
                aligned_sequence=c_global_aa.aligned_query,
                aligned_germline=c_global_aa.aligned_target,
                gapped_germline=complete_c_germline_aa,
                germline_start=ab.c_germline_start // 3,
                is_aa=True,
                ab=ab,
            )
            ab.log("CONSTANT REGION MUTATIONS AA:", ab.c_mutations_aa)
            ab.log("CONSTANT REGION MUTATION COUNT AA:", ab.c_mutation_count_aa)

            # calculate identify (nt and aa)
            ab.c_identity = 1 - ab.c_mutation_count / len(ab.c_germline)
            ab.c_identity_aa = 1 - ab.c_mutation_count_aa / len(ab.c_germline_aa)
            ab.log("CONSTANT REGION IDENTITY:", ab.c_identity)
            ab.log("CONSTANT REGION IDENTITY AA:", ab.c_identity_aa)

            # insertions
            if "-" in c_loc.aligned_target:
                ab.c_insertions = annotate_insertions(
                    aligned_sequence=c_global.aligned_query,
                    aligned_germline=c_global.aligned_target,
                    gapped_germline=complete_c_germline,
                    germline_start=ab.c_germline_start,
                )
                ab.log("CONSTANT REGION INSERTIONS:", ab.c_insertions)

            # deletions
            if "-" in c_loc.aligned_query:
                ab.c_deletions = annotate_deletions(
                    aligned_sequence=c_loc.aligned_query,
                    aligned_germline=c_loc.aligned_target,
                    gapped_germline=complete_c_germline,
                    germline_start=ab.c_germline_start,
                )
                ab.log("CONSTANT REGION DELETIONS:", ab.c_deletions)

    ab.log("\n--------------")
    ab.log(" VDJ ASSEMBLY")
    ab.log("--------------\n")

    # assemble the full V(D)J and germline sequences
    ab.sequence = ab.v_sequence + ab.np1
    ab.germline = ab.v_germline + ab.np1
    ab.sequence_gapped = ab.v_sequence_gapped + ab.np1
    ab.germline_gapped = ab.v_germline_gapped + ab.np1
    if ab.d_sequence is not None:
        ab.sequence += ab.d_sequence
        ab.germline += ab.d_germline
        ab.sequence_gapped += ab.d_sequence
        ab.germline_gapped += ab.d_germline
    if ab.np2 is not None:
        ab.sequence += ab.np2
        ab.germline += ab.np2
        ab.sequence_gapped += ab.np2
        ab.germline_gapped += ab.np2
    ab.sequence += ab.j_sequence
    ab.germline += ab.j_germline
    ab.sequence_gapped += ab.j_sequence
    ab.germline_gapped += ab.j_germline
    if ab.c_sequence is not None:
        # only compute VDJC sequences if there is a constant region
        ab.sequence_vdjc = ab.sequence + ab.c_sequence
        ab.germline_vdjc = ab.germline + ab.c_germline
        ab.sequence_vdjc_gapped = ab.sequence_gapped + ab.c_sequence
        ab.germline_vdjc_gapped = ab.germline_gapped + ab.c_germline

    ab.log("SEQUENCE:", ab.sequence)
    ab.log("GERMLINE:", ab.germline)
    ab.log("GAPPED SEQUENCE:", ab.sequence_gapped)
    ab.log("GAPPED GERMLINE:", ab.germline_gapped)

    ab.log("VDJC SEQUENCE:", ab.sequence_vdjc)
    ab.log("VDJC GERMLINE:", ab.germline_vdjc)
    ab.log("GAPPED VDJC SEQUENCE:", ab.sequence_vdjc_gapped)
    ab.log("GAPPED VDJC GERMLINE:", ab.germline_vdjc_gapped)

    # translated sequences
    ab.sequence_aa = abutils.tl.translate(ab.sequence, frame=ab.frame)
    ab.germline_aa = abutils.tl.translate(ab.germline, frame=ab.frame)
    # ab.sequence_gapped_aa = abutils.tl.translate(
    #     ab.sequence_gapped, frame=ab.frame, allow_dots=True
    # )
    # ab.germline_gapped_aa = abutils.tl.translate(
    #     ab.germline_gapped, frame=ab.frame, allow_dots=True
    # )

    if ab.c_sequence is not None:
        # only compute gapped VDJC sequences if there is a constant region
        ab.sequence_vdjc_aa = abutils.tl.translate(ab.sequence_vdjc, frame=ab.frame)
        ab.germline_vdjc_aa = abutils.tl.translate(ab.germline_vdjc, frame=ab.frame)
        # ab.sequence_vdjc_gapped_aa = abutils.tl.translate(
        #     ab.sequence_vdjc_gapped, frame=ab.frame, allow_dots=True
        # )
        # ab.germline_vdjc_gapped_aa = abutils.tl.translate(
        #     ab.germline_vdjc_gapped, frame=ab.frame, allow_dots=True
        # )
    ab.log("SEQUENCE AA:", ab.sequence_aa)
    ab.log("GERMLINE AA:", ab.germline_aa)
    # ab.log("GAPPED SEQUENCE AA:", ab.sequence_gapped_aa)
    # ab.log("GAPPED GERMLINE AA:", ab.germline_gapped_aa)

    ab.log("VDJC SEQUENCE AA:", ab.sequence_vdjc_aa)
    ab.log("VDJC GERMLINE AA:", ab.germline_vdjc_aa)
    # ab.log("GAPPED VDJC SEQUENCE AA:", ab.sequence_vdjc_gapped_aa)
    # ab.log("GAPPED VDJC GERMLINE AA:", ab.germline_vdjc_gapped_aa)

    # align assembled V(D)J and germline nucleotide sequences
    nt_aln = abutils.tl.global_alignment(
        query=ab.sequence,
        target=ab.germline,
        **ALIGNMENT_PARAMS,
    )
    ab.sequence_alignment = nt_aln.aligned_query
    ab.germline_alignment = nt_aln.aligned_target
    ab.log(
        "SEQUENCE ALIGNMENT:",
        f"     QUERY: {nt_aln.aligned_query}",
        f"            {nt_aln.alignment_midline}",
        f"  GERMLINE: {nt_aln.aligned_target}",
        separator="\n",
    )

    # align assembled V(D)J and germline amino acid sequences
    aa_aln = abutils.tl.global_alignment(
        query=ab.sequence_aa,
        target=ab.germline_aa,
        **ALIGNMENT_PARAMS,
    )
    ab.sequence_alignment_aa = aa_aln.aligned_query
    ab.germline_alignment_aa = aa_aln.aligned_target
    ab.log(
        "SEQUENCE AA ALIGNMENT:",
        f"     QUERY: {aa_aln.aligned_query}",
        f"            {aa_aln.alignment_midline}",
        f"  GERMLINE: {aa_aln.aligned_target}",
        separator="\n",
    )

    # check for complete VDJ
    if ab.v_germline_start == 0 and ab.j_germline_end == len(j_sg.target.sequence):
        ab.complete_vdj = True
    ab.log("COMPLETE VDJ:", ab.complete_vdj)

    ab.log("\n----------")
    ab.log(" JUNCTION")
    ab.log("----------\n")

    # junction start
    # the IMGT start position for FR3 is 196, but 1-indexed so we need to subtract 1 to slice correctly
    # also, since we're trying to identify the junction, we drop the last codon of FR3 (since the last codon is part of the junction)
    # germ_fr3_sequence = ab.v_germline_gapped[196 - 1 : 309].replace(".", "")

    # we include the conserved C codon (which forms the start of the junction) to avoid an edge case
    # in which a sequence has a deletion at the position immediately preceding the conserved C codon
    # because the C is so highly conserved, it's basically impossible for a sequence to have a deletion at that position
    germ_fr3_sequence = ab.v_germline_gapped[196 - 1 : 312].replace(".", "")
    ab.log("FR3 GERMLINE SEQUENCE:", germ_fr3_sequence)

    # edge case where there's an insertion (gaps in the germline when aligned to the query) near the end of the FWR3, which can cause a misalignment.
    # if the region 3' of the indel isn't not long enough to overcome the gap penalty, the sequences can be aligned like so:
    #
    #    QUERY:    ACCGGGACTCTGGGGTCCCAGATAGATTCAGCGGCAGTGGGTCAGGCACTGATTTCACACTGAAAATCAGCAGGGTGGAGGCTGAAGATGTTGTTGGGGTTTATTAC
    #              |||||||||||||||||||||| |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| |||||||  |    |||
    #    GERMLINE: ACCGGGACTCTGGGGTCCCAGACAGATTCAGCGGCAGTGGGTCAGGCACTGATTTCACACTGAAAATCAGCAGGGTGGAGGCTGAGGATGTTGGGGTTTATTAC
    #
    # which should actually be:
    #
    #    QUERY:    ACCGGGACTCTGGGGTCCCAGATAGATTCAGCGGCAGTGGGTCAGGCACTGATTTCACACTGAAAATCAGCAGGGTGGAGGCTGAAGATGTTGTTGGGGTTTATTAC
    #              |||||||||||||||||||||| |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ||   ||||||||||||||||
    #    GERMLINE: ACCGGGACTCTGGGGTCCCAGACAGATTCAGCGGCAGTGGGTCAGGCACTGATTTCACACTGAAAATCAGCAGGGTGGAGGCTGAGGA---TGTTGGGGTTTATTAC
    #
    # to fix, we'll first align the germline region to the pre-aligned germline sequence (from the V-gene semiglobal alignment),
    # which will give us the correct gapped alignment (since the full germline alignment presumably has enough 3' sequence to  overcome the gap penalty)
    # then we can use the re-aligned (and pre-gapped, if necessary) germline for alignment with the query sequence to correctly identify the FWR3 region.

    # since the sequences for alignment shoud be identical (aside from any pre-existing alignment gaps), we can use a smaller gap open penalty
    # this helps us capture gaps very close to the end of the sequence, for which there aren't enough matches on one side of the gap to overcome a high gap penalty.
    REALIGNMENT_PARAMS = ALIGNMENT_PARAMS | {"gap_open": -2}

    # realign germline FWR3 region to the aligned germline sequence and extract the aligned region sequence
    fr3_germ_realign = abutils.tl.semiglobal_alignment(
        query=germ_fr3_sequence,
        target=v_sg.aligned_target,
        **REALIGNMENT_PARAMS,
    )
    # semiglobal alignment is with the full oriented sequence, so there will typically be gaps on both ends of the aligned query (FR3 region)
    realigned_germ_fr3_sequence = fr3_germ_realign.aligned_query.lstrip("-").rstrip("-")
    ab.log("REALIGNED FR3 GERMLINE SEQUENCE:", realigned_germ_fr3_sequence)

    # the same edge case (but in reverse) can happen with a deletion (gaps in the query when aligned to the germline)
    # near the end of the FR3, which causes incorrect identification of the start of the junction.
    # to fix, we'll use the aligned version of the query sequence (from the V-gene semiglobal alignment) and align that to the re-aligned germline FWR3 region.
    # then we need to convert the aligned position of the junction start to the raw (unaligned) position

    # fr3_sg = abutils.tl.semiglobal_alignment(
    #     query=ab.sequence_oriented,
    #     target=realigned_germ_fr3_sequence,
    #     **ALIGNMENT_PARAMS,
    # )
    # ab.junction_start = fr3_sg.query_end + 1

    fr3_sg = abutils.tl.semiglobal_alignment(
        query=v_sg.aligned_query,
        target=realigned_germ_fr3_sequence,
        **ALIGNMENT_PARAMS,  # normal params are fine, since we already put in any gaps we want
    )
    ab.log("FWR3 SG ALIGNMENT")
    ab.log("ALIGNED QUERY:        ", fr3_sg.aligned_query)
    ab.log("                      ", fr3_sg.alignment_midline)
    ab.log("ALIGNED FWR3 GERMLINE:", fr3_sg.aligned_target)
    ab.log("FWR3 ALIGNMENT END:", fr3_sg.query_end)

    ab.junction_start = (
        get_ungapped_position_from_aligned(
            position=fr3_sg.query_end,
            aligned_sequence=v_sg.aligned_query,
        )
        # + 1  # junction start is the position after the end of the FR3 alignment
        - 2  # back up to the start of the conserved C codon (start of the junction)
    )
    ab.log("JUNCTION START:", ab.junction_start)

    # junction end
    # we need to do similar re-alignment gymnastics with the FR4 region as we did with the FR3 above to ensure that
    # FR4 indels near the start of the J gene don't cause misalignment and incorrect identification of the junction end.
    # TODO:

    if ab.locus in ["IGH", "TRA", "TRD"]:
        germ_fr4_sequence = j_sg.target[-34:]
    else:
        germ_fr4_sequence = j_sg.target[-31:]

    fr4_germ_realign = abutils.tl.semiglobal_alignment(
        query=germ_fr4_sequence,
        target=j_sg.aligned_target,
        **REALIGNMENT_PARAMS,
    )
    realigned_germ_fr4_sequence = fr4_germ_realign.aligned_query.lstrip("-").rstrip("-")

    # to prevent possible misalignments, which we can sometimes get when an insertion duplicates a portion of the J-gene,
    # we align the re-aligned germline with just the portion of the query sequence contains the FR4
    fr4_sg = abutils.tl.semiglobal_alignment(
        query=ab.sequence_oriented[ab.junction_start : ab.j_sequence_end],
        # target=germ_fr4_sequence,
        target=realigned_germ_fr4_sequence,
        **ALIGNMENT_PARAMS,
    )
    # junction sequence includes the W/F, so we need to add 3 because the germline FR4 sequence used for alignment also contains the W/F codon
    ab.junction_end = fr4_sg.query_begin + ab.junction_start + 3
    ab.log("FR4 GERMLINE SEQUENCE:", germ_fr4_sequence)
    ab.log("REALIGNED FR4 GERMLINE SEQUENCE:", realigned_germ_fr4_sequence)
    ab.log("FWR4 SG ALIGNMENT")
    ab.log("ALIGNED QUERY        :", fr4_sg.aligned_query)
    ab.log("                      ", fr4_sg.alignment_midline)
    ab.log("ALIGNED FWR4 GERMLINE:", fr4_sg.aligned_target)
    ab.log("FWR4 ALIGNMENT BEGIN:", fr4_sg.query_begin)
    ab.log("JUNCTION END:", ab.junction_end)
    # junction sequence
    ab.junction = ab.sequence_oriented[ab.junction_start : ab.junction_end]
    ab.junction_aa = abutils.tl.translate(ab.junction)
    ab.log("JUNCTION:", ab.junction)
    ab.log("JUNCTION AA:", ab.junction_aa)

    # CDR3 sequence and length
    ab.cdr3 = ab.junction[3:-3]
    ab.cdr3_aa = ab.junction_aa[1:-1]
    ab.cdr3_length = len(ab.cdr3_aa)
    ab.log("CDR3:", ab.cdr3)
    ab.log("CDR3 AA:", ab.cdr3_aa)
    ab.log("CDR3 LENGTH:", ab.cdr3_length)

    ab.log("\n-----------")
    ab.log(" MUTATIONS")
    ab.log("-----------\n")

    # nucleotide mutations
    ab = annotate_v_mutations(
        aligned_sequence=v_global.aligned_query,
        aligned_germline=v_global.aligned_target,
        gapped_germline=ab.v_germline_gapped,
        germline_start=ab.v_germline_start,
        is_aa=False,
        ab=ab,
        debug=debug,
    )
    ab.log("V MUTATIONS:", ab.v_mutations)
    ab.log("V MUTATION COUNT:", ab.v_mutation_count)

    # amino acid mutations
    ab = annotate_v_mutations(
        aligned_sequence=v_global_aa.aligned_query,
        aligned_germline=v_global_aa.aligned_target,
        gapped_germline=ab.v_germline_gapped_aa,
        germline_start=ab.v_germline_start // 3,
        is_aa=True,
        ab=ab,
        debug=debug,
    )
    ab.log("V MUTATIONS AA:", ab.v_mutations_aa)
    ab.log("V MUTATION COUNT AA:", ab.v_mutation_count_aa)

    # calculate identify (nt and aa)
    ab.v_identity = 1 - ab.v_mutation_count / len(ab.v_germline)
    ab.v_identity_aa = 1 - ab.v_mutation_count_aa / len(ab.v_germline_aa)
    ab.log("V IDENTITY:", ab.v_identity)
    ab.log("V IDENTITY AA:", ab.v_identity_aa)

    ab.log("\n---------")
    ab.log(" REGIONS")
    ab.log("---------\n")

    # V regions
    v_regions = ["fwr1", "cdr1", "fwr2", "cdr2", "fwr3"]
    ab.log("ALIGNED SEQUENCE (GLOBAL):", v_global.aligned_query)
    ab.log("ALIGNED GERMLINE (GLOBAL):", v_global.aligned_target)
    ab.log("ALIGNED SEQUENCE (SEMI-GLOBAL):", v_sg.aligned_query)
    ab.log("ALIGNED GERMLINE (SEMI-GLOBAL):", v_sg.aligned_target)
    ab.log("GAPPED GERMLINE:", ab.v_germline_gapped)
    ab.log("GERMLINE START:", ab.v_germline_start + 1)
    ab.log("ALIGNED QUERY AA (GLOBAL):", v_global_aa.aligned_query)
    ab.log("ALIGNED GERMLINE AA (GLOBAL):", v_global_aa.aligned_target)
    ab.log("ALIGNED QUERY AA (SEMI-GLOBAL):", v_sg_aa.aligned_query)
    ab.log("ALIGNED GERMLINE AA (SEMI-GLOBAL):", v_sg_aa.aligned_target)
    ab.log("GAPPED GERMLINE AA:", ab.v_germline_gapped_aa)
    ab.log("GERMLINE START AA:", ab.v_germline_start // 3 + 1)

    for region in v_regions:
        # nucleotide region
        region_start, region_end, region_sequence = get_region_sequence(
            region,
            # aln=v_sg,
            aln=v_global,
            gapped_germline=ab.v_germline_gapped,
            germline_start=ab.v_germline_start + 1,  # needs to be 1-indexed
            ab=ab,
        )
        setattr(ab, f"{region}", region_sequence)
        ab.log(f"{region.upper()} SEQUENCE:", region_sequence)

        # amino acid region
        _, _, region_sequence_aa = get_region_sequence(
            region,
            # aln=v_sg_aa,
            aln=v_global_aa,
            gapped_germline=ab.v_germline_gapped_aa,
            germline_start=ab.v_germline_start // 3 + 1,  # needs to be 1-indexed
            ab=ab,
            aa=True,
            nt_region_start=region_start,
        )
        setattr(ab, f"{region}_aa", region_sequence_aa)
        ab.log(f"{region.upper()} SEQUENCE AA:", region_sequence_aa)

    # J regions
    fwr4_start = fr4_sg.query_begin
    fwr4_end = fr4_sg.query_end + 1  # python end-slicing is exclusive
    ab.fwr4 = fr4_sg.aligned_query[fwr4_start:fwr4_end]
    ab.fwr4_aa = abutils.tl.translate(ab.fwr4)
    ab.log("FR4 SEQUENCE:", ab.fwr4)
    ab.log("FR4 SEQUENCE AA:", ab.fwr4_aa)

    # CDR3 regions
    ab = identify_cdr3_regions(ab)

    ab.log("\n--------------")
    ab.log(" PRODUCTIVITY")
    ab.log("--------------\n")

    ab = assess_productivity(ab)
    ab.log("STOP CODON:", ab.stop_codon)
    ab.log("PRODUCTIVE:", ab.productive)
    ab.log("PRODUCTIVITY ISSUES:", ab.productivity_issues)

    ab.log("\n------------")
    ab.log(" MASKS")
    ab.log("-----------\n")

    ab.cdr_mask = generate_cdr_mask(ab)
    ab.cdr_mask_aa = generate_cdr_mask(ab, aa=True)
    ab.gene_segment_mask = generate_gene_segment_mask(ab)
    ab.gene_segment_mask_aa = generate_gene_segment_mask(ab, aa=True)
    ab.nongermline_mask = generate_nongermline_mask(ab)
    ab.nongermline_mask_aa = generate_nongermline_mask(ab, aa=True)

    ab.log("SEQUENCE:         ", ab.sequence)
    ab.log("CDR MASK:         ", ab.cdr_mask)
    ab.log("GENE SEGMENT MASK:", ab.gene_segment_mask)
    ab.log("NON-GERMLINE MASK:", ab.nongermline_mask)
    ab.log("SEQUENCE AA:         ", ab.sequence_aa)
    ab.log("CDR MASK AA:         ", ab.cdr_mask_aa)
    ab.log("GENE SEGMENT MASK AA:", ab.gene_segment_mask_aa)
    ab.log("NON-GERMLINE MASK AA:", ab.nongermline_mask_aa)

    return ab


ALIGNMENT_PARAMS = {
    "match": 3,
    "mismatch": -2,
    "gap_open": -35,
    "gap_extend": -1,
}
