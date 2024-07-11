# Copyright (c) 2024 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

import abutils

from ..core.antibody import Antibody
from .germline import (
    get_germline,
    process_dgene_alignment,
    process_jgene_alignment,
    process_vgene_alignment,
    realign_germline,
    reassign_dgene,
)
from .indels import annotate_v_deletions, annotate_v_insertions
from .mutations import annotate_mutations
from .regions import get_region_sequence


def annotate(ab: Antibody):
    # log information from the assigner
    ab.log("=" * (len(ab.sequence_id) + 15))
    ab.log(" SEQUENCE ID:", ab.sequence_id)
    ab.log("=" * (len(ab.sequence_id) + 15) + "\n")
    ab.log(f">{ab.sequence_id}\n{ab.sequence_input}\n")
    ab.log(f"{ab.v_call} | {ab.d_call} | {ab.j_call}")
    ab.log("REV COMP:", ab.rev_comp)

    # germline database, locus, and receptor type
    ab.germline_database = "human"  # TODO: replace with the actual germdb variable name
    ab.locus = ab.v_call[:3].upper()
    ab.receptor_type = "bcr" if ab.locus[:2] == "IG" else "tcr"
    ab.log("GERMLINE DATABASE:", ab.germline_database)
    ab.log("LOCUS:", ab.locus)

    # orient the input sequence
    if ab.rev_comp:
        ab.sequence_oriented = abutils.tl.reverse_complement(ab.sequence_input)
    else:
        ab.sequence_oriented = ab.sequence_input
    ab.log("SEQUENCE ORIENTED:", ab.sequence_oriented)

    # check germline calls for species name
    # (mainly for mixed species databases like humouse)
    if len(vsplit := ab.v_call.split("__")) > 1:
        ab.v_call, ab.species = vsplit
    else:
        ab.species = ab.germline_database
    if len(jsplit := ab.j_call.split("__")) > 1:
        ab.j_call, _ = jsplit
    if ab.d_call is not None:
        if len(dsplit := ab.d_call.split("__")) > 1:
            ab.d_call, _ = dsplit

    # get genes from calls
    ab.v_gene = ab.v_call.split("*")[0]
    ab.j_gene = ab.j_call.split("*")[0]
    if ab.d_call is not None:
        ab.d_gene = ab.d_call.split("*")[0]

    ab.log("\n--------")
    ab.log(" V GENE")
    ab.log("--------\n")

    # v-gene realignment
    v_sg, v_loc = realign_germline(
        sequence=ab.sequence_oriented,
        germline_name=ab.v_call,
        germdb_name=ab.germline_database,
        imgt_gapped=False,
        semiglobal_aln_params=ALIGNMENT_PARAMS,
        local_aln_params=ALIGNMENT_PARAMS,
    )
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
    gapped_v_germline = get_germline(
        ab.v_call, ab.germline_database, imgt_gapped=True
    ).sequence
    gapped_v_germline_aa = abutils.tl.translate(gapped_v_germline, allow_dots=True)
    ab.log("GAPPED GERMLINE:", gapped_v_germline)
    ab.log("GAPPED GERMLINE AA:", gapped_v_germline_aa)

    # insertions
    if "-" in v_loc.aligned_target:
        ab = annotate_v_insertions(
            aligned_sequence=v_global.aligned_query,
            aligned_germline=v_global.aligned_target,
            gapped_germline=gapped_v_germline,
            germline_start=ab.v_germline_start,
            ab=ab,
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
        ab = annotate_v_deletions(
            aligned_sequence=v_loc.aligned_query,
            aligned_germline=v_loc.aligned_target,
            gapped_germline=gapped_v_germline,
            germline_start=ab.v_germline_start,
            ab=ab,
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
    j_sg, j_loc = realign_germline(
        sequence=jquery,
        germline_name=ab.j_call,
        germdb_name=ab.germline_database,
        imgt_gapped=False,
        semiglobal_aln_params=ALIGNMENT_PARAMS,
        local_aln_params=ALIGNMENT_PARAMS,
    )
    ab.log("SEMIGLOBAL ALIGNMENT:")
    ab.log(f"     QUERY: {j_sg.aligned_query}")
    ab.log(f"            {j_sg.alignment_midline}")
    ab.log(f"  GERMLINE: {j_sg.aligned_target}")
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
            germline_name=ab.d_call,
            germdb_name=ab.germline_database,
            imgt_gapped=False,
            skip_semiglobal=True,
            local_aln_params=ALIGNMENT_PARAMS,
        )
    # if not, we can try again using local pairwise alignment (for IGH/TRA/TRD chains only)
    elif len(dquery) >= 5 and ab.locus in [
        "IGH",
        "TRA",
        "TRD",
    ]:  # TODO: make the minimum d-gene alignment length a user-controllable parameter?
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

    ab.log("\n--------------")
    ab.log(" VDJ ASSEMBLY")
    ab.log("--------------\n")

    # assemble the full V(D)J and germline sequences
    ab.sequence = ab.v_sequence + ab.np1
    ab.germline = ab.v_germline + ab.np1
    if ab.d_sequence is not None:
        ab.sequence += ab.d_sequence
        ab.germline += ab.d_germline
    if ab.np2 is not None:
        ab.sequence += ab.np2
        ab.germline += ab.np2
    ab.sequence += ab.j_sequence
    ab.germline += ab.j_germline
    ab.log("SEQUENCE:", ab.sequence)
    ab.log("GERMLINE:", ab.germline)

    # translated sequences
    ab.sequence_aa = abutils.tl.translate(ab.sequence, frame=ab.frame)
    ab.germline_aa = abutils.tl.translate(ab.germline, frame=ab.frame)
    ab.log("SEQUENCE AA:", ab.sequence_aa)
    ab.log("GERMLINE AA:", ab.germline_aa)

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
    germ_fr3_sequence = gapped_v_germline[196:309].replace(".", "")
    fr3_sg = abutils.tl.semiglobal_alignment(
        query=ab.sequence_oriented,
        target=germ_fr3_sequence,
        **ALIGNMENT_PARAMS,
    )
    ab.junction_start = fr3_sg.query_end + 1

    # junction end
    germ_fr4_sequence = j_sg.target[-34:]
    fr4_sg = abutils.tl.semiglobal_alignment(
        query=ab.sequence_oriented,
        target=germ_fr4_sequence,
        **ALIGNMENT_PARAMS,
    )
    ab.junction_end = fr4_sg.query_begin + 3  # germ_fr4_sequence includes the W/F

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
    ab = annotate_mutations(
        aligned_sequence=v_global.aligned_query,
        aligned_germline=v_global.aligned_target,
        gapped_germline=gapped_v_germline,
        germline_start=ab.v_germline_start,
        is_aa=False,
        ab=ab,
    )
    ab.log("V MUTATIONS:", ab.v_mutations)
    ab.log("V MUTATION COUNT:", ab.v_mutation_count)

    # amino acid mutations
    ab = annotate_mutations(
        aligned_sequence=v_global_aa.aligned_query,
        aligned_germline=v_global_aa.aligned_target,
        gapped_germline=gapped_v_germline_aa,
        germline_start=ab.v_germline_start // 3,
        is_aa=True,
        ab=ab,
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
    for region in v_regions:
        # nucleotide region
        region_sequence = get_region_sequence(
            region,
            aligned_sequence=v_global.aligned_query,
            aligned_germline=v_global.aligned_target,
            gapped_germline=gapped_v_germline,
            ab=ab,
            germline_start=ab.v_germline_start,
        )
        setattr(ab, f"{region}", region_sequence)
        ab.log(f"{region.upper()} SEQUENCE:", region_sequence)
        # amino acid region
        region_sequence_aa = get_region_sequence(
            region,
            aligned_sequence=v_global_aa.aligned_query,
            aligned_germline=v_global_aa.aligned_target,
            gapped_germline=gapped_v_germline_aa,
            ab=ab,
            germline_start=ab.v_germline_start // 3 + 1,
            aa=True,
        )
        setattr(ab, f"{region}_aa", region_sequence_aa)
        ab.log(f"{region.upper()} SEQUENCE AA:", region_sequence_aa)

    # J regions
    fwr4_start = fr4_sg.query_begin
    fwr4_end = fr4_sg.query_end
    ab.fwr4 = fr4_sg.aligned_query[fwr4_start:fwr4_end]
    ab.fwr4_aa = abutils.tl.translate(ab.fwr4)
    ab.log("FR4 SEQUENCE", ab.fwr4)
    ab.log("FR4 SEQUENCE AA", ab.fwr4_aa)

    ab.log("\n--------------")
    ab.log(" PRODUCTIVITY")
    ab.log("--------------\n")

    # scan for issues
    if "*" in ab.sequence_aa:
        ab.productivity_issues.append("stop codon(s)")
    if "N" in ab.sequence:
        ab.productivity_issues.append("ambiguous nucleotide(s)")

    # flag sequences with issues as non-productive
    if ab.productivity_issues:
        ab.productive = False
    ab.productivity_issues = "|".join(ab.productivity_issues)
    ab.log("PRODUCTIVE:", ab.productive)
    ab.log("PRODUCTIVITY ISSUES:", ab.productivity_issues)

    return ab


ALIGNMENT_PARAMS = {
    "match": 3,
    "mismatch": -2,
    "gap_open": -35,
    "gap_extend": -1,
}
