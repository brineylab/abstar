"""Regression tests for bugs users reported against legacy abstar (GitHub issues)."""

import abstar
import polars as pl

COMP = str.maketrans("ACGTN", "TGCAN")
REGIONS = ["fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4"]

# #82: legacy dropped VRC26.15's first base (EVQLV -> VQLV) and its last codon (TVSS -> TVS).
VRC26_15 = (
    "GAGGTGCAGCTGGTGGAGTCGGGGGGAGGCGTGGTCCAGCCTGGGAGGTCCCTGAGACTCTCCTGTGCAGCCTCTCAATTCACCTTTAGTAATTATGGCA"
    "TGCACTGGGTCCGCCAGGCTCCCGGCAAGGGGCTGGAATGGGTGGCAGGCATCTCATTTGTTGGGATCAGGAAATATTATGCAGATTCCGTGAAGGGCCG"
    "ATTCACCATCTCCAGAGACAACTCCAAGAACACCCTGTATCTGCAAATGAACAGCCTGCGAGCTGCGGACACGGCTCTTTACTTCTGTGCGAGAGATCTG"
    "CGAGAATACGAATGTGAATTGTGGGCGTCTGACTATTACGATTTTGGGAAACCTCAGCCTTGCGAGGACCCACGGGGCGTGGTTGGAACCTCAGATATCT"
    "GGGGCCAAGGGACAATGGTCACCGTCTCTTCA"
)

# #123: reverse-strand IGHV1-2 reads from one clone; Sequence-2 starts inside CDR1. Legacy reported
# stop codons in CDR1/CDR2 for it.
SEQ_1 = (
    "GATGGGGCTGTTGTTTCAGCTGAGGAGACGGTGACCAGGGTTCCCAGGCCCCAGTAGTCAAAGTACCCCCCGTAGCCACTCCAATTCCCTCTCGCACAGT"
    "AATACACGGCCGTGTCGTCAGATCTCAGCCTGCTCAGTTCCATGTAGGCTGTGTTGATGGACGTGTCCCTGGTCATGGTGACCCTGCCCTGAAACTTCTG"
    "TGCATAGTTTGTGCCACCACTGTTAGGGTTGATCCATCCCATCCACTCAAGCCCTTGTCCAGGGGCCTGTCGCACCCAGTGCATATAGTAGCCGGTGAAG"
    "GTGTATCCAGAAGCCTTGCAGGAGACCTTAACTGAGGCCCCAGGCTTCTTCACCTCAGCCCCAGACTGCACCAGCTGCACCTGGGAGTGGGCTCCTGTGG"
    "CTGCTGCCACCAAGAAGAGGATCCTCCAGGTCCAGTCCATGGTGAGGAACGGAGCTCTCAGGGGATTCTCTAGAGGACAGATGTGGTTGCTGGGTGATGC"
    "TCCCCC"
)
SEQ_2 = (
    "GATGGGGCTGTTGTTTCAGCTGAGGAGACGGTGACCAGGGTTCCCAGGCCCCAGTAGTCAAAGTACCCCCCGTAGCCACTCCAATTCCCTCTCGCACAGT"
    "AATACACGGCCGTGTCGTCAGATCTCAGCCTGCTCAGTTCCATGTAGGCTGTGTTGATGGACGTGTCCCTGGTCATGGTGACCCTGCCCTGAAACTTCTG"
    "TGCATAGTTTGTGCCACCACTGTTAGGGTTGATCCATCCCATCCACTCAAGCCCTTGTCCAGGGGCCTGTCGCACCCAGTGCATATAGTAGCCGGTGAAG"
    "CCC"
)

# #112: constant-region-only reads (IgE CH1, reverse strand; kappa constant, reverse strand).
CONSTANT_ONLY = [
    "GTCCCAGGTCACCATCACCGGCTCCGGGAAGTAGCCCGTGGCCAGGCAGCCCAGAGTCACGGAGGTGGCATTGGAGGGAATG",
    "GTCCCAGGTCACCATCACCGGCTCCGGGAAGTAGCCCGTGGCCAGGCAGCCCAGAGTCACGGAGGTGGCATTGGAGGGAATGTTTTT",
    "GCGTTATCCACCTTCCACTGTACTTTGGCCTCTCTGGGATAGAAGTTATTCAGCAGGCACACAACAGAGGCAGTTCCAGATTTCAACTGC",
]


def test_82_vrc26_15_keeps_first_and_last_codons():
    r = abstar.run([("VRC26.15", VRC26_15)]).row(0, named=True)
    assert r["v_sequence_start"] == 1
    assert r["sequence_alignment_aa"].replace(".", "").startswith("EVQLVESGG")
    assert r["fwr4_aa"] == "WGQGTMVTVSS"
    assert r["sequence_vdj_aa"].endswith("TVSS")


def _oriented(df: pl.DataFrame) -> pl.Series:
    return pl.Series([s[::-1].translate(COMP) if rc else s for s, rc in zip(df["sequence"], df["rev_comp"])])


def test_84_region_fields_agree():
    # #84: legacy's region lengths, junction breakdown and junction disagreed. Every region and the
    # junction must be the oriented read's bases at their coordinates, and CDR3 the junction minus
    # its two flanking codons. Reverse-complemented copies check the oriented coordinates.
    gen = abstar.generate(400, "IGH", seed=84, preset="bnab")
    reads = [(f"f{i}", s) for i, s in enumerate(gen["sequence"])]
    reads += [(f"r{i}", s[::-1].translate(COMP)) for i, s in enumerate(gen["sequence"][:200])]
    df = abstar.run(reads).filter(pl.col("junction").is_not_null())
    assert df.height > 550
    seq = _oriented(df)
    for row, s in zip(df.iter_rows(named=True), seq):
        for reg in REGIONS:
            if row[reg]:
                assert s[row[f"{reg}_start"] - 1 : row[f"{reg}_end"]] == row[reg], (row["sequence_id"], reg)
                assert row[f"{reg}_length"] == len(row[f"{reg}_aa"]), (row["sequence_id"], reg)
        j = row["junction"]
        assert j in s and row["junction_length"] == len(j), row["sequence_id"]
        assert row["cdr3"] == j[3:-3], row["sequence_id"]


def test_112_constant_only_reads_are_not_annotated():
    df = abstar.run(CONSTANT_ONLY)
    assert df["failure_reason"].to_list() == ["no_v_hit"] * len(CONSTANT_ONLY)


def test_123_cdr1_starting_read_matches_its_clone_mate():
    for db in ["human", "human+c57bl6"]:
        df = abstar.run([("s1", SEQ_1), ("s2", SEQ_2)], germline_db=db)
        s1, s2 = df.row(0, named=True), df.row(1, named=True)
        for r in (s1, s2):
            assert r["rev_comp"] and r["productive"] and not r["stop_codon"], (db, r["sequence_id"])
            assert r["v_call"] == "IGHV1-2*02" and r["j_call"] == "IGHJ4*02", db
        assert s1["cdr1_aa"] == "GYTFTGYY" and s1["cdr2_aa"] == "INPNSGGT", db
        # Sequence-2 begins three bases before CDR1's TTC; its CDR2 and junction match Sequence-1's.
        assert s2["fwr1"] is None or s2["fwr1"] == "", db
        assert s2["cdr1_aa"].endswith("TGYY") and "*" not in s2["cdr1_aa"], db
        assert (s2["cdr2_aa"], s2["junction_aa"]) == (s1["cdr2_aa"], s1["junction_aa"]), db
