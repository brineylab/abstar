import polars as pl
import pytest

import abstar


def _truth(**over):
    base = {
        "sequence_id": ["a", "b", "c"],
        "v_call": ["IGHV1-18*01", "IGHV1-18*01", "IGHV1-18*01"],
        "v_call_valid": [["IGHV1-18*01", "IGHV1-18*03"]] * 3,
        "unrecoverable_v": [False, False, True],
        "d_call": [None, "D*01", "D*01"],
        "d_call_valid": [[], ["D*01"], ["D*01"]],
        "d_call_null_valid": [True, False, False],
        "unrecoverable_d": [False, False, False],
        "v_sequence_end": [10, 10, 10],
        "v_sequence_end_valid": [[10, 11]] * 3,
        "unrecoverable_j": [False] * 3,
        "high_mutation": [True, False, False],
        "has_indel": [False, True, True],
        "v_germline_start": [1, 1, 1],
        "v_germline_end": [10, 10, 10],
        "v_mutations": [["IMGT3:A>C"], ["IMGT3:A>C", "IMGT10:G>T"], []],
    }
    base.update(over)
    return pl.DataFrame(base)


def _row(df, field, subset="all"):
    return df.filter((pl.col("field") == field) & (pl.col("subset") == subset)).row(0, named=True)


def test_call_in_valid_set():
    g = abstar.grade(_truth(), pl.DataFrame({"sequence_id": ["a", "b", "c"], "v_call": ["IGHV1-18*03", "IGHV1-18*09", "IGHV1-18*09"]}))
    r = _row(g, "v_call")
    assert (r["n"], r["correct"], r["accuracy"]) == (2, 1, 0.5)  # row c is unrecoverable


def test_null_d_when_allowed():
    g = abstar.grade(_truth(), pl.DataFrame({"sequence_id": ["a", "b", "c"], "d_call": [None, None, "D*01"]}))
    r = _row(g, "d_call")
    assert (r["n"], r["correct"]) == (3, 2)  # a: null allowed; b: null not allowed


def test_boundary_lists():
    g = abstar.grade(_truth(), pl.DataFrame({"sequence_id": ["a", "b", "c"], "v_sequence_end": [11, 12, 10]}))
    assert (_row(g, "v_sequence_end")["n"], _row(g, "v_sequence_end")["correct"]) == (2, 1)


def test_shared_span_mutations():
    # The annotator ends V one base later (a valid cut) and reports an extra mutation there,
    # outside truth's span: the shared span is 1..10 and the lists match on it.
    calls = pl.DataFrame({
        "sequence_id": ["a", "b", "c"], "v_call": ["IGHV1-18*01"] * 3, "v_germline_start": [1] * 3, "v_germline_end": [11] * 3,
        "v_sequence_end": [11] * 3,
        "v_mutations": [["IMGT3:A>C", "IMGT11:A>T"], ["IMGT3:A>C", "IMGT10:G>T", "IMGT11:A>T"], ["IMGT11:A>T"]],
        "v_mutation_count": [2, 3, 1],
    })
    g = abstar.grade(_truth(), calls, germline_db="human")
    assert _row(g, "v_mutations")["accuracy"] == 1.0
    assert _row(g, "v_mutation_count")["accuracy"] == 1.0
    # a real disagreement inside the span is wrong
    bad = calls.with_columns(pl.Series("v_mutations", [["IMGT4:A>C", "IMGT11:A>T"], calls["v_mutations"][1].to_list(), ["IMGT11:A>T"]]))
    assert _row(abstar.grade(_truth(), bad), "v_mutations")["correct"] == 1


def test_clipped_span_with_dropped_mutations_is_wrong():
    # V cut 20 bases early, outside every valid cut, with the mutations dropped to match
    calls = pl.DataFrame({
        "sequence_id": ["a", "b", "c"], "v_call": ["IGHV1-18*01"] * 3, "v_germline_start": [1] * 3, "v_germline_end": [4] * 3,
        "v_sequence_end": [4] * 3, "v_mutations": [["IMGT3:A>C"], ["IMGT3:A>C"], []],
    })
    assert _row(abstar.grade(_truth(), calls), "v_mutations")["correct"] == 0


def test_garbage_reported_count_is_wrong():
    calls = pl.DataFrame({
        "sequence_id": ["a", "b", "c"], "v_call": ["IGHV1-18*01"] * 3, "v_germline_start": [1] * 3, "v_germline_end": [10] * 3,
        "v_sequence_end": [10] * 3, "v_mutations": [["IMGT3:A>C"], ["IMGT3:A>C", "IMGT10:G>T"], []],
        "v_mutation_count": [999, 2, 0], "v_mutation_frequency": [0.77, 0.2, 0.0],
    })
    g = abstar.grade(_truth(v_mutation_count=[1, 2, 0]), calls)
    assert _row(g, "v_mutations")["accuracy"] == 1.0
    assert (_row(g, "v_mutation_count")["n"], _row(g, "v_mutation_count")["correct"]) == (2, 1)


def _mask_truth(**over):
    base = {
        "sequence_id": ["a"], "v_call": ["IGHV1-18*01"], "j_call": ["J*01"], "d_call": [None],
        "unrecoverable_v": [False], "unrecoverable_j": [False], "v_sequence_start": [1], "v_germline_start": [1],
        "v_sequence_end": [8], "v_sequence_end_valid": [[7, 8]],
        "gene_segment_mask": ["VVVVVVVVNN"], "gene_segment_mask_aa": ["VVVV"],
    }
    base.update(over)
    return pl.DataFrame(base)


def _mask_calls(mask, aa):
    return pl.DataFrame({"sequence_id": ["a"], "v_call": ["IGHV1-18*01"], "j_call": ["J*01"], "gene_segment_mask": [mask], "gene_segment_mask_aa": [aa]})


def test_mask_difference_at_a_valid_cut():
    # cuts {7, 8} leave base 8 open
    g = abstar.grade(_mask_truth(), _mask_calls("VVVVVVVNNN", "VVVV"))
    assert _row(g, "gene_segment_mask")["correct"] == 1
    # truth outside the valid list: valid [6], true cut 8 -> bases 7-8 open
    g = abstar.grade(_mask_truth(v_sequence_end_valid=[[6]]), _mask_calls("VVVVVVNNNN", "VVVV"))
    assert _row(g, "gene_segment_mask")["correct"] == 1
    g = abstar.grade(_mask_truth(), _mask_calls("VVVVVNNNNN", "VVVV"))
    assert _row(g, "gene_segment_mask")["correct"] == 0


def test_null_d_boundaries():
    truth = pl.DataFrame({
        "sequence_id": ["a", "b"], "d_call_null_valid": [True, False], "unrecoverable_d": [False] * 2,
        "d_sequence_start": [5, 5], "d_sequence_start_valid": [[5, 6]] * 2,
        "d_sequence_end": [9, 9], "d_sequence_end_valid": [[9]] * 2,
    })
    g = abstar.grade(truth, pl.DataFrame({"sequence_id": ["a", "b"], "d_sequence_start": [None, None], "d_sequence_end": [None, None]}))
    assert _row(g, "d_sequence_start")["correct"] == 1 and _row(g, "d_sequence_end")["correct"] == 1


def test_aa_mask_phase():
    # first V base is IMGT 2, so codons start at read base 3: their 3rd bases are reads 5, 8, 11, ...
    truth = _mask_truth(v_germline_start=[2], v_sequence_end=[8], gene_segment_mask_aa=["VVVV"])
    legit = abstar.grade(truth, _mask_calls("VVVVVVVNNN", "VNVV"))  # codon 1 = read 8, the open base
    wrong = abstar.grade(truth, _mask_calls("VVVVVVVNNN", "VVNV"))  # codon 2 = read 11, determined
    assert _row(legit, "gene_segment_mask_aa")["correct"] == 1
    assert _row(wrong, "gene_segment_mask_aa")["correct"] == 0


def test_import_keeps_grade_callable():
    import abstar._grade  # noqa: F401

    assert callable(abstar.grade)


def test_subsets():
    g = abstar.grade(_truth(), pl.DataFrame({"sequence_id": ["a", "b", "c"], "v_call": ["IGHV1-18*01"] * 3}))
    assert [_row(g, "v_call", s)["n"] for s in ("all", "high_mutation", "has_indel")] == [2, 1, 1]
    assert set(g.columns) == {"field", "subset", "n", "correct", "accuracy"}
    assert _row(g, "j_call")["n"] == 0  # absent field reported, not a crash


def test_grade_round_trip():
    truth = abstar.generate(500, "IGH", seed=1)
    # Truth carries the generative allele even when the read holds no base of the segment (empty
    # _valid); there only a null call is correct, so those rows are nulled in the "calls".
    calls = truth
    for seg in "vdjc":
        empty = pl.col(f"{seg}_call_valid").list.len() == 0
        print(f"{seg}: {truth.filter(empty).height} rows with an empty valid set excluded")
        calls = calls.with_columns(pl.when(empty).then(None).otherwise(pl.col(f"{seg}_call")).alias(f"{seg}_call"))
    g = abstar.grade(truth, calls).filter(pl.col("n") > 0)
    boundary = {"v_sequence_end", "d_sequence_start", "d_sequence_end", "j_sequence_start"}
    assert g.height > 100
    for r in g.filter(pl.col("subset") == "all").iter_rows(named=True):
        if r["field"] in boundary:
            print(r["field"], r["accuracy"])
            assert r["accuracy"] >= 0.7
        else:
            assert r["accuracy"] == 1.0, r
    bad = g.filter((pl.col("accuracy") < 1) & ~pl.col("field").is_in(list(boundary)))
    assert bad.height == 0, bad
