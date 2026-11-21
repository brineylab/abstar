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
    # The annotator ended V at 9 (base 10 went to N), so the IMGT10 mutation is outside the shared span.
    calls = pl.DataFrame({
        "sequence_id": ["a", "b", "c"], "v_call": ["IGHV1-18*01"] * 3, "v_germline_start": [1] * 3, "v_germline_end": [9] * 3,
        "v_mutations": [["IMGT3:A>C"], ["IMGT3:A>C"], []], "v_mutation_count": [1, 1, 0],
    })
    g = abstar.grade(_truth(), calls, germline_db="human")
    assert _row(g, "v_mutations")["accuracy"] == 1.0
    assert _row(g, "v_mutation_count")["accuracy"] == 1.0
    # a real disagreement inside the span is wrong
    calls = calls.with_columns(pl.Series("v_mutations", [["IMGT4:A>C"], ["IMGT3:A>C"], []]))
    assert _row(abstar.grade(_truth(), calls), "v_mutations")["correct"] == 1


def test_subsets():
    g = abstar.grade(_truth(), pl.DataFrame({"sequence_id": ["a", "b", "c"], "v_call": ["IGHV1-18*01"] * 3}))
    assert [_row(g, "v_call", s)["n"] for s in ("all", "high_mutation", "has_indel")] == [2, 1, 1]
    assert set(g.columns) == {"field", "subset", "n", "correct", "accuracy"}
    assert _row(g, "j_call")["n"] == 0  # absent field reported, not a crash


def test_grade_round_trip():
    truth = abstar.generate(500, "IGH", seed=1)
    g = abstar.grade(truth, truth).filter(pl.col("n") > 0)
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
