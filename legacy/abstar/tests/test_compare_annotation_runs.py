"""Ordinal-based, exact comparison of corpus reruns."""
import importlib.util
import json
from pathlib import Path

import polars as pl
import pytest

_spec = importlib.util.spec_from_file_location(
    "compare_annotation_runs", Path(__file__).parents[2] / "scripts/compare_annotation_runs.py"
)
compare = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(compare)


@pytest.fixture
def corpus(tmp_path):
    baseline, candidate, fasta = [tmp_path / name for name in ("old", "new", "fasta")]
    fasta.mkdir()
    # The first two records deliberately have identical IDs AND sequences.
    (fasta / "sample.fasta").write_text(">dup\nAAA\n>dup\nAAA\n>10E8\nCCC\n")
    rows = {"sequence_id": ["dup", "dup", "10E8"], "sequence": ["AAA", "AAA", "CCC"],
            "v_call": ["IGHV1", "IGHV2", "IGHV3"], "junction_start": [None, 1, 2],
            "identity": [float("nan"), 0.9, 1.0]}
    frame = pl.DataFrame(rows)
    for root, df in ((baseline, frame[[0, 2]]), (candidate, frame)):
        (root / "parquet").mkdir(parents=True)
        (root / "logs").mkdir()
        df.write_parquet(root / "parquet/sample.parquet")
    header = "input_file\tsample\trow_id\tsequence_id\n"
    (baseline / "logs/failures.tsv").write_text(header + f"{fasta / 'sample.fasta'}\tsample\tabstar_0_1\tdup\n")
    (candidate / "logs/failures.tsv").write_text(header)
    expected = tmp_path / "expected.json"
    expected.write_text(json.dumps({"records": [{"source_file": "sample.fasta", "record_ordinal": 1,
        "sample_ordinal": 0, "expected": {"v_call": "IGHV2", "junction_start": 1}}]}))
    return baseline, candidate, fasta, expected


def test_inserted_recovery_and_duplicate_ids(corpus):
    result = compare.compare_runs(*corpus)
    assert result["ok"]
    assert result["totals"] == {"input_records": 3, "baseline_successes": 2, "candidate_successes": 3, "recovered": 1}


@pytest.mark.parametrize("change", ["biology", "reorder", "missing", "extra", "null", "dtype", "recovery", "sequence"])
def test_rejects_changes(corpus, change):
    path = corpus[1] / "parquet/sample.parquet"
    df = pl.read_parquet(path)
    if change == "biology":
        df = df.with_columns(pl.when(pl.col("v_call") == "IGHV3").then(pl.lit("IGHV4")).otherwise(pl.col("v_call")).alias("v_call"))
    elif change == "reorder":
        df = df[[1, 0, 2]]
    elif change == "missing":
        df = df.head(2)
    elif change == "extra":
        df = pl.concat([df, df.tail(1)])
    elif change == "null":
        df = df.with_columns(pl.col("junction_start").fill_null(0))
    elif change == "dtype":
        df = df.with_columns(pl.col("junction_start").cast(pl.Int32))
    elif change == "recovery":
        df = df.with_columns(pl.when(pl.col("v_call") == "IGHV2").then(pl.lit("IGHV4")).otherwise(pl.col("v_call")).alias("v_call"))
    elif change == "sequence":
        df = df.with_columns(pl.lit("GGG").alias("sequence"))
    df.write_parquet(path)
    assert not compare.compare_runs(*corpus)["ok"]


def test_new_failure_and_incomplete_fixture(corpus):
    (corpus[1] / "logs/failures.tsv").write_text((corpus[0] / "logs/failures.tsv").read_text())
    assert not compare.compare_runs(*corpus)["ok"]
    corpus[3].write_text('{"records": []}')
    assert not compare.compare_runs(*corpus)["ok"]


def test_report_cannot_overwrite_inputs(corpus):
    with pytest.raises(ValueError, match="outside"):
        compare.main(["--baseline", str(corpus[0]), "--candidate", str(corpus[1]),
                      "--fasta-dir", str(corpus[2]), "--expected", str(corpus[3]),
                      "--report", str(corpus[2] / "report.json")])


def test_cli_writes_report_and_returns_failure_status(corpus, tmp_path):
    args = ["--baseline", str(corpus[0]), "--candidate", str(corpus[1]),
            "--fasta-dir", str(corpus[2]), "--expected", str(corpus[3])]
    report = tmp_path / "report.json"
    assert compare.main(args + ["--report", str(report)]) == 0
    assert json.loads(report.read_text())["ok"]
    with pytest.raises(FileExistsError):
        compare.main(args + ["--report", str(report)])
    (corpus[1] / "parquet/sample.parquet").unlink()
    bad_report = tmp_path / "bad.json"
    assert compare.main(args + ["--report", str(bad_report)]) == 1
    assert not json.loads(bad_report.read_text())["ok"]


def test_rejects_duplicate_baseline_failure(corpus):
    path = corpus[0] / "logs/failures.tsv"
    text = path.read_text()
    path.write_text(text + text.splitlines()[1] + "\n")
    with pytest.raises(ValueError, match="Duplicate baseline failure"):
        compare.compare_runs(*corpus)


def test_rejects_baseline_failure_id_mismatch(corpus):
    path = corpus[0] / "logs/failures.tsv"
    path.write_text(path.read_text().replace("\tdup\n", "\tother\n"))
    assert not compare.compare_runs(*corpus)["ok"]
