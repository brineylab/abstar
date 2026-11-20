import gzip
import tomllib
from pathlib import Path

import airr
import polars as pl
import pytest

import abstar
from abstar import _core
from abstar.errors import AbstarError

FIELDS = [f["name"] for f in tomllib.loads(_core.fields_toml())["field"]]


def _fa(path: Path, recs):
    path.parent.mkdir(parents=True, exist_ok=True)
    text = "".join(f">{i}\n{s}\n" for i, s in recs)
    if path.suffix == ".gz":
        path.write_bytes(gzip.compress(text.encode()))
    else:
        path.write_text(text)
    return path


@pytest.fixture
def indir(tmp_path):
    d = tmp_path / "in"
    _fa(d / "x.fasta", [("10E8", "CAGGTGCAG"), ("0012", "ACGT")])
    _fa(d / "sub" / "y.fa.gz", [("r1", "CAGGTGCAG")])
    (d / "notes.txt").write_text("ignore")
    return d


def err_code(fn, *a, **k):
    with pytest.raises(AbstarError) as e:
        fn(*a, **k)
    return e.value.code


def test_single_sequence_dataframe():
    df = abstar.run("CAGGTGCAG")
    assert df.shape[0] == 1
    assert df["sequence_id"][0] == "1"
    assert df["annotation_status"][0] == "unassigned"
    assert df.columns == FIELDS


def test_zero_inputs():
    df = abstar.run([])
    assert df.shape[0] == 0
    assert df.columns == FIELDS


def test_missing_path_not_sequence():
    assert err_code(abstar.run, "sample1.fasta") == "INPUT_NOT_FOUND"
    assert err_code(abstar.run, "data/x") == "INPUT_NOT_FOUND"


def test_tuples_and_generator():
    df = abstar.run((p for p in [("a", "ACGT"), ("b", "")]))
    assert df["sequence_id"].to_list() == ["a", "b"]
    assert df["annotation_status"].to_list() == ["unassigned", "invalid_input"]


def test_directory_input_and_source_file(indir):
    # source_file is the path as joined from the given dir (not resolved)
    df = abstar.run(indir)
    assert df["source_file"].unique(maintain_order=True).to_list() == [
        str(indir / "sub" / "y.fa.gz"),
        str(indir / "x.fasta"),
    ]
    assert df["sequence_id"].to_list() == ["r1", "10E8", "0012"]


def test_files_summary(indir, tmp_path):
    out = tmp_path / "out"
    s = abstar.run(indir, out, formats=["airr", "parquet"])
    assert isinstance(s, abstar.RunSummary)
    names = sorted(p.name for i in s.inputs for p in i.outputs.values())
    assert names == ["x.parquet", "x.tsv", "y.parquet", "y.tsv"]
    assert s.run_json == out / "run.json" and s.run_json.is_file()
    assert set(s.counts) == {"annotated", "unassigned", "invalid_input", "error"}
    assert sum(s.counts.values()) == 3
    assert all(isinstance(p, Path) for i in s.inputs for p in i.outputs.values())


def test_airr_validation(indir, tmp_path):
    out = tmp_path / "out"
    abstar.run(indir, out)
    assert airr.validate_rearrangement(str(out / "x.tsv")) is True


def test_ids_preserved_parquet(indir, tmp_path):
    out = tmp_path / "out"
    abstar.run(indir / "x.fasta", out, formats="parquet")
    df = pl.read_parquet(out / "x.parquet")
    assert df.schema["sequence_id"] == pl.Utf8
    assert df["sequence_id"].to_list() == ["10E8", "0012"]


def test_mixed_inputs_rejected(indir):
    assert err_code(abstar.run, [indir / "x.fasta", "ACGT"]) == "INVALID_ARGUMENT"


@pytest.mark.parametrize("kw", [{"formats": []}, {"formats": "csv"}, {"threads": 0}, {"threads": -1}])
def test_bad_arguments(kw):
    assert err_code(abstar.run, "ACGT", **kw) == "INVALID_ARGUMENT"


def test_sequences_list_ids():
    df = abstar.run(["ACGT", "CAGGTGCAG"])
    assert df["sequence_id"].to_list() == ["1", "2"]


def test_sequences_to_files(tmp_path):
    out = tmp_path / "out"
    s = abstar.run(["ACGT", "GGGG"], out, formats=["airr", "parquet"])
    assert len(s.inputs) == 1
    assert sorted(p.name for p in s.inputs[0].outputs.values()) == ["records.parquet", "records.tsv"]
    assert sum(s.counts.values()) == 2
