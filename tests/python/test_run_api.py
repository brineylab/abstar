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
# IGHV1-2*02 + CC + IGHD3-3*01 + GG + IGHJ4*02 (human): in frame, productive
IGH_READ = (
    "CAGGTGCAGCTGGTGCAGTCTGGGGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGATACACCTTCACCGGCTACTATATGCAC"
    "TGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGATGGATCAACCCTAACAGTGGTGGCACAAACTATGCACAGAAGTTTCAGGGCAGGGTCACCATG"
    "ACCAGGGACACGTCCATCAGCACAGCCTACATGGAGCTGAGCAGGCTGAGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGACCGTATTACGATTTTTGGAG"
    "TGGTTATTATACCGGACTACTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG"
)


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


def test_real_igh_read_is_annotated():
    df = abstar.run([IGH_READ, IGH_READ[::-1].translate(str.maketrans("ACGT", "TGCA"))])
    assert df["annotation_status"].to_list() == ["annotated", "annotated"]
    assert df["failure_reason"].to_list() == [None, None]
    assert df["rev_comp"].to_list() == [False, True]
    assert all(v.startswith("IGHV1-2") for v in df["v_call"])
    assert df["j_call"].to_list() == ["IGHJ4*02", "IGHJ4*02"]
    assert df["junction_aa"].to_list() == ["CARDRITIFGVVIIPDYFDYW"] * 2
    assert df["v_sequence_start"].to_list() == [1, 1] and df["productive"].to_list() == [True, True]


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


def test_long_sequences(tmp_path):
    seq = "ACGT" * 125
    df = abstar.run(seq)
    assert df.shape[0] == 1 and df["annotation_status"][0] == "unassigned"
    assert abstar.run([seq, seq]).shape[0] == 2
    s = abstar.run(seq, tmp_path / "out")
    assert [p.name for p in s.inputs[0].outputs.values()] == ["records.tsv"]


@pytest.mark.parametrize("item", ["sample1", "reads.txt", "ACGT-ACGT"])
def test_non_nucleotide_string_is_not_found(item, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    assert err_code(abstar.run, item) == "INPUT_NOT_FOUND"
    assert err_code(abstar.run, item, tmp_path / "out") == "INPUT_NOT_FOUND"
    assert not (tmp_path / "out").exists()


def test_iupac_and_lowercase_strings_are_sequences():
    assert abstar.run("acgtnRYSWKMBDHV").shape[0] == 1


def test_permission_error_is_unreadable(tmp_path):
    import os

    if os.geteuid() == 0:
        pytest.skip("permissions do not bind root")
    locked = tmp_path / "locked"
    locked.mkdir()
    locked.chmod(0)
    try:
        assert err_code(abstar.run, str(locked / "x.fasta")) == "INPUT_UNREADABLE"
        assert err_code(abstar.run, locked / "x.fasta") == "INPUT_UNREADABLE"
    finally:
        locked.chmod(0o755)


def test_directory_reports_unsupported_compression(indir):
    (indir / "z.fastq.xz").write_bytes(b"\xfd7zXZ\x00")
    with pytest.raises(AbstarError) as e:
        abstar.run(indir)
    assert e.value.code == "INPUT_UNREADABLE" and "use gzip" in e.value.message
    assert "z.fastq.xz" in e.value.message


def test_latin1_header_is_invalid_encoding(tmp_path):
    fa = tmp_path / "l1.fasta"
    fa.write_bytes(b">caf\xe9\nACGT\n>ok\nACGT\n")
    df = abstar.run(fa)
    assert df["sequence_id"].to_list() == ["caf�", "ok"]
    assert df["failure_reason"].to_list() == ["invalid_encoding", "no_v_hit"]
    assert df["annotation_status"].to_list() == ["invalid_input", "unassigned"]
