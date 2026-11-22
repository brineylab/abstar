import csv
import json

import pytest
from typer.testing import CliRunner

from abstar import germline
from abstar.cli import app
from abstar.errors import AbstarError


def v_alleles(db, locus):
    with open(germline.info(db).path / "alleles.tsv", newline="") as f:
        return [r for r in csv.DictReader(f, delimiter="\t") if r["segment"] == "V" and r["locus"] == locus]


def ungapped(r):
    return r["sequence_gapped"].replace(".", "")


def j_record():
    j = next(r for r in csv.DictReader(open(germline.info("human").path / "alleles.tsv"), delimiter="\t") if r["locus"] == "IGH" and r["segment"] == "J")
    return ("IGHJ-test", ungapped(j))


def test_gap_regaps_human_ighv():
    rows = v_alleles("human", "IGH")
    out, report = germline.gap([(r["name"], ungapped(r)) for r in rows], templates=["human"])
    assert [s for _, s in out] == [r["sequence_gapped"] for r in rows]
    assert report.filter(report["status"] == "failed").height == 0


def test_gap_report_columns_and_order():
    r = v_alleles("human", "IGH")[0]
    recs = [(r["name"] + " dup", ungapped(r)), (r["name"] + " dup", ungapped(r)), j_record()]
    out, report = germline.gap(recs, templates=["human"])
    assert report.columns == [
        "name", "status", "template", "template_species", "identity",
        "cdr1_length", "cdr2_length", "warnings", "failure_reason",
    ]
    assert report.height == 3
    assert report["status"].to_list() == ["gapped", "gapped", "failed"]
    assert report["failure_reason"][2] == "no_alignment"
    assert [h for h, _ in out] == [r["name"] + " dup"] * 2


def test_gap_none_gapped_raises():
    with pytest.raises(AbstarError) as e:
        germline.gap([j_record()], templates=["human"])
    assert e.value.code == "NO_ALLELES_GAPPED"


def write_fasta(tmp_path):
    rows = v_alleles("human", "IGH")[:2]
    p = tmp_path / "in.fasta"
    p.write_text("".join(f">{r['name']} some description\n{ungapped(r)}\n" for r in rows))
    return p, rows


def test_cli_gap_stdout_and_report(tmp_path):
    p, rows = write_fasta(tmp_path)
    res = CliRunner().invoke(app, ["germline", "gap", str(p), "--templates", "human", "--report", str(tmp_path / "r.tsv")])
    assert res.exit_code == 0, res.output
    lines = res.stdout.splitlines()
    assert [l for l in lines if l.startswith(">")] == [f">{r['name']} some description" for r in rows]
    assert lines[1] == rows[0]["sequence_gapped"]
    assert len((tmp_path / "r.tsv").read_text().splitlines()) == 3


def test_cli_gap_json_summary(tmp_path):
    p, _ = write_fasta(tmp_path)
    res = CliRunner().invoke(app, ["germline", "gap", str(p), "--templates", "human", "--output", str(tmp_path / "o.fasta"), "--json"])
    assert res.exit_code == 0, res.output
    s = json.loads(res.stdout)
    assert s["n_input"] == 2 and s["n_failed"] == 0 and s["templates"] == ["human"]


def test_gap_templates_replace_default():
    r = v_alleles("human", "IGK")[0]
    _, report = germline.gap([(r["name"], ungapped(r))], templates=["c57bl6"])
    assert report["template_species"][0] == "mouse"
