"""FWR4 must terminate at the final retained J boundary, including repeats."""
import csv
import json
from pathlib import Path

from Bio.Seq import Seq
import pytest

from abstar.annotation.antibody import Antibody
from abstar.annotation.annotator import annotate_single_sequence


CASES = json.loads(
    (Path(__file__).parents[1] / "test_data/fwr4_endpoints.json").read_text()
)["records"]
REGIONS = ("fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4")


def assert_expected(row, expected):
    for field, value in expected.items():
        if field in ("cdr_mask_length", "cdr_mask_aa_length"):
            assert len(row[field.removesuffix("_length")]) == value, field
        else:
            assert row[field] == value, field


@pytest.mark.parametrize("case", CASES, ids=lambda c: c["name"])
@pytest.mark.parametrize("reverse", [False, True])
def test_fwr4_uses_final_j_endpoint_without_changing_gene_evidence(case, reverse):
    assignment = dict(case["assignment"])
    if reverse:
        assignment["sequence_input"] = str(Seq(assignment["sequence_input"]).reverse_complement())
        assignment["rev_comp"] = True
    ab = annotate_single_sequence(Antibody(**assignment), "human")
    assert_expected(ab.to_dict(), case["expected"])
    assert ab.rev_comp is reverse
    assert ab.fwr4_start == ab.junction_end - 3 == ab.cdr3_end
    assert ab.fwr4_end == ab.j_sequence_end
    vdj = ab.sequence_oriented[ab.v_sequence_start:ab.j_sequence_end]
    assert "".join(getattr(ab, r) for r in REGIONS) == vdj == ab.sequence
    assert ab.sequence_alignment.replace("-", "") == vdj
    if ab.productive:
        coding = vdj[ab.frame - 1:]
        protein = str(Seq(coding[:len(coding) // 3 * 3]).translate())
        assert "".join(getattr(ab, r + "_aa") for r in REGIONS) == protein
        assert len(ab.cdr_mask_aa) == len(protein)


@pytest.mark.parametrize("name,j_start,j_end,nt,aa", [
    ("heavy_productive", 462, 485, "TGGATGGGCCAGGGAACCCTGGTCACCG", "WMGQGTLVT"),
    ("kappa", 393, 407, "TTCGGCCAA", "FGQ"),
])
def test_truncated_j_preserves_retained_endpoint(name, j_start, j_end, nt, aa):
    case = next(c for c in CASES if c["name"] == name)
    assignment = dict(case["assignment"])
    assignment["sequence_input"] = assignment["sequence_input"][:case["expected"]["j_sequence_end"] - 6]
    assignment["c_call"] = None
    ab = annotate_single_sequence(Antibody(**assignment), "human")
    assert (ab.j_sequence_start, ab.j_sequence_end) == (j_start, j_end)
    assert ab.fwr4_start == case["expected"]["fwr4_start"]
    assert ab.fwr4_end == j_end
    assert (ab.fwr4, ab.fwr4_aa) == (nt, aa)
    assert ab.junction == case["expected"]["junction"]
    assert ab.productive is True


@pytest.mark.e2e
@pytest.mark.parametrize("entrypoint", ["api", "cli"])
def test_final_fwr4_endpoint_survives_public_serialization(tmp_path, entrypoint):
    import abstar
    import polars as pl
    from click.testing import CliRunner
    from abstar.scripts.abstar import cli

    source = tmp_path / "endpoints.fasta"
    source.write_text("".join(
        f">{c['assignment']['sequence_id']}\n{c['assignment']['sequence_input']}\n" for c in CASES
    ))
    project = tmp_path / "project"
    if entrypoint == "api":
        abstar.run(str(source), project_path=str(project), output_format=["airr", "parquet"],
                   n_processes=1, chunksize=1, mmseqs_threads=1, strict=True)
    else:
        result = CliRunner().invoke(cli, [
            "run", str(source), str(project), "-o", "airr", "-o", "parquet",
            "--n_processes", "2", "--chunksize", "3", "--mmseqs_threads", "1",
            "--strict", "--quiet",
        ])
        assert result.exit_code == 0, (result.output, result.exception)
    rows = pl.read_parquet(project / "parquet/endpoints.parquet").to_dicts()
    with (project / "airr/endpoints.tsv").open() as handle:
        airr = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == len(airr) == len(CASES)
    for row, text, case in zip(rows, airr, CASES):
        assert_expected(row, case["expected"])
        assert row["sequence"] == case["assignment"]["sequence_input"]
        assert int(text["fwr4_start"]) == row["fwr4_start"] + 1
        assert int(text["fwr4_end"]) == int(text["j_sequence_end"]) == row["j_sequence_end"]
        assert text["fwr4"] == row["fwr4"]
        assert text["fwr4_aa"] == row["fwr4_aa"]
        assert text["productive"] == ("T" if row["productive"] else "F")
    with (project / "logs/failures.tsv").open() as handle:
        assert list(csv.DictReader(handle, delimiter="\t")) == []
    assert not list((project / "tmp").glob("*.parquet"))
