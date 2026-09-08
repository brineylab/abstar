"""Retained V evidence can end before an otherwise present FWR3 endpoint."""
import json
from pathlib import Path

import abutils
import pytest
from Bio.Seq import Seq

from abstar.annotation.antibody import Antibody
from abstar.annotation.annotator import annotate_single_sequence


CASES = json.loads(
    (Path(__file__).parents[1] / "test_data/missing_fwr3.json").read_text()
)["records"]
UNRECOVERABLE = json.loads(
    (Path(__file__).parents[1] / "test_data/missing_fwr3.json").read_text()
)["unrecoverable"]
REGIONS = ("fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4")


@pytest.mark.parametrize("removed,remaining", [(320, 94), (370, 44)])
def test_read_beginning_inside_fwr3_remains_a_valid_partial_region(removed, remaining):
    assignment = dict(CASES[0]["assignment"])
    assignment["sequence_input"] = assignment["sequence_input"][removed:]
    ab = annotate_single_sequence(Antibody(**assignment), "human")
    assert (ab.fwr3_start, ab.fwr3_end) == (0, remaining)
    assert ab.fwr3 == CASES[0]["expected"]["fwr3"][removed - 300:]
    assert ab.junction == CASES[0]["expected"]["junction"]
    assert ab.productive is True


@pytest.mark.parametrize("case", UNRECOVERABLE, ids=lambda c: c["assignment"]["sequence_id"])
def test_unsupported_fwr3_is_an_explicit_record_error(case):
    with pytest.raises(ValueError, match=case["expected_error"]):
        annotate_single_sequence(Antibody(**case["assignment"]), "human")


@pytest.mark.e2e
def test_unsupported_fwr3_diagnostics_preserve_other_records(tmp_path):
    import csv
    import abstar
    import polars as pl

    inputs = [case["assignment"] for case in [*UNRECOVERABLE, CASES[0]]]
    source = tmp_path / "fwr3.fasta"
    source.write_text("".join(f">{r['sequence_id']}\n{r['sequence_input']}\n" for r in inputs))
    project = tmp_path / "project"
    with pytest.warns(RuntimeWarning, match="1 sequence.*failed annotation"):
        abstar.run(str(source), project_path=str(project), output_format="parquet",
                   n_processes=1, mmseqs_threads=1)
    frame = pl.read_parquet(project / "parquet/fwr3.parquet")
    assert frame["sequence_id"].to_list() == [CASES[0]["expected"]["sequence_id"]]
    assert frame["fwr3"].to_list() == [CASES[0]["expected"]["fwr3"]]
    with (project / "logs/failures.tsv").open() as handle:
        failures = list(csv.DictReader(handle, delimiter="\t"))
    assert [r["sequence_id"] for r in failures] == [r["sequence_id"] for r in inputs[:len(UNRECOVERABLE)]]
    for failure, case in zip(failures, UNRECOVERABLE):
        diagnostic = (project / "logs" / failure["diagnostic_path"]).read_text()
        assert case["expected_error"] in diagnostic
        assert case["assignment"]["sequence_input"] in diagnostic


@pytest.mark.parametrize("deleted_start", [298, 299])
def test_recovered_fwr3_preserves_codon_deletion_boundary(deleted_start):
    assignment = dict(CASES[0]["assignment"])
    sequence = assignment["sequence_input"]
    assignment["sequence_input"] = sequence[:deleted_start] + sequence[deleted_start + 3:]
    ab = annotate_single_sequence(Antibody(**assignment), "human")
    assert (ab.cdr2_end, ab.fwr3_start, ab.fwr3_end) == (300, 300, 411)
    assert ab.fwr3 == (
        "TATGCACAGAAGTTCCAGGGCAGAGTCACCATGACCAGGAACACCTCCATAAGCACAGCCTAC"
        "ATGGAGCTGAGCAGCCTGAGATCTGAGGACACGGCCGTGTATTACTGC"
    )
    assert ab.fwr3_aa == "YAQKFQGRVTMTRNTSISTAYMELSSLRSEDTAVYYC"
    assert ab.junction == CASES[0]["expected"]["junction"]
    assert ab.productive is True
    vdj = ab.sequence_oriented[126:484]
    assert "".join(getattr(ab, r) for r in REGIONS) == vdj
    assert len(ab.cdr_mask) == 358
    assert "".join(getattr(ab, r + "_aa") for r in REGIONS) == ab.sequence_aa


@pytest.mark.parametrize("case", CASES, ids=lambda c: c["name"])
@pytest.mark.parametrize("reverse", [False, True])
def test_missing_fwr3_recovered_from_oriented_query(case, reverse):
    assignment = dict(case["assignment"])
    if reverse:
        assignment["sequence_input"] = str(Seq(assignment["sequence_input"]).reverse_complement())
        assignment["rev_comp"] = True
    ab = annotate_single_sequence(Antibody(**assignment), "human")
    for field, expected in case["expected"].items():
        assert getattr(ab, field) == expected, field
    assert ab.rev_comp is reverse
    assert "FWR3 REGION RECOVERY" in ab.format_log()
    # This reference is independent of both gene and FWR/CDR assembly.
    vdj = ab.sequence_oriented[ab.v_sequence_start:ab.j_sequence_end]
    assert ab.sequence == ab.sequence_alignment.replace("-", "") == vdj
    assert "".join(getattr(ab, region) for region in REGIONS) == vdj
    assert len(ab.cdr_mask) == len(vdj)
    if case["expected"]["productive"]:
        coding = vdj[ab.frame - 1:]
        protein = str(Seq(coding[:len(coding) // 3 * 3]).translate())
        assert "".join(getattr(ab, region + "_aa") for region in REGIONS) == protein
        assert len(ab.cdr_mask_aa) == len(protein)


@pytest.mark.e2e
@pytest.mark.parametrize("entrypoint", ["api", "cli"])
def test_missing_fwr3_public_entrypoints(tmp_path, entrypoint):
    case = CASES[0]
    assignment = case["assignment"]
    if entrypoint == "api":
        import abstar

        result = abstar.run(
            abutils.Sequence(assignment["sequence_input"], id=assignment["sequence_id"]),
            n_processes=1, mmseqs_threads=1, strict=True,
        )
        row = result
    else:
        import polars as pl
        from click.testing import CliRunner
        from abstar.scripts.abstar import cli

        source = tmp_path / "missing.fasta"
        source.write_text(f">{assignment['sequence_id']}\n{assignment['sequence_input']}\n")
        project = tmp_path / "project"
        result = CliRunner().invoke(cli, [
            "run", str(source), str(project), "-o", "parquet", "--strict",
            "--n_processes", "2", "--mmseqs_threads", "1", "--quiet",
        ])
        assert result.exit_code == 0, (result.output, result.exception)
        frame = pl.read_parquet(project / "parquet/missing.parquet")
        assert frame.height == 1
        row = frame.row(0, named=True)
        assert not list((project / "tmp").glob("*.parquet"))
    for field, expected in case["expected"].items():
        assert row[field] == expected, field
    assert len(row["cdr_mask"]) == 361
    assert len(row["cdr_mask_aa"]) == 120
