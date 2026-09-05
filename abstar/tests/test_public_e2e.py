"""Public input, ordering, return-shape and CLI contracts on accepted BCRs."""

from pathlib import Path

import abstar
from abutils import Sequence
from click.testing import CliRunner
import polars as pl
import pytest

from abstar.annotation.schema import OUTPUT_SCHEMA
from abstar.scripts.abstar import cli
from abstar.tests import helpers


ANNOTATION_FIELDS = (
    "sequence_id", "sequence_input", "sequence_oriented", "sequence",
    "annotation_status", "failure_reason", "germline_database", "species",
    "locus", "rev_comp", "v_call", "d_call", "j_call", "c_call",
    "v_score", "d_score", "j_score", "c_score",
    "v_sequence_start", "v_sequence_end", "v_germline_start", "v_germline_end",
    "j_sequence_start", "j_sequence_end", "j_germline_start", "j_germline_end",
    "c_sequence_start", "c_sequence_end", "c_germline_start", "c_germline_end",
    "junction", "junction_aa", "cdr3", "cdr3_aa", "cdr3_length",
    "productive", "productivity_issues", "vj_in_frame", "stop_codon",
    "v_insertions", "v_deletions",
)
INPUT_FORMS = (
    "list", "iterator", "generator", "fasta", "fastq",
    "flat_directory", "nested_directory",
)
WORKER_CHUNKS = ((1, 1), (1, 2), (2, 1), (2, 3))


@pytest.fixture(scope="module")
def public_bcr_baseline(public_bcr_cases):
    """Serial list baseline, independently checked against accepted evidence."""
    result = abstar.run([case.as_sequence() for case in public_bcr_cases],
                        n_processes=1, chunksize=1, mmseqs_threads=1)
    assert isinstance(result, list) and len(result) == 3
    assert all(isinstance(record, Sequence) for record in result)
    assert [record.id for record in result] == [case.sequence_id for case in public_bcr_cases]
    rows = [dict(record.annotations) for record in result]
    assert [row["locus"] for row in rows] == ["IGH", "IGK", "IGL"]
    for row, case in zip(rows, public_bcr_cases):
        assert row["sequence_input"] == row["sequence_oriented"] == case.sequence
        assert row["annotation_status"] == case.expected["status"] == "annotated"
        assert row["failure_reason"] is None
        assert row["germline_database"] == "human"
        assert "row_id" not in row
        for field, expected in case.expected.items():
            if field in ("status", "junction_start", "junction_end"):
                continue  # Junction boundaries are internal, not public fields.
            if field.endswith("_call") and isinstance(expected, tuple):
                calls = {call.split("*")[0] for call in row[field].split(",")}
                assert calls and calls <= set(expected), (case.sequence_id, field)
            elif field == "productivity_issues":
                assert (tuple(row[field].split("|")) if row[field] else ()) == expected
            else:
                assert row[field] == expected, (case.sequence_id, field)
        assert row["junction"] == case.sequence[
            case.expected["junction_start"]:case.expected["junction_end"]]
    return rows


def test_shared_comparator_checks_count_order_fields_and_value_types():
    compare = getattr(helpers, "assert_same_annotations", None)
    assert callable(compare), "the shared ordered annotation comparator is missing"
    baseline = [{"sequence_id": "00123", "productive": False},
                {"sequence_id": "10E8", "productive": True}]
    fields = ("sequence_id", "productive")
    compare(baseline, [dict(row) for row in baseline], fields)
    for changed in (baseline[:1], baseline[::-1],
                    [{"sequence_id": 123, "productive": False}, baseline[1]],
                    [{"sequence_id": "00123", "productive": 0}, baseline[1]],
                    [{"sequence_id": "00123"}, baseline[1]]):
        with pytest.raises(AssertionError):
            compare(baseline, changed, fields)


@pytest.mark.e2e
@pytest.mark.parametrize("input_form", INPUT_FORMS)
@pytest.mark.parametrize("n_processes,chunksize", WORKER_CHUNKS)
def test_api_input_and_worker_matrix_matches_authenticated_serial_baseline(
    public_bcr_inputs, public_bcr_baseline, input_form, n_processes, chunksize,
):
    result = abstar.run(public_bcr_inputs[input_form](), n_processes=n_processes,
                        chunksize=chunksize, mmseqs_threads=1)
    assert isinstance(result, list)
    assert all(isinstance(record, Sequence) for record in result)
    assert [record.id for record in result] == [row["sequence_id"] for row in public_bcr_baseline]
    assert all("row_id" not in record.annotations for record in result)
    helpers.assert_same_annotations(public_bcr_baseline,
                                    [record.annotations for record in result], ANNOTATION_FIELDS)


@pytest.mark.e2e
def test_api_arbitrary_iterable_is_consumed_once(public_bcr_cases, public_bcr_baseline):
    class OnceOnly:
        def __init__(self):
            self.consumed = False

        def __iter__(self):
            assert not self.consumed, "input iterable was consumed twice"
            self.consumed = True
            return (case.as_sequence() for case in public_bcr_cases)

    result = abstar.run(OnceOnly(), n_processes=1, mmseqs_threads=1)
    helpers.assert_same_annotations(public_bcr_baseline,
                                    [record.annotations for record in result], ANNOTATION_FIELDS)


@pytest.mark.e2e
@pytest.mark.parametrize("as_dataframe", (False, True))
def test_api_one_record_return_shape(public_bcr_cases, public_bcr_baseline, as_dataframe):
    result = abstar.run(public_bcr_cases[0].as_sequence(), as_dataframe=as_dataframe,
                        n_processes=1, mmseqs_threads=1)
    if as_dataframe:
        assert isinstance(result, pl.DataFrame)
        rows = result.to_dicts()
    else:
        assert isinstance(result, Sequence)
        assert result.id == public_bcr_cases[0].sequence_id
        rows = [result.annotations]
    helpers.assert_same_annotations(public_bcr_baseline[:1], rows, ANNOTATION_FIELDS)


@pytest.mark.e2e
@pytest.mark.parametrize("n_processes,chunksize", WORKER_CHUNKS)
def test_api_mixed_outcomes_keep_multiple_input_return_shape(
    public_bcr_cases, public_bcr_baseline, n_processes, chunksize,
):
    result = abstar.run([public_bcr_cases[0].as_sequence(), Sequence("N", id="00123")],
                        n_processes=n_processes, chunksize=chunksize, mmseqs_threads=1)
    assert isinstance(result, list) and all(isinstance(record, Sequence) for record in result)
    assert [record.id for record in result] == [public_bcr_cases[0].sequence_id, "00123"]
    helpers.assert_same_annotations(public_bcr_baseline[:1],
                                    [result[0].annotations], ANNOTATION_FIELDS)
    assert result[1]["annotation_status"] == "unassigned"
    assert result[1]["failure_reason"] == "no compatible V gene assignment"
    assert result[1]["sequence_input"] == "N"
    assert result[1]["v_call"] is result[1]["j_call"] is None
    assert result[1]["productive"] is result[1]["productivity_issues"] is None
    assert "row_id" not in result[1].annotations


@pytest.mark.e2e
@pytest.mark.parametrize("input_form", ("list", "iterator", "generator", "directory"))
def test_api_empty_inputs_raise_before_creating_project(tmp_path, input_form):
    empty = tmp_path / "empty"
    empty.mkdir()
    inputs = {"list": [], "iterator": iter(()), "generator": (item for item in ()),
              "directory": str(empty)}
    project = tmp_path / "must-not-exist"
    with pytest.raises(ValueError, match="empty|No supported FASTA or FASTQ"):
        abstar.run(inputs[input_form], project_path=str(project), n_processes=1)
    assert not project.exists()


@pytest.mark.e2e
@pytest.mark.parametrize("parameter,value", [
    (parameter, value) for parameter in ("n_processes", "chunksize", "mmseqs_chunksize")
    for value in (0, -1, 1.5, "2", True)
] + [("output_format", value) for value in ("csv", [], ["airr", "csv"], [None], None)])
def test_api_invalid_options_raise_before_creating_project(
    public_bcr_cases, tmp_path, parameter, value,
):
    project = tmp_path / "must-not-exist"
    with pytest.raises(ValueError, match="positive integer|output.format|Unsupported output format"):
        abstar.run(public_bcr_cases[0].as_sequence(), project_path=str(project),
                   **{parameter: value})
    assert not project.exists()


def test_documented_namespaces_resolve_to_active_implementations():
    from abstar import gl, pp, tl
    from abstar.annotation.germline import get_germline, get_germline_database_path
    from abstar.annotation.umi import parse_umis
    from abstar.core.germline import build_germline_database
    from abstar.preprocess.merging import merge_fastqs

    assert gl.get_germline is get_germline
    assert gl.get_germline_database_path is get_germline_database_path
    assert pp.merge_fastqs is merge_fastqs
    assert tl.parse_umis is parse_umis
    assert tl.build_germline_database is build_germline_database


@pytest.mark.e2e
def test_api_project_returns_none_and_preserves_nested_input_paths(
    public_bcr_inputs, public_bcr_baseline, tmp_path,
):
    source = Path(public_bcr_inputs["nested_directory"]())
    project = tmp_path / "project"
    result = abstar.run(str(source), project_path=str(project), output_format=["airr", "parquet"],
                        n_processes=1, mmseqs_threads=1, copy_inputs_to_project=True)
    assert result is None
    rows = []
    for ordinal in (1, 2, 10):
        relative = Path("batch") / f"sample{ordinal}" / "reads.fasta"
        assert (project / "input" / relative).read_bytes() == (source / relative).read_bytes()
        name = f"sample{ordinal}__reads"
        parquet = pl.read_parquet(project / "parquet" / f"{name}.parquet")
        airr = pl.read_csv(project / "airr" / f"{name}.tsv", separator="\t",
                           schema_overrides=OUTPUT_SCHEMA)
        helpers.assert_same_annotations(parquet.to_dicts(), airr.to_dicts(), ANNOTATION_FIELDS)
        assert "row_id" not in parquet.columns and "row_id" not in airr.columns
        rows.extend(parquet.to_dicts())
    assert len(list((project / "input").rglob("*.fasta"))) == 3
    assert len(list((project / "parquet").glob("*.parquet"))) == 3
    assert len(list((project / "airr").glob("*.tsv"))) == 3
    helpers.assert_same_annotations(public_bcr_baseline, rows, ANNOTATION_FIELDS)


@pytest.mark.e2e
@pytest.mark.parametrize("arguments,expected", [
    (["--help"], ("run", "build_germline_database")),
    (["run", "--help"], ("INPUT_PATH PROJECT_PATH", "--receptor", "--output_format",
                             "--n_processes", "--chunksize", "--copy-inputs")),
])
def test_cli_public_help(arguments, expected):
    result = CliRunner().invoke(cli, arguments)
    assert result.exit_code == 0, result.output
    assert all(text in result.output for text in expected)


@pytest.mark.e2e
def test_cli_real_airr_run_matches_api(public_bcr_inputs, public_bcr_baseline, tmp_path):
    project = tmp_path / "cli-project"
    result = CliRunner().invoke(cli, ["run", public_bcr_inputs["fastq"](), str(project),
                                    "--output_format", "airr", "--n_processes", "2",
                                    "--chunksize", "2", "--mmseqs_threads", "1", "--quiet"])
    assert result.exit_code == 0, (result.output, result.exception)
    rows = pl.read_csv(project / "airr" / "controls.tsv", separator="\t",
                       schema_overrides=OUTPUT_SCHEMA).to_dicts()
    helpers.assert_same_annotations(public_bcr_baseline, rows, ANNOTATION_FIELDS)
    assert all("row_id" not in row for row in rows)
    assert (project / "input" / "controls.fastq").read_bytes() == Path(
        public_bcr_inputs["fastq"]()).read_bytes()


def _crash_annotation_worker(*args, **kwargs):
    raise RuntimeError("forced public worker failure")


@pytest.mark.e2e
@pytest.mark.parametrize("boundary", ("assignment", "annotation"))
def test_cli_structured_failure_reports_surviving_artifact(
    public_bcr_inputs, public_bcr_cases, monkeypatch, tmp_path, boundary,
):
    if boundary == "assignment":
        from abstar.assigners.mmseqs import AssignmentExternalToolError

        def fail_assignment(*args, **kwargs):
            raise AssignmentExternalToolError("forced public tool failure")

        monkeypatch.setattr("abstar.core.abstar.MMseqs.__call__", fail_assignment)
        summary = "assignment/external_tool=1"
        diagnostic = "forced public tool failure"
    else:
        monkeypatch.setattr("abstar.core.abstar.annotate", _crash_annotation_worker)
        summary = "annotation/internal_error=3"
        diagnostic = "forced public worker failure"
    project = tmp_path / "failed-cli-project"
    result = CliRunner().invoke(cli, ["run", public_bcr_inputs["fasta"](), str(project),
                                    "--n_processes", "1", "--chunksize", "2",
                                    "--mmseqs_threads", "1", "--quiet"])
    assert result.exit_code != 0
    artifact = project / "logs" / "controls.failed"
    assert artifact.is_file() and diagnostic in artifact.read_text()
    assert summary in result.output
    assert str(artifact) in result.output
    error = result.exception
    while error is not None and not isinstance(error, abstar.AnnotationRunError):
        error = error.__cause__ or error.__context__
    assert isinstance(error, abstar.AnnotationRunError)
    assert all(failure.stage == boundary for failure in error.failures)
    if boundary == "annotation":
        assert len(error.failures) == 3
        assert [failure.sequence_id for failure in error.failures] == [
            case.sequence_id for case in public_bcr_cases]
        assert [failure.row_id for failure in error.failures] == [
            "abstar_0_0", "abstar_0_1", "abstar_0_2"]
        assert error.partial_output_paths
        assert all(Path(path).is_file() for path in error.partial_output_paths)
