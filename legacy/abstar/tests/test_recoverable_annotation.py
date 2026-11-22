"""Record failures must survive without discarding other records or samples."""
import csv
import json
import subprocess
import sys
from pathlib import Path

import abstar
from abutils import Sequence
import polars as pl
import pytest

from abstar.tests.recoverable_annotation_helpers import (
    FAILURE_MESSAGE,
    INJECTED_FAILURE_SEQUENCE,
    annotate_with_record_failure,
)


@pytest.fixture
def inject_record_failure(monkeypatch):
    monkeypatch.setattr("abstar.core.abstar.annotate", annotate_with_record_failure)


@pytest.mark.e2e
@pytest.mark.parametrize("entrypoint", ["api", "cli"])
def test_failed_records_do_not_discard_chunk_or_later_samples(tmp_path, public_bcr_cases, entrypoint, inject_record_failure):
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    good = public_bcr_cases[0]
    (inputs / "first.fasta").write_text(
        f">duplicate\n{INJECTED_FAILURE_SEQUENCE}\n>{good.sequence_id}\n{good.sequence}\n"
        f">duplicate\n{INJECTED_FAILURE_SEQUENCE}\n>unassigned\nN\n"
    )
    (inputs / "second.fasta").write_text(f">duplicate\n{INJECTED_FAILURE_SEQUENCE}\n")
    (inputs / "third.fasta").write_text(f">{good.sequence_id}\n{good.sequence}\n")
    project = tmp_path / "project"
    if entrypoint == "api":
        with pytest.warns(RuntimeWarning, match="3 sequence.*failed annotation"):
            abstar.run(str(inputs), project_path=str(project), output_format=["airr", "parquet"],
                       chunksize=2, n_processes=1, mmseqs_threads=1)
    else:
        result = subprocess.run(
            [sys.executable, "-c", "import importlib; "
             "from abstar.tests.recoverable_annotation_helpers import annotate_with_record_failure; "
             "importlib.import_module('abstar.core.abstar').annotate = annotate_with_record_failure; "
             "from abstar.scripts.abstar import cli; cli()",
             "run", str(inputs), str(project), "-o", "parquet",
             "--chunksize", "2", "--n_processes", "2", "--mmseqs_threads", "1"],
            capture_output=True, text=True, check=False,
        )
        assert result.returncode == 0, result.stdout + result.stderr
        assert "3 sequences failed annotation" in result.stdout + result.stderr
    for sample in ("first", "third"):
        frame = pl.read_parquet(project / "parquet" / f"{sample}.parquet")
        assert frame["sequence_id"].to_list() == ([good.sequence_id, "unassigned"] if sample == "first" else [good.sequence_id])
        assert frame.filter(pl.col("annotation_status") == "annotated")["junction"].to_list() == [good.expected["junction"]]
        if sample == "first":
            unassigned = frame.row(1, named=True)
            assert unassigned["annotation_status"] == "unassigned"
            assert unassigned["productive"] is None
            assert unassigned["failure_reason"] == "no compatible V gene assignment"
    assert pl.read_parquet(project / "parquet" / "second.parquet").is_empty()
    with (project / "logs" / "failures.tsv").open() as handle:
        failures = list(csv.DictReader(handle, delimiter="\t"))
    assert len(failures) == 3
    assert len({f["row_id"] for f in failures}) == 3
    assert [f["sample"] for f in failures] == ["first", "first", "second"]
    assert all(f["sequence_id"] == "duplicate" for f in failures)
    assert all(f["exception_type"] == "RuntimeError" for f in failures)
    paths = [project / "logs" / f["diagnostic_path"] for f in failures]
    assert len(set(paths)) == 3
    assert [p.parent.name for p in paths] == ["first", "first", "second"]
    for path in paths:
        diagnostic = path.read_text()
        assert INJECTED_FAILURE_SEQUENCE in diagnostic
        assert "Traceback" in diagnostic and FAILURE_MESSAGE in diagnostic
        assert "IGLV2-14" in diagnostic
    assert not list((project / "tmp").glob("*.parquet"))
    metadata = json.loads((project / "logs" / "run.json").read_text())
    assert metadata["parameters"]["strict"] is False
    assert metadata["versions"]["abutils"]


@pytest.mark.e2e
@pytest.mark.parametrize("as_dataframe", [False, True])
def test_no_project_all_failed_retains_diagnostics(tmp_path, monkeypatch, as_dataframe, inject_record_failure):
    import tempfile
    monkeypatch.setattr(tempfile, "tempdir", str(tmp_path))
    with pytest.warns(RuntimeWarning, match="1 sequence.*failed annotation") as warnings:
        result = abstar.run(INJECTED_FAILURE_SEQUENCE, n_processes=1, mmseqs_threads=1, as_dataframe=as_dataframe)
    assert len(result) == 0
    indexes = list(tmp_path.glob("abstar-*/logs/failures.tsv"))
    assert len(indexes) == 1
    assert str(indexes[0]) in str(warnings[0].message)
    assert list(indexes[0].parent.glob("sequences/*.failed"))
    assert not list(indexes[0].parent.parent.glob("tmp/*.parquet"))


@pytest.mark.e2e
def test_strict_mode_still_raises_for_record_error(tmp_path, inject_record_failure):
    with pytest.raises(abstar.AnnotationRunError) as caught:
        abstar.run(INJECTED_FAILURE_SEQUENCE, project_path=str(tmp_path / "project"), strict=True,
                   n_processes=1, mmseqs_threads=1)
    assert caught.value.failures[0].sequence_id
    assert caught.value.failures[0].message == FAILURE_MESSAGE
    assert any(Path(p).suffix == ".failed" for p in caught.value.partial_output_paths)


@pytest.mark.e2e
def test_cli_strict_mode_keeps_nonzero_exit_and_diagnostics(tmp_path, inject_record_failure):
    from click.testing import CliRunner
    from abstar.scripts.abstar import cli
    source = tmp_path / "sample.fasta"
    source.write_text(f">failed\n{INJECTED_FAILURE_SEQUENCE}\n")
    project = tmp_path / "project"
    result = CliRunner().invoke(cli, ["run", str(source), str(project), "--strict",
                                     "--n_processes", "1", "--mmseqs_threads", "1", "--quiet"])
    assert result.exit_code != 0
    assert "annotation/internal_error=1" in result.output
    paths = list((project / "logs" / "sample").glob("*.failed"))
    assert len(paths) == 1 and str(paths[0]) in result.output
    assert not list((project / "airr").glob("*.tsv"))


def test_record_diagnostics_encode_ids_and_preserve_duplicate_records(tmp_path):
    directory = tmp_path / "sample"
    directory.mkdir()
    identifiers = ["../../escape", "a/b", "a%2Fb", "x" * 500, "duplicate", "duplicate"]
    records = [dict(row_id=f"abstar_0_{i}", sequence_id=identifier,
                    sequence_input=INJECTED_FAILURE_SEQUENCE, v_call="IGLV2-14*01",
                    j_call="IGLJ2*01", rev_comp=False)
               for i, identifier in enumerate(identifiers)]
    source = tmp_path / "assigned.parquet"
    pl.DataFrame(records).write_parquet(source)
    result = annotate_with_record_failure(str(source), output_directory=str(tmp_path),
                                          germline_database="human", failure_directory=str(directory))
    assert len(result.failures) == len(identifiers)
    assert all(failure.message == FAILURE_MESSAGE for failure in result.failures)
    files = list(directory.glob("*.failed"))
    assert len(files) == len(identifiers)
    saved = [json.loads(path.read_text().splitlines()[1]) for path in files]
    assert sorted(r["sequence_id"] for r in saved) == sorted(identifiers)
    assert all(r["sequence_input"] == INJECTED_FAILURE_SEQUENCE for r in saved)
    assert all(len(path.name.encode()) < 255 for path in files)
    assert not (tmp_path.parent / "escape.failed").exists()


def test_failure_log_write_error_escapes_record_recovery(tmp_path, monkeypatch):
    source = tmp_path / "assigned.parquet"
    pl.DataFrame([dict(row_id="abstar_0_0", sequence_id="failed",
                       sequence_input=INJECTED_FAILURE_SEQUENCE, v_call="IGLV2-14*01",
                       j_call="IGLJ2*01", rev_comp=False)]).write_parquet(source)
    original = Path.open

    def fail_write(path, *args, **kwargs):
        if path.suffix == ".failed":
            raise PermissionError("diagnostics unavailable")
        return original(path, *args, **kwargs)

    monkeypatch.setattr(Path, "open", fail_write)
    with pytest.raises(PermissionError, match="diagnostics unavailable"):
        annotate_with_record_failure(str(source), output_directory=str(tmp_path),
                                     germline_database="human", failure_directory=str(tmp_path))
    assert not (tmp_path / "assigned.parquet_annotated.parquet").exists()


def test_diagnostic_reruns_preserve_prior_index_and_record(tmp_path):
    from abstar.core.diagnostics import initialize_diagnostics, sample_directory
    initialize_diagnostics(tmp_path, {"strict": False}, ["source/sample.fasta"])
    first = Path(sample_directory(tmp_path, "sample"))
    (first / "record.failed").write_text("original failure")
    (tmp_path / "failures.tsv").write_text("original index")
    initialize_diagnostics(tmp_path, {"strict": True}, ["source/sample.fasta"])
    second = Path(sample_directory(tmp_path, "sample"))
    assert first.name == "sample" and second != first
    assert (first / "record.failed").read_text() == "original failure"
    assert [p.read_text() for p in tmp_path.glob("failures.*.tsv")] == ["original index"]
    assert json.loads((tmp_path / "run.json").read_text())["parameters"]["strict"] is True


@pytest.mark.parametrize("error", [OSError("germline read failed"), MemoryError("allocation failed")])
def test_infrastructure_errors_are_not_record_failures(tmp_path, monkeypatch, error):
    from abstar.annotation.annotator import annotate
    source = tmp_path / "assigned.parquet"
    pl.DataFrame([dict(row_id="abstar_0_0", sequence_id="failed", sequence_input="ACGT",
                       v_call="IGLV2-14*01", j_call="IGLJ2*01")]).write_parquet(source)

    def fail(**kwargs):
        raise error

    monkeypatch.setattr("abstar.annotation.annotator.annotate_single_sequence", fail)
    with pytest.raises(type(error), match=str(error)):
        annotate(str(source), str(tmp_path), "human", failure_directory=str(tmp_path))


def test_promoted_warning_retains_no_project_diagnostics(tmp_path, monkeypatch):
    import tempfile
    import warnings
    from abstar.core.abstar import _project_workspace
    monkeypatch.setattr(tempfile, "tempdir", str(tmp_path))
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        with pytest.raises(RuntimeWarning, match="record failed"):
            with _project_workspace(None, False, []) as workspace:
                logs = Path(workspace) / "logs"
                logs.mkdir()
                (logs / "failures.tsv").write_text("header\nfailed record\n")
                (logs / "record.failed").write_text("diagnostic")
                warnings.warn("record failed", RuntimeWarning)
    assert (logs / "record.failed").read_text() == "diagnostic"


@pytest.mark.e2e
def test_index_storage_failure_is_fatal(tmp_path, monkeypatch, public_bcr_cases, inject_record_failure):
    import csv
    original = csv.DictWriter.writerow

    def fail_index(writer, row):
        if row.get("category") == "internal_error":
            raise OSError("failure index unavailable")
        return original(writer, row)

    monkeypatch.setattr(csv.DictWriter, "writerow", fail_index)
    project = tmp_path / "project"
    with pytest.raises(abstar.AnnotationRunError) as caught:
        abstar.run([public_bcr_cases[0].as_sequence(), Sequence(INJECTED_FAILURE_SEQUENCE)],
                   project_path=str(project), n_processes=1, mmseqs_threads=1)
    assert caught.value.failures[0].stage == "output"
    assert "failure index unavailable" in caught.value.failures[0].message
    assert not list((project / "airr").glob("*.tsv"))
    assert list((project / "logs" / "sequences").glob("*.failed"))


@pytest.mark.parametrize("boundary", ["initialize_diagnostics", "sample_directory"])
def test_diagnostic_setup_failure_aborts_before_output(tmp_path, monkeypatch, public_bcr_cases, boundary):
    def fail(*args, **kwargs):
        raise PermissionError("diagnostic storage unavailable")

    monkeypatch.setattr(f"abstar.core.abstar.{boundary}", fail)
    project = tmp_path / "project"
    with pytest.raises(abstar.AnnotationRunError) as caught:
        abstar.run(public_bcr_cases[0].as_sequence(), project_path=str(project),
                   n_processes=1, mmseqs_threads=1)
    assert caught.value.failures[0].stage == "output"
    assert caught.value.failures[0].message == "diagnostic storage unavailable"
    assert any("diagnostic storage unavailable" in Path(path).read_text()
               for path in caught.value.partial_output_paths if path.endswith(".failed"))
    assert not list((project / "airr").glob("*.tsv"))


@pytest.mark.parametrize("strict", ["false", 1, None])
def test_invalid_strict_option_is_rejected_before_project_creation(tmp_path, strict):
    project = tmp_path / "project"
    with pytest.raises(ValueError, match="strict must be a boolean"):
        abstar.run("ACGT", project_path=str(project), strict=strict)
    assert not project.exists()
