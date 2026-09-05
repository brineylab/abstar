# Copyright (c) 2026 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

import pickle
import shlex
import subprocess
import tempfile
from concurrent.futures import Future
from dataclasses import FrozenInstanceError
from pathlib import Path

import polars as pl
import pytest
import abstar
import abutils
from click.testing import CliRunner

from ..annotation.annotator import annotate, annotate_single_sequence
from ..core.results import AnnotationChunkResult, AnnotationRunError, RecordFailure
from ..core.abstar import _assert_record_conservation
from ..scripts.abstar import cli


def test_record_failure_is_immutable_and_pickleable():
    failure = RecordFailure(
        row_id="abstar_0_1",
        sequence_id="sample-001",
        stage="annotation",
        category="internal_error",
        message="annotation failed",
        traceback_text="private traceback",
    )

    assert pickle.loads(pickle.dumps(failure)) == failure
    with pytest.raises(FrozenInstanceError):
        failure.message = "changed"


def test_annotation_chunk_result_is_immutable_and_pickleable():
    failure = RecordFailure(
        row_id="abstar_0_1",
        sequence_id="sample-001",
        stage="assignment",
        category="unassigned",
        message="no compatible gene call",
    )
    result = AnnotationChunkResult(
        output_path="chunk.parquet",
        failures=(failure,),
        failed_log_path="chunk.failed",
        succeeded_log_path=None,
    )

    assert pickle.loads(pickle.dumps(result)) == result
    with pytest.raises(FrozenInstanceError):
        result.output_path = "changed.parquet"


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("row_id", []),
        ("sequence_id", object()),
        ("stage", []),
        ("category", object()),
        ("message", ["annotation failed"]),
        ("traceback_text", object()),
    ],
)
def test_record_failure_rejects_non_string_diagnostic_fields(field, value):
    values = {
        "row_id": "abstar_0_1",
        "sequence_id": "sample-001",
        "stage": "annotation",
        "category": "internal_error",
        "message": "annotation failed",
        "traceback_text": None,
    }
    values[field] = value

    with pytest.raises(TypeError, match=field):
        RecordFailure(**values)


@pytest.mark.parametrize(
    ("field", "value"),
    [("stage", "alignment"), ("category", "unexpected")],
)
def test_record_failure_rejects_unknown_failure_literals(field, value):
    values = {
        "row_id": "abstar_0_1",
        "sequence_id": "sample-001",
        "stage": "annotation",
        "category": "internal_error",
        "message": "annotation failed",
    }
    values[field] = value

    with pytest.raises(ValueError, match=field):
        RecordFailure(**values)


def test_annotation_chunk_result_copies_caller_owned_failure_list():
    failure = RecordFailure(
        row_id="abstar_0_1",
        sequence_id="sample-001",
        stage="assignment",
        category="unassigned",
        message="no compatible gene call",
    )
    failures = [failure]

    result = AnnotationChunkResult(
        output_path="chunk.parquet",
        failures=failures,
        failed_log_path=None,
        succeeded_log_path=None,
    )
    failures.clear()

    assert result.failures == (failure,)


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("output_path", object()),
        ("failures", ["not a record failure"]),
        ("failed_log_path", []),
        ("succeeded_log_path", object()),
    ],
)
def test_annotation_chunk_result_rejects_invalid_fields(field, value):
    values = {
        "output_path": "chunk.parquet",
        "failures": (),
        "failed_log_path": None,
        "succeeded_log_path": None,
    }
    values[field] = value

    with pytest.raises(TypeError, match=field):
        AnnotationChunkResult(**values)


@pytest.mark.parametrize(
    ("record_type", "values", "invalid_field"),
    [
        (
            RecordFailure,
            {
                "row_id": "abstar_0_1",
                "sequence_id": "sample-001",
                "stage": "annotation",
                "category": "internal_error",
                "message": ["invalid serialized message"],
                "traceback_text": None,
            },
            "message",
        ),
        (
            AnnotationChunkResult,
            {
                "output_path": "chunk.parquet",
                "failures": ("invalid serialized failure",),
                "failed_log_path": None,
                "succeeded_log_path": None,
            },
            "failures",
        ),
    ],
)
def test_dataclass_pickle_reconstructs_through_validation(
    record_type, values, invalid_field
):
    invalid_record = object.__new__(record_type)
    for field, value in values.items():
        object.__setattr__(invalid_record, field, value)

    payload = pickle.dumps(invalid_record)

    with pytest.raises(TypeError, match=invalid_field):
        pickle.loads(payload)


def test_annotation_run_error_summarizes_counts_without_record_details():
    failures = (
        RecordFailure(
            row_id="internal-row-2",
            sequence_id="private-sequence-2",
            stage="annotation",
            category="internal_error",
            message="second private message",
            traceback_text="second private traceback",
        ),
        RecordFailure(
            row_id="internal-row-1",
            sequence_id="private-sequence-1",
            stage="assignment",
            category="unassigned",
            message="first private message",
            traceback_text="first private traceback",
        ),
        RecordFailure(
            row_id="internal-row-3",
            sequence_id="private-sequence-3",
            stage="annotation",
            category="internal_error",
            message="third private message",
        ),
    )

    error = AnnotationRunError(
        (failure for failure in failures),
        (path for path in ("/private/chunk-1.parquet",)),
    )

    assert str(error) == (
        "abstar run failed: annotation/internal_error=2, assignment/unassigned=1"
    )
    assert error.failures == failures
    assert error.partial_output_paths == ("/private/chunk-1.parquet",)
    for private_text in (
        "internal-row",
        "private-sequence",
        "private message",
        "private traceback",
        "/private/",
    ):
        assert private_text not in str(error)


def test_annotation_run_error_has_clear_empty_summary():
    assert str(AnnotationRunError([])) == (
        "abstar run failed: no record failures reported"
    )


def test_annotation_run_error_copies_and_validates_partial_output_paths():
    partial_output_paths = ["chunk.parquet"]
    error = AnnotationRunError([], partial_output_paths)
    partial_output_paths.clear()

    assert error.partial_output_paths == ("chunk.parquet",)
    for invalid_paths in ([object()], "chunk.parquet"):
        with pytest.raises(TypeError, match="partial_output_paths"):
            AnnotationRunError([], invalid_paths)


def test_only_run_error_is_exported_from_package_root():
    import abstar

    assert abstar.AnnotationRunError is AnnotationRunError
    assert not hasattr(abstar, "RecordFailure")
    assert not hasattr(abstar, "AnnotationChunkResult")


def test_record_conservation_accepts_exactly_accounted_records():
    failure = RecordFailure(
        row_id="abstar_0_4",
        sequence_id="broken",
        stage="annotation",
        category="internal_error",
        message="annotation failed",
    )

    _assert_record_conservation(5, 3, 1, [failure])


def test_record_conservation_raises_structured_run_error_on_mismatch():
    with pytest.raises(AnnotationRunError) as exc_info:
        _assert_record_conservation(5, 3, 1, [])

    assert len(exc_info.value.failures) == 1
    failure = exc_info.value.failures[0]
    assert failure == RecordFailure(
        row_id="run",
        sequence_id="<run>",
        stage="output",
        category="internal_error",
        message="record conservation failed: input=5, accounted=4",
    )


def test_annotate_returns_structured_failure_after_writing_diagnostics(
    tmp_path, monkeypatch
):
    input_path = tmp_path / "assigned.parquet"
    output_directory = tmp_path / "output"
    log_directory = tmp_path / "logs"
    output_directory.mkdir()
    log_directory.mkdir()
    pl.DataFrame(
        {
            "row_id": ["abstar_0_0"],
            "sequence_id": ["duplicate"],
            "sequence_input": ["ACGT"],
            "quality": [""],
            "rev_comp": [False],
            "v_call": ["IGHV3-23*01"],
            "v_support": [1e-20],
            "d_call": [None],
            "d_support": [None],
            "j_call": ["IGHJ4*02"],
            "j_support": [1e-10],
            "c_call": [None],
            "c_support": [None],
        }
    ).write_parquet(input_path)

    def fail_annotation(**kwargs):
        raise RuntimeError("deliberate annotation failure")

    monkeypatch.setattr(
        "abstar.annotation.annotator.annotate_single_sequence", fail_annotation
    )

    result = annotate(
        str(input_path),
        str(output_directory),
        "human",
        log_directory=str(log_directory),
    )

    assert pl.read_parquet(result.output_path).is_empty()
    assert len(result.failures) == 1
    assert result.failures[0].row_id == "abstar_0_0"
    assert result.failures[0].sequence_id == "duplicate"
    assert result.failures[0].stage == "annotation"
    assert result.failures[0].category == "internal_error"
    assert "deliberate annotation failure" in result.failures[0].message
    assert result.failed_log_path is not None
    assert "deliberate annotation failure" in Path(result.failed_log_path).read_text()


def test_unassigned_reverse_strand_preserves_oriented_sequence(tmp_path):
    input_path = tmp_path / "reverse-unassigned.parquet"
    output_directory = tmp_path / "output"
    output_directory.mkdir()
    pl.DataFrame(
        {
            "row_id": ["abstar_0_0"],
            "sequence_id": ["reverse-unassigned"],
            "sequence_input": ["AACCGT"],
            "quality": [""],
            "rev_comp": [True],
            "v_call": ["IGHV3-23*01"],
            "v_support": [1e-20],
            "d_call": [None],
            "d_support": [None],
            "j_call": [None],
            "j_support": [None],
            "c_call": [None],
            "c_support": [None],
        }
    ).write_parquet(input_path)

    result = annotate(str(input_path), str(output_directory), "human")
    record = pl.read_parquet(result.output_path).row(0, named=True)

    assert record["annotation_status"] == "unassigned"
    assert record["sequence_oriented"] == "ACGGTT"
    assert record["sequence"] == "ACGGTT"


def _assert_assignment_failure_category(
    monkeypatch, tmp_path, error, stage, category
):
    def fail_assignment(*args, **kwargs):
        raise error

    monkeypatch.setattr("abstar.core.abstar.MMseqs.__call__", fail_assignment)
    project_path = tmp_path / category
    with pytest.raises(AnnotationRunError) as exc_info:
        from ..core.abstar import run

        run("ACGT", project_path=str(project_path), n_processes=1)

    assert len(exc_info.value.failures) == 1
    assert exc_info.value.failures[0].stage == stage
    assert exc_info.value.failures[0].category == category
    failure_log = project_path / "logs" / "sequences.failed"
    assert failure_log.is_file()
    assert str(error) in failure_log.read_text()


def test_controller_classifies_assignment_programming_error_as_internal(
    monkeypatch, tmp_path
):
    _assert_assignment_failure_category(
        monkeypatch,
        tmp_path,
        ValueError("dataframe cardinality bug"),
        "assignment",
        "internal_error",
    )


def test_controller_classifies_input_error_at_input_boundary(monkeypatch, tmp_path):
    from ..assigners.mmseqs import AssignmentInputError

    _assert_assignment_failure_category(
        monkeypatch,
        tmp_path,
        AssignmentInputError("malformed FASTQ"),
        "preprocess",
        "invalid_input",
    )


def test_controller_classifies_tool_error_at_external_boundary(monkeypatch, tmp_path):
    from ..assigners.mmseqs import AssignmentExternalToolError

    _assert_assignment_failure_category(
        monkeypatch,
        tmp_path,
        AssignmentExternalToolError("MMseqs exited nonzero"),
        "assignment",
        "external_tool",
    )


def test_controller_worker_failure_keeps_diagnostics_and_partial_work(
    monkeypatch, tmp_path, single_hc_sequence
):
    class FailingExecutor:
        def __init__(self, *args, **kwargs):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        def submit(self, *args, **kwargs):
            future = Future()
            future.set_exception(RuntimeError("worker process crashed"))
            return future

    monkeypatch.setattr("abstar.core.abstar.ProcessPoolExecutor", FailingExecutor)
    project_path = tmp_path / "worker-failure"

    with pytest.raises(AnnotationRunError) as exc_info:
        from ..core.abstar import run

        run(
            single_hc_sequence,
            project_path=str(project_path),
            n_processes=1,
        )

    assert len(exc_info.value.failures) == 1
    assert exc_info.value.failures[0].stage == "annotation"
    assert exc_info.value.failures[0].category == "internal_error"
    assert all(Path(path).is_file() for path in exc_info.value.partial_output_paths)
    failure_log = project_path / "logs" / "sequences.failed"
    assert failure_log.is_file()
    assert "worker process crashed" in failure_log.read_text()


def _assert_no_final_outputs(project):
    assert not list((project / "airr").glob("*.tsv"))
    assert not list((project / "parquet").glob("*.parquet"))


@pytest.mark.parametrize("project_output", (False, True))
def test_public_missing_translate_fails_actionably_before_output(
    monkeypatch, tmp_path, single_hc_sequence, project_output,
):
    project = tmp_path / "missing-capability"
    with monkeypatch.context() as patch:
        patch.delattr(abutils.tl, "translate")
        with pytest.raises(AnnotationRunError) as captured:
            abstar.run(
                single_hc_sequence,
                project_path=str(project) if project_output else None,
                output_format=["airr", "parquet"], as_dataframe=True,
                n_processes=1, mmseqs_threads=1,
            )
    assert callable(abutils.tl.translate)
    assert len(captured.value.failures) == 1
    failure = captured.value.failures[0]
    assert (failure.stage, failure.category) == ("preprocess", "internal_error")
    assert "abutils" in failure.message and "translate" in failure.message
    assert "install" in failure.message.lower()
    assert captured.value.partial_output_paths == ()
    assert not project.exists()


def test_cli_checked_mmseqs_failure_retains_command_and_streams(
    monkeypatch, tmp_path, small_fasta_file,
):
    owned_paths = []
    commands = []

    def fail_process(command, **kwargs):
        commands.append(command)
        owned = Path(command[5])
        owned_paths.append(owned)
        assert owned.is_dir()
        (owned / "partial-search").write_text("interrupted work")
        raise subprocess.CalledProcessError(
            37, command, output="TASK17-STDOUT", stderr="TASK17-STDERR",
        )

    monkeypatch.setattr(subprocess, "run", fail_process)
    monkeypatch.setenv("TASK17_UNRELATED_SECRET", "must-never-appear-in-diagnostics")
    project = tmp_path / "quoted ' project; $(false)"
    result = CliRunner().invoke(cli, [
        "run", small_fasta_file, str(project), "--n_processes", "1",
        "--mmseqs_threads", "1", "--quiet",
    ])
    assert result.exit_code != 0, result.output
    assert "assignment/external_tool=1" in result.output
    diagnostic_path = project / "logs" / "test_sequences.failed"
    diagnostic = diagnostic_path.read_text()
    assert str(diagnostic_path) in result.output
    assert "TASK17-STDOUT" in diagnostic
    assert "TASK17-STDERR" in diagnostic
    assert "37" in diagnostic
    assert shlex.join(commands[0]) in diagnostic
    assert "must-never-appear-in-diagnostics" not in diagnostic
    assert owned_paths and all(not path.exists() for path in owned_paths)
    _assert_no_final_outputs(project)


def test_checked_mmseqs_success_uses_argv_and_cleans_owned_scratch(monkeypatch, tmp_path):
    from ..assigners.mmseqs import MMseqs

    query = tmp_path / "query ' ; $(false).fasta"
    target = tmp_path / "target with spaces"
    output = tmp_path / "hits.tsv"
    scratch_paths = []

    def completed_process(command, *, check, capture_output, text):
        assert isinstance(command, list)
        assert command[1:5] == ["easy-search", str(query), str(target), str(output)]
        assert command[6:] == [
            "--search-type", "3", "-s", "7.5", "--max-seqs", "25",
            "-e", "1000.0", "--format-mode", "4", "--format-output",
            "query,target,bits", "--threads", "1", "--min-aln-len", "12",
            "-k", "5", "--alignment-mode", "3",
        ]
        assert check and capture_output and text
        scratch = Path(command[5])
        assert scratch.is_dir()
        scratch_paths.append(scratch)
        (scratch / "search-work").write_text("owned work")
        output.write_text("query\ttarget\tbits\nabstar_0_0\tIGHV3-23*01\t100\n")
        return subprocess.CompletedProcess(command, 0, "search completed", "")

    monkeypatch.setattr(subprocess, "run", completed_process)
    log = tmp_path / "mmseqs.log"
    MMseqs._run_mmseqs_search(
        query=str(query), target=str(target), output_path=str(output), search_type=3,
        sensitivity=7.5, max_seqs=25, max_evalue=1000.0, format_mode=4,
        format_output="query,target,bits", threads=1,
        additional_cli_args="--min-aln-len 12 -k 5 --alignment-mode 3", log_to=str(log),
    )
    assert output.read_text().endswith("abstar_0_0\tIGHV3-23*01\t100\n")
    assert "search completed" in log.read_text()
    assert scratch_paths and all(not path.exists() for path in scratch_paths)


@pytest.mark.parametrize("content", ("", " \t\n"))
def test_empty_fasta_is_rejected_before_project_creation(tmp_path, content):
    source = tmp_path / "empty.fasta"
    source.write_text(content)
    project = tmp_path / "absent-project"
    with pytest.raises(ValueError, match="empty"):
        abstar.run(str(source), project_path=str(project), n_processes=1)
    assert not project.exists()


@pytest.mark.parametrize("suffix,content,diagnostic", (
    ("fastq", "@00123\nACGT\n+\nIII\n", "parse input"),
    ("fastq", "@10E8\nACGT\n+\nIIII\n@00123\nACGT\n+\nIII\n", "parse input"),
    ("fasta", ">10E8\nACGTX\n", "10E8 contains non-IUPAC characters: X"),
))
def test_invalid_sequence_content_retains_distinct_input_diagnostic(
    tmp_path, suffix, content, diagnostic,
):
    source = tmp_path / f"invalid.{suffix}"
    source.write_text(content)
    project = tmp_path / "invalid-project"
    with pytest.raises(AnnotationRunError) as captured:
        abstar.run(str(source), project_path=str(project),
                   output_format=["airr", "parquet"], n_processes=1)
    assert len(captured.value.failures) == 1
    failure = captured.value.failures[0]
    assert (failure.stage, failure.category) == ("preprocess", "invalid_input")
    assert diagnostic in failure.message
    assert diagnostic in (project / "logs" / "invalid.failed").read_text()
    _assert_no_final_outputs(project)


@pytest.mark.parametrize("missing", ("ungapped/v.fasta", "imgt_gapped/v.fasta", "mmseqs/v.dbtype"))
def test_incomplete_database_fails_at_assignment_with_context(
    monkeypatch, tmp_path, small_fasta_file, missing,
):
    from ..annotation.germline import get_germline_database_path

    source = Path(get_germline_database_path("human", "tcr"))
    database = tmp_path / "germline_dbs" / "tcr" / "incomplete"
    database.mkdir(parents=True)
    for directory in ("ungapped", "imgt_gapped", "mmseqs"):
        (database / directory).mkdir()
        for path in (source / directory).iterdir():
            if str(Path(directory) / path.name) != missing:
                (database / directory / path.name).symlink_to(path)
    monkeypatch.setattr(
        "abstar.assigners.assigner.get_germline_database_path",
        lambda germdb_name, receptor: str(database),
    )
    project = tmp_path / "invalid-database"
    with pytest.raises(AnnotationRunError) as captured:
        abstar.run(str(small_fasta_file), project_path=str(project),
                   germline_database="incomplete", receptor="tcr", n_processes=1)
    assert len(captured.value.failures) == 1
    failure = captured.value.failures[0]
    assert (failure.stage, failure.category) == ("assignment", "invalid_input")
    for text in ("incomplete", "tcr", missing):
        assert text in failure.message
        assert text in (project / "logs" / "test_sequences.failed").read_text()
    _assert_no_final_outputs(project)


def test_public_valid_multiline_fastq_retains_annotation(public_bcr_cases, tmp_path):
    case = public_bcr_cases[0]
    source = tmp_path / "multiline.fastq"
    source.write_text(
        f"@00123\n{case.sequence[:90]}\n{case.sequence[90:]}\n+\n"
        f"{'I' * 90}\n{'I' * (len(case.sequence) - 90)}\n"
    )
    row = abstar.run(str(source), as_dataframe=True, n_processes=1, mmseqs_threads=1).row(0, named=True)
    assert row["sequence_id"] == "00123"
    assert row["sequence_input"] == case.sequence
    assert row["annotation_status"] == "annotated"
    for segment in ("v", "j"):
        assert row[f"{segment}_gene"] in case.expected[f"{segment}_call"]
    for field in ("junction", "productive"):
        assert row[field] == case.expected[field]


def _crash_chunk_with_second_record(input_file, **kwargs):
    frame = pl.read_parquet(input_file)
    if "abstar_0_1" in frame["row_id"].to_list():
        raise RuntimeError("TASK17-WORKER-CRASH")
    return annotate(input_file, **kwargs)


def _fail_second_annotation(**kwargs):
    if kwargs["ab"].row_id == "abstar_0_1":
        raise RuntimeError("TASK17-RECORD-CRASH")
    return annotate_single_sequence(**kwargs)


def _annotate_with_second_record_failure(input_file, **kwargs):
    # This top-level worker is spawn-pickleable. Only the injected per-record
    # operation changes; input loading, accounting, writing and futures are real.
    import importlib
    module = importlib.import_module("abstar.annotation.annotator")
    original = module.annotate_single_sequence
    module.annotate_single_sequence = _fail_second_annotation
    try:
        return annotate(input_file, **kwargs)
    finally:
        module.annotate_single_sequence = original


@pytest.mark.parametrize("n_processes,chunksize", ((1, 1), (2, 1), (1, 2), (2, 2)))
@pytest.mark.parametrize("crash_scope", ("worker", "record"))
def test_public_mixed_worker_failures_account_for_every_row(
    monkeypatch, tmp_path, public_bcr_cases, n_processes, chunksize, crash_scope,
):
    records = [abutils.Sequence(case.sequence, id=identifier)
               for case, identifier in zip(public_bcr_cases[:2], ("10E8", '00123"opaque'))]
    replacement = (_crash_chunk_with_second_record if crash_scope == "worker"
                   else _annotate_with_second_record_failure)
    monkeypatch.setattr("abstar.core.abstar.annotate", replacement)
    project = tmp_path / "mixed-failure"
    with pytest.raises(AnnotationRunError) as captured:
        abstar.run(records, project_path=str(project), output_format=["airr", "parquet"],
                   n_processes=n_processes, chunksize=chunksize, mmseqs_threads=1)
    error = captured.value
    failed_indices = [0, 1] if crash_scope == "worker" and chunksize == 2 else [1]
    assert [(f.row_id, f.sequence_id) for f in error.failures] == [
        (f"abstar_0_{i}", records[i].id) for i in failed_indices]
    diagnostic = (project / "logs" / "sequences.failed").read_text()
    for failure in error.failures:
        assert (failure.stage, failure.category) == ("annotation", "internal_error")
        assert f"TASK17-{crash_scope.upper()}-CRASH" in failure.message
        assert failure.row_id in diagnostic
        assert failure.sequence_id in diagnostic
    assert error.partial_output_paths
    survivors = pl.concat([pl.read_parquet(path) for path in error.partial_output_paths])
    assert survivors["row_id"].to_list() == [
        f"abstar_0_{i}" for i in range(2) if i not in failed_indices]
    assert survivors.height + len(error.failures) == 2
    if survivors.height:
        survivor = survivors.row(0, named=True)
        assert survivor["sequence_id"] == "10E8"
        assert survivor["locus"] == "IGH"
        assert survivor["germline_database"] == "human"
        assert survivor["junction"] == public_bcr_cases[0].expected["junction"]
        assert survivor["productive"] is True
    _assert_no_final_outputs(project)
    assert not list((project / "tmp").glob("mmseqs-*"))


@pytest.mark.parametrize("writer", ("airr", "parquet"))
def test_public_output_failure_retains_diagnostics_without_final_success(
    monkeypatch, tmp_path, public_bcr_cases, writer,
):
    project = tmp_path / "writer-failure"
    if writer == "airr":
        def fail_airr(frame, path):
            Path(path).write_text("sequence_id\tsequence\npartial")
            raise PermissionError("TASK17-AIRR-UNWRITABLE")

        monkeypatch.setattr("abstar.core.abstar.write_airr_tsv", fail_airr)
    else:
        original = pl.DataFrame.write_parquet

        def fail_final_parquet(frame, path, *args, **kwargs):
            # Internal work Parquets include row_id; only the public writer
            # fails, after leaving partial bytes as a real I/O failure can.
            if "row_id" not in frame.columns:
                Path(path).write_bytes(b"PAR1partial")
                raise PermissionError("TASK17-PARQUET-UNWRITABLE")
            return original(frame, path, *args, **kwargs)

        monkeypatch.setattr(pl.DataFrame, "write_parquet", fail_final_parquet)
    with pytest.raises(AnnotationRunError) as captured:
        abstar.run(public_bcr_cases[0].as_sequence(), project_path=str(project),
                   output_format=["airr", "parquet"], n_processes=1, mmseqs_threads=1)
    error = captured.value
    assert len(error.failures) == 1
    failure = error.failures[0]
    assert (failure.stage, failure.category) == ("output", "internal_error")
    assert f"TASK17-{writer.upper()}-UNWRITABLE" in failure.message
    assert failure.message in (project / "logs" / "sequences.failed").read_text()
    assert error.partial_output_paths
    assert all(Path(path).is_file() for path in error.partial_output_paths)
    _assert_no_final_outputs(project)
    assert not list((project / "tmp").glob("output-*"))
    assert not list((project / "tmp").glob("mmseqs-*"))


def test_output_directory_creation_error_is_structured_and_retains_diagnostic(
    tmp_path, small_fasta_file,
):
    project = tmp_path / "blocked-output"
    project.mkdir()
    blocker = project / "airr"
    blocker.write_text("caller-owned file")
    with pytest.raises(AnnotationRunError) as captured:
        abstar.run(small_fasta_file, project_path=str(project), n_processes=1)
    failure, = captured.value.failures
    assert (failure.stage, failure.category) == ("output", "internal_error")
    assert str(blocker) in failure.message
    assert blocker.read_text() == "caller-owned file"
    assert any(failure.message in path.read_text() for path in (project / "logs").glob("*.failed"))
    _assert_no_final_outputs(project)


@pytest.fixture
def owned_temporary_paths(monkeypatch, tmp_path):
    paths = []
    original = tempfile.TemporaryDirectory
    monkeypatch.setattr(tempfile, "tempdir", str(tmp_path))

    def track_temporary_directory(*args, **kwargs):
        directory = original(*args, **kwargs)
        paths.append(Path(directory.name))
        return directory

    monkeypatch.setattr(tempfile, "TemporaryDirectory", track_temporary_directory)
    return paths


@pytest.mark.parametrize("outcome", ("success", "assignment", "worker"))
def test_api_owned_workspace_cleanup_preserves_failure_artifacts(
    monkeypatch, public_bcr_cases, owned_temporary_paths, outcome,
):
    records = [case.as_sequence() for case in public_bcr_cases[:2]]
    if outcome == "assignment":
        def fail_search(*args, **kwargs):
            raise subprocess.CalledProcessError(
                41, args[0], output="TASK17-API-STDOUT", stderr="TASK17-API-STDERR",
            )
        monkeypatch.setattr(subprocess, "run", fail_search)
    elif outcome == "worker":
        monkeypatch.setattr("abstar.core.abstar.annotate", _crash_chunk_with_second_record)
    if outcome == "success":
        result = abstar.run(records, as_dataframe=True, n_processes=2, chunksize=1, mmseqs_threads=1)
        assert result["sequence_id"].to_list() == [record.id for record in records]
        assert result["annotation_status"].to_list() == ["annotated", "annotated"]
    else:
        with pytest.raises(AnnotationRunError) as captured:
            abstar.run(records, as_dataframe=True, n_processes=2, chunksize=1, mmseqs_threads=1)
        artifacts = [Path(path) for path in captured.value.partial_output_paths]
        assert artifacts and all(path.is_file() for path in artifacts)
        logs = [path for path in artifacts if path.suffix == ".failed"]
        assert logs
        expected = "TASK17-API-STDERR" if outcome == "assignment" else "TASK17-WORKER-CRASH"
        assert any(expected in path.read_text() for path in logs)
        assert all(not any(path.is_relative_to(owned) for owned in owned_temporary_paths)
                   for path in artifacts)
    assert owned_temporary_paths
    assert all(not path.exists() for path in owned_temporary_paths)


def test_cli_output_failure_reports_surviving_diagnostic(monkeypatch, tmp_path, small_fasta_file):
    def fail_airr(frame, path):
        raise PermissionError("TASK17-CLI-OUTPUT")

    monkeypatch.setattr("abstar.core.abstar.write_airr_tsv", fail_airr)
    project = tmp_path / "cli-output-failure"
    result = CliRunner().invoke(cli, [
        "run", small_fasta_file, str(project), "--n_processes", "1",
        "--mmseqs_threads", "1", "--quiet",
    ])
    diagnostic = project / "logs" / "test_sequences.failed"
    assert result.exit_code != 0
    assert "output/internal_error=1" in result.output
    assert str(diagnostic) in result.output
    assert "TASK17-CLI-OUTPUT" in diagnostic.read_text()
    _assert_no_final_outputs(project)


def test_real_mmseqs_accepts_quoted_project_paths(tmp_path, public_bcr_cases):
    project = tmp_path / "quotes ' and spaces; literal $(false)"
    abstar.run(public_bcr_cases[0].as_sequence(), project_path=str(project),
               output_format=["airr", "parquet"], n_processes=1, mmseqs_threads=1)
    row = pl.read_parquet(project / "parquet" / "sequences.parquet").row(0, named=True)
    assert row["sequence_id"] == public_bcr_cases[0].sequence_id
    assert row["junction"] == public_bcr_cases[0].expected["junction"]
    assert row["productive"] is True
    assert not list((project / "tmp").glob("mmseqs-*"))


def _crash_second_sample(input_file, **kwargs):
    if pl.read_parquet(input_file)["row_id"][0] == "abstar_1_0":
        raise RuntimeError("TASK17-SECOND-SAMPLE")
    return annotate(input_file, **kwargs)


def test_later_sample_failure_reports_earlier_files_as_partial(
    monkeypatch, tmp_path, public_bcr_cases,
):
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    for i, case in enumerate(public_bcr_cases[:2]):
        (inputs / f"sample{i}.fasta").write_text(f">{case.sequence_id}\n{case.sequence}\n")
    monkeypatch.setattr("abstar.core.abstar.annotate", _crash_second_sample)
    project = tmp_path / "partial-run"
    with pytest.raises(AnnotationRunError) as captured:
        abstar.run(str(inputs), project_path=str(project), output_format=["airr", "parquet"],
                   n_processes=1, mmseqs_threads=1)
    failure, = captured.value.failures
    assert failure.row_id == "abstar_1_0"
    assert failure.sequence_id == public_bcr_cases[1].sequence_id
    for path in (project / "airr" / "sample0.tsv", project / "parquet" / "sample0.parquet"):
        assert path.is_file()
        assert str(path) in captured.value.partial_output_paths
    assert not (project / "airr" / "sample1.tsv").exists()
    assert not (project / "parquet" / "sample1.parquet").exists()


@pytest.mark.parametrize("segment", ("d", "c"))
@pytest.mark.parametrize("remaining", ("other-components", "gapped-only", "sidecar-only"))
def test_partial_optional_database_segment_is_invalid_input(
    monkeypatch, tmp_path, small_fasta_file, segment, remaining,
):
    from ..annotation.germline import get_germline_database_path

    source = Path(get_germline_database_path("human", "bcr"))
    database = tmp_path / "incomplete"
    kept = {"gapped-only": f"imgt_gapped/{segment}.fasta",
            "sidecar-only": f"mmseqs/{segment}.lookup"}
    for directory in ("ungapped", "imgt_gapped", "mmseqs"):
        (database / directory).mkdir(parents=True)
        for path in (source / directory).iterdir():
            relative = f"{directory}/{path.name}"
            optional = path.name == segment or path.name.startswith((f"{segment}.", f"{segment}_"))
            if optional:
                if remaining == "other-components":
                    if relative in (f"ungapped/{segment}.fasta", f"mmseqs/{segment}"):
                        continue
                elif relative != kept[remaining]:
                    continue
            (database / directory / path.name).symlink_to(path)
    monkeypatch.setattr(
        "abstar.assigners.assigner.get_germline_database_path",
        lambda germdb_name, receptor: str(database),
    )
    project = tmp_path / "failed-project"
    with pytest.raises(AnnotationRunError) as captured:
        abstar.run(small_fasta_file, project_path=str(project), germline_database="incomplete",
                   receptor="bcr", n_processes=1, mmseqs_threads=1)
    failure, = captured.value.failures
    assert (failure.stage, failure.category) == ("assignment", "invalid_input")
    for detail in ("incomplete", "bcr", f"ungapped/{segment}.fasta", f"mmseqs/{segment}"):
        assert detail in failure.message
        assert detail in (project / "logs" / "test_sequences.failed").read_text()
    _assert_no_final_outputs(project)


@pytest.mark.parametrize("receptor,locus", (("bcr", "IGK"), ("bcr", "IGL"), ("tcr", "TRA"), ("tcr", "TRG")))
def test_vj_only_database_without_optional_components_annotates_vj_chains(
    monkeypatch, tmp_path, public_bcr_cases, receptor, locus,
):
    from ..annotation.germline import get_germline_database_path

    source = Path(get_germline_database_path("human", receptor))
    database = tmp_path / "vj-only"
    for directory in ("ungapped", "imgt_gapped", "mmseqs"):
        (database / directory).mkdir(parents=True)
        for path in (source / directory).iterdir():
            if path.name[0] in ("v", "j"):
                (database / directory / path.name).symlink_to(path)
    monkeypatch.setattr(
        "abstar.assigners.assigner.get_germline_database_path",
        lambda germdb_name, receptor: str(database),
    )
    if receptor == "bcr":
        case = public_bcr_cases[1 if locus == "IGK" else 2]
        expected = case.expected
    else:
        from .test_tcr_e2e import load_tcr_cases
        case = next(case for case in load_tcr_cases() if case.locus == locus)
        expected = {"junction": case.junction, "productive": case.productive}
    # The selected assignment sources are the unchanged human V/J genes; the
    # spawn worker uses their same packaged gapped counterparts for realignment.
    row = abstar.run(case.as_sequence(), as_dataframe=True, germline_database="human",
                     receptor=receptor, n_processes=1, mmseqs_threads=1).row(0, named=True)
    assert row["sequence_id"] == case.sequence_id
    assert row["locus"] == locus
    assert row["d_call"] is None and row["c_call"] is None
    assert row["annotation_status"] == "annotated"
    assert row["germline_database"] == "human"
    for field in ("junction", "productive"):
        assert row[field] == expected[field]


@pytest.mark.parametrize("blocked", ("project", "logs", "tmp"))
@pytest.mark.parametrize("entrypoint", ("api-file", "api-sequence", "cli"))
def test_initial_project_storage_failure_is_structured_and_preserves_caller_file(
    monkeypatch, tmp_path, small_fasta_file, single_hc_sequence, blocked, entrypoint,
):
    monkeypatch.setattr(tempfile, "tempdir", str(tmp_path))
    project = tmp_path / "blocked-project"
    if blocked == "project":
        blocker = project
    else:
        project.mkdir()
        blocker = project / blocked
    original_bytes = b"caller-owned\x00file\nunchanged"
    blocker.write_bytes(original_bytes)
    if entrypoint == "cli":
        result = CliRunner().invoke(cli, [
            "run", small_fasta_file, str(project), "--n_processes", "1", "--quiet",
        ])
        assert result.exit_code != 0
        assert "output/internal_error=1" in result.output
        error = result.exception
        while error is not None and not isinstance(error, AnnotationRunError):
            error = error.__cause__ or error.__context__
        assert isinstance(error, AnnotationRunError)
    else:
        with pytest.raises(AnnotationRunError) as captured:
            abstar.run(small_fasta_file if entrypoint == "api-file" else single_hc_sequence,
                       project_path=str(project), n_processes=1)
        error = captured.value
    failure, = error.failures
    assert (failure.stage, failure.category) == ("output", "internal_error")
    assert str(blocker) in failure.message
    assert isinstance(error.__cause__, OSError)
    assert blocker.read_bytes() == original_bytes
    artifacts = [Path(path) for path in error.partial_output_paths]
    if blocked in ("project", "logs"):
        assert artifacts
        assert all(not path.is_relative_to(project) for path in artifacts)
    else:
        artifacts.extend((project / "logs").glob("*.failed"))
    assert artifacts and all(path.is_file() for path in artifacts)
    assert any(failure.message in path.read_text() for path in artifacts)
    if entrypoint == "cli":
        assert any(str(path) in result.output for path in artifacts)
    _assert_no_final_outputs(project)


@pytest.mark.parametrize("failure_point", ("allocate", "write"))
def test_fallback_diagnostic_failure_preserves_original_output_error(
    monkeypatch, tmp_path, small_fasta_file, failure_point,
):
    import builtins

    monkeypatch.setattr(tempfile, "tempdir", str(tmp_path))
    project = tmp_path / "blocked-project"
    project.write_bytes(b"caller file")
    if failure_point == "allocate":
        def fail_allocation(*args, **kwargs):
            raise PermissionError("TASK17-FALLBACK-ALLOCATION")

        monkeypatch.setattr(tempfile, "mkdtemp", fail_allocation)
        sentinel = "TASK17-FALLBACK-ALLOCATION"
    else:
        original_open = builtins.open

        def fail_fallback_write(path, *args, **kwargs):
            if Path(path).parent.name.startswith("abstar-failed-"):
                with original_open(path, "w") as handle:
                    handle.write("partial diagnostic")
                raise PermissionError("TASK17-FALLBACK-WRITE")
            return original_open(path, *args, **kwargs)

        monkeypatch.setattr(builtins, "open", fail_fallback_write)
        sentinel = "TASK17-FALLBACK-WRITE"
    with pytest.raises(AnnotationRunError) as captured:
        abstar.run(small_fasta_file, project_path=str(project), n_processes=1)
    error = captured.value
    failure, = error.failures
    assert (failure.stage, failure.category) == ("output", "internal_error")
    assert isinstance(error.__cause__, NotADirectoryError)
    assert str(project) in failure.message
    assert sentinel in failure.message
    assert str(project) in failure.traceback_text
    assert "NotADirectoryError" in failure.traceback_text
    assert error.partial_output_paths == ()
    assert not list(tmp_path.glob("abstar-failed-*"))
    assert project.read_bytes() == b"caller file"
