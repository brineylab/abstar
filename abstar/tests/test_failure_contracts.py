# Copyright (c) 2026 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

import pickle
from dataclasses import FrozenInstanceError

import pytest

from ..core.results import AnnotationChunkResult, AnnotationRunError, RecordFailure


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


def test_only_run_error_is_exported_from_package_root():
    import abstar

    assert abstar.AnnotationRunError is AnnotationRunError
    assert not hasattr(abstar, "RecordFailure")
    assert not hasattr(abstar, "AnnotationChunkResult")
