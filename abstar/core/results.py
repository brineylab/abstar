# Copyright (c) 2026 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

from collections import Counter
from dataclasses import dataclass
from typing import Iterable, Literal, get_args


FailureStage = Literal["preprocess", "assignment", "annotation", "output"]
FailureCategory = Literal[
    "invalid_input", "unassigned", "external_tool", "internal_error"
]


def _validate_string(value: object, field_name: str, *, optional: bool = False):
    if optional and value is None:
        return
    if not isinstance(value, str):
        expected = "a string or None" if optional else "a string"
        raise TypeError(f"{field_name} must be {expected}")


def _normalize_failures(failures: Iterable["RecordFailure"]):
    try:
        normalized = tuple(failures)
    except TypeError:
        raise TypeError(
            "failures must be an iterable of RecordFailure objects"
        ) from None
    if not all(isinstance(failure, RecordFailure) for failure in normalized):
        raise TypeError("failures must contain only RecordFailure objects")
    return normalized


def _normalize_paths(paths: Iterable[str], field_name: str):
    if isinstance(paths, str):
        raise TypeError(f"{field_name} must be an iterable of strings")
    try:
        normalized = tuple(paths)
    except TypeError:
        raise TypeError(f"{field_name} must be an iterable of strings") from None
    if not all(isinstance(path, str) for path in normalized):
        raise TypeError(f"{field_name} must contain only strings")
    return normalized


@dataclass(frozen=True, slots=True)
class RecordFailure:
    row_id: str
    sequence_id: str
    stage: FailureStage
    category: FailureCategory
    message: str
    traceback_text: str | None = None

    def __post_init__(self):
        _validate_string(self.row_id, "row_id")
        _validate_string(self.sequence_id, "sequence_id")
        _validate_string(self.stage, "stage")
        _validate_string(self.category, "category")
        _validate_string(self.message, "message")
        _validate_string(self.traceback_text, "traceback_text", optional=True)
        if self.stage not in get_args(FailureStage):
            raise ValueError(f"invalid failure stage: {self.stage!r}")
        if self.category not in get_args(FailureCategory):
            raise ValueError(f"invalid failure category: {self.category!r}")

    def __reduce__(self):
        return (
            type(self),
            (
                self.row_id,
                self.sequence_id,
                self.stage,
                self.category,
                self.message,
                self.traceback_text,
            ),
        )


@dataclass(frozen=True, slots=True)
class AnnotationChunkResult:
    output_path: str
    failures: tuple[RecordFailure, ...]
    failed_log_path: str | None
    succeeded_log_path: str | None

    def __post_init__(self):
        _validate_string(self.output_path, "output_path")
        failures = _normalize_failures(self.failures)
        _validate_string(self.failed_log_path, "failed_log_path", optional=True)
        _validate_string(
            self.succeeded_log_path, "succeeded_log_path", optional=True
        )
        object.__setattr__(self, "failures", failures)

    def __reduce__(self):
        return (
            type(self),
            (
                self.output_path,
                self.failures,
                self.failed_log_path,
                self.succeeded_log_path,
            ),
        )


class AnnotationRunError(RuntimeError):
    def __init__(
        self,
        failures: Iterable[RecordFailure],
        partial_output_paths: Iterable[str] = (),
    ):
        self.failures = _normalize_failures(failures)
        self.partial_output_paths = _normalize_paths(
            partial_output_paths, "partial_output_paths"
        )
        counts = Counter(
            (failure.stage, failure.category) for failure in self.failures
        )
        summary = ", ".join(
            f"{stage}/{category}={count}"
            for (stage, category), count in sorted(counts.items())
        )
        if not summary:
            summary = "no record failures reported"
        super().__init__(f"abstar run failed: {summary}")
