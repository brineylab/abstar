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


@dataclass(frozen=True, slots=True)
class RecordFailure:
    row_id: str
    sequence_id: str
    stage: FailureStage
    category: FailureCategory
    message: str
    traceback_text: str | None = None

    def __post_init__(self):
        if self.stage not in get_args(FailureStage):
            raise ValueError(f"invalid failure stage: {self.stage!r}")
        if self.category not in get_args(FailureCategory):
            raise ValueError(f"invalid failure category: {self.category!r}")


@dataclass(frozen=True, slots=True)
class AnnotationChunkResult:
    output_path: str
    failures: tuple[RecordFailure, ...]
    failed_log_path: str | None
    succeeded_log_path: str | None


class AnnotationRunError(RuntimeError):
    def __init__(
        self,
        failures: Iterable[RecordFailure],
        partial_output_paths: Iterable[str] = (),
    ):
        self.failures = tuple(failures)
        self.partial_output_paths = tuple(partial_output_paths)
        counts = Counter(
            (failure.stage, failure.category) for failure in self.failures
        )
        summary = ", ".join(
            f"{stage}/{category}={count}"
            for (stage, category), count in sorted(counts.items())
        )
        super().__init__(f"abstar run failed: {summary}")
