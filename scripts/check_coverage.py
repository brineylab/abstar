#!/usr/bin/env python3
"""Enforce per-file combined statement and branch coverage floors."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any


def _load_json(path: Path, label: str) -> tuple[Any | None, str | None]:
    try:
        with path.open(encoding="utf-8") as handle:
            return json.load(handle), None
    except FileNotFoundError:
        return None, f"{label} not found: {path}"
    except json.JSONDecodeError:
        return None, f"{label} is not valid JSON: {path}"
    except (OSError, UnicodeError) as error:
        return None, f"could not read {label}: {path}: {error}"


def _normalized_key(value: str) -> str:
    key = value.replace("\\", "/")
    while key.startswith("./"):
        key = key[2:]
    return key


def _normalize_mapping(
    values: dict[str, Any], label: str
) -> tuple[dict[str, Any] | None, str | None]:
    normalized: dict[str, Any] = {}
    for key, value in values.items():
        if not isinstance(key, str) or not _normalized_key(key):
            return None, f"{label} file keys must be nonempty strings"
        normalized_key = _normalized_key(key)
        if normalized_key in normalized:
            return None, f"{label} contains duplicate normalized path: {normalized_key}"
        normalized[normalized_key] = value
    return normalized, None


def _coverage_files(report: Any) -> tuple[dict[str, Any] | None, str | None]:
    if not isinstance(report, dict):
        return None, "coverage report must be an object"
    meta = report.get("meta")
    if not isinstance(meta, dict) or meta.get("branch_coverage") is not True:
        return None, "coverage report must contain branch coverage"
    files = report.get("files")
    if not isinstance(files, dict):
        return None, "coverage report field 'files' must be an object"
    return _normalize_mapping(files, "coverage report")


def _coverage_floors(floors: Any) -> tuple[dict[str, Any] | None, str | None]:
    if not isinstance(floors, dict):
        return None, "coverage floors must be an object"
    files = floors.get("files")
    if not isinstance(files, dict):
        return None, "coverage floors field 'files' must be an object"
    normalized, error = _normalize_mapping(files, "coverage floors")
    if error is not None:
        return None, error
    assert normalized is not None
    for name in sorted(normalized):
        floor = normalized[name]
        if isinstance(floor, bool) or not isinstance(floor, int) or not 0 <= floor <= 100:
            return (
                None,
                f"coverage floor for {name} must be an integer from 0 through 100",
            )
    return normalized, None


def check_coverage(report_path, floors_path) -> list[str]:
    """Return sorted diagnostics for critical files below their coverage floors."""
    report_path = Path(report_path)
    floors_path = Path(floors_path)

    report, error = _load_json(report_path, "coverage report")
    if error is not None:
        return [error]
    floors, error = _load_json(floors_path, "coverage floors")
    if error is not None:
        return [error]

    report_files, error = _coverage_files(report)
    if error is not None:
        return [error]
    floor_files, error = _coverage_floors(floors)
    if error is not None:
        return [error]
    assert report_files is not None
    assert floor_files is not None

    diagnostics = []
    for name, floor in floor_files.items():
        if name not in report_files:
            diagnostics.append(
                f"{name}: missing from coverage report (floor {floor}%)"
            )
            continue
        file_report = report_files[name]
        if not isinstance(file_report, dict) or not isinstance(
            file_report.get("summary"), dict
        ):
            return [f"coverage report entry for {name} must contain a summary object"]
        percent = file_report["summary"].get("percent_covered")
        if (
            isinstance(percent, bool)
            or not isinstance(percent, (int, float))
            or not math.isfinite(percent)
        ):
            return [f"coverage report percentage for {name} must be a number"]
        if percent < floor:
            diagnostics.append(f"{name}: {percent}% is below {floor}% floor")
    return sorted(diagnostics)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Check coverage.py JSON against per-file coverage floors."
    )
    parser.add_argument("report", help="coverage.py JSON report")
    parser.add_argument("floors", help="per-file coverage floor JSON")
    args = parser.parse_args(argv)

    diagnostics = check_coverage(args.report, args.floors)
    if diagnostics:
        print(*diagnostics, sep="\n")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
