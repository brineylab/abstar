import json
import subprocess
import sys
from pathlib import Path

import pytest

from scripts.check_coverage import check_coverage


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
CHECKER = REPOSITORY_ROOT / "scripts" / "check_coverage.py"


def _write_json(path, value):
    path.write_text(json.dumps(value), encoding="utf-8")
    return path


def _report(files):
    return {
        "meta": {"branch_coverage": True},
        "files": {
            name: {"summary": {"percent_covered": percent}}
            for name, percent in files.items()
        },
    }


def _floors(files):
    return {"files": files}


def test_check_coverage_accepts_modules_at_or_above_floor(tmp_path):
    report = _write_json(
        tmp_path / "coverage.json",
        _report(
            {
                "abstar/a.py": 80,
                "abstar/b.py": 91.125,
                "abstar/zero.py": 0,
                "abstar/full.py": 100,
            }
        ),
    )
    floors = _write_json(
        tmp_path / "floors.json",
        _floors(
            {
                "abstar/a.py": 80,
                "abstar/b.py": 91,
                "abstar/zero.py": 0,
                "abstar/full.py": 100,
            }
        ),
    )

    assert check_coverage(report, floors) == []


def test_check_coverage_sorts_below_floor_and_missing_diagnostics(tmp_path):
    report = _write_json(
        tmp_path / "coverage.json",
        _report({"abstar/z.py": 79.125}),
    )
    floors = _write_json(
        tmp_path / "floors.json",
        _floors({"abstar/z.py": 80, "abstar/a.py": 75}),
    )

    assert check_coverage(report, floors) == [
        "abstar/a.py: missing from coverage report (floor 75%)",
        "abstar/z.py: 79.125% is below 80% floor",
    ]


def test_check_coverage_normalizes_file_keys_across_platforms(tmp_path):
    report = _write_json(
        tmp_path / "coverage.json",
        _report({r".\abstar\annotation\umi.py": 80.05115089514067}),
    )
    floors = _write_json(
        tmp_path / "floors.json",
        _floors({"abstar/annotation/umi.py": 80}),
    )

    assert check_coverage(report, floors) == []


@pytest.mark.parametrize("source", ["report", "floors"])
def test_check_coverage_rejects_normalized_path_aliases(tmp_path, source):
    report = _write_json(
        tmp_path / "coverage.json",
        _report(
            {
                r"abstar\annotation\umi.py": 80,
                "abstar/annotation/umi.py": 81,
            }
        ),
    )
    floors = _write_json(
        tmp_path / "floors.json",
        _floors({"abstar/annotation/umi.py": 80}),
    )
    if source == "floors":
        report = _write_json(
            tmp_path / "coverage.json",
            _report({"abstar/annotation/umi.py": 80}),
        )
        floors = _write_json(
            tmp_path / "floors.json",
            _floors(
                {
                    r"abstar\annotation\umi.py": 80,
                    "abstar/annotation/umi.py": 81,
                }
            ),
        )

    assert check_coverage(report, floors) == [
        f"coverage {source} contains duplicate normalized path: "
        "abstar/annotation/umi.py"
    ]


def test_check_coverage_reports_missing_and_malformed_inputs(tmp_path):
    missing = tmp_path / "missing.json"
    floors = _write_json(tmp_path / "floors.json", _floors({"abstar/a.py": 80}))
    malformed = tmp_path / "malformed.json"
    malformed.write_text("{", encoding="utf-8")
    bad_schema = _write_json(tmp_path / "bad-schema.json", {"files": []})

    assert check_coverage(missing, floors) == [
        f"coverage report not found: {missing}"
    ]
    assert check_coverage(malformed, floors) == [
        f"coverage report is not valid JSON: {malformed}"
    ]
    assert check_coverage(_write_json(tmp_path / "report.json", _report({})), bad_schema) == [
        "coverage floors field 'files' must be an object"
    ]


def test_check_coverage_rejects_invalid_floor_and_report_percentage(tmp_path):
    report = _write_json(
        tmp_path / "coverage.json",
        _report({"abstar/a.py": "eighty"}),
    )
    floors = _write_json(
        tmp_path / "floors.json",
        _floors({"abstar/a.py": True}),
    )

    assert check_coverage(report, floors) == [
        "coverage floor for abstar/a.py must be an integer from 0 through 100"
    ]

    floors = _write_json(
        tmp_path / "floors.json",
        _floors({"abstar/a.py": 80}),
    )
    assert check_coverage(report, floors) == [
        "coverage report percentage for abstar/a.py must be a number from 0 through 100"
    ]


@pytest.mark.parametrize(
    ("source", "value", "expected"),
    [
        (
            "report",
            True,
            "coverage report percentage for abstar/a.py must be a number from 0 through 100",
        ),
        (
            "report",
            -0.01,
            "coverage report percentage for abstar/a.py must be a number from 0 through 100",
        ),
        (
            "report",
            100.01,
            "coverage report percentage for abstar/a.py must be a number from 0 through 100",
        ),
        (
            "floors",
            True,
            "coverage floor for abstar/a.py must be an integer from 0 through 100",
        ),
        (
            "floors",
            -1,
            "coverage floor for abstar/a.py must be an integer from 0 through 100",
        ),
        (
            "floors",
            101,
            "coverage floor for abstar/a.py must be an integer from 0 through 100",
        ),
    ],
)
def test_check_coverage_rejects_boolean_and_out_of_range_values(
    tmp_path, source, value, expected
):
    report_value = value if source == "report" else 80
    floor_value = value if source == "floors" else 80
    report = _write_json(
        tmp_path / "coverage.json", _report({"abstar/a.py": report_value})
    )
    floors = _write_json(
        tmp_path / "floors.json", _floors({"abstar/a.py": floor_value})
    )

    assert check_coverage(report, floors) == [expected]


@pytest.mark.parametrize("source", ["report", "floors"])
@pytest.mark.parametrize("constant", ["NaN", "Infinity", "-Infinity"])
def test_check_coverage_rejects_nonfinite_json_numbers(tmp_path, source, constant):
    report = tmp_path / "coverage.json"
    report.write_text(
        '{"meta":{"branch_coverage":true},"files":{"abstar/a.py":'
        '{"summary":{"percent_covered":80}}}}',
        encoding="utf-8",
    )
    floors = tmp_path / "floors.json"
    floors.write_text('{"files":{"abstar/a.py":80}}', encoding="utf-8")
    target = report if source == "report" else floors
    target.write_text(
        target.read_text(encoding="utf-8").replace("80", constant),
        encoding="utf-8",
    )

    assert check_coverage(report, floors) == [
        f"coverage {source} contains invalid JSON number: {constant}"
    ]


@pytest.mark.parametrize(
    ("source", "content", "duplicate"),
    [
        (
            "report",
            '{"meta":{"branch_coverage":true},"files":{},"files":{}}',
            "files",
        ),
        (
            "report",
            '{"meta":{"branch_coverage":true},"files":{"abstar/a.py":'
            '{"summary":{"percent_covered":80,"percent_covered":81}}}}',
            "percent_covered",
        ),
        ("floors", '{"files":{"abstar/a.py":80,"abstar/a.py":81}}', "abstar/a.py"),
    ],
)
def test_check_coverage_rejects_literal_duplicate_json_keys(
    tmp_path, source, content, duplicate
):
    report = _write_json(tmp_path / "coverage.json", _report({"abstar/a.py": 80}))
    floors = _write_json(tmp_path / "floors.json", _floors({"abstar/a.py": 80}))
    target = report if source == "report" else floors
    target.write_text(content, encoding="utf-8")

    expected = f"coverage {source} contains duplicate JSON object key: {duplicate}"
    assert check_coverage(report, floors) == [expected]

    failed = subprocess.run(
        [sys.executable, str(CHECKER), str(report), str(floors)],
        check=False,
        capture_output=True,
        text=True,
    )
    assert failed.returncode == 1
    assert failed.stdout == expected + "\n"
    assert failed.stderr == ""


def test_check_coverage_cli_prints_sorted_failures_and_sets_status(tmp_path):
    report = _write_json(
        tmp_path / "coverage.json",
        _report({"abstar/z.py": 79.125}),
    )
    floors = _write_json(
        tmp_path / "floors.json",
        _floors({"abstar/z.py": 80, "abstar/a.py": 75}),
    )

    failed = subprocess.run(
        [sys.executable, str(CHECKER), str(report), str(floors)],
        check=False,
        capture_output=True,
        text=True,
    )
    assert failed.returncode == 1
    assert failed.stdout.splitlines() == [
        "abstar/a.py: missing from coverage report (floor 75%)",
        "abstar/z.py: 79.125% is below 80% floor",
    ]
    assert failed.stderr == ""

    passing_report = _write_json(
        tmp_path / "passing-coverage.json",
        _report({"abstar/z.py": 80, "abstar/a.py": 75}),
    )
    passed = subprocess.run(
        [sys.executable, str(CHECKER), str(passing_report), str(floors)],
        check=False,
        capture_output=True,
        text=True,
    )
    assert passed.returncode == 0
    assert passed.stdout == ""
    assert passed.stderr == ""
