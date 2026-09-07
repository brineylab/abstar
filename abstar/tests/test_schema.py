# Copyright (c) 2025 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

from pathlib import Path
import os
import subprocess
import sys

import polars as pl
import pytest

from ..annotation.schema import ANNOTATION_WORK_SCHEMA, OUTPUT_SCHEMA, NoneDict


def test_noneddict_missing_key_returns_none():
    nd = NoneDict({"a": 1})
    assert nd["a"] == 1
    assert nd["missing"] is None


def test_output_schema_is_polars_schema_like():
    # Create a DataFrame with a subset of the schema to ensure it can be used as polars schema
    subset = {k: OUTPUT_SCHEMA[k] for k in ["sequence_id", "sequence", "productive"]}
    df = pl.DataFrame(
        [{"sequence_id": "x", "sequence": "ACGT", "productive": True}],
        schema=NoneDict(subset),
    )
    assert df.shape == (1, 3)


def test_output_schema_exposes_annotation_outcome_without_internal_row_id():
    assert OUTPUT_SCHEMA["annotation_status"] == pl.String
    assert OUTPUT_SCHEMA["failure_reason"] == pl.String
    assert "row_id" not in OUTPUT_SCHEMA


def test_annotation_work_schema_adds_internal_row_id_to_public_output_schema():
    assert list(ANNOTATION_WORK_SCHEMA) == ["row_id", *OUTPUT_SCHEMA]
    assert ANNOTATION_WORK_SCHEMA["row_id"] == pl.String
    assert {
        field: dtype
        for field, dtype in ANNOTATION_WORK_SCHEMA.items()
        if field != "row_id"
    } == OUTPUT_SCHEMA


def test_pytest_configuration_registers_test_scopes(pytestconfig):
    assert Path(pytestconfig.getini("testpaths")[0]) == Path("abstar/tests")
    markers = "\n".join(pytestconfig.getini("markers"))
    for name in ("integration", "e2e", "slow"):
        assert f"{name}:" in markers


def test_declared_pytest_floor_supports_strict_configuration():
    from packaging.requirements import Requirement

    requirements = Path(__file__).resolve().parents[2] / "requirements-test.txt"
    requirement = next(Requirement(line) for line in requirements.read_text().splitlines()
                       if line.startswith("pytest>="))
    assert "9.0.0" in requirement.specifier
    assert "8.4.2" not in requirement.specifier


@pytest.mark.parametrize("invalid", ["config", "marker"])
def test_project_pytest_configuration_rejects_unknown_settings(tmp_path, invalid):
    config = Path(__file__).resolve().parents[2] / "pyproject.toml"
    project = config.read_text()
    if invalid == "config":
        project += '\nunknown_abstar_option = true\n'
        code = 'def test_control():\n    assert True\n'
        diagnostic = "Unknown config option: unknown_abstar_option"
    else:
        code = 'import pytest\n@pytest.mark.unknown_abstar_marker\ndef test_control():\n    assert True\n'
        diagnostic = "unknown_abstar_marker"
    (tmp_path / "pyproject.toml").write_text(project)
    (tmp_path / "test_contract.py").write_text(code)
    result = subprocess.run(
        [sys.executable, "-m", "pytest", "test_contract.py", "-q"],
        cwd=tmp_path, env={**os.environ, "PYTEST_ADDOPTS": ""}, capture_output=True, text=True,
    )
    assert result.returncode != 0
    assert diagnostic in result.stdout + result.stderr
