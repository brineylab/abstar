# Copyright (c) 2025 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

from pathlib import Path

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
