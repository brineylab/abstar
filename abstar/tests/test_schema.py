# Copyright (c) 2025 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

import polars as pl

from ..annotation.schema import OUTPUT_SCHEMA, NoneDict


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


def test_output_schema_includes_indel_aware_identity_fields():
    assert OUTPUT_SCHEMA["v_identity_indel"] == pl.Float64
    assert OUTPUT_SCHEMA["v_identity_aa_indel"] == pl.Float64
    assert OUTPUT_SCHEMA["c_identity_indel"] == pl.Float64
    assert OUTPUT_SCHEMA["c_identity_aa_indel"] == pl.Float64
