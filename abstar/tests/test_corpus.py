# Copyright (c) 2025 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

"""Deterministic discovery contracts using only in-memory and temporary data."""

import csv
import hashlib
import random
from collections import Counter
from dataclasses import FrozenInstanceError, replace
from itertools import product

import pytest

from .corpus import (
    CorpusRecord,
    load_cellranger_annotations,
    normalize_gene,
    select_sweep,
    selection_metadata,
    stable_rank,
)


@pytest.fixture
def records():
    return [
        CorpusRecord(dataset, f"contig-{i}", "ACGT", donor, flow, chain)
        for dataset in ("00123", "10E8")
        for i, (donor, flow, chain, _) in enumerate(
            product(("donor-1", "donor-2"), ("IgG", "IgM"), ("IGH", "IGK", "IGL"), range(25))
        )
    ]


def test_stable_rank_has_versioned_identity_only_payload():
    record = CorpusRecord("dataset", "10E8", "ACGT", "donor", "IgG", "IGH")
    assert record.row_key == "dataset\0" + "10E8"
    assert stable_rank(record, seed="abstar-real-bcr-v1") == (
        "768b40148810d5a7ad981775c49b6f4275e41e199f923a61ed85794db48894da"
    )
    assert stable_rank(record, seed="x") == stable_rank(
        replace(record, sequence="TGCA", chain="IGK"), seed="x"
    )
    assert stable_rank(record, seed="x") != stable_rank(record, seed="y")
    assert stable_rank(record, seed="x") != stable_rank(record, seed="x", algorithm_version=2)
    assert stable_rank(record, seed="x") != stable_rank(replace(record, dataset="other"), seed="x")
    with pytest.raises(FrozenInstanceError):
        record.sequence_id = "changed"


def test_sweep_order_quota_and_all_strata(records):
    selected = select_sweep(records)
    shuffled = records.copy()
    random.Random(918).shuffle(shuffled)
    assert selected == select_sweep(reversed(records)) == select_sweep(shuffled)
    assert Counter(record.dataset for record in selected) == {"00123": 200, "10E8": 200}
    assert len({record.row_key for record in selected}) == 400
    assert { (r.dataset, r.donor, r.flow_class, r.chain) for r in selected } == {
        (r.dataset, r.donor, r.flow_class, r.chain) for r in records
    }

    # Independent payload hashing checks representatives and fill priority.
    rank = lambda r: hashlib.sha256(
        ("1\0abstar-real-bcr-v1\0" + r.dataset + "\0" + r.sequence_id).encode()
    ).hexdigest()
    for dataset in ("00123", "10E8"):
        candidates = [r for r in records if r.dataset == dataset]
        chosen = [r for r in selected if r.dataset == dataset]
        representatives = {
            min((r for r in candidates if (r.donor, r.flow_class, r.chain) == stratum), key=rank)
            for stratum in product(("donor-1", "donor-2"), ("IgG", "IgM"), ("IGH", "IGK", "IGL"))
        }
        expected = representatives | set(sorted(set(candidates) - representatives, key=rank)[:188])
        assert chosen == sorted(expected, key=rank)


def test_sweep_small_quota_uses_sorted_strata(records):
    selected = select_sweep(records, per_dataset=2)
    assert selected == select_sweep(reversed(records), per_dataset=2)
    assert Counter(r.dataset for r in selected) == {"00123": 2, "10E8": 2}
    assert {(r.donor, r.flow_class, r.chain) for r in selected} == {
        ("donor-1", "IgG", "IGH"), ("donor-1", "IgG", "IGK")
    }
    assert select_sweep(records, per_dataset=0) == []


def test_sweep_materializes_once_and_handles_unknown_chain():
    class Once:
        def __iter__(self):
            assert not getattr(self, "consumed", False), "input consumed twice"
            self.consumed = True
            yield CorpusRecord("d", "00123", "ACGT", "donor", "IgG", "IGH")
            yield CorpusRecord("d", "10E8", "ACGT", "donor", "IgG", None)

    assert [r.sequence_id for r in select_sweep(Once(), per_dataset=1)] == ["10E8"]
    assert len(select_sweep(Once())) == 2
    assert select_sweep(iter(())) == []


@pytest.mark.parametrize("quota", [-1, 1.5, True])
def test_sweep_rejects_invalid_quota(quota):
    with pytest.raises(ValueError, match="per_dataset"):
        select_sweep([], per_dataset=quota)


def test_sweep_rejects_duplicate_identity(records):
    with pytest.raises(ValueError, match="cardinality"):
        select_sweep([records[0], replace(records[0], sequence="TGCA")])


def test_selection_report_metadata_replays_configuration(records):
    metadata = selection_metadata(per_dataset=3, seed="custom")
    assert metadata == {"algorithm_version": 1, "seed": "custom", "per_dataset": 3}
    assert selection_metadata() == {
        "algorithm_version": 1, "seed": "abstar-real-bcr-v1", "per_dataset": 200
    }
    selected = select_sweep(records, per_dataset=metadata["per_dataset"], seed=metadata["seed"])
    assert len(selected) == 6
    assert selected != select_sweep(records, per_dataset=3)


@pytest.mark.parametrize("raw, expected", [
    (None, ()), ("", ()), ("None", ()),
    ("IGHV1-2*01", ("IGHV1-2",)),
    ("IGHV3-23*02, IGHV1-2*01,IGHV3-23*01", ("IGHV1-2", "IGHV3-23")),
    ("IGHV3-30-3*01", ("IGHV3-30-3",)),
    ("IGHA1", ("IGHA1",)),
])
def test_normalize_gene_preserves_distinct_genes_and_removes_only_alleles(raw, expected):
    assert normalize_gene(raw) == expected


def write_csv(path, fields, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


@pytest.fixture
def annotations(tmp_path):
    manifest = tmp_path / "sample_manifest.csv"
    write_csv(manifest, ["dataset", "donor", "flow_class"], [
        {"dataset": dataset, "donor": "001", "flow_class": "IgG"}
        for dataset in ("00123", "10E8")
    ])
    fields = ["contig_id", "chain", "v_gene", "d_gene", "j_gene", "c_gene",
              "productive", "cdr3", "cdr3_nt", "reads", "umis", "barcode"]
    rows = [
        dict(zip(fields, [identifier, "IGH", "IGHV1-2*01,IGHV1-2*02", "None",
                          "IGHJ4*02", "IGHG1", "True", "CAR", "TGTGCTCGT", "00123", "10E8", "00123"]))
        for identifier in ("10E8", "00123")
    ]
    paths = {dataset: tmp_path / f"{dataset}.csv" for dataset in ("00123", "10E8")}
    for path in paths.values():
        write_csv(path, fields, rows)
    return manifest, paths, fields, rows


def test_loader_preserves_source_strings_raw_evidence_and_cross_dataset_ids(annotations):
    manifest, paths, _, _ = annotations
    evidence = load_cellranger_annotations(manifest, paths)
    assert list(evidence) == ["00123\0" + "00123", "00123\0" + "10E8", "10E8\0" + "00123", "10E8\0" + "10E8"]
    item = evidence["00123\0" + "10E8"]
    assert (item.dataset, item.sequence_id, item.donor, item.flow_class, item.chain) == (
        "00123", "10E8", "001", "IgG", "IGH"
    )
    assert (item.v_gene, item.d_gene, item.j_gene, item.c_gene) == (
        "IGHV1-2*01,IGHV1-2*02", "None", "IGHJ4*02", "IGHG1"
    )
    assert item.normalized_genes == {"v": ("IGHV1-2",), "d": (), "j": ("IGHJ4",), "c": ("IGHG1",)}
    assert (item.productive, item.cdr3, item.cdr3_nt, item.reads, item.umis) == (
        "True", "CAR", "TGTGCTCGT", "00123", "10E8"
    )
    assert item.row_key == "00123\0" + "10E8"
    with pytest.raises(FrozenInstanceError):
        item.v_gene = "IGHV3-23"


def test_loader_rejects_duplicate_contig_within_dataset(annotations):
    manifest, paths, fields, rows = annotations
    write_csv(paths["00123"], fields, [rows[0], rows[0]])
    with pytest.raises(ValueError, match="cardinality.*00123.*10E8"):
        load_cellranger_annotations(manifest, paths)


def test_loader_rejects_duplicate_manifest_dataset(annotations):
    manifest, paths, _, _ = annotations
    with manifest.open("a") as handle:
        handle.write("00123,001,IgG\n")
    with pytest.raises(ValueError, match="cardinality.*00123"):
        load_cellranger_annotations(manifest, paths)


def test_loader_requires_matching_manifest_and_annotation_datasets(annotations):
    manifest, paths, _, _ = annotations
    with pytest.raises(ValueError, match="dataset"):
        load_cellranger_annotations(manifest, {"00123": paths["00123"]})
    with pytest.raises(ValueError, match="dataset"):
        load_cellranger_annotations(manifest, {**paths, "extra": paths["00123"]})


def test_loader_supports_explicit_manifest_column_mapping(annotations):
    manifest, paths, _, _ = annotations
    manifest.write_text(manifest.read_text().replace("dataset,donor,flow_class", "sample,subject,sort"))
    evidence = load_cellranger_annotations(
        manifest, paths, dataset_column="sample", donor_column="subject", flow_class_column="sort"
    )
    assert len(evidence) == 4
    assert evidence["00123\0" + "10E8"].donor == "001"


def test_loader_rejects_missing_evidence_columns(annotations):
    manifest, paths, _, _ = annotations
    paths["00123"].write_text("contig_id,chain\n10E8,IGH\n")
    with pytest.raises(ValueError, match="columns"):
        load_cellranger_annotations(manifest, paths)


@pytest.mark.parametrize("identifier", ["", "bad\0id"])
def test_loader_rejects_invalid_source_identity(annotations, identifier):
    manifest, paths, fields, rows = annotations
    write_csv(paths["00123"], fields, [{**rows[0], "contig_id": identifier}])
    with pytest.raises(ValueError, match="contig_id"):
        load_cellranger_annotations(manifest, paths)
