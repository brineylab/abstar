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


@pytest.fixture
def discovery_module():
    import importlib.util
    from pathlib import Path

    path = Path(__file__).resolve().parents[2] / "scripts" / "discover_bcr_cases.py"
    assert path.is_file(), "optional discovery command is missing"
    spec = importlib.util.spec_from_file_location("discover_bcr_cases", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def discovery_corpus(tmp_path, annotations):
    manifest, paths, fields, rows = annotations
    fasta_dir = tmp_path / "fastas"
    fasta_dir.mkdir()
    root = tmp_path / "raw"
    for dataset in paths:
        target = root / "test_data" / dataset / "outs" / "per_sample_outs" / dataset / "vdj_b"
        target.mkdir(parents=True)
        write_csv(target / "filtered_contig_annotations.csv", fields, rows)
        (fasta_dir / f"{dataset}.fasta").write_text(
            ">10E8\nACGT\n>00123\nTGCA\n", encoding="utf-8"
        )
    return fasta_dir, manifest, root


@pytest.mark.integration
def test_discovery_cli_is_path_free_and_discovery_order_invariant(
    discovery_module, discovery_corpus, tmp_path, monkeypatch
):
    """Real file/CSV/report boundary; annotation is deliberately injected."""
    import json
    module = discovery_module
    fasta_dir, manifest, root = discovery_corpus
    report = tmp_path / "candidates.jsonl"

    def annotate(records, **kwargs):
        return {r.row_key: {"status": "annotated", "exception_category": None,
                "failure_scope": None, "failure_group": None,
                "annotation": {"v_call": "IGHV1-2*02,IGHV1-2*01", "j_call": "IGHJ4*02",
                               "productive": True, "junction_aa": "CAR", "junction": "TGTGCTCGT",
                               "cdr3_aa": "A", "cdr3": "GCT", "v_insertions": "3:AAA"}}
                for r in records}

    monkeypatch.setattr(module, "annotate_records", annotate)
    args = ["--fasta-dir", str(fasta_dir), "--manifest", str(manifest),
            "--cellranger-root", str(root), "--output", str(report),
            "--per-dataset", "2", "--seed", "fixture-seed", "--n-processes", "1"]
    assert module.main(args) == 0
    first = report.read_bytes()
    scan = module.discover_annotation_files
    monkeypatch.setattr(module, "discover_annotation_files", lambda root: list(reversed(scan(root))))
    assert module.main(args) == 0
    assert report.read_bytes() == first
    assert b"/home/" not in first and str(tmp_path).encode() not in first
    header, *candidates = [json.loads(line) for line in first.splitlines()]
    assert header["schema_version"] == 2
    assert header["selection"] == {"algorithm_version": 1, "seed": "fixture-seed", "per_dataset": 2}
    assert header["receptor"] == "bcr" and header["database"] == "human"
    assert header["abstar"]["version"] and header["abstar"]["git_revision"]
    assert all(header["dependencies"][key] for key in ("abutils", "polars", "pyarrow", "parasail"))
    assert header["external_tools"]["mmseqs"]["version"]
    assert len(header["external_tools"]["mmseqs"]["sha256"]) == 64
    assert len(header["germline_manifests"]) >= 1
    assert all(len(item["sha256"]) == 64 for item in header["germline_manifests"])
    assert len(candidates) == 4
    assert {(r["source"]["dataset"], r["source"]["sequence_id"]) for r in candidates} == {
        (d, s) for d in ("00123", "10E8") for s in ("00123", "10E8")}
    for row in candidates:
        assert row["status"] == "annotated"
        assert row["failure_scope"] is None and row["failure_group"] is None
        assert row["source"]["reads"] == "00123"
        assert row["source"]["umis"] == "10E8"
        assert row["normalized_calls"]["abstar"]["v"] == ["IGHV1-2"]
        assert row["comparisons"]["junction_aa"] is True
        assert row["comparisons"]["cdr3_aa"] is False
        assert row["comparisons"]["productive"] is True
        assert row["no_d"] == {"abstar": True, "cellranger": True}
        assert row["ties"]["abstar"]["v"] == ["IGHV1-2*01", "IGHV1-2*02"]
        assert row["indels"]["v_insertions"] == "3:AAA"
        assert row["selection_reasons"]
        assert row["lengths"]["input_nt"] == 4
        assert row["sequence_sha256"] == hashlib.sha256(row["sequence"].encode()).hexdigest()


@pytest.mark.parametrize("mode, expected", [
    ("missing", ["annotated", "exception"]),
    ("raise", ["exception", "exception"]),
    ("duplicate", ["exception", "annotated"]),
    ("unknown", ["exception", "exception"]),
    ("unassigned", ["unassigned", "unassigned"]),
])
def test_discovery_annotation_accounts_for_every_selected_identity(
    discovery_module, tmp_path, monkeypatch, mode, expected
):
    """Unit test of missing/exception/cardinality accounting at abstar.run boundary."""
    import polars as pl
    module = discovery_module
    records = [CorpusRecord("d", identifier, "ACGT", "donor", "IgG", "IGH")
               for identifier in ("10E8", "00123")]
    projects = []

    def run(sequences, project_path, **kwargs):
        from pathlib import Path
        projects.append(Path(project_path))
        assert kwargs["receptor"] == "bcr" and kwargs["germline_database"] == "human"
        assert kwargs["n_processes"] == 1 and kwargs["output_format"] == "parquet"
        ids = [s.id for s in sequences]
        assert len(set(ids)) == 2 and all(i not in ("10E8", "00123") for i in ids)
        if mode == "raise":
            raise RuntimeError("/home/private/input must never be serialized")
        out = Path(project_path) / "parquet"
        out.mkdir()
        result_ids = [ids[0]] if mode == "missing" else [ids[0], ids[0], ids[1]] if mode == "duplicate" else ids if mode == "unassigned" else ["unknown"]
        pl.DataFrame({"sequence_id": result_ids, "v_call": [None if mode == "unassigned" else "IGHV1-2*01"] * len(result_ids),
                      "j_call": [None if mode == "unassigned" else "IGHJ4*02"] * len(result_ids)}).write_parquet(out / "sequences.parquet")

    monkeypatch.setattr(module.abstar, "run", run)
    outcomes = module.annotate_records(records, n_processes=1, germline_database="human")
    assert list(outcomes) == [r.row_key for r in records]
    assert [v["status"] for v in outcomes.values()] == expected
    categories = [v["exception_category"] for v in outcomes.values()]
    assert categories == {"missing": [None, "missing_output"], "raise": ["RuntimeError"] * 2,
                          "duplicate": ["output_cardinality", None],
                          "unknown": ["output_cardinality"] * 2, "unassigned": [None, None]}[mode]
    assert [v["failure_scope"] for v in outcomes.values()] == {
        "missing": [None, "record"], "raise": ["batch", "batch"],
        "duplicate": ["record", None], "unknown": ["batch", "batch"],
        "unassigned": [None, None]}[mode]
    groups = [v["failure_group"] for v in outcomes.values()]
    if mode in ("raise", "unknown"):
        assert groups[0] == groups[1] and groups[0].startswith("batch-")
    else:
        assert groups == [None, None]
    from .corpus import CandidateEvidence
    for record in records:
        outcome = outcomes[record.row_key]
        source = CandidateEvidence(record.dataset, record.sequence_id, record.donor, record.flow_class,
                                   record.chain, None, None, None, None, None, None, None, None, None)
        candidate = module._candidate(record, source, outcome, "fixture-seed")
        assert candidate["failure_scope"] == outcome["failure_scope"]
        assert candidate["failure_group"] == outcome["failure_group"]
        if outcome["failure_scope"]:
            assert f"annotation_{outcome['failure_scope']}_failure" in candidate["selection_reasons"]
    assert all(not p.exists() for p in projects)


@pytest.mark.parametrize("problem", ["missing_fasta", "ambiguous_csv", "duplicate_fasta_id", "output_in_source", "unmatched_fasta_id", "invalid_processes"])
@pytest.mark.integration
def test_discovery_rejects_invalid_inputs_before_annotation(
    discovery_module, discovery_corpus, tmp_path, monkeypatch, problem
):
    module = discovery_module
    fasta_dir, manifest, root = discovery_corpus
    report = tmp_path / "report.jsonl"
    n_processes = "1"
    if problem == "missing_fasta":
        (fasta_dir / "00123.fasta").unlink()
    elif problem == "ambiguous_csv":
        target = root / "00123"
        target.mkdir()
        (target / "filtered_contig_annotations.csv").write_text("contig_id\n")
    elif problem == "duplicate_fasta_id":
        with (fasta_dir / "00123.fasta").open("a") as handle:
            handle.write(">10E8\nAAAA\n")
    elif problem == "unmatched_fasta_id":
        with (fasta_dir / "00123.fasta").open("a") as handle:
            handle.write(">missing\nAAAA\n")
    elif problem == "output_in_source":
        report = fasta_dir / "report.jsonl"
    else:
        n_processes = "0"
    def unexpected(*args, **kwargs):
        pytest.fail("invalid inputs reached annotation")
    monkeypatch.setattr(module, "annotate_records", unexpected)
    with pytest.raises((ValueError, FileNotFoundError)):
        module.main(["--fasta-dir", str(fasta_dir), "--manifest", str(manifest),
                     "--cellranger-root", str(root), "--output", str(report),
                     "--n-processes", n_processes])
    assert not report.exists()


@pytest.mark.e2e
def test_discovery_real_annotation_preserves_source_identity(discovery_module, single_hc_sequence):
    """The real abstar.run boundary must return an accounted, nonempty annotation."""
    record = CorpusRecord("00123", "10E8", single_hc_sequence.sequence, "001", "IgG", "IGH")
    outcomes = discovery_module.annotate_records([record], n_processes=1, germline_database="human")
    assert list(outcomes) == ["00123\0" + "10E8"]
    outcome = outcomes[record.row_key]
    assert outcome["status"] == "annotated"
    assert outcome["exception_category"] is None
    assert normalize_gene(outcome["annotation"]["v_call"]) == ("IGHV3-15",)
    assert normalize_gene(outcome["annotation"]["j_call"]) == ("IGHJ1",)
    assert outcome["annotation"]["junction_aa"] == "CARTGKYYDFWSGYPPGEEYFQDW"


def test_discovery_hash_supports_python_310(discovery_module, tmp_path, monkeypatch):
    """Python 3.10 has no hashlib.file_digest; file hashing must stay portable."""
    path = tmp_path / "source"
    path.write_bytes(b"ACGT\n")
    monkeypatch.delattr(hashlib, "file_digest", raising=False)
    assert discovery_module.sha256(path) == hashlib.sha256(b"ACGT\n").hexdigest()


@pytest.mark.integration
@pytest.mark.parametrize("mode", ["valid", "missing", "failed", "empty"])
def test_discovery_records_checked_mmseqs_provenance(discovery_module, tmp_path, monkeypatch, mode):
    """Resolve/hash an executable and check its actual version subprocess."""
    import subprocess
    executable = tmp_path / "mmseqs"
    contents = {"valid": b'#!/bin/sh\nprintf "fixture-mmseqs-1\\n"\n',
                "failed": b'#!/bin/sh\nexit 7\n', "empty": b'#!/bin/sh\nexit 0\n'}
    if mode != "missing":
        executable.write_bytes(contents[mode])
        executable.chmod(0o700)
    monkeypatch.setattr(discovery_module, "get_binary_path", lambda name: str(executable))
    if mode == "valid":
        assert discovery_module.mmseqs_metadata() == {
            "version": "fixture-mmseqs-1", "sha256": hashlib.sha256(contents[mode]).hexdigest()}
    else:
        exception = {"missing": FileNotFoundError, "failed": subprocess.CalledProcessError,
                     "empty": ValueError}[mode]
        with pytest.raises(exception):
            discovery_module.mmseqs_metadata()


@pytest.mark.integration
@pytest.mark.parametrize("database", ["/home/private/database", "../human", "human/../other", "human/child", "human\\child", ".", ".."])
def test_discovery_rejects_database_paths(discovery_module, discovery_corpus, tmp_path, monkeypatch, database):
    """A database argument must never become a serialized filesystem path."""
    fasta_dir, manifest, root = discovery_corpus
    output = tmp_path / "report.jsonl"
    def unexpected(*args, **kwargs):
        pytest.fail("database path reached discovery/annotation")
    monkeypatch.setattr(discovery_module, "_select", unexpected)
    with pytest.raises(ValueError, match="logical database name"):
        discovery_module.main(["--fasta-dir", str(fasta_dir), "--manifest", str(manifest),
                               "--cellranger-root", str(root), "--output", str(output),
                               "--germline-database", database])
    assert not output.exists()


def test_discovery_mixed_batch_failure_does_not_blame_valid_record(discovery_module, monkeypatch):
    """A valid record is affected by a batch exception, not labeled its cause."""
    import json
    import polars as pl
    from pathlib import Path
    from .corpus import CandidateEvidence
    module = discovery_module
    records = [CorpusRecord("d", "valid", "ACGT", "donor", "IgG", "IGH"),
               CorpusRecord("d", "bad", "NNNN", "donor", "IgG", "IGH")]
    def run(sequences, project_path, **kwargs):
        good = [s for s in sequences if s.sequence == "ACGT"]
        out = Path(project_path) / "parquet"
        out.mkdir()
        pl.DataFrame({"sequence_id": [s.id for s in good],
                      "v_call": ["IGHV1-2*01"] * len(good),
                      "j_call": ["IGHJ4*02"] * len(good)}).write_parquet(out / "sequences.parquet")
        if any(s.sequence == "NNNN" for s in sequences):
            raise RuntimeError("failure on second input /home/private/corpus")
    monkeypatch.setattr(module.abstar, "run", run)
    alone = module.annotate_records(records[:1], n_processes=1, germline_database="human")
    assert alone[records[0].row_key]["status"] == "annotated"
    combined = module.annotate_records(records, n_processes=1, germline_database="human")
    reversed_outcomes = module.annotate_records(reversed(records), n_processes=1, germline_database="human")
    for record in records:
        outcome = combined[record.row_key]
        assert outcome["failure_scope"] == "batch"
        assert outcome["failure_group"] == reversed_outcomes[record.row_key]["failure_group"]
        assert outcome["failure_group"] == combined[records[0].row_key]["failure_group"]
        assert outcome["exception_category"] == "RuntimeError"
        source = CandidateEvidence(record.dataset, record.sequence_id, record.donor, record.flow_class,
                                   record.chain, "IGHV1-2", None, "IGHJ4", None, "True",
                                   "CAR", "TGTGCTCGT", "1", "1")
        candidate = module._candidate(record, source, outcome, "fixture-seed")
        assert candidate["failure_scope"] == "batch"
        assert candidate["failure_group"] == outcome["failure_group"]
        assert "annotation_batch_failure" in candidate["selection_reasons"]
        assert "annotation_exception" not in candidate["selection_reasons"]
        assert "/home/" not in json.dumps(candidate, sort_keys=True)
