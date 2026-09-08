"""Small corpus harness contracts; the committed large cohort runs only in CI."""

import csv
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path

import polars as pl
import pytest


ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture
def runner():
    path = ROOT / "scripts/run_corpus.py"
    assert path.exists(), "The committed corpus runner has not been implemented"
    spec = importlib.util.spec_from_file_location("ci_corpus_runner", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def make_corpus(path, sequences=(("10E8", "ACGT"), ("10E8", "TGCA"))):
    path.mkdir()
    with gzip.open(path / "sequences.fasta.gz", "wt") as handle:
        for name, sequence in sequences:
            handle.write(f">{name}\n{sequence}\n")
    rows = [dict(corpus_ordinal=i, source_file="original.fasta", record_ordinal=i,
                 sequence_id=name, sequence_sha256=hashlib.sha256(seq.encode()).hexdigest(),
                 donor="donor", flow_class="IgM", locus="IGH", panel="representative",
                 selection_reasons=["representative"])
            for i, (name, seq) in enumerate(sequences)]
    pl.DataFrame(rows).write_parquet(path / "records.parquet")
    manifest = dict(schema_version=1, corpus_version="bcr-ci-v1", record_count=len(rows),
                    files={}, sources=[], selection={}, expected_failures=[],
                    parameters=dict(receptor="bcr", germline_database="human",
                                    n_processes=2, mmseqs_threads=2, chunksize=500))
    (path / "manifest.json").write_text(json.dumps(manifest))
    refresh_hashes(path)
    return path


def refresh_hashes(path):
    manifest = json.loads((path / "manifest.json").read_text())
    manifest["files"] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                         for p in path.iterdir() if p.name != "manifest.json"}
    (path / "manifest.json").write_text(json.dumps(manifest))


def test_load_preserves_duplicate_ids_and_fixed_ordinals(runner, tmp_path):
    corpus = make_corpus(tmp_path / "corpus")
    _, records, sequences = runner.load_corpus(corpus, record_baseline=True)
    assert records["sequence_id"].to_list() == ["10E8", "10E8"]
    assert sequences == [("10E8", "ACGT"), ("10E8", "TGCA")]


@pytest.mark.parametrize("damage", ["hash", "metadata_order", "duplicate_identity", "numeric_id", "sequence", "count", "unsafe_file"])
def test_load_rejects_corrupt_inputs(runner, tmp_path, damage):
    corpus = make_corpus(tmp_path / "corpus")
    records = pl.read_parquet(corpus / "records.parquet")
    if damage == "hash":
        (corpus / "sequences.fasta.gz").write_bytes(b"broken")
    elif damage in {"count", "unsafe_file"}:
        manifest = json.loads((corpus / "manifest.json").read_text())
        if damage == "count":
            manifest["record_count"] += 1
        else:
            manifest["files"]["../secret"] = "0" * 64
        (corpus / "manifest.json").write_text(json.dumps(manifest))
    else:
        if damage == "metadata_order":
            records = records.reverse()
        elif damage == "duplicate_identity":
            records = records.with_columns(pl.lit(0).alias("record_ordinal"))
        elif damage == "numeric_id":
            records = records.with_columns(pl.lit(10).alias("sequence_id"))
        elif damage == "sequence":
            records = records.with_columns(pl.lit("0" * 64).alias("sequence_sha256"))
        records.write_parquet(corpus / "records.parquet")
        refresh_hashes(corpus)
    with pytest.raises(ValueError):
        runner.load_corpus(corpus, record_baseline=True)


def test_output_mapping_never_joins_external_ids(runner, tmp_path):
    corpus = make_corpus(tmp_path / "corpus")
    _, records, sequences = runner.load_corpus(corpus, record_baseline=True)
    actual = pl.DataFrame({"sequence_id": ["10E8", "10E8"], "sequence": ["ACGT", "TGCA"]})
    mapped = runner.map_output(actual, records, sequences, [], tmp_path / "sequences.fasta")
    assert mapped["corpus_ordinal"].to_list() == [0, 1]
    for malformed in (actual.head(1), pl.concat([actual, actual.head(1)]), actual.reverse(),
                      pl.concat([actual.head(1), actual.head(1)])):
        with pytest.raises(ValueError):
            runner.map_output(malformed, records, sequences, [], tmp_path / "sequences.fasta")


def test_failure_mapping_and_individual_expectations(runner, tmp_path):
    corpus = make_corpus(tmp_path / "corpus")
    _, records, sequences = runner.load_corpus(corpus, record_baseline=True)
    source = tmp_path / "sequences.fasta"
    failure = dict(input_file=str(source), sample="sequences", row_id="abstar_0_1",
                   sequence_id="10E8", stage="annotation", category="record_error", exception_type="ValueError")
    actual = pl.DataFrame({"sequence_id": ["10E8"], "sequence": ["ACGT"]})
    mapped = runner.map_output(actual, records, sequences, [failure], source)
    assert mapped["corpus_ordinal"].to_list() == [0]
    expected = [dict(corpus_ordinal=1, stage="annotation", category="record_error", exception_type="ValueError")]
    assert runner.compare_failures([failure], expected, records) == []
    assert runner.compare_failures([], expected, records)[0]["kind"] == "record_outcome"
    assert runner.compare_failures([failure], [], records)[0]["kind"] == "record_outcome"
    assert runner.compare_failures([dict(failure, category="other")], expected, records)
    for damaged in ([failure, failure], [dict(failure, row_id="abstar_1_1")],
                    [dict(failure, row_id="abstar_0_2")], [dict(failure, sequence_id="wrong")]):
        with pytest.raises(ValueError):
            runner.map_output(actual, records, sequences, damaged, source)


def test_exact_comparison_reports_identity_and_every_changed_field(runner, tmp_path):
    corpus = make_corpus(tmp_path / "corpus")
    _, records, _ = runner.load_corpus(corpus, record_baseline=True)
    expected = pl.DataFrame(dict(corpus_ordinal=[0, 1], sequence_id=["10E8", "10E8"],
                                 v_call=["IGHV1-2*02"] * 2, v_support=[1.0, 2.0]))
    actual = expected.with_columns(pl.when(pl.col("corpus_ordinal") == 1)
        .then(pl.lit("IGHV1-2*04")).otherwise(pl.col("v_call")).alias("v_call"),
        (pl.col("v_support") + 1e-12).alias("v_support"))
    differences = runner.compare_frames(expected, actual, records)
    assert len(differences) == 3
    assert {d["field"] for d in differences} == {"v_call", "v_support"}
    assert all(d["source_file"] == "original.fasta" for d in differences)
    assert runner.compare_frames(expected, expected, records) == []
    assert runner.compare_frames(expected, expected.with_columns(pl.col("corpus_ordinal").cast(pl.Int32)), records)[0]["kind"] == "schema"


def test_validation_exception_retains_report_and_rejects_existing_output(runner, tmp_path):
    output = tmp_path / "run"
    assert runner.main(["--corpus", str(tmp_path / "missing"), "--output", str(output)]) == 1
    report = json.loads((output / "report.json").read_text())
    assert report["status"] == "failed"
    assert report["error"]["type"] == "FileNotFoundError"
    assert "abstar/tests/README.md" in report["debugging_guide"]
    with pytest.raises(FileExistsError):
        runner.run_corpus(tmp_path / "missing", output)
    with pytest.raises(ValueError, match="outside"):
        runner.run_corpus(tmp_path / "missing", ROOT / "forbidden-corpus-run")


@pytest.mark.parametrize("failure", [None, {}, {"corpus_ordinal": 5, "stage": "annotation", "category": "error", "exception_type": "ValueError"}])
def test_expected_failures_are_explicit_validated_records(runner, tmp_path, failure):
    corpus = make_corpus(tmp_path / "corpus")
    manifest = json.loads((corpus / "manifest.json").read_text())
    manifest["expected_failures"] = failure if failure is None else [failure]
    (corpus / "manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(ValueError):
        runner.load_corpus(corpus, record_baseline=True)


def fake_public_run(monkeypatch, *, failed_ordinal=None, inconsistent=False, crash=False):
    """Substitute only the costly external boundary to inject pipeline outcomes."""
    import abstar
    from abstar.annotation.schema import OUTPUT_SCHEMA
    from abstar.core.diagnostics import FIELDS

    def run(sequences, project_path, **kwargs):
        output = Path(project_path)
        (output / "parquet").mkdir()
        (output / "logs").mkdir()
        if crash:
            raise RuntimeError("MMseqs failed with saved stderr")
        lines = Path(sequences).read_text().splitlines()
        rows, failures = [], []
        for ordinal, (header, sequence) in enumerate(zip(lines[::2], lines[1::2])):
            if ordinal == failed_ordinal:
                failures.append(dict(input_file=sequences, sample="sequences", row_id=f"abstar_0_{ordinal}",
                                     sequence_id=header[1:], stage="annotation", category="record_error",
                                     exception_type="ValueError"))
            else:
                row = dict(sequence_id=header[1:], sequence=sequence,
                           annotation_status="annotated" if inconsistent else "unassigned")
                rows.append(row)
        pl.DataFrame(rows, schema=OUTPUT_SCHEMA).write_parquet(output / "parquet/sequences.parquet")
        with (output / "logs/failures.tsv").open("w") as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
            writer.writeheader()
            writer.writerows(failures)

    monkeypatch.setattr(abstar, "run", run)


@pytest.mark.parametrize("scenario", ["unexpected_failure", "expected_failure", "wrong_category", "inconsistent", "tool_crash"])
def test_bootstrap_never_accepts_unreviewed_failures_or_inconsistency(runner, tmp_path, monkeypatch, scenario):
    corpus = make_corpus(tmp_path / "corpus")
    if scenario in {"expected_failure", "wrong_category"}:
        manifest = json.loads((corpus / "manifest.json").read_text())
        manifest["expected_failures"] = [dict(corpus_ordinal=1, stage="annotation",
            category="other" if scenario == "wrong_category" else "record_error", exception_type="ValueError")]
        (corpus / "manifest.json").write_text(json.dumps(manifest))
    monkeypatch.setattr(runner, "collect_environment", lambda: {"reference_sha256": {}})
    fake_public_run(monkeypatch, failed_ordinal=1 if "failure" in scenario or scenario == "wrong_category" else None,
                    inconsistent=scenario == "inconsistent", crash=scenario == "tool_crash")
    output = tmp_path / "proposed"
    report = runner.run_corpus(corpus, output, record_baseline=True)
    if scenario == "expected_failure":
        assert report["status"] == "baseline_recorded", report
        assert pl.read_parquet(output / "baseline.parquet")["corpus_ordinal"].to_list() == [0]
        assert report["counts"]["failures"] == 1
    else:
        assert report["status"] == "failed"
        assert not (output / "baseline.parquet").exists()
        assert not (output / "baseline-metadata.json").exists()
        if scenario == "inconsistent":
            assert report["counts"]["inconsistent"] == 2
            assert report["differences"][0]["kind"] == "consistency"
        elif scenario == "tool_crash":
            assert report["error"]["type"] == "RuntimeError"
        else:
            assert report["differences"][0]["kind"] == "record_outcome"
    assert json.loads((output / "report.json").read_text())["status"] == report["status"]


@pytest.mark.parametrize("damage", ["parameters", "reference", "corpus_hash"])
def test_baseline_metadata_drift_blocks_annotation(runner, tmp_path, monkeypatch, damage):
    corpus = make_corpus(tmp_path / "corpus")
    monkeypatch.setattr(runner, "collect_environment", lambda: {"reference_sha256": {}})
    fake_public_run(monkeypatch)
    proposed = tmp_path / "proposed"
    assert runner.run_corpus(corpus, proposed, record_baseline=True)["status"] == "baseline_recorded"
    metadata = json.loads((proposed / "baseline-metadata.json").read_text())
    if damage == "parameters":
        metadata["parameters"]["chunksize"] = 1
    elif damage == "reference":
        metadata["environment"]["reference_sha256"] = {"v.fasta": "changed"}
    else:
        metadata["corpus_files_sha256"]["records.parquet"] = "changed"
    (corpus / "baseline-metadata.json").write_text(json.dumps(metadata))
    (corpus / "baseline.parquet").write_bytes((proposed / "baseline.parquet").read_bytes())
    refresh_hashes(corpus)
    report = runner.run_corpus(corpus, tmp_path / "check")
    assert report["status"] == "failed"
    assert "drift" in report["error"]["message"]
    assert not (tmp_path / "check/parquet").exists()


def test_environment_rejects_import_from_another_checkout(runner, tmp_path, monkeypatch):
    import abstar
    from abstar.annotation import germline
    monkeypatch.setattr(abstar, "__file__", str(tmp_path / "installed/abstar/__init__.py"))
    monkeypatch.setattr(germline, "get_germline_database_path", lambda *args:
                        str(tmp_path / "installed/abstar/germline_dbs/bcr/human"))
    with pytest.raises(ValueError, match="imported abstar.*checkout"):
        runner.collect_environment()


@pytest.mark.integration
def test_real_public_api_baseline_check_and_changed_gene(runner, tmp_path, monkeypatch):
    from abstar.tests.corpus import load_real_bcr_cases
    from abstar.annotation import germline
    original_expanduser = germline.os.path.expanduser
    monkeypatch.setattr(germline.os.path, "expanduser", lambda value:
                        str(tmp_path / "user-db") if str(value).startswith("~/.abstar/") else original_expanduser(value))
    case = next(c for c in load_real_bcr_cases() if c.sequence_id == "GCTGCGAGTCCTGCTT-1_contig_1")
    corpus = make_corpus(tmp_path / "corpus", [("0001", case.sequence), ("0001", case.sequence)])
    proposed = tmp_path / "proposed"
    report = runner.run_corpus(corpus, proposed, record_baseline=True)
    assert report["status"] == "baseline_recorded", report
    assert "POLARS_MAX_THREADS" in report["environment"]["thread_environment"]
    assert len(report["environment"]["runner_sha256"]) == 64
    assert {"pandas", "click", "pyfastx"}.issubset(report["environment"]["packages"])
    baseline = pl.read_parquet(proposed / "baseline.parquet")
    assert baseline["sequence_id"].to_list() == ["0001", "0001"]
    assert baseline["junction_aa"].to_list() == ["CARYHPVLRNGFDVW"] * 2
    assert baseline["productive"].to_list() == [True, True]
    for field in ("v_sequence_start", "v_sequence_end", "j_sequence_start", "j_sequence_end", "rev_comp", "cdr3"):
        assert baseline[field].to_list() == [case.expected[field]] * 2
    for name in ("baseline.parquet", "baseline-metadata.json"):
        (corpus / name).write_bytes((proposed / name).read_bytes())
    refresh_hashes(corpus)
    before = {p.name: p.read_bytes() for p in corpus.iterdir()}
    assert runner.run_corpus(corpus, tmp_path / "check")["status"] == "passed"
    assert before == {p.name: p.read_bytes() for p in corpus.iterdir()}
    baseline.with_columns(pl.lit("IGHV9-99*99").alias("v_call")).write_parquet(corpus / "baseline.parquet")
    refresh_hashes(corpus)
    report = runner.run_corpus(corpus, tmp_path / "changed")
    assert report["status"] == "failed"
    assert [d["field"] for d in report["differences"]] == ["v_call", "v_call"]
