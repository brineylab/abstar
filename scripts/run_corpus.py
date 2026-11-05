#!/usr/bin/env python3
"""Run the fixed committed BCR cohort and compare all native output fields exactly.

The output directory must be new and outside the checkout and corpus trees.
Baseline recording proposes files in that output directory only. It never
accepts unexpected failures or inconsistencies, or updates the committed corpus.
"""

import argparse
from collections import Counter
import csv
import gzip
import hashlib
import importlib.metadata
import importlib.util
import json
import math
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys
import time
import traceback

import polars as pl


ROOT = Path(__file__).resolve().parents[1]
PARAMETERS = {"receptor": "bcr", "germline_database": "human",
              "n_processes": 2, "mmseqs_threads": 2, "chunksize": 500}
# Specify all execution-sensitive options, including assignment batch size.
RUN_OPTIONS = {**PARAMETERS, "mmseqs_chunksize": 1_000_000,
               "output_format": "parquet", "debug": False, "strict": False,
               "merge": False, "interleaved_fastq": False,
               "umi_pattern": None, "umi_length": None,
               "copy_inputs_to_project": False, "verbose": False,
               "concise_logging": True, "as_dataframe": False,
               "started_from_cli": False}
POLICY = "Exact ordered schema and every public field; no float tolerances or gene normalization."
GUIDE = "abstar/tests/README.md"
RECORD_COLUMNS = {"corpus_ordinal", "source_file", "record_ordinal", "sequence_id",
                  "sequence_sha256", "donor", "flow_class", "locus", "panel",
                  "selection_reasons"}
FAILURE_FIELDS = ("stage", "category", "exception_type")


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_write(path, data):
    Path(path).write_text(json.dumps(data, indent=2, ensure_ascii=False, default=str) + "\n")


def load_corpus(corpus, record_baseline=False):
    """Authenticate bytes, metadata types, physical order, and original identities."""
    corpus = Path(corpus).resolve()
    manifest = json.loads((corpus / "manifest.json").read_text())
    if manifest.get("schema_version") != 1 or manifest.get("corpus_version") != "bcr-ci-v1":
        raise ValueError("Unsupported corpus schema/version")
    if manifest.get("parameters") != PARAMETERS:
        raise ValueError(f"Corpus parameters must be exactly {PARAMETERS}")
    required = {"sequences.fasta.gz", "records.parquet"}
    if not record_baseline:
        required.update({"baseline.parquet", "baseline-metadata.json"})
    hashes = manifest.get("files", {})
    if not isinstance(hashes, dict) or not required.issubset(hashes):
        raise ValueError(f"Manifest files must authenticate {sorted(required)}")
    for name, expected_hash in hashes.items():
        path = corpus / name
        if Path(name).name != name or not path.resolve().is_relative_to(corpus):
            raise ValueError(f"Unsafe manifest filename: {name}")
        if not isinstance(expected_hash, str) or not re.fullmatch(r"[a-f0-9]{64}", expected_hash):
            raise ValueError(f"Invalid SHA256 for {name}")
        if sha256(path) != expected_hash:
            raise ValueError(f"SHA256 mismatch: {name}")
    records = pl.read_parquet(corpus / "records.parquet")
    if set(records.columns) != RECORD_COLUMNS:
        raise ValueError("Unexpected records.parquet schema columns")
    for column, dtype in records.schema.items():
        if column in {"corpus_ordinal", "record_ordinal"}:
            valid = dtype.is_integer()
        elif column == "selection_reasons":
            valid = dtype == pl.List(pl.String)
        else:
            valid = dtype == pl.String
        if not valid or records[column].null_count():
            raise ValueError(f"Invalid records metadata schema/nulls: {column}")
    count = manifest.get("record_count")
    if type(count) is not int or count <= 0 or records.height != count:
        raise ValueError("Corpus record count must be positive and match metadata")
    if records["corpus_ordinal"].to_list() != list(range(count)):
        raise ValueError("Metadata must be in contiguous corpus ordinal order")
    if records["record_ordinal"].min() < 0:
        raise ValueError("Source record ordinals must be nonnegative")
    if records.select("source_file", "record_ordinal").n_unique() != count:
        raise ValueError("Duplicate original source/record ordinal identity")
    sequences = []
    # Parse the committed two-line format strictly. Blank, truncated, malformed,
    # or description-bearing headers must not be silently skipped by a parser.
    with gzip.open(corpus / "sequences.fasta.gz", "rt") as handle:
        for row in records.iter_rows(named=True):
            header, sequence = handle.readline().rstrip("\r\n"), handle.readline().rstrip("\r\n")
            if not header.startswith(">") or not header[1:] or any(c.isspace() for c in header[1:]):
                raise ValueError(f"Invalid FASTA header at corpus ordinal {row['corpus_ordinal']}")
            if not re.fullmatch(r"[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+", sequence):
                raise ValueError(f"Invalid FASTA sequence at corpus ordinal {row['corpus_ordinal']}")
            if header[1:] != row["sequence_id"] or hashlib.sha256(sequence.encode()).hexdigest() != row["sequence_sha256"]:
                raise ValueError(f"FASTA identity/sequence mismatch at corpus ordinal {row['corpus_ordinal']}")
            sequences.append((header[1:], sequence))
        if handle.read():
            raise ValueError("Extra FASTA records/content beyond metadata count")
    failures = manifest.get("expected_failures")
    if not isinstance(failures, list):
        raise ValueError("Manifest must list expected_failures individually")
    seen = set()
    for failure in failures:
        if not isinstance(failure, dict) or set(failure) != {"corpus_ordinal", *FAILURE_FIELDS}:
            raise ValueError("Malformed expected failure")
        ordinal = failure["corpus_ordinal"]
        if type(ordinal) is not int or not 0 <= ordinal < count or ordinal in seen:
            raise ValueError("Invalid or duplicate expected failure ordinal")
        if any(not isinstance(failure[k], str) or not failure[k] for k in FAILURE_FIELDS):
            raise ValueError("Expected failures require explicit stage/category/exception_type")
        seen.add(ordinal)
    return manifest, records, sequences


def _identity(records, ordinal):
    row = records.row(ordinal, named=True)
    return {k: row[k] for k in ("corpus_ordinal", "source_file", "record_ordinal", "sequence_id", "panel")}


def _failure_ordinal(failure):
    match = re.fullmatch(r"abstar_0_(\d+)", failure.get("row_id", ""))
    if match is None:
        raise ValueError(f"Invalid failure row identity: {failure}")
    return int(match[1])


def read_failures(output):
    with (Path(output) / "logs/failures.tsv").open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"input_file", "sample", "row_id", "sequence_id", *FAILURE_FIELDS}.issubset(reader.fieldnames or []):
            raise ValueError("Malformed failure index")
        return list(reader)


def map_output(frame, records, sequences, failures, input_path):
    """Recover ordinals from ordered rows and failure row keys, never external IDs."""
    failed = set()
    for failure in failures:
        ordinal = _failure_ordinal(failure)
        if ordinal in failed or not 0 <= ordinal < records.height:
            raise ValueError(f"Duplicate or out-of-range failure ordinal: {ordinal}")
        if (Path(failure["input_file"]).resolve() != Path(input_path).resolve()
                or failure["sample"] != Path(input_path).stem
                or failure["sequence_id"] != sequences[ordinal][0]):
            raise ValueError(f"Failure source/sequence identity mismatch: {ordinal}")
        failed.add(ordinal)
    ordinals = [i for i in range(records.height) if i not in failed]
    if frame.height != len(ordinals):
        raise ValueError(f"Record conservation failed: input={records.height}, output={frame.height}, failures={len(failed)}")
    if not {"sequence_id", "sequence"}.issubset(frame.columns) or frame.schema["sequence_id"] != pl.String:
        raise ValueError("Output must preserve string sequence_id and original sequence")
    if "corpus_ordinal" in frame.columns:
        raise ValueError("Unexpected internal corpus ordinal in public output")
    for ordinal, row in zip(ordinals, frame.select("sequence_id", "sequence").iter_rows()):
        if row != sequences[ordinal]:
            raise ValueError(f"Output order/sequence mismatch at corpus ordinal {ordinal}: {_identity(records, ordinal)}")
    return frame.with_columns(pl.Series("corpus_ordinal", ordinals, dtype=pl.Int64))


def compare_failures(failures, expected, records):
    actual = {_failure_ordinal(f): {k: f[k] for k in FAILURE_FIELDS} for f in failures}
    baseline = {f["corpus_ordinal"]: {k: f[k] for k in FAILURE_FIELDS} for f in expected}
    return [{"kind": "record_outcome", **_identity(records, ordinal),
             "expected": baseline.get(ordinal, "output row"),
             "actual": actual.get(ordinal, "output row")}
            for ordinal in sorted(actual.keys() | baseline.keys())
            if actual.get(ordinal) != baseline.get(ordinal)]


def _same(left, right):
    if isinstance(left, float) and isinstance(right, float) and math.isnan(left) and math.isnan(right):
        return True
    return left == right


def compare_frames(expected, actual, records):
    """Exact all-field comparison; identify changed records by immutable ordinals."""
    differences = []
    if list(expected.schema.items()) != list(actual.schema.items()):
        differences.append({"kind": "schema", "expected": {k: str(v) for k, v in expected.schema.items()},
                            "actual": {k: str(v) for k, v in actual.schema.items()}})
    for label, frame in (("expected", expected), ("actual", actual)):
        if "corpus_ordinal" not in frame.columns or not frame.schema["corpus_ordinal"].is_integer():
            raise ValueError(f"{label} baseline/output has no integer corpus ordinal")
        ordinals = frame["corpus_ordinal"].to_list()
        if any(type(i) is not int or not 0 <= i < records.height for i in ordinals):
            raise ValueError(f"{label} baseline/output ordinal out of range")
        if ordinals != sorted(set(ordinals)):
            raise ValueError(f"{label} baseline/output ordinals duplicated or reordered")
    old_ordinals, new_ordinals = set(expected["corpus_ordinal"]), set(actual["corpus_ordinal"])
    for ordinal in sorted(old_ordinals ^ new_ordinals):
        differences.append({"kind": "record_outcome", **_identity(records, ordinal),
                            "expected": "output row" if ordinal in old_ordinals else "failure",
                            "actual": "output row" if ordinal in new_ordinals else "failure"})
    common = sorted(old_ordinals & new_ordinals)
    old = expected.filter(pl.col("corpus_ordinal").is_in(common))
    new = actual.filter(pl.col("corpus_ordinal").is_in(common))
    # Column equality is cheap for unchanged columns, avoiding materializing
    # millions of Python values in the usual successful large-corpus check.
    for field in expected.columns:
        if field not in actual.columns or old[field].equals(new[field], check_dtypes=True):
            continue
        for ordinal, left, right in zip(common, old[field], new[field]):
            if not _same(left, right):
                differences.append({"kind": "field", **_identity(records, ordinal),
                                    "field": field, "expected": left, "actual": right})
    return differences


def collect_environment():
    import abutils
    import abstar
    from abstar.annotation.germline import get_germline_database_path

    imported_package = Path(abstar.__file__).resolve()
    if imported_package != (ROOT / "abstar/__init__.py").resolve():
        raise ValueError(f"The imported abstar package is outside this checkout: {imported_package}; install this checkout with python -m pip install -e .")
    resolved = Path(get_germline_database_path("human", "bcr")).resolve()
    packaged = Path(abstar.__file__).resolve().parent / "germline_dbs/bcr/human"
    if resolved != packaged.resolve():
        raise ValueError(f"User database shadows packaged BCR human reference: {resolved}; use an environment without this override")
    binary = Path(abutils.bin.get_path("mmseqs")).resolve()
    version = subprocess.run([str(binary), "version"], check=True, capture_output=True, text=True)
    packages = {}
    for package in ("abstar", "abutils", "biopython", "polars", "pyarrow", "parasail",
                    "numpy", "pandas", "click", "pyfastx", "pytest"):
        try:
            packages[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            packages[package] = "unavailable"
    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True, check=True)
    status = subprocess.run(["git", "status", "--short"], cwd=ROOT, capture_output=True, text=True, check=True)
    return {"python": sys.version, "executable": sys.executable, "platform": platform.platform(),
            "packages": packages, "git_sha": commit.stdout.strip(), "git_status": status.stdout,
            "runner_sha256": sha256(__file__),
            "thread_environment": {name: os.environ.get(name) for name in
                                   ("POLARS_MAX_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")},
            "cpu_count": os.cpu_count(),
            "cpu_affinity": sorted(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else None,
            "source_sha256": {str(p.relative_to(ROOT)): sha256(p) for p in sorted((ROOT / "abstar").rglob("*.py"))},
            "mmseqs": {"path": str(binary), "version": version.stdout.strip(), "sha256": sha256(binary)},
            "reference_path": str(resolved),
            "reference_sha256": {str(p.relative_to(resolved)): sha256(p) for p in sorted(resolved.rglob("*")) if p.is_file()}}


def _audit(frame):
    spec = importlib.util.spec_from_file_location("corpus_consistency_auditor", ROOT / "scripts/audit_annotation_consistency.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.audit_frame(frame)


def run_corpus(corpus, output, record_baseline=False):
    """Run once, retain native output/diagnostics, and always attempt a report."""
    corpus, output = Path(corpus).resolve(), Path(output).resolve()
    if any(output.is_relative_to(p) for p in (ROOT, corpus)):
        raise ValueError("Output must be outside checkout and corpus/input trees")
    output.mkdir(parents=True, exist_ok=False)
    started = time.perf_counter()
    report = {"status": "failed", "mode": "record_baseline" if record_baseline else "check",
              "corpus": str(corpus), "output": str(output), "parameters": RUN_OPTIONS,
              "comparison_policy": POLICY, "debugging_guide": GUIDE,
              "differences": [], "timing_seconds": {}}
    try:
        manifest, records, sequences = load_corpus(corpus, record_baseline)
        report.update(corpus_version=manifest["corpus_version"], manifest_sha256=sha256(corpus / "manifest.json"),
                      input_files_sha256=manifest["files"], counts={"input": records.height},
                      panels=dict(Counter(records["panel"].to_list())))
        report["environment"] = environment = collect_environment()
        baseline_metadata = None
        if not record_baseline:
            baseline_metadata = json.loads((corpus / "baseline-metadata.json").read_text())
            for field, expected in (("parameters", RUN_OPTIONS), ("comparison_policy", POLICY),
                                    ("corpus_files_sha256", {k: manifest["files"][k] for k in ("records.parquet", "sequences.fasta.gz")})):
                if baseline_metadata.get(field) != expected:
                    raise ValueError(f"Baseline metadata drift: {field}")
            for field in ("reference_sha256",):
                if baseline_metadata.get("environment", {}).get(field) != environment[field]:
                    raise ValueError(f"Baseline reference drift: {field}")
            report["environment_differences"] = [
                {"field": field, "expected": baseline_metadata["environment"].get(field), "actual": environment[field]}
                for field in ("python", "platform", "packages", "mmseqs", "thread_environment")
                if baseline_metadata["environment"].get(field) != environment[field]]
        input_path = output / "sequences.fasta"
        with gzip.open(corpus / "sequences.fasta.gz", "rb") as source, input_path.open("xb") as destination:
            shutil.copyfileobj(source, destination)
        report["timing_seconds"]["validation"] = time.perf_counter() - started
        import abstar
        annotation_started = time.perf_counter()
        abstar.run(str(input_path), project_path=str(output), **RUN_OPTIONS)
        report["timing_seconds"]["annotation"] = time.perf_counter() - annotation_started
        paths = sorted((output / "parquet").glob("*.parquet"))
        if paths != [output / "parquet/sequences.parquet"]:
            raise ValueError(f"Unexpected native output files: {paths}")
        frame = pl.read_parquet(paths[0])
        failures = read_failures(output)
        actual = map_output(frame, records, sequences, failures, input_path)
        report["counts"].update(output=frame.height, failures=len(failures),
                                statuses=dict(Counter(frame["annotation_status"].to_list())))
        report["differences"].extend(compare_failures(failures, manifest["expected_failures"], records))
        audit_started = time.perf_counter()
        findings = _audit(frame)
        report["counts"]["inconsistent"] = findings.height
        for finding in findings.iter_rows(named=True):
            ordinal = actual["corpus_ordinal"][finding["parquet_row"]]
            report["differences"].append({"kind": "consistency", **_identity(records, ordinal), **finding})
        report["timing_seconds"]["audit"] = time.perf_counter() - audit_started
        if not record_baseline:
            comparison_started = time.perf_counter()
            expected = pl.read_parquet(corpus / "baseline.parquet")
            report["differences"].extend(compare_frames(expected, actual, records))
            report["timing_seconds"]["comparison"] = time.perf_counter() - comparison_started
        report["changed_panels"] = dict(Counter(d["panel"] for d in report["differences"] if "panel" in d))
        if not report["differences"]:
            if record_baseline:
                actual.write_parquet(output / "baseline.parquet", compression="zstd")
                _json_write(output / "baseline-metadata.json", {
                    "schema_version": 1, "parameters": RUN_OPTIONS, "comparison_policy": POLICY,
                    "corpus_files_sha256": {k: manifest["files"][k] for k in ("records.parquet", "sequences.fasta.gz")},
                    "environment": environment, "counts": report["counts"],
                    "debugging_guide": GUIDE,
                })
                report["status"] = "baseline_recorded"
            else:
                report["status"] = "passed"
    except Exception as error:
        # This boundary turns any failure into a failing report/exit status. It
        # never converts a record, tool, programming, or storage error to success.
        report["error"] = {"type": type(error).__name__, "message": str(error),
                           "traceback": traceback.format_exc()}
    finally:
        report["timing_seconds"]["total"] = time.perf_counter() - started
        _json_write(output / "report.json", report)
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--corpus", type=Path, default=ROOT / "test_data/bcr_corpus")
    parser.add_argument("--output", type=Path, required=True, help="New run directory outside checkout and input trees")
    parser.add_argument("--record-baseline", action="store_true", help="Propose a baseline in the output directory; never modify corpus inputs")
    args = parser.parse_args(argv)
    report = run_corpus(args.corpus, args.output, args.record_baseline)
    print(json.dumps({k: report[k] for k in ("status", "output", "timing_seconds")}, indent=2))
    if "error" in report:
        print(f"{report['error']['type']}: {report['error']['message']}", file=sys.stderr)
    return 0 if report["status"] in {"passed", "baseline_recorded"} else 1


if __name__ == "__main__":
    raise SystemExit(main())
