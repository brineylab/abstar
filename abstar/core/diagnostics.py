"""Persistent, record-addressable diagnostics for recoverable annotation errors."""

import csv
import hashlib
import json
import platform
import traceback
import uuid
from datetime import datetime, timezone
from importlib.metadata import version, PackageNotFoundError
from pathlib import Path
from urllib.parse import quote


FIELDS = (
    "input_file", "sample", "row_id", "sequence_id", "stage", "category",
    "exception_type", "message", "traceback_location", "diagnostic_path",
)


def safe_component(value):
    """Encode path syntax, bounding bytes while retaining a collision-resistant suffix."""
    encoded = quote(str(value), safe="-_")
    if encoded in ("", ".", ".."):
        encoded = "sequence-" + hashlib.sha256(str(value).encode()).hexdigest()[:16]
    if len(encoded) > 120:
        encoded = encoded[:90] + "-" + hashlib.sha256(str(value).encode()).hexdigest()[:24]
    return encoded


def failure_path(directory, sequence_id, row_id):
    return Path(directory) / f"{safe_component(sequence_id)}__{safe_component(row_id)}.failed"


def sample_directory(log_directory, sample):
    """Reserve a fresh namespace, including for reruns and duplicate input stems."""
    stem = safe_component(sample)
    index = 1
    while True:
        path = Path(log_directory) / (stem if index == 1 else f"{stem}__{index}")
        try:
            path.mkdir()
            return str(path)
        except FileExistsError:
            index += 1


def initialize_diagnostics(log_directory, parameters, inputs):
    log_directory = Path(log_directory)
    # Keep previous run indexes and metadata alongside their existing diagnostics.
    previous = uuid.uuid4().hex
    for name in ("failures.tsv", "run.json"):
        path = log_directory / name
        if path.exists():
            path.rename(path.with_name(f"{path.stem}.{previous}{path.suffix}"))
    with (log_directory / "failures.tsv").open("w", newline="") as handle:
        csv.DictWriter(handle, FIELDS, delimiter="\t").writeheader()
    versions = {"python": platform.python_version()}
    for package in ("abstar", "abutils", "biopython", "polars", "pyarrow", "parasail"):
        try:
            versions[package] = version(package)
        except PackageNotFoundError:
            versions[package] = "unavailable"
    metadata = {
        "started_at": datetime.now(timezone.utc).isoformat(),
        "versions": versions,
        "parameters": parameters,
        "inputs": inputs,
        # Editable checkouts can share a package version but differ in implementation.
        "source_sha256": {
            str(path.relative_to(Path(__file__).parents[1])): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted(Path(__file__).parents[1].rglob("*.py"))
        },
    }
    (log_directory / "run.json").write_text(json.dumps(metadata, indent=2) + "\n")


def write_record_diagnostic(directory, record, antibody, error):
    path = failure_path(directory, record["sequence_id"], record["row_id"])
    with path.open("x") as handle:
        handle.write("ASSIGNMENT RECORD\n" + json.dumps(record, ensure_ascii=False) + "\n\n")
        handle.write(f"ROW ID: {record['row_id']}\nSEQUENCE ID: {record['sequence_id']}\n\n")
        handle.write(antibody.format_log())
        handle.write("\n" + "".join(traceback.format_exception(error)))


def index_failures(log_directory, directory, input_file, sample, failures):
    paths = []
    with (Path(log_directory) / "failures.tsv").open("a", newline="") as handle:
        writer = csv.DictWriter(handle, FIELDS, delimiter="\t")
        for failure in failures:
            path = failure_path(directory, failure.sequence_id, failure.row_id)
            # Missing diagnostics are a storage/pipeline failure, never partial success.
            if not path.is_file():
                raise FileNotFoundError(f"Missing record diagnostic: {path}")
            lines = (failure.traceback_text or "").strip().splitlines()
            locations = [line.strip() for line in lines if line.lstrip().startswith('File "')]
            writer.writerow(dict(
                input_file=str(input_file), sample=sample, row_id=failure.row_id,
                sequence_id=failure.sequence_id, stage=failure.stage, category=failure.category,
                exception_type=lines[-1].split(":", 1)[0] if lines else "",
                message=failure.message, traceback_location=locations[-1] if locations else "",
                diagnostic_path=str(path.relative_to(log_directory)),
            ))
            paths.append(str(path))
    return paths
