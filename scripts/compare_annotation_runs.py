#!/usr/bin/env python
"""Compare corpus reruns by immutable original FASTA record ordinal.

Inputs are explicit run roots containing parquet/<FASTA-stem>.parquet and
logs/failures.tsv. FASTAs are directly inside --fasta-dir. Expected JSON has a
`records` list of {source_file, record_ordinal, expected: {field: value}};
ordinals are zero-based. Optional sample_ordinal validates failure row keys.
Every baseline failure must have exactly one expected record. Other fixture
metadata is ignored. Comparisons include column order, dtypes, nulls, NaNs and
all public values, with no floating-point tolerance. Identical public rows
cannot reveal a permutation of physically indistinguishable input records.

Memory is bounded by one sample's Parquet frames and 8192 streamed FASTA
records. Reports must be outside the repository and all input/run trees.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import re

from Bio.SeqIO.FastaIO import SimpleFastaParser
import polars as pl
import pyarrow.parquet as pq


def _failures(root, required=False):
    path = root / "logs/failures.tsv"
    if not path.exists() and not required:
        return []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"input_file", "sample", "row_id", "sequence_id"}.issubset(reader.fieldnames or []):
            raise ValueError(f"Invalid failure index: {path}")
        return list(reader)


def _same(a, b):
    if isinstance(a, float) and isinstance(b, float) and math.isnan(a) and math.isnan(b):
        return True
    return type(a) is type(b) and a == b


def compare_runs(baseline, candidate, fasta_dir, expected):
    baseline, candidate, fasta_dir, expected = map(Path, (baseline, candidate, fasta_dir, expected))
    report = {"ok": False, "totals": dict.fromkeys(
        ("input_records", "baseline_successes", "candidate_successes", "recovered"), 0),
        "samples": [], "discrepancy_count": 0, "discrepancies": []}

    def error(message):
        report["discrepancy_count"] += 1
        if len(report["discrepancies"]) < 100:
            report["discrepancies"].append(message)

    files = sorted(p for p in fasta_dir.iterdir() if p.suffix.lower() in (".fasta", ".fa", ".fna"))
    if not files or len({p.stem for p in files}) != len(files):
        raise ValueError("FASTA inputs must be nonempty with unique stems")
    names = {p.name for p in files}
    failures = {}
    for failure in _failures(baseline, required=True):
        source = Path(failure["input_file"])
        match = re.fullmatch(r"abstar_(\d+)_(\d+)", failure["row_id"])
        if source.name not in names or failure["sample"] != source.stem or match is None:
            raise ValueError(f"Invalid baseline failure identity: {failure}")
        key = (source.name, int(match[2]))
        if key in failures:
            raise ValueError(f"Duplicate baseline failure: {key}")
        failures[key] = (failure["sequence_id"], int(match[1]))
    if _failures(candidate):
        error("Candidate has record failures")
    fixtures = {}
    for record in json.loads(expected.read_text())["records"]:
        key = (record["source_file"], record["record_ordinal"])
        if type(key[1]) is not int or key[1] < 0 or key in fixtures or not record["expected"]:
            raise ValueError(f"Invalid or duplicate expected record: {key}")
        fixtures[key] = record
    if failures.keys() != fixtures.keys():
        error("Expected recovery keys do not exactly match baseline failure keys")
    expected_outputs = {f"{p.stem}.parquet" for p in files}
    for label, root in (("baseline", baseline), ("candidate", candidate)):
        actual = {p.name for p in (root / "parquet").glob("*.parquet")}
        if actual != expected_outputs:
            error(f"{label}: missing or extra sample outputs: {sorted(actual ^ expected_outputs)}")
    for source in files:
        old_path = baseline / "parquet" / f"{source.stem}.parquet"
        new_path = candidate / "parquet" / f"{source.stem}.parquet"
        if not old_path.exists() or not new_path.exists():
            continue
        old, new = pl.read_parquet(old_path), pl.read_parquet(new_path)
        failed = {ordinal: value for (name, ordinal), value in failures.items() if name == source.name}
        count = 0
        old_offset = 0
        chunk = []

        def check_chunk(records, offset):
            nonlocal old_offset
            ids, sequences, old_ids, old_sequences = [], [], [], []
            for relative, (identifier, sequence) in enumerate(records):
                ordinal = offset + relative
                ids.append(identifier)
                sequences.append(sequence)
                if ordinal in failed:
                    if failed[ordinal][0] != identifier:
                        error(f"{source.name}:{ordinal}: baseline failure ID differs from FASTA")
                else:
                    old_ids.append(identifier)
                    old_sequences.append(sequence)
            for label, frame, start, expected_ids, expected_sequences in (
                ("candidate", new, offset, ids, sequences),
                ("baseline", old, old_offset, old_ids, old_sequences),
            ):
                for field, values in (("sequence_id", expected_ids), ("sequence", expected_sequences)):
                    if field not in frame.columns or frame[field].slice(start, len(values)).to_list() != values:
                        error(f"{source.name}: {label} {field} differs from FASTA at chunk {offset}")
            old_offset += len(old_ids)

        with source.open() as handle:
            for title, sequence in SimpleFastaParser(handle):
                if not title.split():
                    raise ValueError(f"Empty FASTA ID in {source}")
                chunk.append((title.split()[0], sequence.upper()))
                count += 1
                if len(chunk) == 8192:
                    check_chunk(chunk, count - len(chunk))
                    chunk = []
            if chunk:
                check_chunk(chunk, count - len(chunk))
        if any(ordinal >= count for ordinal in failed):
            error(f"{source.name}: baseline failure ordinal outside FASTA")
        if new.height != count or old.height != count - len(failed):
            error(f"{source.name}: record conservation failed: input={count}, baseline={old.height}, candidate={new.height}, old_failures={len(failed)}")
        schema_equal = (
            list(old.schema.items()) == list(new.schema.items())
            and pq.read_schema(old_path).equals(pq.read_schema(new_path), check_metadata=False)
        )
        if not schema_equal:
            error(f"{source.name}: public schema, column order or dtype changed")
        else:
            # Filter by ordinal, never by user ID or sequence content.
            kept = new.filter(~pl.int_range(0, pl.len()).is_in(sorted(failed)))
            if not old.equals(kept, null_equal=True):
                changed = [name for name in old.columns if not old[name].equals(kept[name], null_equal=True)]
                error(f"{source.name}: old successful rows changed or reordered; fields={changed}")
        recovered = 0
        for ordinal, (_, sample_ordinal) in sorted(failed.items()):
            fixture = fixtures.get((source.name, ordinal))
            if fixture is None or ordinal >= new.height:
                continue
            if "sample_ordinal" in fixture and fixture["sample_ordinal"] != sample_ordinal:
                error(f"{source.name}:{ordinal}: fixture sample ordinal mismatch")
            row = new.row(ordinal, named=True)
            for field, value in fixture["expected"].items():
                if field not in row or not _same(row[field], value):
                    error(f"{source.name}:{ordinal}: recovered {field} differs: expected={value!r}, actual={row.get(field)!r}")
            recovered += 1
        stats = {"input_records": count, "baseline_successes": old.height,
                 "candidate_successes": new.height, "recovered": recovered}
        report["samples"].append({"source_file": source.name, **stats})
        for key, value in stats.items():
            report["totals"][key] += value
    report["ok"] = report["discrepancy_count"] == 0
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("baseline", "candidate", "fasta-dir", "expected", "report"):
        parser.add_argument(f"--{name}", required=True, type=Path)
    args = parser.parse_args(argv)
    destination = args.report.resolve()
    protected = [Path(__file__).resolve().parents[1], args.baseline.resolve(),
                 args.candidate.resolve(), args.fasta_dir.resolve()]
    if destination == args.expected.resolve() or any(destination.is_relative_to(root) for root in protected):
        raise ValueError("Report must be outside repository, input and run trees")
    try:
        result = compare_runs(args.baseline, args.candidate, args.fasta_dir, args.expected)
    except (ValueError, OSError, KeyError, TypeError, pl.exceptions.PolarsError) as exc:
        result = {"ok": False, "error": f"{type(exc).__name__}: {exc}"}
    destination.parent.mkdir(parents=True, exist_ok=True)
    # Never overwrite a report or other preexisting artifact.
    with destination.open("x") as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write("\n")
    print(json.dumps({"ok": result["ok"], "report": str(destination), "totals": result.get("totals")}))
    return 0 if result["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
