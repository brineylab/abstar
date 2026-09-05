#!/usr/bin/env python
# Copyright (c) 2025 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

"""Optionally discover BCR regression candidates from a read-only local corpus.

Run with --help for explicit input paths. FASTAs are <dataset>.fasta directly
inside --fasta-dir. Exactly one filtered_contig_annotations.csv must occur
below a directory component equal to each manifest dataset. Cell Ranger calls
are comparison evidence, never expected answers. Output is external JSONL:
one reproducibility header and one outcome per selected original contig.
No timestamp, runtime, source path or exception message enters the report.
"""

import argparse
import csv
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import tempfile
from collections.abc import Sequence
from dataclasses import asdict

import abstar
import polars as pl
from abutils import Sequence as AbSequence
from abutils.tools.search import get_binary_path
from Bio import SeqIO

from abstar.annotation.germline import get_germline_database_path
from abstar.tests.corpus import (
    CorpusRecord, load_cellranger_annotations, normalize_gene, select_sweep,
    selection_metadata, stable_rank,
)


SEGMENTS = ("v", "d", "j", "c")
BATCH_SIZE = 2000


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def discover_annotation_files(root):
    return list(root.rglob("filtered_contig_annotations.csv"))


def _paths(args):
    fasta_dir = Path(args.fasta_dir).resolve(strict=True)
    manifest = Path(args.manifest).resolve(strict=True)
    root = Path(args.cellranger_root).resolve(strict=True)
    output = Path(args.output).resolve()
    if not fasta_dir.is_dir() or not root.is_dir() or not manifest.is_file():
        raise ValueError("inputs must be FASTA/Cell Ranger directories and a manifest file")
    if not output.parent.is_dir() or output.is_dir():
        raise ValueError("output must have an existing parent directory")
    source_tree = Path(__file__).resolve().parents[1]
    if any(output.is_relative_to(p) for p in (fasta_dir, root, source_tree)) or output == manifest:
        raise ValueError("output must be outside input directories and the source tree")
    if args.n_processes <= 0:
        raise ValueError("n_processes must be positive")
    selection_metadata(per_dataset=args.per_dataset, seed=args.seed)
    with manifest.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not {"dataset", "donor", "flow_class"} <= set(reader.fieldnames or ()):
            raise ValueError("manifest requires dataset, donor, flow_class columns")
        samples = list(reader)
    datasets = [r["dataset"] for r in samples]
    if not datasets or len(set(datasets)) != len(datasets):
        raise ValueError("manifest must contain unique datasets")
    for row in samples:
        for field in ("dataset", "donor", "flow_class"):
            value = row[field]
            if not value or "\0" in value:
                raise ValueError(f"invalid manifest {field}")
        if Path(row["dataset"]).name != row["dataset"] or row["dataset"] in (".", ".."):
            raise ValueError("dataset must be a single directory/file component")
    files = discover_annotation_files(root)
    paths = {}
    for dataset in sorted(datasets):
        matches = [p.resolve(strict=True) for p in files if dataset in p.relative_to(root).parts[:-1]]
        if len(matches) != 1:
            raise ValueError(f"dataset {dataset!r} requires exactly one annotation CSV; found {len(matches)}")
        fasta = (fasta_dir / f"{dataset}.fasta").resolve(strict=True)
        if not fasta.is_file() or not fasta.is_relative_to(fasta_dir) or not matches[0].is_relative_to(root):
            raise ValueError("corpus files must remain inside their explicit input roots")
        if output in (fasta, matches[0]):
            raise ValueError("output cannot replace an input")
        paths[dataset] = (fasta, matches[0])
    if len({p[1] for p in paths.values()}) != len(paths):
        raise ValueError("annotation CSV cannot represent multiple datasets")
    return manifest, output, sorted(samples, key=lambda r: r["dataset"]), paths


def _select(samples, paths, *, per_dataset, seed):
    selected, evidence, sources = [], {}, []
    # Task 5's checked loader requires exact manifest coverage. Give it a
    # temporary one-dataset manifest to bound memory without weakening checks.
    with tempfile.TemporaryDirectory(prefix="abstar-discovery-") as directory:
        manifest = Path(directory) / "manifest.csv"
        for sample in samples:
            dataset = sample["dataset"]
            fasta, annotation = paths[dataset]
            with manifest.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=("dataset", "donor", "flow_class"))
                writer.writeheader()
                writer.writerow({k: sample[k] for k in writer.fieldnames})
            source = load_cellranger_annotations(manifest, {dataset: annotation})
            records = []
            with fasta.open(encoding="utf-8") as handle:
                for sequence in SeqIO.parse(handle, "fasta"):
                    key = dataset + "\0" + sequence.id
                    if key not in source:
                        raise ValueError(f"FASTA contig has no Cell Ranger evidence: {key!r}")
                    item = source[key]
                    if not sequence.seq:
                        raise ValueError(f"empty FASTA sequence: {key!r}")
                    records.append(CorpusRecord(dataset, sequence.id, str(sequence.seq),
                                                item.donor, item.flow_class, item.chain))
            if not records:
                raise ValueError(f"empty FASTA dataset: {dataset!r}")
            chosen = select_sweep(records, per_dataset=per_dataset, seed=seed)
            selected.extend(chosen)
            evidence.update({r.row_key: source[r.row_key] for r in chosen})
            sources.append({"dataset": dataset, "fasta_sha256": sha256(fasta),
                            "cellranger_sha256": sha256(annotation),
                            "input_records": len(records), "selected_records": len(chosen),
                            "cellranger_records": len(source)})
            print(f"Selected {len(chosen)} of {len(records)} from {dataset}", file=sys.stderr)
    return selected, evidence, sources


def annotate_records(records, *, n_processes, germline_database):
    """Account for each immutable source key across the actual public API.

    Replace external IDs only at the annotation boundary with unique ASCII row
    tokens. Restore source identity using the explicit token map, never inferred
    dataframe keys. A batch exception/cardinality error invalidates its entire
    batch; a silent omission is explicitly missing_output, not non-assignment.
    Exception messages may contain paths, so only class/category is retained.
    """
    records = list(records)
    if len({r.row_key for r in records}) != len(records):
        raise ValueError("selected row identity must be unique")
    outcomes = {}
    for start in range(0, len(records), BATCH_SIZE):
        batch = records[start:start + BATCH_SIZE]
        tokens = {f"discovery{start + i:09d}": r.row_key for i, r in enumerate(batch)}
        queries = [AbSequence(r.sequence, id=token) for token, r in zip(tokens, batch)]
        rows, category = {}, None
        with tempfile.TemporaryDirectory(prefix="abstar-discovery-annotation-") as project:
            try:
                abstar.run(queries, project_path=project, receptor="bcr",
                           germline_database=germline_database, output_format="parquet",
                           n_processes=n_processes, mmseqs_threads=n_processes)
                for path in sorted((Path(project) / "parquet").glob("*.parquet")):
                    for row in pl.read_parquet(path).iter_rows(named=True):
                        token = row.get("sequence_id")
                        if not isinstance(token, str) or token not in tokens or token in rows:
                            category = "output_cardinality"
                            break
                        rows[token] = row
                    if category:
                        break
            except Exception as error:
                # This optional command explicitly reports failures; it never
                # turns an exception into an ordinary successful empty result.
                category = type(error).__name__
        for token, key in tokens.items():
            row = rows.get(token)
            exception = category or ("missing_output" if row is None else None)
            outcomes[key] = {
                "status": "exception" if exception else "annotated" if row.get("v_call") and row.get("j_call") else "unassigned",
                "exception_category": exception,
                "annotation": None if exception else row,
            }
        print(f"Accounted for {min(start + BATCH_SIZE, len(records))}/{len(records)} selected records", file=sys.stderr)
    return outcomes


def mmseqs_metadata():
    """Identify the executable using the resolver used by mmseqs_search."""
    executable = Path(get_binary_path("mmseqs")).resolve(strict=True)
    version = subprocess.run([str(executable), "version"], check=True,
                             capture_output=True, text=True).stdout.strip()
    if not version:
        raise ValueError("MMseqs executable returned no version")
    return {"version": version, "sha256": sha256(executable)}


def _metadata(args, manifest, sources):
    database = Path(get_germline_database_path(args.germline_database, receptor="bcr"))
    manifests = sorted(database.rglob("manifest.txt"))
    if not manifests:
        raise ValueError("germline database must provide provenance manifests")
    repository = Path(abstar.__file__).resolve().parents[1]
    git = subprocess.run(["git", "-C", str(repository), "rev-parse", "HEAD"],
                         check=True, capture_output=True, text=True)
    # Capture the resolved installed dependency closure, including transitive
    # MMseqs wrapper/scientific dependencies, without serializing install paths.
    dependencies = {}
    pending = ["abstar"]
    from packaging.requirements import Requirement
    while pending:
        name = pending.pop()
        distribution = importlib.metadata.distribution(name)
        canonical = distribution.metadata["Name"].lower().replace("_", "-")
        if canonical in dependencies:
            continue
        dependencies[canonical] = distribution.version
        for raw in distribution.requires or ():
            requirement = Requirement(raw)
            if requirement.marker is None or requirement.marker.evaluate():
                pending.append(requirement.name)
    return {
        "type": "metadata", "schema_version": 1,
        "selection": selection_metadata(per_dataset=args.per_dataset, seed=args.seed),
        "abstar": {"version": importlib.metadata.version("abstar"), "git_revision": git.stdout.strip()},
        "python": platform.python_version(), "dependencies": dependencies,
        "external_tools": {"mmseqs": mmseqs_metadata()},
        "receptor": "bcr", "database": args.germline_database,
        "germline_manifests": [{"sha256": sha256(p)} for p in manifests],
        "manifest_sha256": sha256(manifest), "sources": sources,
        "n_processes": args.n_processes, "annotation_batch_size": BATCH_SIZE,
        "comparison_evidence": "Cell Ranger; not an oracle",
        "discovery_code_sha256": sha256(__file__),
    }


def _equal(left, right):
    return None if left is None or right is None else left == right


def _productive(raw):
    if raw is None or raw == "None" or raw == "":
        return None
    if raw.lower() not in ("true", "false"):
        raise ValueError("Cell Ranger productive must be True, False or missing")
    return raw.lower() == "true"


def _candidate(record, source, outcome, seed):
    annotation = outcome["annotation"] or {}
    raw = {"cellranger": {s: getattr(source, f"{s}_gene") for s in SEGMENTS},
           "abstar": {s: annotation.get(f"{s}_call") for s in SEGMENTS}}
    normalized = {name: {s: normalize_gene(call) for s, call in calls.items()} for name, calls in raw.items()}
    ties = {name: {s: sorted(set(c.strip() for c in (call or "").split(",") if c.strip() not in ("", "None")))
                   for s, call in calls.items()} for name, calls in raw.items()}
    comparisons = {s: _equal(normalized["abstar"][s], normalized["cellranger"][s])
                   if outcome["annotation"] is not None else None for s in SEGMENTS}
    comparisons.update({field: _equal(annotation.get(field), getattr(source, target))
                        for field, target in (("junction", "cdr3_nt"), ("junction_aa", "cdr3"),
                                              ("cdr3", "cdr3_nt"), ("cdr3_aa", "cdr3"))})
    comparisons["productive"] = _equal(annotation.get("productive"), _productive(source.productive))
    reasons = ["deterministic_stratified_sweep"]
    if outcome["exception_category"]:
        reasons.append("annotation_exception")
    reasons.extend(f"{key}_disagreement" for key, value in comparisons.items() if value is False)
    return {
        "type": "candidate", "row_key": record.row_key, "source": asdict(source),
        "sequence": record.sequence, "sequence_sha256": hashlib.sha256(record.sequence.encode()).hexdigest(),
        "status": outcome["status"], "exception_category": outcome["exception_category"],
        "raw_calls": raw, "normalized_calls": normalized, "comparisons": comparisons,
        "junction": {"abstar_nt": annotation.get("junction"), "abstar_aa": annotation.get("junction_aa"),
                     "cellranger_nt": source.cdr3_nt, "cellranger_aa": source.cdr3},
        "cdr3": {"abstar_nt": annotation.get("cdr3"), "abstar_aa": annotation.get("cdr3_aa")},
        "productivity": {"abstar": annotation.get("productive"), "cellranger": _productive(source.productive),
                         "abstar_issues": annotation.get("productivity_issues")},
        "ties": ties,
        "indels": {f"{s}_{kind}": annotation.get(f"{s}_{kind}") for s in ("v", "c") for kind in ("insertions", "deletions")},
        "no_d": {"abstar": not normalized["abstar"]["d"] if outcome["annotation"] is not None else None,
                 "cellranger": not normalized["cellranger"]["d"]},
        "lengths": {"input_nt": len(record.sequence), "abstar_cdr3_aa": annotation.get("cdr3_length"),
                    "cellranger_cdr3_aa": len(source.cdr3) if source.cdr3 is not None else None},
        "selection_rank": stable_rank(record, seed=seed), "selection_reasons": sorted(reasons),
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("fasta-dir", "manifest", "cellranger-root", "output"):
        parser.add_argument(f"--{name}", required=True)
    parser.add_argument("--per-dataset", type=int, default=200)
    parser.add_argument("--seed", default="abstar-real-bcr-v1")
    parser.add_argument("--n-processes", type=int, default=1)
    parser.add_argument("--germline-database", default="human")
    args = parser.parse_args(argv)
    manifest, output, samples, paths = _paths(args)
    records, evidence, sources = _select(samples, paths, per_dataset=args.per_dataset, seed=args.seed)
    header = _metadata(args, manifest, sources)
    outcomes = annotate_records(records, n_processes=args.n_processes, germline_database=args.germline_database)
    if set(outcomes) != {r.row_key for r in records}:
        raise ValueError("annotation outcomes must cover exactly the selected records")
    # Stage externally: malformed comparison data cannot leave a partial report.
    with tempfile.TemporaryDirectory(prefix="abstar-discovery-report-", dir=output.parent) as directory:
        staged = Path(directory) / "report.jsonl"
        with staged.open("w", encoding="utf-8", newline="\n") as handle:
            handle.write(json.dumps(header, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n")
            for record in records:
                row = _candidate(record, evidence[record.row_key], outcomes[record.row_key], args.seed)
                handle.write(json.dumps(row, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n")
        os.replace(staged, output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
