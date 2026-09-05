# Copyright (c) 2025 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

"""Source evidence and deterministic selection for optional BCR discovery.

These helpers never annotate sequences or derive expected biological answers.
Cell Ranger values are comparison evidence, retained without type inference.
"""

import hashlib
from collections import defaultdict
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from os import PathLike

import polars as pl


ALGORITHM_VERSION = 1
DEFAULT_SEED = "abstar-real-bcr-v1"
DEFAULT_PER_DATASET = 200


@dataclass(frozen=True, slots=True)
class CorpusRecord:
    dataset: str
    sequence_id: str
    sequence: str
    donor: str
    flow_class: str
    chain: str | None

    @property
    def row_key(self) -> str:
        return f"{self.dataset}\0{self.sequence_id}"


@dataclass(frozen=True, slots=True)
class CandidateEvidence:
    """Immutable source identity and raw Cell Ranger comparison fields.

    Empty CSV cells are None; other cells, including counts and productivity,
    remain strings. A consumer can explicitly parse these when comparing them.
    Sequence data and its checksum belong to the joined CorpusRecord/report.
    """

    dataset: str
    sequence_id: str
    donor: str
    flow_class: str
    chain: str | None
    v_gene: str | None
    d_gene: str | None
    j_gene: str | None
    c_gene: str | None
    productive: str | None
    cdr3: str | None
    cdr3_nt: str | None
    reads: str | None
    umis: str | None

    @property
    def row_key(self) -> str:
        return f"{self.dataset}\0{self.sequence_id}"

    @property
    def normalized_genes(self) -> dict[str, tuple[str, ...]]:
        """Return a fresh comparison view without changing raw calls."""
        return {
            segment: normalize_gene(getattr(self, f"{segment}_gene"))
            for segment in ("v", "d", "j", "c")
        }


def stable_rank(
    record: CorpusRecord, *, seed: str, algorithm_version: int = 1
) -> str:
    """Hash immutable source identity using a versioned, UTF-8 payload."""
    payload = f"{algorithm_version}\0{seed}\0{record.row_key}".encode()
    return hashlib.sha256(payload).hexdigest()


def selection_metadata(
    *, per_dataset: int = DEFAULT_PER_DATASET, seed: str = DEFAULT_SEED
) -> dict[str, str | int]:
    """Return JSON-ready configuration for a selection report header."""
    if (
        isinstance(per_dataset, bool)
        or not isinstance(per_dataset, int)
        or per_dataset < 0
    ):
        raise ValueError("per_dataset must be a nonnegative integer")
    return {
        "algorithm_version": ALGORITHM_VERSION,
        "seed": seed,
        "per_dataset": per_dataset,
    }


def select_sweep(
    records: Iterable[CorpusRecord],
    *,
    per_dataset: int = DEFAULT_PER_DATASET,
    seed: str = DEFAULT_SEED,
) -> list[CorpusRecord]:
    """Choose sorted stratum representatives, then fill quotas by stable hash.

    Datasets are ordered lexically; within each dataset results are hash-ordered.
    When the quota is smaller than the stratum count, the first sorted strata
    receive slots. None sorts before a known chain. Duplicate source identities
    are cardinality errors, even if their sequences match. Quota zero is valid.
    Use selection_metadata with the same arguments to record configuration.
    """
    selection_metadata(per_dataset=per_dataset, seed=seed)
    materialized = list(records)
    by_dataset: dict[str, list[CorpusRecord]] = defaultdict(list)
    seen: set[str] = set()
    for record in materialized:
        _validate_identity(record.dataset, "dataset")
        _validate_identity(record.sequence_id, "sequence_id")
        if record.row_key in seen:
            raise ValueError(f"corpus cardinality error: duplicate {record.row_key!r}")
        seen.add(record.row_key)
        by_dataset[record.dataset].append(record)

    selected: list[CorpusRecord] = []
    for dataset in sorted(by_dataset):
        candidates = by_dataset[dataset]
        quota = min(per_dataset, len(candidates))
        by_stratum: dict[tuple[str, str, str | None], list[CorpusRecord]] = defaultdict(list)
        for record in candidates:
            by_stratum[(record.donor, record.flow_class, record.chain)].append(record)

        ranks = {
            r.row_key: stable_rank(r, seed=seed, algorithm_version=ALGORITHM_VERSION)
            for r in candidates
        }
        # The row key also breaks a hypothetical digest collision deterministically.
        def rank(record):
            return ranks[record.row_key], record.row_key

        strata = sorted(
            by_stratum,
            key=lambda value: (
                tuple("" if item is None else item for item in value),
                tuple(item is not None for item in value),
            ),
        )
        chosen: dict[str, CorpusRecord] = {}
        for stratum in strata[:quota]:
            record = min(by_stratum[stratum], key=rank)
            chosen[record.row_key] = record
        for record in sorted(candidates, key=rank):
            if len(chosen) >= quota:
                break
            chosen.setdefault(record.row_key, record)
        selected.extend(sorted(chosen.values(), key=rank))
    return selected


def normalize_gene(raw: str | None) -> tuple[str, ...]:
    """Normalize comma-separated allele ties for comparison/ranking only.

    Trim surrounding whitespace and allele suffixes; deduplicate and sort names.
    The Cell Ranger missing-call sentinel 'None' and empty calls become ().
    Distinct genes and constant-region names are never collapsed heuristically.
    """
    if raw is None:
        return ()
    genes = {
        call.strip().split("*", 1)[0]
        for call in raw.split(",")
        if call.strip() not in ("", "None")
    }
    return tuple(sorted(genes))


def _validate_identity(value: str | None, field: str) -> None:
    if not isinstance(value, str) or not value or "\0" in value:
        raise ValueError(f"{field} must be a nonempty string without NUL: {value!r}")


def _read_strings(path: str | PathLike, columns: Iterable[str]) -> pl.DataFrame:
    columns = tuple(columns)
    frame = pl.read_csv(
        path,
        schema_overrides={column: pl.String for column in columns},
        infer_schema=False,
    )
    missing = set(columns) - set(frame.columns)
    if missing:
        raise ValueError(f"missing required columns in {path}: {sorted(missing)}")
    return frame.select(columns)


def load_cellranger_annotations(
    manifest_path: str | PathLike,
    annotation_paths: Mapping[str, str | PathLike],
    *,
    dataset_column: str = "dataset",
    donor_column: str = "donor",
    flow_class_column: str = "flow_class",
) -> dict[str, CandidateEvidence]:
    """Read explicitly mapped CSV files with a one-row-per-dataset manifest.

    annotation_paths must cover exactly the manifest datasets. The caller owns
    filesystem discovery; no directory layout, home path or corpus is assumed.
    All selected CSV columns have explicit string schemas and inference is
    disabled for additional columns. Repeated contig IDs across datasets are
    valid; duplicate manifest datasets or within-dataset contigs raise ValueError
    identifying the cardinality violation. Output is sorted by dataset/contig ID.
    """
    manifest = _read_strings(
        manifest_path, (dataset_column, donor_column, flow_class_column)
    )
    samples: dict[str, tuple[str, str]] = {}
    for dataset, donor, flow_class in manifest.iter_rows():
        _validate_identity(dataset, dataset_column)
        _validate_identity(donor, donor_column)
        _validate_identity(flow_class, flow_class_column)
        if dataset in samples:
            raise ValueError(f"manifest cardinality error: duplicate dataset {dataset!r}")
        samples[dataset] = (donor, flow_class)
    if set(samples) != set(annotation_paths):
        raise ValueError("annotation datasets must exactly match manifest datasets")

    evidence: dict[str, CandidateEvidence] = {}
    fields = (
        "chain", "v_gene", "d_gene", "j_gene", "c_gene", "productive",
        "cdr3", "cdr3_nt", "reads", "umis",
    )
    for dataset in sorted(samples):
        donor, flow_class = samples[dataset]
        frame = _read_strings(annotation_paths[dataset], ("contig_id", *fields))
        for row in frame.iter_rows(named=True):
            sequence_id = row.pop("contig_id")
            _validate_identity(sequence_id, "contig_id")
            item = CandidateEvidence(dataset, sequence_id, donor, flow_class, **row)
            if item.row_key in evidence:
                raise ValueError(
                    f"annotation cardinality error: dataset {dataset!r}, "
                    f"duplicate contig_id {sequence_id!r}"
                )
            evidence[item.row_key] = item
    return dict(sorted(evidence.items()))
