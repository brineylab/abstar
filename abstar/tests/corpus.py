# Copyright (c) 2025 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

"""Source evidence and deterministic selection for optional BCR discovery.

These helpers never annotate sequences or derive expected biological answers.
Cell Ranger values are comparison evidence, retained without type inference.
"""

import hashlib
import json
import re
from collections import defaultdict
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from os import PathLike
from pathlib import Path
from types import MappingProxyType

from abutils import Sequence
from Bio.Seq import Seq

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


REAL_BCR_DIRECTORY = Path(__file__).parent / "data" / "real_bcr"
PILOT_LOSS_IDS = (
    "ACGATACCAGGTTTCA-1_contig_2", "CAAGATCAGAGCTTCT-1_contig_2",
    "CCATGTCCAGTCTTCC-1_contig_1", "CTAAGACAGCAATCTC-1_contig_2",
    "CTGTTTACAGGTGCCT-1_contig_1", "GCTGCTTAGAAACGAG-1_contig_2",
    "GGTGCGTAGGGAAACA-1_contig_1", "GTAACTGAGTATTGGA-1_contig_2",
)
_REQUIRED_EXPECTED = frozenset({
    "locus", "rev_comp", "status", "v_call", "j_call", "v_sequence_start",
    "v_sequence_end", "j_sequence_start", "j_sequence_end", "junction_start",
    "junction_end", "junction", "junction_aa", "cdr3", "cdr3_aa",
    "productive", "productivity_issues",
})
_EXPECTED_FIELDS = _REQUIRED_EXPECTED | {
    "d_call", "c_call", "d_sequence_start", "d_sequence_end",
    "c_sequence_start", "c_sequence_end", "v_insertions", "v_deletions",
}


def _freeze(value):
    """Copy every container so callers cannot mutate shared fixture state."""
    if isinstance(value, Mapping):
        return MappingProxyType({key: _freeze(item) for key, item in value.items()})
    if isinstance(value, (list, tuple)):
        return tuple(_freeze(item) for item in value)
    return value


@dataclass(frozen=True, slots=True)
class RealBCRCase:
    dataset: str
    sequence_id: str
    sequence: str
    sequence_sha256: str
    selection_reasons: tuple[str, ...]
    source: Mapping[str, object]
    expected: Mapping[str, object]
    evidence: tuple[str, ...]

    def __post_init__(self):
        for name in ("selection_reasons", "source", "expected", "evidence"):
            object.__setattr__(self, name, _freeze(getattr(self, name)))

    def as_sequence(self) -> Sequence:
        return Sequence(self.sequence, id=self.sequence_id)


def _nonempty_strings(value, field, *, allow_empty=False):
    if not isinstance(value, list) or (not value and not allow_empty) or any(
        not isinstance(item, str) or not item.strip() for item in value
    ):
        raise ValueError(f"{field} must contain nonempty strings")


def _validate_expected(expected, sequence):
    if not isinstance(expected, dict) or not _REQUIRED_EXPECTED <= expected.keys():
        raise ValueError("missing required expected fields")
    if expected.keys() - _EXPECTED_FIELDS:
        raise ValueError("unknown expected field")
    if expected["locus"] not in ("IGH", "IGK", "IGL"):
        raise ValueError("unsupported BCR locus")
    if expected["status"] != "annotated":
        raise ValueError("curated biological expectations require annotated status")
    for field in ("rev_comp", "productive"):
        if type(expected[field]) is not bool:
            raise ValueError(f"{field} must be Boolean")
    for field in ("v_call", "d_call", "j_call", "c_call"):
        if field not in expected:
            continue
        call = expected[field]
        if call is None and field in ("d_call", "c_call"):
            continue
        if isinstance(call, str):
            valid = re.fullmatch(r"IG[HKL][A-Za-z0-9/.-]+\*[A-Za-z0-9_]+", call)
        else:
            valid = (isinstance(call, list) and bool(call)
                     and all(isinstance(gene, str) and re.fullmatch(
                         r"IG[HKL][A-Za-z0-9/.-]+", gene) for gene in call)
                     and call == sorted(set(call)))
        if not valid:
            raise ValueError(f"{field} requires an exact allele or sorted gene-level allowed set")
    for prefix in ("v_sequence", "d_sequence", "j_sequence", "c_sequence", "junction"):
        start, end = expected.get(prefix + "_start"), expected.get(prefix + "_end")
        if start is None and end is None and prefix in ("d_sequence", "c_sequence"):
            continue
        if type(start) is not int or type(end) is not int or not 0 <= start < end <= len(sequence):
            raise ValueError(f"{prefix} requires zero-based half-open coordinates")
    for field in ("v_insertions", "v_deletions"):
        if field not in expected:
            continue
        indels = expected[field]
        if not isinstance(indels, list):
            raise ValueError(f"{field} must be an ordered indel list")
        for indel in indels:
            if not isinstance(indel, dict) or indel.keys() != {"query_start", "query_end", "sequence"}:
                raise ValueError("indels require query boundaries and observed/removed bases")
            start, end, bases = (indel[k] for k in ("query_start", "query_end", "sequence"))
            if (type(start) is not int or type(end) is not int or not 0 <= start <= end <= len(sequence)
                    or not isinstance(bases, str) or not re.fullmatch(r"[ACGT]+", bases)):
                raise ValueError("invalid indel coordinates or sequence")
            if field == "v_deletions" and start != end:
                raise ValueError("a deletion has zero width in query coordinates")
            if field == "v_insertions" and end - start != len(bases):
                raise ValueError("insertion interval must span inserted bases")
    oriented = sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1] if expected["rev_comp"] else sequence
    junction = expected["junction"]
    if junction != oriented[expected["junction_start"]:expected["junction_end"]]:
        raise ValueError("junction does not match oriented query interval")
    for field in ("junction", "junction_aa", "cdr3", "cdr3_aa"):
        if not isinstance(expected[field], str) or not expected[field]:
            raise ValueError(f"{field} must be a nonempty string")
    if len(junction) % 3 or str(Seq(junction).translate()) != expected["junction_aa"]:
        raise ValueError("junction translation does not match expected amino acids")
    if expected["cdr3"] != junction[3:-3] or expected["cdr3_aa"] != expected["junction_aa"][1:-1]:
        raise ValueError("CDR3 must exclude junction anchor codons")
    _nonempty_strings(expected["productivity_issues"], "productivity_issues", allow_empty=True)
    if expected["productive"] == bool(expected["productivity_issues"]):
        raise ValueError("productivity and issue evidence disagree")


def load_real_bcr_cases(directory=None) -> tuple[RealBCRCase, ...]:
    """Load checked, recursively immutable cases without annotation or MMseqs.

    FASTA and JSON order jointly identify records, allowing unchanged external
    contig IDs to recur across datasets. Coordinates are in the oriented full
    input query (zero-based, half-open), never the trimmed VDJ or AIRR space.
    Every call reparses files and returns fresh cases and immutable containers.
    """
    directory = REAL_BCR_DIRECTORY if directory is None else Path(directory)
    raw_cases = json.loads((directory / "cases.json").read_text(encoding="utf-8"))
    if not isinstance(raw_cases, list) or not raw_cases:
        raise ValueError("cases must be a nonempty ordered JSON list")
    records = []
    for line in (directory / "sequences.fasta").read_text(encoding="ascii").splitlines():
        if line.startswith(">"):
            records.append([line[1:], ""])
        elif line:
            if not records or not re.fullmatch(r"[ACGTN]+", line):
                raise ValueError("invalid fixture FASTA sequence")
            records[-1][1] += line
    if len(records) != len(raw_cases):
        raise ValueError("FASTA/JSON record count mismatch")
    cases, seen = [], set()
    fields = set(RealBCRCase.__dataclass_fields__) - {"sequence"}
    for raw, (identifier, sequence) in zip(raw_cases, records):
        if not isinstance(raw, dict) or raw.keys() != fields:
            raise ValueError("invalid case schema")
        for field in ("dataset", "sequence_id"):
            _validate_identity(raw[field], field)
        key = (raw["dataset"], raw["sequence_id"])
        if key in seen:
            raise ValueError("duplicate dataset/sequence_id")
        seen.add(key)
        if identifier != raw["sequence_id"] or not sequence:
            raise ValueError("FASTA/JSON order or identity mismatch")
        if hashlib.sha256(sequence.encode("ascii")).hexdigest() != raw["sequence_sha256"]:
            raise ValueError("sequence SHA-256 mismatch")
        _nonempty_strings(raw["selection_reasons"], "selection_reasons")
        _nonempty_strings(raw["evidence"], "evidence")
        source = raw["source"]
        if not isinstance(source, dict) or not all(isinstance(source.get(k), dict) and source[k]
                for k in ("publication", "cellranger", "alignment")):
            raise ValueError("publication, Cell Ranger and direct alignment evidence are required")
        _validate_expected(raw["expected"], sequence)
        cases.append(RealBCRCase(sequence=sequence, **raw))
    return tuple(cases)
