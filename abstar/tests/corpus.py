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
from collections import Counter, defaultdict
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from importlib.resources import files
from os import PathLike
from pathlib import Path
from types import MappingProxyType

from abutils import Sequence
from Bio import SeqIO
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


def _gene_pattern(locus, segment):
    if segment == "v":
        return locus + r"V[1-9][0-9]*D?-(?:[1-9][0-9]*(?:-[1-9][0-9]*)?D?|NL[1-9][0-9]*)"
    if segment == "d":
        return r"IGHD[1-9][0-9]*-[1-9][0-9]*" if locus == "IGH" else r"(?!)"
    if segment == "j":
        return locus + {"IGH": r"J[1-6]", "IGK": r"J[1-5]", "IGL": r"J[1-7]"}[locus]
    return {"IGH": r"IGH(?:M|D|E|A[12]|G[1-4]|G4A)", "IGK": r"IGKC", "IGL": r"IGLC[1-7]"}[locus]


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
        pattern = _gene_pattern(expected["locus"], field[0])
        if isinstance(call, str):
            valid = re.fullmatch(pattern + r"\*(?:[0-9]{2,}|i[0-9]{2,})(?:_[acgt][0-9]+[acgt])*", call)
        else:
            valid = (isinstance(call, list) and bool(call)
                     and all(isinstance(gene, str) and re.fullmatch(
                         pattern, gene) for gene in call)
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


def _packaged_bcr_references():
    """Read fixed package resources, never the configurable user database path."""
    root = files("abstar").joinpath("germline_dbs").joinpath("bcr").joinpath("human")
    references = {}
    for kind, segment in (("ungapped", "v"), ("ungapped", "j"), ("imgt_gapped", "v")):
        entries = {}
        with root.joinpath(kind).joinpath(segment + ".fasta").open("r", encoding="ascii") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                name = record.id.removesuffix("__homo_sapiens")
                if name in entries:
                    raise ValueError("duplicate named packaged germline reference")
                entries[name] = str(record.seq)
        references[kind, segment] = entries
    return references


def _trace_evidence(trace, sequence, references):
    """Authenticate retained bases and recover the reference-to-query mapping."""
    if not isinstance(trace, dict):
        raise ValueError("missing retained segment trace")
    starts = (trace.get("query_start"), trace.get("germline_start"))
    ends = (trace.get("query_end"), trace.get("germline_end"))
    if any(type(x) is not int for x in (*starts, *ends)) or not (
        0 <= starts[0] < ends[0] <= len(sequence) and 0 <= starts[1] < ends[1]
    ):
        raise ValueError("invalid retained trace coordinates")
    query, germline = trace.get("query_aligned"), trace.get("germline_aligned")
    if (not isinstance(query, str) or not isinstance(germline, str)
            or not re.fullmatch(r"[ACGTN-]+", query) or not re.fullmatch(r"[ACGTN-]+", germline)
            or len(query) != len(germline)
            or query.replace("-", "") != sequence[starts[0]:ends[0]]
            or len(germline.replace("-", "")) != ends[1] - starts[1]
            or any(q == g == "-" for q, g in zip(query, germline))):
        raise ValueError("retained trace does not match query/reference spans")
    reference = trace.get("reference")
    if (not isinstance(reference, str) or reference not in references
            or ends[1] > len(references[reference])
            or germline.replace("-", "") != references[reference][starts[1]:ends[1]]):
        raise ValueError("retained reference trace does not match its named packaged germline")
    qpos, gpos = starts
    mapping = {}
    for q, g in zip(query, germline):
        if g != "-":
            mapping[gpos] = qpos if q != "-" else None
            gpos += 1
        if q != "-":
            qpos += 1
    indels = {"v_insertions": [], "v_deletions": []}
    for field, gapped, bases in (("v_insertions", germline, query),
                                  ("v_deletions", query, germline)):
        for match in re.finditer("-+", gapped):
            start = starts[0] + len(query[:match.start()].replace("-", ""))
            end = start + len(match[0]) if field == "v_insertions" else start
            indels[field].append({"query_start": start, "query_end": end,
                                  "sequence": bases[match.start():match.end()]})
    return mapping, indels


def _mapped_anchor(trace, mapping, offset, *, project=False):
    if type(offset) is not int:
        raise ValueError("anchor germline offset must be an integer")
    if offset in mapping and mapping[offset] is not None:
        return mapping[offset]
    distance = trace["germline_start"] - offset
    if project and distance in (1, 2) and trace.get("anchor_projection"):
        return trace["query_start"] - distance
    raise ValueError("anchor is not supported by the retained trace")


def _validate_alignment_evidence(expected, source, sequence, references):
    oriented = str(Seq(sequence).reverse_complement()) if expected["rev_comp"] else sequence
    alignment = source["alignment"]
    if (alignment.get("germline_receptor"), alignment.get("germline_database")) != ("bcr", "human"):
        raise ValueError("curated references require the packaged human BCR database")
    maps, indels = {}, {}
    for segment in ("v", "j"):
        trace = alignment.get(segment)
        maps[segment], indels[segment] = _trace_evidence(trace, oriented, references["ungapped", segment])
        reference = trace.get("reference")
        call = expected[segment + "_call"]
        if not isinstance(reference, str) or not (
            reference == call if isinstance(call, str) else reference.split("*", 1)[0] in call
        ):
            raise ValueError("expected call does not include the retained reference")
    v, j = alignment["v"], alignment["j"]
    if ((expected["v_sequence_start"], expected["v_sequence_end"])
            != (v["query_start"], v["query_end"])
            or (expected["j_sequence_start"], expected["j_sequence_end"])
            != (max(v["query_end"], j["query_start"]), j["query_end"])):
        raise ValueError("expected segment spans disagree with retained traces and V-priority ownership")
    if not expected["v_sequence_end"] <= expected["j_sequence_start"] < expected["j_sequence_end"]:
        raise ValueError("expected V/J spans are out of order")
    for field in ("v_insertions", "v_deletions"):
        if field in expected and expected[field] != indels["v"][field]:
            raise ValueError("expected indel boundaries/bases disagree with retained reference trace")
    voffset, joffset = alignment.get("v_imgt104_ungapped_offset"), alignment.get("j_anchor_germline_offset")
    if type(voffset) is not int or type(joffset) is not int:
        raise ValueError("anchor germline offsets must be integers")
    gapped_v = references["imgt_gapped", "v"].get(v["reference"], "")
    if (gapped_v.replace(".", "") != references["ungapped", "v"][v["reference"]]
            or gapped_v[309:312] not in ("TGT", "TGC")
            or voffset != len(gapped_v[:309].replace(".", ""))):
        raise ValueError("V anchor must be IMGT104 in the named packaged gapped reference")
    # The short packaged J references have one locus-compatible W/F-G-X-G motif.
    # Search nucleotide offsets directly, including J sequences starting mid-codon.
    jreference = references["ungapped", "j"][j["reference"]]
    anchor_pattern = ("TGG" if expected["locus"] == "IGH" else "TT[TC]") + r"GG[ACGT][ACGT]{3}GG[ACGT]"
    joffsets = [i for i in range(len(jreference)) if re.match(anchor_pattern, jreference[i:])]
    if (len(joffsets) != 1 or joffset != joffsets[0]
            or alignment.get("j_anchor_germline_codon") != jreference[joffset:joffset + 3]):
        raise ValueError("J anchor must be the W/F-G-X-G motif in its named packaged reference")
    vanchor = _mapped_anchor(v, maps["v"], voffset)
    janchor = _mapped_anchor(j, maps["j"], joffset, project=True)
    if (vanchor != expected["junction_start"] or vanchor != alignment.get("v_imgt104_query_start")
            or janchor != expected["junction_end"] - 3 or janchor != alignment.get("j_anchor_query_start")):
        raise ValueError("expected junction anchors disagree with retained mapping")
    motif = "W" if expected["locus"] == "IGH" else "F"
    scope = alignment.get("coding_scope")
    if scope == "through_secondary_j_repeat":
        endpoint = alignment.get("j_secondary_repeat")
        _trace_evidence(endpoint, oriented, references["ungapped", "j"])
    elif scope == "through_primary_j":
        endpoint = j
    else:
        raise ValueError("coding scope must name the primary or retained secondary J trace")
    start, end = alignment.get("coding_start"), alignment.get("coding_end")
    if (type(start) is not int or type(end) is not int
            or start != v["query_start"] + (-v["germline_start"] % 3)
            or end != endpoint["query_end"] - (endpoint["query_end"] - start) % 3
            or not 0 <= start <= vanchor < janchor + 3 <= end <= len(oriented)
            or (vanchor - start) % 3
            or str(Seq(oriented[start:end]).translate()) != alignment.get("coding_translation")):
        raise ValueError("coding origin, endpoint, translation or junction frame lacks trace support")
    issues = []
    if "*" in alignment["coding_translation"]:
        issues.append("stop codon(s)")
    if set(oriented) - set("ACGT"):
        issues.append("ambiguous nucleotide(s)")
    if expected["junction_aa"][0] != "C":
        issues.append("junction does not start with conserved C")
    if expected["junction_aa"][-1] != motif:
        issues.append("junction does not end with conserved " + motif)
    if expected["productivity_issues"] != issues or expected["productive"] != (not issues):
        raise ValueError("productivity/reason codes disagree with the retained ORF and anchor evidence")


PILOT_JUNCTION_IDS = frozenset({
    "ATCATCTTCAGCAACT-1_contig_1", "GTTACAGCACATAACC-1_contig_2",
})
_REQUIRED_BUCKETS = {
    "pilot_loss": 8, "pilot_junction_disagreement": 2,
    "concordant_IGH": 2, "concordant_IGK": 2, "concordant_IGL": 2,
    "productivity_disagreement_IGH": 2, "productivity_disagreement_IGK": 2,
    "productivity_disagreement_IGL": 2, "tied_call_IGH": 2, "tied_call_IGK": 2,
    "tied_call_IGL": 2, "insertion": 2, "deletion": 2, "no_d_IGH": 2,
    "shortest_junction": 1, "longest_junction": 1,
}


def _validate_nomination(case):
    """Verify bucket evidence locally; labels alone cannot nominate a record."""
    if len(case.selection_reasons) != 1 or case.selection_reasons[0] not in _REQUIRED_BUCKETS:
        raise ValueError("each case requires one supported nomination bucket")
    reason = case.selection_reasons[0]
    expected, source = case.expected, case.source
    cr, observed = source["cellranger"], source.get("abstar_evidence", {})
    if not isinstance(observed, Mapping):
        raise ValueError("reported nomination evidence must be a mapping")
    if cr.get("chain") != expected["locus"]:
        raise ValueError("nomination locus disagrees with source and expected locus")
    if reason == "pilot_loss":
        valid = (case.dataset == "1279068" and case.sequence_id in PILOT_LOSS_IDS
                 and observed.get("status") == "missing_output" and observed.get("scope") == "record")
    elif reason == "pilot_junction_disagreement":
        valid = (case.dataset == "1279068" and case.sequence_id in PILOT_JUNCTION_IDS
                 and isinstance(observed.get("junction"), str) and bool(cr.get("cdr3_nt"))
                 and observed["junction"] != cr["cdr3_nt"]
                 and observed["junction"] == expected["junction"])
    else:
        raw_calls, junction_evidence, productivity = (
            observed.get(field) for field in ("raw_calls", "junction", "productivity")
        )
        if not all(isinstance(value, Mapping) for value in (raw_calls, junction_evidence, productivity)):
            raise ValueError("calls, junction and productivity nomination evidence must be mappings")
        calls = raw_calls.get("abstar")
        if not isinstance(calls, Mapping):
            raise ValueError("abstar nomination calls must be a mapping")
        junction = junction_evidence.get("abstar_nt")
        reported_productive = productivity.get("abstar")
        raw_productive = cr.get("productive")
        if (observed.get("status") != "annotated" or type(reported_productive) is not bool
                or raw_productive not in ("true", "false")
                or productivity.get("cellranger") is not (raw_productive == "true")):
            raise ValueError("missing or inconsistent reported nomination evidence")
        cr_productive = raw_productive == "true"
        junction_agrees = bool(junction) and junction == cr.get("cdr3_nt") == expected["junction"]
        if reason.startswith("concordant_"):
            valid = (reason == "concordant_" + expected["locus"] and junction_agrees
                     and reported_productive == cr_productive == expected["productive"]
                     and all(normalize_gene(calls.get(s)) == normalize_gene(cr.get(s + "_gene"))
                             and bool(normalize_gene(calls.get(s))) for s in ("v", "j")))
        elif reason.startswith("productivity_disagreement_"):
            valid = (reason == "productivity_disagreement_" + expected["locus"]
                     and reported_productive != cr_productive)
        elif reason.startswith("tied_call_"):
            valid = (reason == "tied_call_" + expected["locus"]
                     and any(isinstance(calls.get(s), str)
                             and len({x.strip() for x in calls[s].split(",") if x.strip()}) > 1
                             for s in ("v", "j")))
        elif reason in ("insertion", "deletion"):
            field = "v_" + reason + "s"
            reported_indels = observed.get("indels")
            if not isinstance(reported_indels, Mapping):
                raise ValueError("reported nomination indels must be a mapping")
            raw = reported_indels.get(field)
            parsed = re.fullmatch(r"[0-9]+(?:-[0-9]+)?:([1-9][0-9]*)>([ACGT]+)", raw or "")
            retained = expected.get(field, ())
            valid = (parsed is not None and bool(retained)
                     and int(parsed[1]) == len(parsed[2]) == sum(len(x["sequence"]) for x in retained))
        elif reason == "no_d_IGH":
            residual = source["alignment"].get("d_residual", {})
            if not isinstance(residual, Mapping):
                raise ValueError("no-D residual evidence must be a mapping")
            start, end = expected["v_sequence_end"], expected["j_sequence_start"]
            oriented = str(Seq(case.sequence).reverse_complement()) if expected["rev_comp"] else case.sequence
            valid = (expected["locus"] == "IGH" and "d_call" in expected and expected["d_call"] is None
                     and not normalize_gene(calls.get("d")) and not normalize_gene(cr.get("d_gene"))
                     and 0 < end - start <= 4 and residual.get("query_start") == start
                     and residual.get("query_end") == end and residual.get("sequence") == oriented[start:end])
        else:
            positions = source.get("selection", {}).get("bucket_positions", {})
            valid = (junction_agrees and type(positions.get(reason)) is int and positions[reason] > 0)
    if not valid:
        raise ValueError(f"unsupported {reason} nomination for {case.dataset}/{case.sequence_id}")


def validate_real_bcr_cohort(cases):
    """Check all required buckets and original pilot identities without annotation."""
    cases = tuple(cases)
    keys = {(case.dataset, case.sequence_id) for case in cases}
    if len(keys) != len(cases):
        raise ValueError("cohort contains duplicated source records")
    for case in cases:
        _validate_nomination(case)
    counts = Counter(reason for case in cases for reason in case.selection_reasons)
    if counts != _REQUIRED_BUCKETS:
        raise ValueError("cohort does not contain every required bucket count")
    for reason, ids in (("pilot_loss", PILOT_LOSS_IDS), ("pilot_junction_disagreement", PILOT_JUNCTION_IDS)):
        actual = {(case.dataset, case.sequence_id) for case in cases if reason in case.selection_reasons}
        if actual != {("1279068", identifier) for identifier in ids}:
            raise ValueError("cohort does not contain the exact original pilot identities")
    lengths = [len(case.expected["junction"]) for case in cases]
    for case in cases:
        reason = case.selection_reasons[0]
        if reason in ("shortest_junction", "longest_junction"):
            endpoint = min(lengths) if reason == "shortest_junction" else max(lengths)
            if len(case.expected["junction"]) != endpoint:
                raise ValueError("selected junction extreme disagrees with the retained cohort")


def load_real_bcr_cases(directory=None) -> tuple[RealBCRCase, ...]:
    """Load checked, recursively immutable cases without annotation or MMseqs.

    FASTA and JSON order jointly identify records, allowing unchanged external
    contig IDs to recur across datasets. Coordinates are in the oriented full
    input query (zero-based, half-open), never the trimmed VDJ or AIRR space.
    Every call reparses files and returns fresh cases and immutable containers.
    """
    complete_cohort = directory is None
    directory = REAL_BCR_DIRECTORY if complete_cohort else Path(directory)
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
    references = _packaged_bcr_references()
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
        _validate_alignment_evidence(raw["expected"], source, sequence, references)
        case = RealBCRCase(sequence=sequence, **raw)
        _validate_nomination(case)
        cases.append(case)
    if complete_cohort:
        validate_real_bcr_cohort(cases)
    return tuple(cases)
