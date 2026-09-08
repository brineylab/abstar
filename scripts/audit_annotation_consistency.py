#!/usr/bin/env python3
"""Audit final/native Parquets against a direct oriented-input V(D)J slice.

This is a read-only discovery tool, not a biological adjudicator. Coordinates
must be zero-based and half-open; AIRR TSV coordinates are not accepted.
Exit 1 means inconsistent records were found; 0 means all tested checks passed.
"""
import argparse
from collections import Counter
from functools import lru_cache
from itertools import product
import json
from pathlib import Path

from Bio.Seq import translate
import polars as pl


REGIONS = ("fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4")
SEGMENTS = ("v_sequence", "np1", "d_sequence", "np2", "j_sequence")
MASKS = tuple(m + suffix for m in ("cdr_mask", "gene_segment_mask", "nongermline_mask")
              for suffix in ("", "_aa"))
REQUIRED = {
    "sequence_id", "annotation_status", "sequence_oriented", "sequence_alignment",
    "frame", "v_sequence_start", "j_sequence_end", *SEGMENTS, *MASKS,
    *(r + suffix for r in REGIONS for suffix in ("", "_aa", "_start", "_end")),
}


@lru_cache(maxsize=1)
def _codons():
    codons = ("".join(c) for c in product("ACGTRYSWKMBDHVN", repeat=3))
    return {codon: translate(codon) for codon in codons}


def audit_frame(frame):
    """Return candidate rows; preserve physical Parquet ordinals and exact IDs."""
    missing = REQUIRED - set(frame.columns)
    if missing:
        raise ValueError(f"Missing annotation columns: {', '.join(sorted(missing))}")
    if frame.schema["sequence_id"] != pl.String:
        raise ValueError("sequence_id must have a string schema")
    integer_columns = {"frame", "v_sequence_start", "j_sequence_end",
                       *(r + suffix for r in REGIONS for suffix in ("_start", "_end"))}
    for column in REQUIRED:
        dtype = frame.schema[column]
        if dtype == pl.Null:
            frame = frame.with_columns(pl.col(column).cast(pl.Int64 if column in integer_columns else pl.String))
        elif column in integer_columns:
            if not dtype.is_integer():
                raise ValueError(f"{column} must have an integer schema")
        elif dtype != pl.String:
            raise ValueError(f"{column} must have a string schema")
    statuses = set(frame["annotation_status"].unique().to_list())
    if statuses - {"annotated", "unassigned"}:
        raise ValueError(f"Unexpected annotation_status values: {statuses}")
    frame = frame.with_row_index("parquet_row").filter(pl.col("annotation_status") == "annotated")
    start, end, phase = pl.col("v_sequence_start"), pl.col("j_sequence_end"), pl.col("frame")
    valid = (
        start.is_not_null() & end.is_not_null() & phase.is_not_null()
        & (start >= 0) & (end > start)
        & (end <= pl.col("sequence_oriented").str.len_chars()) & phase.is_in([1, 2, 3])
    ).fill_null(False)
    frame = frame.with_columns(valid.alias("reference_valid"))
    # Invalid coordinates are reported explicitly, never accepted via Python's
    # forgiving slice bounds or negative indexing.
    frame = frame.with_columns(
        pl.when(pl.col("reference_valid"))
        .then(pl.col("sequence_oriented").str.slice(start, (end - start).clip(lower_bound=0)))
        .otherwise(None).alias("vdj_nt"),
        pl.concat_str([pl.col(c).fill_null("") for c in SEGMENTS]).alias("segments_nt"),
        pl.concat_str([pl.col(r) for r in REGIONS]).alias("regions_nt"),
        pl.concat_str([pl.col(r + "_aa") for r in REGIONS]).alias("regions_aa"),
    )
    # Codon translation is independent of regional AA annotation and of the
    # public AA alignment, whose columns can contain codon-spanning gaps.
    frame = frame.with_columns(
        pl.col("vdj_nt").str.slice((phase - 1).clip(lower_bound=0))
        .str.to_uppercase().str.extract_all("...")
        .list.eval(pl.element().replace_strict(_codons(), default="?"))
        .list.join("").alias("vdj_aa"),
    ).with_columns(
        pl.col("vdj_nt").str.len_chars().alias("vdj_nt_length"),
        pl.col("vdj_aa").str.len_chars().alias("vdj_aa_length"),
    )
    checks = {
        "vdj_reference_invalid": ~pl.col("reference_valid"),
        "vdj_alphabet_invalid": ~pl.col("vdj_nt").str.contains("(?i)^[ACGTRYSWKMBDHVN]+$"),
        "segments_nt_content": pl.col("segments_nt") != pl.col("vdj_nt"),
        "alignment_nt_content": pl.col("sequence_alignment").is_null() | (pl.col("sequence_alignment").str.replace_all("-", "", literal=True) != pl.col("vdj_nt")),
        "regions_nt_content": pl.col("regions_nt") != pl.col("vdj_nt"),
        "regions_aa_content": pl.col("regions_aa") != pl.col("vdj_aa"),
    }
    for mask in MASKS:
        length = "vdj_aa_length" if mask.endswith("_aa") else "vdj_nt_length"
        checks[mask + "_null"] = pl.col(mask).is_null()
        checks[mask + "_length"] = pl.col(mask).str.len_chars() != pl.col(length)
    for region in REGIONS:
        rstart, rend = pl.col(region + "_start"), pl.col(region + "_end")
        present = pl.col(region).str.len_chars() > 0
        coordinates = ((rstart >= start) & (rend <= end) & (rend > rstart)).fill_null(False)
        content = pl.col("sequence_oriented").str.slice(rstart, (rend - rstart).clip(lower_bound=0))
        checks[region + "_coordinates"] = present & (
            ~coordinates | (pl.col(region) != content)
        )
        checks[region + "_null"] = pl.col(region).is_null() | pl.col(region + "_aa").is_null()
    for left, right in zip(REGIONS, REGIONS[1:]):
        checks[left + "_" + right + "_boundary"] = pl.col(left + "_end") != pl.col(right + "_start")
    frame = frame.with_columns(
        pl.concat_list([
            pl.when(expr.fill_null(False)).then(pl.lit(name)).otherwise(None)
            for name, expr in checks.items()
        ]).list.drop_nulls().alias("checks")
    )
    columns = ["parquet_row", "sequence_id", "vdj_nt_length", "vdj_aa_length", "checks"]
    columns += [c for c in ("productive", "locus", "v_call", "j_call") if c in frame.columns]
    return frame.filter(pl.col("checks").list.len() > 0).select(columns)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("parquet", type=Path, help="Parquet file or directory of Parquet files")
    parser.add_argument("--report", required=True, type=Path, help="New JSON report outside input and repository trees")
    args = parser.parse_args(argv)
    source, destination = args.parquet.resolve(), args.report.resolve()
    protected = [Path(__file__).resolve().parents[1], source if source.is_dir() else source.parent]
    if any(destination.is_relative_to(root) for root in protected):
        raise ValueError("Report must be outside repository and input trees")
    if destination.exists():
        raise FileExistsError(destination)
    files = sorted(source.glob("*.parquet")) if source.is_dir() else [source]
    if not files:
        raise ValueError("No Parquet files found")
    totals = Counter(rows=0, annotated=0, unassigned=0, inconsistent=0)
    check_counts, samples, candidates = Counter(), [], []
    for path in files:
        frame = pl.read_parquet(path)
        found = audit_frame(frame).to_dicts()
        counts = dict(rows=frame.height,
                      annotated=frame.filter(pl.col("annotation_status") == "annotated").height,
                      unassigned=frame.filter(pl.col("annotation_status") == "unassigned").height,
                      inconsistent=len(found))
        totals.update(counts)
        samples.append(dict(source_parquet=str(path), **counts))
        for candidate in found:
            candidate["source_parquet"] = str(path)
            candidates.append(candidate)
            check_counts.update(candidate["checks"])
    report = dict(reference="sequence_oriented[v_sequence_start:j_sequence_end]; frame is one-based",
                  totals=dict(totals), check_counts=dict(check_counts), samples=samples, candidates=candidates)
    with destination.open("x") as handle:
        json.dump(report, handle, indent=2)
        handle.write("\n")
    print(json.dumps(dict(totals)))
    return int(bool(candidates))


if __name__ == "__main__":
    raise SystemExit(main())
