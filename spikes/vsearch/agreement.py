"""Gene-level agreement of spike V calls with legacy abstar v_call.

usage: python agreement.py <preds.tsv> <baseline.parquet>
Allele suffixes are stripped; a record agrees if the spike's tied best genes and
the legacy (comma-separated, possibly tied) genes share any gene.
"""
import sys

import polars as pl

genes = lambda col: (
    pl.col(col).str.split(",").list.eval(pl.element().str.split("*").list.first()).list.unique()
)
p = pl.read_csv(sys.argv[1], separator="\t", schema_overrides={"best_genes": pl.Utf8})
b = pl.read_parquet(sys.argv[2], columns=["corpus_ordinal", "v_call", "locus"])
d = (
    b.join(p, on="corpus_ordinal", how="left")
    .filter(pl.col("v_call").is_not_null())
    .with_columns(legacy=genes("v_call"), spike=genes("best_genes"))
    .with_columns(agree=pl.col("legacy").list.set_intersection("spike").list.len() > 0)
)
print(f"records compared: {d.height} (baseline rows {b.height}, null v_call {b['v_call'].null_count()})")
print(f"agreement: {d['agree'].mean():.4%}")
print(d.group_by("locus").agg(n=pl.len(), agreement=pl.col("agree").mean()).sort("locus"))
