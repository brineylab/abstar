"""Short-J (ruling B1) calibration: false J calls on random tails after a V, and J accuracy on 3'-truncated reads
by J nt left. Measures the installed build; not a CI dependency.

    .venv/bin/python scripts/benchmark/short_j_null.py [--n 100000]

Null: for each of five V alleles, the full allele plus random tails (seeded, so reruns use identical reads), first
100,000 tails of 250 nt (record 0029's null), then 100,000 of 1-150 nt. Reported per 100,000 reads.
Truncated: the benchmark's igh_truncated set (5,000 reads, seed 20261205), J call accuracy (primary in
j_call_valid) by J nt left in the read.
"""
import argparse
import csv
import random

import polars as pl

import abstar

VS = ("IGHV3-23*01", "IGKV3-20*01", "IGLV2-14*01", "TRAV1-1*01", "TRBV20-1*01")
BUCKETS = ((0, 10), (10, 15), (15, 20), (20, 10**6))


def null(n):
    path = abstar.germline.resolve("human")[0] / "alleles.tsv"
    seqs = {r["name"]: r["sequence_gapped"].replace(".", "").upper() for r in csv.DictReader(open(path), delimiter="\t")}
    rng = random.Random(7)
    out = {}
    for v in VS:
        for label, length in (("250", lambda: 250), ("1-150", lambda: rng.randint(1, 150))):
            reads = [(f"n{i}", seqs[v] + "".join(rng.choice("ACGT") for _ in range(length()))) for i in range(n)]
            out[(v[:4], label)] = abstar.run(reads)["j_call"].is_not_null().sum() * 100_000 / n
    return out


def truncated():
    t = abstar.generate(5000, "IGH", seed=20261205, preset="memory", truncate_5p="uniform:0-60", truncate_3p="uniform:0-40")
    c = abstar.run(list(zip(t["sequence_id"], t["sequence"])))
    m = t.join(c.select("sequence_id", pl.col("j_call").alias("call")), on="sequence_id")
    m = m.filter(~pl.col("unrecoverable_j"), pl.col("j_call_valid").list.len() > 0)
    m = m.with_columns(left=pl.col("j_sequence_end") - pl.col("j_sequence_start") + 1,
                       ok=pl.struct("call", "j_call_valid").map_elements(lambda r: r["call"] in r["j_call_valid"], return_dtype=pl.Boolean))
    rows = [(f"{lo}-{hi - 1}" if hi < 10**6 else f">={lo}", x.height, x["ok"].mean())
            for lo, hi in BUCKETS for x in [m.filter((pl.col("left") >= lo) & (pl.col("left") < hi))]]
    return rows, m["ok"].mean()


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--n", type=int, default=100_000)
    a = p.parse_args()
    rows, overall = truncated()
    print("igh_truncated J call by J nt left: " + "; ".join(f"{b} (n={n}) {100 * acc:.1f}%" for b, n, acc in rows)
          + f"; all {100 * overall:.2f}%")
    res = null(a.n)
    print("false J per 100k, tails 250 nt / 1-150 nt: "
          + "; ".join(f"{v} {res[(v, '250')]:.0f} / {res[(v, '1-150')]:.0f}" for v in dict.fromkeys(k[0] for k in res)))


if __name__ == "__main__":
    main()
