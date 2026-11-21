"""Score annotator output against `abstar generate` truth (decisions 0025, 0026)."""
from __future__ import annotations

import csv
import math
import re
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import polars as pl

_FLAGS = ("high_mutation", "long_cdr3", "has_indel", "truncated_5p", "truncated_3p", "ambiguous_v", "ambiguous_d", "ambiguous_j")
_EXACT = ("junction", "junction_aa", "productive", "vj_in_frame", "stop_codon")
_REGIONS = ("fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4")
_REGION_PARTS = ("", "_aa", "_start", "_end", "_length")
_MASKS = ("cdr_mask", "gene_segment_mask", "nongermline_mask", "cdr_mask_aa", "gene_segment_mask_aa", "nongermline_mask_aa")
_MUT_PARTS = ("mutations", "mutations_aa", "insertions", "deletions", "mutation_count", "mutation_count_aa",
              "mutation_frequency", "mutation_frequency_aa")
_POS = re.compile(r"^(?:IMGT)?(\d+):")


def _pos(s):
    return int(_POS.match(s).group(1))


def _ceil3(x):
    return -(-x // 3)


def _null_eq(t, c):
    return t is None and c is None


def _gapped_maps(germline_db):
    """V allele name -> IMGT position of each ungapped base (index 0 is base 1)."""
    from . import germline

    path = germline.resolve(germline_db)[0] / "alleles.tsv"
    with open(path, newline="") as f:
        return {
            r["name"]: [i + 1 for i, ch in enumerate(r["sequence_gapped"]) if ch != "."]
            for r in csv.DictReader(f, delimiter="\t")
            if r["segment"] == "V"
        }


def grade(truth: "pl.DataFrame", calls: "pl.DataFrame", germline_db: str = "human") -> "pl.DataFrame":
    """Grade `calls` against `truth`: one row per (field, subset) with n, correct, accuracy.

    Subsets are "all" plus each truth difficulty flag that is true for a row. A field absent
    from `calls` (or `truth`) is reported with n == 0. `germline_db` supplies the IMGT gapping
    that V mutation positions need for the shared-span rule.
    """
    import polars as pl

    tcols, ccols = set(truth.columns), set(calls.columns)
    crows = {r["sequence_id"]: r for r in calls.iter_rows(named=True)}
    gap = {}

    def gmap(allele):
        if not gap:
            gap.update(_gapped_maps(germline_db))
        return gap.get(allele)

    def ok(t, seg):  # recoverable
        return not t.get(f"unrecoverable_{seg}")

    def same_allele(t, c, seg):
        return t.get(f"{seg}_call") is not None and c.get(f"{seg}_call") == t[f"{seg}_call"]

    def call_fn(seg):
        def f(t, c):
            if not ok(t, seg):
                return None
            v = c.get(f"{seg}_call")
            if v is None:
                return bool(t.get("d_call_null_valid")) if seg == "d" else False
            valid = t.get(f"{seg}_call_valid")
            # An empty valid list means the segment left no base in the read; only the true allele is then known.
            return v in valid if valid else v == t.get(f"{seg}_call")
        return f

    def boundary_fn(col, seg):
        def f(t, c):
            if not ok(t, seg):
                return None
            valid, v = t.get(f"{col}_valid") or [], c.get(col)
            return v is None if not valid else v in valid
        return f

    def exact_fn(col):
        return lambda t, c: (t[col] == c.get(col)) if ok(t, "v") and ok(t, "j") else None

    def region_fn(col, segs):
        def f(t, c):
            if not all(ok(t, s) and same_allele(t, c, s) for s in segs):
                return None
            return t[col] == c.get(col)
        return f

    def mask_segs(t):
        return ("v", "j") + (("d",) if t.get("d_call") is not None else ())

    def uncertain(t):
        """Read positions (1-based) whose segment assignment is not fixed by the _valid lists."""
        out = set()
        for col, lo_off, hi_off in (("v_sequence_end", 1, 0), ("d_sequence_start", 0, -1), ("d_sequence_end", 1, 0), ("j_sequence_start", 0, -1)):
            v = t.get(f"{col}_valid") or []
            if len(v) > 1:
                out.update(range(min(v) + lo_off, max(v) + hi_off + 1))
        return out

    def mask_fn(col):
        gcol = "gene_segment_mask_aa" if col.endswith("_aa") else "gene_segment_mask"

        def f(t, c):
            if not all(ok(t, s) and same_allele(t, c, s) for s in mask_segs(t)):
                return None
            tm, cm = t[col], c.get(col)
            if tm is None or cm is None:
                return tm == cm
            if len(tm) != len(cm):
                return False
            tg, cg = t[gcol], c.get(gcol)
            if tg is None or cg is None or len(tg) != len(tm) or len(cg) != len(cm):
                return False
            if col == gcol:  # differences are allowed only where a boundary is ambiguous
                amb, start = uncertain(t), t["v_sequence_start"]
                step, off = (3, 2) if col.endswith("_aa") else (1, 0)
                return all(a == b or (start + step * i + off) in amb for i, (a, b) in enumerate(zip(tm, cm)))
            return all(a == b for a, b, x, y in zip(tm, cm, tg, cg) if x == y)
        return f

    def span(t, c, seg):
        """None (no span), or (ulo, uhi, lo, hi, nlo, nhi): ungapped span, nt filter span, codon span."""
        ts, te, cs, ce = t.get(f"{seg}_germline_start"), t.get(f"{seg}_germline_end"), c.get(f"{seg}_germline_start"), c.get(f"{seg}_germline_end")
        ulo, uhi = max(ts, cs), min(te, ce)
        if seg == "v":
            m = gmap(t["v_call"])
            if m is None or uhi > len(m):  # ponytail: allele missing from germline_db, V span left unrestricted
                return ulo, uhi, 0, 10**9, 0, 10**9
            lo, hi = (m[ulo - 1], m[uhi - 1]) if ulo <= uhi else (1, 0)
            return ulo, uhi, lo, hi, _ceil3(lo + 2), hi // 3
        return ulo, uhi, ulo, uhi, _ceil3(ulo + 2), (uhi - 2) // 3

    def mut_fn(seg, part):
        lists = {
            "mutations": ("mutations",), "mutations_aa": ("mutations_aa",), "insertions": ("insertions",), "deletions": ("deletions",),
            "mutation_count": ("mutations",), "mutation_count_aa": ("mutations_aa",),
            "mutation_frequency": ("mutations", "deletions"), "mutation_frequency_aa": ("mutations_aa",),
        }[part]

        def f(t, c):
            if not (ok(t, seg) and same_allele(t, c, seg)):
                return None
            cols = [f"{seg}_{x}" for x in lists]
            if any(t[x] is None or c.get(x) is None for x in cols) or t[f"{seg}_germline_start"] is None or c.get(f"{seg}_germline_start") is None:
                return all(_null_eq(t[x], c.get(x)) for x in cols) and _null_eq(t[f"{seg}_germline_start"], c.get(f"{seg}_germline_start"))
            ulo, uhi, lo, hi, nlo, nhi = span(t, c, seg)

            def val(row):
                aa = part.endswith("_aa")
                muts = sorted(m for m in row[f"{seg}_{'mutations_aa' if aa else 'mutations'}"] if (nlo <= _pos(m) <= nhi if aa else lo <= _pos(m) <= hi))
                if part in ("mutations", "mutations_aa"):
                    return muts
                if part in ("insertions", "deletions"):
                    return sorted(x for x in row[f"{seg}_{part}"] if lo <= _pos(x) <= hi)
                if part.startswith("mutation_count"):
                    return len(muts)
                if aa:
                    denom = max(0, nhi - nlo + 1)
                else:
                    deleted = sum(len(x.split(":")[2]) for x in row[f"{seg}_deletions"] if lo <= _pos(x) <= hi)
                    denom = max(0, uhi - ulo + 1) - deleted
                return len(muts) / denom if denom > 0 else None

            a, b = val(t), val(c)
            if isinstance(a, float) and isinstance(b, float):
                return math.isclose(a, b, abs_tol=1e-9)
            return a == b
        return f

    fields: dict = {}
    for seg in "vdjc":
        fields[f"{seg}_call"] = ((f"{seg}_call",), call_fn(seg) if seg != "c" else lambda t, c: (c.get("c_call") in t["c_call_valid"]) if t["c_call_valid"] else None)
    for col, seg in (("v_sequence_end", "v"), ("d_sequence_start", "d"), ("d_sequence_end", "d"), ("j_sequence_start", "j")):
        fields[col] = ((col,), boundary_fn(col, seg))
    for col in _EXACT:
        fields[col] = ((col,), exact_fn(col))
    for r in _REGIONS:
        for p in _REGION_PARTS:
            fields[r + p] = ((r + p,), region_fn(r + p, ("v", "j") if r in ("cdr3", "fwr4") else ("v",)))
    for col in _MASKS:
        gcol = "gene_segment_mask_aa" if col.endswith("_aa") else "gene_segment_mask"
        fields[col] = ((col, gcol), mask_fn(col))
    for seg in "vdjc":
        for p in _MUT_PARTS:
            need = {f"{seg}_{p}", f"{seg}_germline_start", f"{seg}_germline_end", f"{seg}_call"}
            if p.startswith("mutation_"):
                need.add(f"{seg}_mutations_aa" if p.endswith("_aa") else f"{seg}_mutations")
            if p == "mutation_frequency":
                need.add(f"{seg}_deletions")
            fields[f"{seg}_{p}"] = (tuple(sorted(need)), mut_fn(seg, p))

    subsets = ["all"] + [f for f in _FLAGS if f in tcols]
    tally = {(f, s): [0, 0] for f in fields for s in subsets}
    # Count and frequency are recomputed from the lists, so truth need not hold the column itself.
    derived = lambda f: f[2:].startswith("mutation_")  # noqa: E731
    live = {
        f: all(x in ccols and (x in tcols or (x == f and derived(f))) for x in cols) and (f in tcols or derived(f))
        for f, (cols, _) in fields.items()
    }
    for t in truth.iter_rows(named=True):
        c = crows.get(t["sequence_id"], {})
        sub = ["all"] + [f for f in subsets[1:] if t.get(f)]
        for name, (_, fn) in fields.items():
            if not live[name]:
                continue
            r = fn(t, c)
            if r is None:
                continue
            for s in sub:
                tally[(name, s)][0] += 1
                tally[(name, s)][1] += bool(r)
    return pl.DataFrame(
        [(f, s, n, k, k / n if n else None) for (f, s), (n, k) in tally.items()],
        schema={"field": pl.String, "subset": pl.String, "n": pl.Int64, "correct": pl.Int64, "accuracy": pl.Float64},
        orient="row",
    )
