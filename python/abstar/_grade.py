"""Score annotator output against `abstar generate` truth (decisions 0025, 0026)."""
from __future__ import annotations

import csv
import math
import re
import warnings
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


def _alleles(germline_db):
    """Allele name -> (gapped sequence, anchor or None)."""
    from . import germline

    path = germline.resolve(germline_db)[0] / "alleles.tsv"
    with open(path, newline="") as f:
        return {r["name"]: (r["sequence_gapped"], int(r["anchor"]) if r["anchor"] else None) for r in csv.DictReader(f, delimiter="\t")}


def grade(truth: "pl.DataFrame", calls: "pl.DataFrame", germline_db: str = "human") -> "pl.DataFrame":
    """Grade `calls` against `truth`: one row per (field, subset) with n, correct, accuracy.

    Subsets are "all" plus each truth difficulty flag that is true for a row. A field absent
    from `calls` (or `truth`) is reported with n == 0. `germline_db` supplies the IMGT gapping
    that V mutation positions need for the shared-span rule (a V allele missing from it leaves
    the V span unrestricted). An annotator span edge must match truth, or (for the ambiguous
    cuts) be a valid cut with a consistent germline shift, else the segment's mutation fields are
    wrong. Reported counts and frequencies are always checked.

    Rules worth knowing: an empty `<seg>_call_valid` makes only a null call correct;
    `gene_segment_mask` may differ from truth only at read positions the `_valid` boundary lists
    leave open, and the other masks are compared only where both gene masks agree; aa mutation
    codons count as the generator counts them: V non-gap IMGT codons wholly inside the span, D/J/C
    codons in the generator's frames (J: anchor % 3, D: the V phase, C: continuing J), minus codons
    with a deleted base.
    """
    import polars as pl

    tcols, ccols = set(truth.columns), set(calls.columns)
    crows = {}
    for r in calls.iter_rows(named=True):
        if r["sequence_id"] in crows:
            warnings.warn(f"duplicate sequence_id {r['sequence_id']!r} in calls; the last row is graded", stacklevel=2)
        crows[r["sequence_id"]] = r
    allele_db = {}

    def allele(name):
        if not allele_db:
            allele_db.update(_alleles(germline_db))
        return allele_db.get(name)

    def gmap(name):
        """IMGT position of each ungapped base of a V allele (index 0 is base 1), or None."""
        a = allele(name)
        return None if a is None else [i + 1 for i, ch in enumerate(a[0]) if ch != "."]

    def ok(t, seg):  # recoverable
        return not t.get(f"unrecoverable_{seg}")

    def same_allele(t, c, seg):
        return t.get(f"{seg}_call") is not None and c.get(f"{seg}_call") == t[f"{seg}_call"]

    def call_fn(seg):
        def f(t, c):
            if not ok(t, seg):
                return None
            v, valid = c.get(f"{seg}_call"), t.get(f"{seg}_call_valid")
            if v is None and seg == "d" and t.get("d_call_null_valid"):
                return True
            # An empty valid set means the read holds no base of the segment: only a null call is correct.
            return v is None if not valid else v in valid
        return f

    def boundary_fn(col, seg):
        def f(t, c):
            if not ok(t, seg):
                return None
            valid, v = t.get(f"{col}_valid") or [], c.get(col)
            if v is None and seg == "d" and t.get("d_call_null_valid"):
                return True  # a valid null D has no boundaries
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
        """Read positions (1-based) whose segment assignment the _valid lists (plus the true cut) leave open."""
        out = set()
        for col, lo_off, hi_off in (("v_sequence_end", 1, 0), ("d_sequence_start", 0, -1), ("d_sequence_end", 1, 0), ("j_sequence_start", 0, -1)):
            cuts = set(t.get(f"{col}_valid") or [])
            if t.get(col) is not None:
                cuts.add(t[col])
            if cuts:
                out.update(range(min(cuts) + lo_off, max(cuts) + hi_off + 1))
        return out

    def codon_phase(t):
        """Index (0-based from the first V base) of the first complete codon: codons are in phase with IMGT numbering."""
        m = gmap(t["v_call"])
        if m is None or not 1 <= t["v_germline_start"] <= len(m):
            return 0
        return (3 - (m[t["v_germline_start"] - 1] - 1) % 3) % 3

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
                aa = col.endswith("_aa")
                step, off = (3, codon_phase(t) + 2) if aa else (1, 0)
                return all(a == b or (start + step * i + off) in amb for i, (a, b) in enumerate(zip(tm, cm)))
            return all(a == b for a, b, x, y in zip(tm, cm, tg, cg) if x == y)
        return f

    def frame(t, seg):
        """0-based germline index of the segment's codon 1 (the generator's frames), or None."""
        if seg == "v":
            return None
        if seg == "d":
            try:
                phase = (t["d_sequence_start"] - t["v_sequence_end"] + t["v_germline_end"] - 1) % 3
                return (t["d_germline_start"] - 1 - phase) % 3
            except TypeError:
                return None
        a = allele(t["j_call"]) if t.get("j_call") else None
        if a is None or a[1] is None:
            return None
        return a[1] % 3 if seg == "j" else (a[1] - len(a[0].replace(".", ""))) % 3

    def mkspan(seg, t, ulo, uhi):
        """(ulo, uhi, lo, hi, nlo, nhi): ungapped span, nt filter span, codon filter span (codons wholly inside)."""
        if seg == "v":
            m = gmap(t["v_call"])
            if m is None or uhi > len(m):  # ponytail: allele missing from germline_db, V span left unrestricted
                return ulo, uhi, 0, 10**9, 0, 10**9
            lo, hi = (m[ulo - 1], m[uhi - 1]) if ulo <= uhi else (1, 0)
            return ulo, uhi, lo, hi, _ceil3(lo + 2), hi // 3
        s0 = frame(t, seg)
        if s0 is None:  # no frame known: codons unrestricted
            return ulo, uhi, ulo, uhi, 0, 10**9
        return ulo, uhi, ulo, uhi, max(1, _ceil3(ulo - s0 + 2)), (uhi - s0) // 3

    def aa_denom(row, seg, t, sp):
        """Codons counted for the aa frequency: those wholly inside the span, minus gap codons and codons with a deleted base."""
        nlo, nhi = sp[4], sp[5]
        if nhi >= 10**9:
            return 10**9
        codons = set(range(max(nlo, 1), nhi + 1))
        if seg == "v":
            gapped = allele(t["v_call"])[0]
            codons = {n for n in codons if "." not in gapped[3 * n - 3 : 3 * n]}
            m = gmap(t["v_call"])
            inv = {p: u for u, p in enumerate(m, 1)}
            for x in row["v_deletions"]:
                u = inv.get(_pos(x))
                for i in range(len(x.split(":")[2])) if u else ():
                    if u - 1 + i < len(m):
                        codons.discard(_ceil3(m[u - 1 + i]))
        else:
            s0 = frame(t, seg)
            for x in row[f"{seg}_deletions"]:
                for p in range(_pos(x), _pos(x) + len(x.split(":")[2])):
                    codons.discard((p - s0 - 1) // 3 + 1)
        return len(codons)

    # read column of the ends whose cut is ambiguous; every other end must match truth.
    cuts = {"v": {"end": "v_sequence_end"}, "d": {"start": "d_sequence_start", "end": "d_sequence_end"}, "j": {"start": "j_sequence_start"}}

    def span_consistent(t, c, seg):
        for end in ("start", "end"):
            g = f"{seg}_germline_{end}"
            col = cuts.get(seg, {}).get(end)
            if col is None:
                if t[g] != c[g]:
                    return False
                continue
            sc, st = c.get(col), t.get(col)
            if sc is None or st is None or not (sc in (t.get(f"{col}_valid") or []) or sc == st):
                return False
            # ponytail: assumes no indel between the two cuts, so germline and read shift together
            if c[g] - t[g] != sc - st:
                return False
        return True

    def mut_fn(seg, part):
        lists = {
            "mutations": ("mutations",), "mutations_aa": ("mutations_aa",), "insertions": ("insertions",), "deletions": ("deletions",),
            "mutation_count": ("mutations",), "mutation_count_aa": ("mutations_aa",),
            "mutation_frequency": ("mutations", "deletions"), "mutation_frequency_aa": ("mutations_aa",),
        }[part]
        derived_part = part.startswith("mutation_")

        def val(row, sp, t):
            ulo, uhi, lo, hi, nlo, nhi = sp
            aa = part.endswith("_aa")
            muts = sorted(m for m in row[f"{seg}_{'mutations_aa' if aa else 'mutations'}"] if (nlo <= _pos(m) <= nhi if aa else lo <= _pos(m) <= hi))
            if part in ("mutations", "mutations_aa"):
                return muts
            if part in ("insertions", "deletions"):
                return sorted(x for x in row[f"{seg}_{part}"] if lo <= _pos(x) <= hi)
            if part.startswith("mutation_count"):
                return len(muts)
            if aa:
                denom = aa_denom(row, seg, t, sp)
            else:
                deleted = sum(len(x.split(":")[2]) for x in row[f"{seg}_deletions"] if lo <= _pos(x) <= hi)
                denom = max(0, uhi - ulo + 1) - deleted
            return len(muts) / denom if denom > 0 else None

        def eq(a, b):
            if isinstance(a, float) and isinstance(b, float):
                return math.isclose(a, b, abs_tol=1e-9)
            return a == b

        def f(t, c):
            if not (ok(t, seg) and same_allele(t, c, seg)):
                return None
            cols = [f"{seg}_{x}" for x in lists]
            g = [f"{seg}_germline_start", f"{seg}_germline_end"]
            if any(t[x] is None or c.get(x) is None for x in cols) or any(t[x] is None or c.get(x) is None for x in g):
                return all(_null_eq(t[x], c.get(x)) for x in cols + g)
            if not span_consistent(t, c, seg):
                return False
            tspan, cspan = mkspan(seg, t, t[g[0]], t[g[1]]), mkspan(seg, t, c[g[0]], c[g[1]])
            shared = mkspan(seg, t, max(t[g[0]], c[g[0]]), min(t[g[1]], c[g[1]]))
            if not derived_part:
                return eq(val(t, shared, t), val(c, shared, t))
            field, reported = f"{seg}_{part}", c.get(f"{seg}_{part}")
            if tspan == cspan:  # identical spans: the reported value is compared with truth directly
                return eq(reported, t[field] if t.get(field) is not None else val(t, tspan, t))
            # differing (valid) spans: the reported value must match the annotator's own list over its own span,
            # and both lists must agree over the shared span
            return eq(reported, val(c, cspan, t)) and eq(val(t, shared, t), val(c, shared, t))
        return f

    fields: dict = {}
    for seg in "vdjc":
        fields[f"{seg}_call"] = ((f"{seg}_call",), call_fn(seg))
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
            if p in ("mutation_frequency", "mutation_frequency_aa"):
                need.add(f"{seg}_deletions")
            need |= {col for col in cuts.get(seg, {}).values()}
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
