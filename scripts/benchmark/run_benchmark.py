"""Benchmark new abstar against IgBLAST and legacy abstar (generator truth) and legacy's real-data baseline.

Not a CI dependency. Needs the project venv (abstar.generate/run/grade) and the abstar-bench env
(scripts/benchmark/setup_env.sh) plus IgBLAST databases (scripts/benchmark/make_igblast_db.py).

    .venv/bin/python scripts/benchmark/run_benchmark.py --work WORK [--igblast-db WORK/igblast_db]

Raw outputs go to WORK (reused when present; delete a tool's directory to rerun it). Summary tables go
to --out (default docs/superpowers/benchmarks/).
"""
import argparse
import gzip
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

import polars as pl

import abstar

ROOT = Path(__file__).resolve().parents[2]
N = int(os.environ.get("BENCH_N", 5000))  # smaller N for a smoke run only
THREADS = os.cpu_count()
TARGET = 0.995
# name: (locus, seed, generate CLI options)
SETS = {
    "igh_naive": ("IGH", 20261201, ["--preset", "naive"]),
    "igh_memory": ("IGH", 20261202, ["--preset", "memory"]),
    "igh_bnab": ("IGH", 20261203, ["--preset", "bnab"]),
    "igh_indels": ("IGH", 20261204, ["--preset", "memory", "--indels", "0.5"]),
    "igh_truncated": ("IGH", 20261205, ["--preset", "memory", "--truncate-5p", "uniform:0-60", "--truncate-3p", "uniform:0-40"]),
    "igk": ("IGK", 20261206, ["--preset", "naive"]),
    "igl": ("IGL", 20261207, ["--preset", "naive"]),
    "tra": ("TRA", 20261208, ["--preset", "tcr"]),
    "trb": ("TRB", 20261209, ["--preset", "tcr"]),
}
TOOLS = ("new", "igblast", "legacy")
PRIMARY = ("v_call", "d_call", "j_call", "v_sequence_end", "d_sequence_start", "d_sequence_end", "j_sequence_start",
           "junction", "cdr3", "productive")
CALLS = ("v_call", "d_call", "j_call", "c_call")


def sh(cmd, log, env=None):
    t0 = time.perf_counter()
    with open(log, "w") as f:
        subprocess.run(cmd, check=True, stdout=f, stderr=subprocess.STDOUT, env=env)
    return time.perf_counter() - t0


class Runner:
    def __init__(self, a):
        self.venv = ROOT / ".venv" / "bin"
        self.mm = [a.micromamba, "run", "-n", a.env]
        self.env_prefix = Path(subprocess.run(self.mm + ["sh", "-c", "echo $CONDA_PREFIX"], capture_output=True, text=True,
                                              check=True).stdout.strip())
        self.db = Path(a.igblast_db)

    def new(self, fasta, out):
        return sh([self.venv / "abstar", "run", fasta, out, "--format", "parquet", "--threads", str(THREADS), "-q"],
                  out.with_suffix(".log")), out / f"{fasta.stem}.parquet"

    def igblast(self, fasta, out, tcr):
        out.mkdir(parents=True, exist_ok=True)
        kind = "tr" if tcr else "ig"
        env = {**os.environ, "IGDATA": str(self.env_prefix / "share" / "igblast")}
        tsv = out / "igblast.tsv"
        cmd = [self.env_prefix / "bin" / "igblastn", "-query", fasta, "-organism", "human", "-domain_system", "imgt",
               "-ig_seqtype", "TCR" if tcr else "Ig", "-outfmt", "19", "-num_threads", str(THREADS), "-out", tsv,
               "-custom_internal_data", self.db / "abstar.ndm", "-auxiliary_data", self.db / "abstar.aux"]
        for seg in "VDJ":
            cmd += [f"-germline_db_{seg}", self.db / f"{kind}_{seg}"]
        return sh(cmd, out / "run.log", env), tsv

    def legacy(self, fasta, out, tcr):
        cmd = self.mm + ["abstar", "run", fasta, out, "-o", "parquet", "--receptor", "tcr" if tcr else "bcr",
                         "--mmseqs_threads", str(THREADS), "--n_processes", str(THREADS)]
        return sh(cmd, out.with_suffix(".log")), out / "parquet" / f"{fasta.stem}.parquet"


def cached(work, key, fn):
    """Run fn() -> (seconds, path) once; the timing is kept in WORK/timings.json."""
    tf = work / "timings.json"
    times = json.loads(tf.read_text()) if tf.exists() else {}
    if key not in times:
        sec, path = fn()
        times[key] = {"seconds": sec, "path": str(path)}
        tf.write_text(json.dumps(times, indent=1))
    return times[key]["seconds"], Path(times[key]["path"])


# legacy list formats -> abstar's: 70:A>C, 77-85:9>SEQ (deletion), 118:9>SEQ (insertion)
def _legacy_item(col, x):
    pos = re.match(r"\d+", x).group(0)
    if col.endswith("_deletions"):
        return f"IMGT{pos}:del:{x.split('>', 1)[1]}"
    if col.endswith("_insertions"):
        return f"IMGT{pos}:ins:{x.split('>', 1)[1]}"
    return f"IMGT{x}"


def load(tool, path, truth):
    """Tool output -> frame with truth's column names and dtypes; calls reduced to their first allele."""
    df = pl.read_csv(path, separator="\t", infer_schema_length=0, quote_char=None) if path.suffix == ".tsv" else pl.read_parquet(path)
    cols = [pl.col("sequence_id").cast(pl.String)]
    for c in df.columns:
        if c not in truth.columns or c in ("sequence_id", "sequence"):
            continue
        dt, s = truth.schema[c], pl.col(c)
        if df.schema[c] == pl.String:
            s = pl.when(s == "").then(None).otherwise(s)
        if c in CALLS:
            s = s.str.split(",").list.first()
        elif dt == pl.List(pl.String) and df.schema[c] == pl.String:  # legacy '|'-joined lists
            seg_call = pl.col(f"{c[0]}_call")
            s = (pl.when(s.is_null() & seg_call.is_not_null()).then(pl.lit([], pl.List(pl.String)))
                 .otherwise(s.str.split("|").list.eval(pl.element().map_elements(lambda x, c=c: _legacy_item(c, x), return_dtype=pl.String))))
        elif dt == pl.Boolean and df.schema[c] == pl.String:
            s = s.replace_strict({"T": True, "F": False, "true": True, "false": False, "True": True, "False": False}, default=None,
                                 return_dtype=pl.Boolean)
        elif dt != df.schema[c]:
            s = s.cast(dt, strict=False)
        cols.append(s.alias(c))
    return df.select(cols)


def fails(truth, calls):
    """Per-read failures on the primary fields (the grader's rules), for examples. Counts are checked against grade()."""
    crows = {r["sequence_id"]: r for r in calls.iter_rows(named=True)}
    out = []
    for t in truth.iter_rows(named=True):
        c = crows.get(t["sequence_id"], {})
        for f in PRIMARY:
            if f not in calls.columns:
                continue
            seg = f[0] if f[1] == "_" else None
            if seg:  # call or boundary
                if t.get(f"unrecoverable_{seg}"):
                    continue
                valid = t.get(f"{f}_valid")
                v = c.get(f)
                ok = (v is None and seg == "d" and t["d_call_null_valid"]) or (v is None if not valid else v in valid)
            elif f in ("junction", "productive"):
                if t["unrecoverable_v"] or t["unrecoverable_j"]:
                    continue
                ok = t[f] == c.get(f)
            else:  # cdr3: graded only when both V and J alleles agree
                if any(t[f"unrecoverable_{s}"] or t[f"{s}_call"] is None or c.get(f"{s}_call") != t[f"{s}_call"] for s in "vj"):
                    continue
                ok = t[f] == c.get(f)
            if not ok:
                out.append({"sequence_id": t["sequence_id"], "field": f, "truth": str(t.get(f)),
                            "valid": str(t.get(f"{f}_valid")), "call": str(c.get(f)),
                            "flags": ",".join(k for k in ("high_mutation", "long_cdr3", "has_indel", "truncated_5p", "truncated_3p",
                                                          "ambiguous_v", "ambiguous_d", "ambiguous_j") if t.get(k))})
    return pl.DataFrame(out, schema={k: pl.String for k in ("sequence_id", "field", "truth", "valid", "call", "flags")})


def db_names(db="human"):
    import csv

    with open(abstar.germline.resolve(db)[0] / "alleles.tsv", newline="") as f:
        return {r["name"] for r in csv.DictReader(f, delimiter="\t")}


def generated(a, run, work):
    names, seqs = db_names(), allele_seqs()
    acc, timing, naming, short = [], [], [], []
    for name, (locus, seed, opts) in SETS.items():
        d = work / name
        gen = d / "gen"
        if not (gen / "generated.truth.parquet").exists():
            subprocess.run([run.venv / "abstar", "generate", "--n", str(N), "--locus", locus, "--seed", str(seed), *opts,
                            "--format", "parquet", "-q", gen], check=True)
        truth = pl.read_parquet(gen / "generated.truth.parquet")
        fasta = gen / "generated.fasta"
        tcr = locus.startswith("TR")
        outs = {
            "new": cached(work, f"{name}/new", lambda: run.new(fasta, d / "new")),
            "igblast": cached(work, f"{name}/igblast", lambda: run.igblast(fasta, d / "igblast", tcr)),
            "legacy": cached(work, f"{name}/legacy", lambda: run.legacy(fasta, d / "legacy", tcr)),
        }
        for tool, (sec, path) in outs.items():
            calls = load(tool, path, truth)
            g = abstar.grade(truth, calls)
            live = g.filter(pl.col("n") > 0)
            acc.append(live.with_columns(set=pl.lit(name), tool=pl.lit(tool)))
            fl = fails(truth, calls)
            fl.write_csv(d / f"{tool}.fails.tsv", separator="\t")
            for f in PRIMARY:  # the example extractor must agree with the grader
                row = g.filter((pl.col("field") == f) & (pl.col("subset") == "all"))
                if row["n"][0]:
                    assert row["n"][0] - row["correct"][0] == fl.filter(pl.col("field") == f).height, (name, tool, f)
            timing.append({"set": name, "tool": tool, "reads": N, "returned": calls.height, "seconds": round(sec, 2),
                           "reads_per_s": round(N / sec, 1), "threads": THREADS})
            if tool == "legacy":
                naming += legacy_naming(name, truth, path, names)
            if tool == "new":
                short.append(shortfalls(name, truth, path, seqs))
        print(f"{name}: done", file=sys.stderr)
    return (pl.concat(acc).select("set", "tool", "field", "subset", "n", "correct", "accuracy"), pl.DataFrame(timing),
            pl.DataFrame(naming), pl.concat(short).sort("set", "field", "count", descending=[False, False, True]))


def shortfalls(name, truth, path, seqs):
    """Classify each new-abstar miss on a primary field (see the report's shortfall section)."""
    calls = pl.read_parquet(path)
    crow = {r["sequence_id"]: r for r in calls.iter_rows(named=True)}
    trow = {r["sequence_id"]: r for r in truth.iter_rows(named=True)}
    fl = pl.read_csv(Path(path).parents[1] / "new.fails.tsv", separator="\t", infer_schema_length=0)
    failed = {}
    for r in fl.iter_rows(named=True):
        failed.setdefault(r["sequence_id"], set()).add(r["field"])

    def score(region, allele, both=False):
        g = seqs.get(allele or "")
        if not g or not region:
            return None
        return max(best_local(region, g), best_local(region, g[::-1].translate(_RC))) if both else best_local(region, g)

    out = []
    for sid, fields in failed.items():
        t, c = trow[sid], crow.get(sid, {})
        read = t["sequence"]
        for f in fields:
            seg = f[0] if f[1] == "_" else None
            if c.get("failure_reason") == "no_j_hit" and f not in ("v_call", "v_sequence_end"):
                jl = (t["j_sequence_end"] or 0) - (t["j_sequence_start"] or 1) + 1
                cls = ("no J hit: J below the score floor (< 20 J nt left)" if jl < 20 else
                       "no J hit: heavy mutation" if t["high_mutation"] else "no J hit: other")
            elif seg == "d":
                ve = [x for x in (t["v_sequence_end"], c.get("v_sequence_end")) if x]
                js = [x for x in (t["j_sequence_start"], c.get("j_sequence_start")) if x]
                region = read[min(ve) if ve else 0:(max(js) - 1) if js else len(read)]
                st, sc = score(region, t["d_call"], True), score(region, c.get("d_call"), True)
                if f != "d_call" and "d_call" in fields:
                    cls = "D boundary: follows the D call miss"
                elif f != "d_call":
                    cls = "D boundary: right D, cut outside _valid"
                elif c.get("d_call") is None:
                    cls = f"D not called (true D best score {'<=6' if (st or 0) <= 6 else '7-9' if st <= 9 else '>=10'})"
                elif not t["d_call_valid"]:
                    cls = "D called, but the read holds no D base (only null is correct, A11): chance placement over N"
                else:
                    cls = ("D call: read scores the call above the true D (chance placement, A26)" if sc > st else
                           "D call: ties the true D (read can't decide)" if sc == st else
                           "D call: true D scores higher (annotator D error)")
                if t["high_mutation"]:
                    cls += " [high_mutation]"
            elif f in ("v_call", "j_call"):
                st, sc = score(read, t[f]), score(read, c.get(f))
                cls = (f"{f}: null call" if c.get(f) is None else
                       f"{f}: read scores the call above the true allele" if sc > st else
                       f"{f}: ties the true allele (a tie the _valid set misses)" if sc == st else
                       f"{f}: true allele scores higher (annotator error)")
                if t["high_mutation"]:
                    cls += " [high_mutation]"
            elif f in ("v_sequence_end", "j_sequence_start"):
                cls = f"{f}: follows the {seg} call miss" if f"{seg}_call" in fields else f"{f}: right allele, cut outside _valid"
            else:  # junction, cdr3, productive
                if fields & {"v_call", "j_call"}:
                    cls = "follows a V/J call miss"
                elif fields & {"v_sequence_end", "j_sequence_start"}:
                    cls = "follows a V-end/J-start cut outside _valid"
                elif c.get("junction") is None and t["junction"] is not None:
                    vmax = max(t["v_sequence_end_valid"] or [0])
                    cls = ("annotator junction null: Cys104 beyond every valid V end (truth places it from origin)"
                           if t["cdr3_start"] and vmax < t["cdr3_start"] - 3 else "annotator junction null: other")
                elif t["junction"] is None:
                    cls = "truth junction null, annotator reports one"
                else:
                    cls = "junction differs: other"
            out.append({"set": name, "field": f, "class": cls, "sequence_id": sid})
    df = pl.DataFrame(out, schema={"set": pl.String, "field": pl.String, "class": pl.String, "sequence_id": pl.String})
    return df.group_by("set", "field", "class").agg(pl.len().alias("count"), pl.col("sequence_id").sort().head(3).str.join(" ").alias("examples"))


def legacy_naming(name, truth, path, names):
    """How much of legacy's call error is naming: calls absent from abstar's DB, and accuracy with '__species' stripped."""
    raw = pl.read_parquet(path).select("sequence_id", *[pl.col(f"{s}_call").str.split(",").list.first().alias(f"{s}_raw") for s in "vdj"])
    m = truth.join(raw, on="sequence_id", how="left")
    rows = []
    for s in "vdj":
        valid = f"{s}_call_valid"
        sub = m.filter(~pl.col(f"unrecoverable_{s}"))
        n_called = sub.filter(pl.col(f"{s}_raw").is_not_null()).height
        not_in_db = sub.filter(pl.col(f"{s}_raw").is_not_null() & ~pl.col(f"{s}_raw").is_in(list(names))).height

        stripped = sub.with_columns(pl.col(f"{s}_raw").str.replace(r"__.*$", ""))
        wrong_not_in_db = sum(
            (x is not None and x not in names) and not (vl and x in vl)
            for x, vl in zip(sub[f"{s}_raw"].to_list(), sub[valid].to_list()))
        rows.append({"set": name, "segment": s, "n": sub.height, "called": n_called, "not_in_abstar_db": not_in_db,
                     "wrong_and_not_in_db": wrong_not_in_db, "accuracy_as_is": round(_acc(sub, s), 4),
                     "accuracy_suffix_stripped": round(_acc(stripped, s), 4)})
    return rows


def _acc(sub, s):
    """Primary-call accuracy under the grader's call rules."""
    v = sub[f"{s}_raw"].to_list()
    ok = sum((x is None and s == "d" and nv) or (x is None if not vl else x in vl)
             for x, vl, nv in zip(v, sub[f"{s}_call_valid"].to_list(), sub["d_call_null_valid"].to_list()))
    return ok / sub.height


def real_data(a, run, work):
    """New abstar (and IgBLAST) gene calls vs legacy's committed corpus baseline, joined on corpus_ordinal."""
    corpus = ROOT / "legacy" / "test_data" / "bcr_corpus"
    d = work / "corpus"
    d.mkdir(parents=True, exist_ok=True)
    fasta = d / "corpus.fasta"
    if not fasta.exists():
        with gzip.open(corpus / "sequences.fasta.gz", "rt") as f:
            seqs = [l.strip() for l in f if not l.startswith(">")]
        fasta.write_text("".join(f">c{i}\n{s}\n" for i, s in enumerate(seqs)))
    n = sum(1 for l in fasta.open() if l.startswith(">"))
    outs = {
        "new": cached(work, "corpus/new", lambda: run.new(fasta, d / "new")),
        "igblast": cached(work, "corpus/igblast", lambda: run.igblast(fasta, d / "igblast", False)),
        "legacy": cached(work, "corpus/legacy", lambda: run.legacy(fasta, d / "legacy", False)),
    }
    timing = [{"set": "corpus_real_bcr", "tool": t, "reads": n, "seconds": round(s, 2), "reads_per_s": round(n / s, 1),
               "threads": THREADS} for t, (s, _) in outs.items()]
    cols = ("locus", "productive", "junction_aa", "v_call", "d_call", "j_call", "v_sequence_end", "j_sequence_start")
    base = pl.read_parquet(corpus / "baseline.parquet").select("corpus_ordinal", *cols)

    def tool_frame(tool, path):
        df = pl.read_csv(path, separator="\t", infer_schema_length=0, quote_char=None) if path.suffix == ".tsv" else pl.read_parquet(path)
        keep = [c for c in (*cols, "d_inverted", "v_call_ties", "d_call_ties", "j_call_ties") if c in df.columns]
        df = df.select(pl.col("sequence_id").str.slice(1).cast(pl.Int64).alias("corpus_ordinal"), *keep)
        df = df.with_columns([pl.when(pl.col(c) == "").then(None).otherwise(pl.col(c)).alias(c) for c in keep if df.schema[c] == pl.String])
        return df.with_columns([pl.col(c).cast(pl.Int64) for c in ("v_sequence_end", "j_sequence_start") if c in keep])

    frames = {t: tool_frame(t, p) for t, (_, p) in outs.items()}
    names, seqs = db_names(), allele_seqs()
    reads = [l.strip() for l in fasta.open() if not l.startswith(">")]
    rows, examples = [], []
    for tool, ref in (("new", "baseline"), ("igblast", "baseline"), ("new", "igblast"), ("legacy", "baseline")):
        left = frames[tool]
        right = base if ref == "baseline" else frames[ref].drop([c for c in frames[ref].columns if c.endswith("_ties")])
        m = left.join(right, on="corpus_ordinal", how="full", suffix="_ref", coalesce=True)
        evidence = tool == "new" and ref == "baseline"
        for s in "vdj":
            ties = f"{s}_call_ties" if f"{s}_call_ties" in m.columns else None
            for r in m.iter_rows(named=True):
                cls = classify(r.get(f"{s}_call"), r.get(f"{s}_call_ref"), r.get(ties) if ties else None, names)
                row = {"comparison": f"{tool} vs {ref}", "segment": s, "class": cls, "read_evidence": None}
                if evidence and not cls.startswith(("same allele", "both null")):
                    ev = read_evidence(s, r, reads[r["corpus_ordinal"]], seqs)
                    row["read_evidence"] = ev.pop("verdict")
                    examples.append({"segment": s, "class": cls, "read_evidence": row["read_evidence"], "corpus_ordinal": r["corpus_ordinal"],
                                     "locus": r.get("locus") or r.get("locus_ref"), "new": r.get(f"{s}_call"),
                                     "new_ties": ",".join(r.get(ties) or []) if ties else None, "baseline": r.get(f"{s}_call_ref"), **ev})
                rows.append(row)
        if evidence:
            for r in m.iter_rows(named=True):
                rows.append({"comparison": "new vs baseline", "segment": "junction_aa", "class": junction_class(r), "read_evidence": None})
                rows.append({"comparison": "new vs baseline", "segment": "locus", "read_evidence": None,
                             "class": "agree" if r["locus"] == r["locus_ref"] else f"{r['locus_ref']}->{r['locus']}"})
    pl.DataFrame(examples).write_csv(d / "new_vs_baseline.disagreements.tsv", separator="\t")
    summary = (pl.DataFrame(rows).group_by("comparison", "segment", "class", "read_evidence").len("count")
               .sort("comparison", "segment", "count", descending=[False, False, True]))
    return summary, pl.DataFrame(timing)


def junction_class(r):
    a, b = r["junction_aa"], r["junction_aa_ref"]
    if a == b:
        return "agree" if a else "both null"
    if a is None or b is None:
        return "only baseline has a junction" if a is None else "only new has a junction"
    gene = lambda x: re.sub(r"__.*$", "", x or "").split(",")[0].split("*")[0]  # noqa: E731
    if any(gene(r[f"{s}_call"]) != gene(r[f"{s}_call_ref"]) for s in "vj"):
        return "differ: V or J gene differs"
    return "differ: length" if len(a) != len(b) else "differ: same length"


def allele_seqs(db="human"):
    import csv

    with open(abstar.germline.resolve(db)[0] / "alleles.tsv", newline="") as f:
        return {r["name"]: r["sequence_gapped"].replace(".", "").upper() for r in csv.DictReader(f, delimiter="\t")}


def best_local(read, allele):
    """Best gapless local score (+1 match, -1 mismatch) of allele on read, either orientation for D."""
    import numpy as np

    r, a = np.frombuffer(read.encode(), np.uint8), np.frombuffer(allele.encode(), np.uint8)
    best = 0
    for off in range(-len(a) + 1, len(r)):
        lo, hi = max(0, off), min(len(r), off + len(a))
        if hi - lo <= best:
            continue
        c = np.concatenate(([0], np.where(r[lo:hi] == a[lo - off:hi - off], 1, -1).cumsum()))
        best = max(best, int((c - np.minimum.accumulate(c)).max()))
    return best


_RC = str.maketrans("ACGTN", "TGCAN")


def read_evidence(s, r, read, seqs):
    """Score both tools' alleles against the read with one rule: V and J over the whole read, D (both orientations)
    over the junction span both tools leave between V and J. A gapless score is approximate on reads with an indel."""
    region = read
    if s == "d":
        ve = [x for x in (r.get("v_sequence_end"), r.get("v_sequence_end_ref")) if x]
        js = [x for x in (r.get("j_sequence_start"), r.get("j_sequence_start_ref")) if x]
        region = read[min(ve) if ve else 0:(max(js) - 1) if js else len(read)]

    def score(call):
        name = re.sub(r"__.*$", "", call.split(",")[0]) if call else None
        if name is None or name not in seqs:
            return None
        g = seqs[name]
        return max(best_local(region, g), best_local(region, g[::-1].translate(_RC))) if s == "d" else best_local(region, g)

    sn, sb = score(r.get(f"{s}_call")), score(r.get(f"{s}_call_ref"))
    if sn is None or sb is None:
        x = sn if sb is None else sb
        verdict = "one null" if s != "d" or x is None else f"one null (called D scores {'<=6' if x <= 6 else '7-9' if x <= 9 else '>=10'})"
    else:
        verdict = "read favors new" if sn > sb else "tie (read can't decide)" if sn == sb else "read favors baseline"
    return {"verdict": verdict, "score_new": sn, "score_baseline": sb}


def classify(call, ref, ties, names):
    """Agreement class of one primary call against the reference (first allele of the reference's comma list)."""
    tagged = "__" in (ref or "").split(",")[0] or "__" in (call or "").split(",")[0]
    strip = lambda x: re.sub(r"__.*$", "", x)  # noqa: E731  legacy tags some alleles with __homo_sapiens
    refs = [strip(x) for x in (ref or "").split(",") if x]
    call = strip(call.split(",")[0]) if call else None
    if call is None and not refs:
        return "both null"
    if call is None:
        return "only reference called"
    if not refs:
        return "only tool called"
    if call == refs[0]:
        return "same allele (suffix stripped)" if tagged else "same allele"
    tie = set(ties or []) | {call}
    if call in refs or refs[0] in tie:
        return "tie overlap"
    gene = lambda x: x.split("*")[0]  # noqa: E731
    if gene(call) == gene(refs[0]):
        return "same gene, other allele" + ("" if refs[0] in names else " (reference allele not in abstar db)")
    return "different gene" + ("" if refs[0] in names else " (reference allele not in abstar db)")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--work", required=True, type=Path)
    p.add_argument("--igblast-db", default=None)
    p.add_argument("--out", type=Path, default=ROOT / "docs" / "superpowers" / "benchmarks")
    p.add_argument("--micromamba", default=str(Path.home() / ".local" / "bin" / "micromamba"))
    p.add_argument("--env", default="abstar-bench")
    a = p.parse_args()
    a.work = a.work.resolve()
    a.igblast_db = a.igblast_db or str(a.work / "igblast_db")
    a.work.mkdir(parents=True, exist_ok=True)
    a.out.mkdir(parents=True, exist_ok=True)
    run = Runner(a)
    acc, timing, naming, short = generated(a, run, a.work)
    real, rtiming = real_data(a, run, a.work)
    stem = a.out / "2026-11-20-annotator-v1"
    acc.sort("set", "tool", "field", "subset").write_csv(f"{stem}-accuracy.tsv", separator="\t", float_precision=4)
    pl.concat([timing, rtiming.with_columns(returned=pl.lit(None, pl.Int64)).select(timing.columns)]).write_csv(f"{stem}-speed.tsv", separator="\t")
    naming.write_csv(f"{stem}-legacy-naming.tsv", separator="\t")
    short.write_csv(f"{stem}-shortfalls.tsv", separator="\t")
    real.write_csv(f"{stem}-realdata.tsv", separator="\t")
    print(f"wrote {stem}-{{accuracy,speed,legacy-naming,shortfalls,realdata}}.tsv", file=sys.stderr)


if __name__ == "__main__":
    main()
