"""Accuracy gate: generate fixed-seed sets, annotate, grade against truth, hold each cell to a floor.

`ABSTAR_WRITE_FLOORS=1` rewrites data/accuracy_floors.toml (nucleotide) and data/accuracy_floors_aa.toml
(amino acid: truth `sequence_aa` of in-frame reads, `sequence_type="aa"`) as measured - 0.5 pp (capped at
the 99.5% target). Cells graded on fewer than `MIN_N` reads are reported but carry no floor.
"""
import os
import random
import tomllib
from pathlib import Path

import airr
import pytest

import abstar

N, TARGET, MIN_N = 2000, 0.995, 5
FLOORS = Path(__file__).parent / "data" / "accuracy_floors.toml"
FLOORS_AA = Path(__file__).parent / "data" / "accuracy_floors_aa.toml"
SETS = {  # name: (locus, generate options)
    "igh_naive": ("IGH", {"preset": "naive"}),
    "igh_memory": ("IGH", {"preset": "memory"}),
    "igh_bnab": ("IGH", {"preset": "bnab", "productive_only": True}),
    "igk": ("IGK", {"preset": "naive"}),
    "igl": ("IGL", {"preset": "naive"}),
    "tra": ("TRA", {"preset": "tcr"}),
    "trb": ("TRB", {"preset": "tcr"}),
    # with 60 nt of constant region, so C calls and C mutation fields are graded
    "igh_memory_c": ("IGH", {"preset": "memory", "constant": 60}),
    "igh_bnab_c": ("IGH", {"preset": "bnab", "productive_only": True, "constant": 60}),
    "igk_c": ("IGK", {"preset": "naive", "constant": 60}),
}


def _measure(truth):
    calls = abstar.run(list(zip(truth["sequence_id"], truth["sequence"])))
    g = abstar.grade(truth, calls).filter(_not_null())
    return {f"{r['field']}/{r['subset']}": (r["accuracy"], r["n"]) for r in g.iter_rows(named=True)}


def _not_null():
    import polars as pl

    return pl.col("accuracy").is_not_null()


@pytest.fixture(scope="module")
def truths():
    return {n: abstar.generate(N, locus, seed=20261121, **opts) for n, (locus, opts) in SETS.items()}


@pytest.fixture(scope="module")
def measured(truths):
    return {n: _measure(t) for n, t in truths.items()}


def test_accuracy_floors(measured):
    _check_floors(measured, FLOORS)


V_REGIONS = ("fwr1", "cdr1", "fwr2", "cdr2", "fwr3")


def _imgt_1_100(muts):
    return sorted(m for m in muts or [] if 1 <= int(m[4 : m.index(":")]) <= 100)


def _grade_aa(truth, calls):
    """{field: (accuracy, n)}; region fields only on reads whose V (FWR1-3), J (FWR4) or both (CDR3) are right."""
    hits = {}
    rows = truth.join(calls, on="sequence_id", how="left", suffix="_c")
    for r in rows.iter_rows(named=True):
        ties = {s: r[f"{s}_call_ties"] or [] for s in "vj"}
        ok = {s: r[f"{s}_call"] in ties[s] for s in "vj"}
        cells = {f"{s}_call": ok[s] for s in "vj"}
        cells |= {f"{s}_gene": r[f"{s}_gene"] in {t.split("*")[0] for t in ties[s]} for s in "vj"}
        cells |= {f: r[f] == r[f"{f}_c"] for f in ("locus", "productive", "stop_codon")}
        gated = {**{reg: ok["v"] for reg in V_REGIONS}, "fwr4": ok["j"], "cdr3": ok["v"] and ok["j"]}
        cells |= {f"{reg}{x}": r[f"{reg}{x}"] == r[f"{reg}{x}_c"] for reg, g in gated.items() if g for x in ("_aa", "_length")}
        if ok["v"] and ok["j"]:
            cells["junction_aa"] = r["junction_aa"] == r["junction_aa_c"]
        if ok["v"]:
            cells["v_mutations_aa"] = _imgt_1_100(r["v_mutations_aa"]) == _imgt_1_100(r["v_mutations_aa_c"])
        for f, hit in cells.items():
            hits.setdefault(f, []).append(hit)
    return {f: (sum(h) / len(h), len(h)) for f, h in hits.items()}


def _run_aa(truth):
    t = truth.filter(truth["vj_in_frame"])
    return t, abstar.run(list(zip(t["sequence_id"], t["sequence_aa"])), sequence_type="aa")


def test_accuracy_floors_aa(truths):
    _check_floors({n: _grade_aa(*_run_aa(t)) for n, t in truths.items()}, FLOORS_AA)


def test_aa_negative_controls(truths):
    ubiquitin = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
    gfp = (
        "MSKGEELFTGVVPVLVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFF"
        "KDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLS"
        "TQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
    )
    lysozyme = (
        "KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDG"
        "NGMNAWVAWRNRCKGTDVQAWIRGCRL"
    )
    rng = random.Random(20261121)
    rand = ["".join(rng.choices("ACDEFGHIKLMNPQRSTVWY", k=120)) for _ in range(20)]
    seqs = [ubiquitin, gfp, truths["igh_memory"]["sequence"][0], lysozyme, *rand]
    calls = abstar.run([(f"neg{i}", s) for i, s in enumerate(seqs)], sequence_type="aa")
    assert calls["failure_reason"].to_list() == ["no_v_hit"] * len(seqs), calls.select("sequence_id", "v_call", "v_score", "failure_reason")


def _check_floors(measured, path):
    if os.environ.get("ABSTAR_WRITE_FLOORS"):
        lines = []
        for s, cells in measured.items():
            lines.append(f"[{s}]")
            lines += [f'"{k}" = {max(0.0, min(TARGET, a - 0.005)):.4f}' for k, (a, n) in cells.items() if n >= MIN_N]
            lines.append("")
        path.write_text("\n".join(lines))
    floors = tomllib.loads(path.read_text())
    bad = []
    print(f"\n{'set':<11} {'field/subset':<44} {'n':>5} {'measured':>9} {'target':>7} {'floor':>7}")
    for s, cells in measured.items():
        for k, (a, n) in cells.items():
            fl = floors.get(s, {}).get(k)
            if a < TARGET or os.environ.get("ABSTAR_FULL_TABLE"):
                print(f"{s:<11} {k:<44} {n:>5} {a:>9.4f} {TARGET:>7.3f} {fl or 0:>7.4f}{' <target' if a < TARGET else ''}")
            if n >= MIN_N and (fl is None or a < fl):
                bad.append(f"{s} {k}: {a:.4f} < floor {fl}")
    # every floor still has its cell, and every cell graded on MIN_N reads or more has a floor
    want = {(s, k) for s, cells in measured.items() for k, (_, n) in cells.items() if n >= MIN_N}
    have = {(s, k) for s, cells in floors.items() for k in cells}
    bad += [f"{s} {k}: floor without a measured cell" for s, k in sorted(have - want)]
    n_below = sum(a < TARGET for c in measured.values() for a, _ in c.values())
    print(f"{n_below} of {sum(len(c) for c in measured.values())} cells below the {TARGET} target")
    assert not bad, "\n".join(bad)


def test_airr_validates_2000_reads(tmp_path):
    truth = abstar.generate(N, "IGH", seed=20261121, preset="memory")
    fa = tmp_path / "reads.fasta"
    fa.write_text("".join(f">{i}\n{s}\n" for i, s in zip(truth["sequence_id"], truth["sequence"])))
    out = tmp_path / "out"
    abstar.run(fa, out)
    assert airr.validate_rearrangement(str(out / "reads.tsv")) is True
    assert sum(1 for _ in airr.read_rearrangement(str(out / "reads.tsv"))) == N
