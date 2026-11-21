"""Accuracy gate: generate fixed-seed sets, annotate, grade against truth, hold each cell to a floor.

`ABSTAR_WRITE_FLOORS=1` rewrites data/accuracy_floors.toml as measured - 0.5 pp (capped at the 99.5% target).
"""
import os
import tomllib
from pathlib import Path

import airr
import pytest

import abstar

N, TARGET = 2000, 0.995
FLOORS = Path(__file__).parent / "data" / "accuracy_floors.toml"
SETS = {  # name: (locus, generate options)
    "igh_naive": ("IGH", {"preset": "naive"}),
    "igh_memory": ("IGH", {"preset": "memory"}),
    "igh_bnab": ("IGH", {"preset": "bnab", "productive_only": True}),
    "igk": ("IGK", {"preset": "naive"}),
    "igl": ("IGL", {"preset": "naive"}),
    "tra": ("TRA", {"preset": "tcr"}),
    "trb": ("TRB", {"preset": "tcr"}),
}


def _measure(name):
    locus, opts = SETS[name]
    truth = abstar.generate(N, locus, seed=20261121, **opts)
    calls = abstar.run(list(zip(truth["sequence_id"], truth["sequence"])))
    g = abstar.grade(truth, calls).filter(_not_null())
    return {f"{r['field']}/{r['subset']}": (r["accuracy"], r["n"]) for r in g.iter_rows(named=True)}


def _not_null():
    import polars as pl

    return pl.col("accuracy").is_not_null()


@pytest.fixture(scope="module")
def measured():
    return {n: _measure(n) for n in SETS}


def test_accuracy_floors(measured):
    if os.environ.get("ABSTAR_WRITE_FLOORS"):
        lines = []
        for s, cells in measured.items():
            lines.append(f"[{s}]")
            lines += [f'"{k}" = {max(0.0, min(TARGET, a - 0.005)):.4f}' for k, (a, _) in cells.items()]
            lines.append("")
        FLOORS.write_text("\n".join(lines))
    floors = tomllib.loads(FLOORS.read_text())
    bad = []
    print(f"\n{'set':<11} {'field/subset':<44} {'n':>5} {'measured':>9} {'target':>7} {'floor':>7}")
    for s, cells in measured.items():
        for k, (a, n) in cells.items():
            fl = floors.get(s, {}).get(k)
            if a < TARGET or os.environ.get("ABSTAR_FULL_TABLE"):
                print(f"{s:<11} {k:<44} {n:>5} {a:>9.4f} {TARGET:>7.3f} {fl or 0:>7.4f}{' <target' if a < TARGET else ''}")
            if fl is None or a < fl:
                bad.append(f"{s} {k}: {a:.4f} < floor {fl}")
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
