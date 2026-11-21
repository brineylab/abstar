"""Build IgBLAST databases from abstar's built-in germlines so allele names match the truth sets.

Writes, under OUT (default: igblast_db/ next to the cwd):
  {ig,tr}_{V,D,J}  BLAST nucleotide DBs (makeblastdb -parse_seqids -dbtype nucl)
  abstar.ndm       custom internal data (-custom_internal_data): FWR/CDR ends per V allele from its IMGT gaps
  abstar.aux       auxiliary data (-auxiliary_data): J coding frame and CDR3 end from abstar's J anchors

Run with the project venv (it reads the database through abstar.germline) and makeblastdb on PATH
or in MAKEBLASTDB: .venv/bin/python scripts/benchmark/make_igblast_db.py [OUT] [DB]
"""
import csv
import os
import subprocess
import sys
from pathlib import Path

from abstar import germline

# IMGT nt positions (gapped, 1-based) of the V regions; FWR3 includes Cys104 (310-312).
REGIONS = ((1, 78), (79, 114), (115, 165), (166, 195), (196, 312))
CHAIN = {"IGH": "H", "IGK": "K", "IGL": "L", "TRA": "A", "TRB": "B", "TRG": "G", "TRD": "D"}


def ndm_row(name, locus, gapped):
    """FWR1..FWR3 start/stop on the ungapped allele (1-based, -1 if absent), chain type, coding frame (0-based)."""
    pos = [i + 1 for i, ch in enumerate(gapped) if ch != "."]  # IMGT position of each ungapped base
    cols = []
    for lo, hi in REGIONS:
        idx = [u for u, p in enumerate(pos, 1) if lo <= p <= hi]
        cols += [idx[0], idx[-1]] if idx else [-1, -1]
    frame = (3 - (pos[0] - 1) % 3) % 3  # first base of the first complete IMGT codon
    return [name, *cols, "V" + CHAIN[locus], frame]


def main(out="igblast_db", db="human"):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    with open(germline.resolve(db)[0] / "alleles.tsv", newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    fasta = {}
    ndm, aux = [], []
    for r in rows:
        seg, locus = r["segment"], r["locus"]
        if seg not in "VDJ" or locus not in CHAIN:
            continue
        seq = r["sequence_gapped"].replace(".", "").upper()
        kind = "ig" if locus.startswith("IG") else "tr"
        fasta.setdefault(f"{kind}_{seg}", []).append(f">{r['name']}\n{seq}\n")
        if seg == "V":
            ndm.append(ndm_row(r["name"], locus, r["sequence_gapped"]))
        elif seg == "J" and r["anchor"]:
            a = int(r["anchor"])  # 0-based start of the Trp/Phe118 codon
            aux.append([r["name"], a % 3, "J" + CHAIN[locus], a - 1, (len(seq) - a % 3) % 3])
    mk = os.environ.get("MAKEBLASTDB", "makeblastdb")
    for name, recs in fasta.items():
        fa = out / f"{name}.fasta"
        fa.write_text("".join(recs))
        subprocess.run([mk, "-parse_seqids", "-dbtype", "nucl", "-in", str(fa), "-out", str(out / name)],
                       check=True, stdout=subprocess.DEVNULL)
    (out / "abstar.ndm").write_text("".join("\t".join(map(str, r)) + "\n" for r in ndm))
    (out / "abstar.aux").write_text("".join("\t".join(map(str, r)) + "\n" for r in aux))
    print(f"{out}: " + ", ".join(f"{k} {len(v)}" for k, v in sorted(fasta.items())) + f"; ndm {len(ndm)}, aux {len(aux)}")


if __name__ == "__main__":
    main(*sys.argv[1:])
