"""Rebuild the built-in germline databases from germline-sources/recipes.

Usage: python scripts/build_builtin_germlines.py [NAME ...]   (default: all)
"""
import sys
from pathlib import Path

from abstar import germline

ROOT = Path(__file__).resolve().parents[1]
RECIPES = ROOT / "germline-sources" / "recipes"
OUT = ROOT / "python" / "abstar" / "germline_dbs"


def main(names: list[str]) -> None:
    # macaque is not shipped until #137 regaps it to IMGT (decision 0036); the recipe stays for that
    for name in names or sorted(p.stem for p in RECIPES.glob("*.toml") if p.stem != "macaque"):
        r = germline.build(RECIPES / f"{name}.toml", output_dir=OUT / name, overwrite=True)
        print(f"{r.name}: {r.n_alleles} alleles, {len(r.warnings)} warnings -> {r.path}")
        for w in r.warnings:
            print(f"  {w['allele']}: {w['message']}")


if __name__ == "__main__":
    main(sys.argv[1:])
