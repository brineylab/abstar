import csv
import re
from collections import defaultdict
from pathlib import Path

import pytest

from abstar import _core, germline

LEGACY = Path(__file__).parents[2] / "legacy" / "abstar" / "germline_dbs"
DATABASES = ["balbc", "c57bl6", "human", "human+c57bl6", "macaque"]
LEGACY_SOURCES = {
    "human": ["bcr/human", "tcr/human"],
    "c57bl6": ["bcr/c57bl6"],
    "balbc": ["bcr/balbc"],
    "human+c57bl6": ["bcr/human+c57bl6"],
    "macaque": ["bcr/macaque"],
}
# legacy macaque V alleles abstar germline gap cannot gap (framework insertion, decision 0037)
MACAQUE_DROPOUTS = {
    n: "framework_insertion"
    for n in [
        "IGKV2-39*01", "IGKV2-87*01", "IGKV6-67*01", "IGLV3-18*01",
        "IGLV3-18*03", "IGLV4-17*01", "IGLV4-17*02",
    ]
}
SUFFIX_SPECIES = {"homo_sapiens": "human", "mus_musculus": "mouse"}


def alleles(db):
    with open(germline.info(db).path / "alleles.tsv", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def read_fasta(path):
    out, name = {}, None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].strip()
            out[name] = ""
        elif name is not None:
            out[name] += line.strip()
    return out


def test_builtin_names():
    assert germline.builtin_names() == DATABASES


@pytest.mark.parametrize(
    "db,expected",
    [
        ("human", {"V": 588, "D": 37, "J": 109, "C": 124}),
        ("c57bl6", {"V": 196, "D": 8, "J": 22, "C": 35}),
        ("balbc", {"V": 265, "D": 10, "J": 22, "C": 35}),
        ("human+c57bl6", {"V": 538, "D": 39, "J": 45, "C": 138}),
        ("macaque", {"V": 1064, "D": 52, "J": 28, "C": 46}),
    ],
)
def test_segment_counts(db, expected):
    per_segment = defaultdict(int)
    for key, n in germline.info(db).counts.items():
        per_segment[key.split("/")[1]] += n
    assert dict(per_segment) == expected


@pytest.mark.parametrize("db", DATABASES)
@pytest.mark.parametrize("seg", ["V", "D", "J", "C"])
def test_sequences_match_legacy(db, seg):
    mixed = db == "human+c57bl6"
    legacy = set()
    for src in LEGACY_SOURCES[db]:
        for header, seq in read_fasta(LEGACY / src / "ungapped" / f"{seg.lower()}.fasta").items():
            m = re.fullmatch(r"(.+?)(?:__([a-z_]+))?", header)
            name, suffix = m.group(1), m.group(2)
            legacy.add((SUFFIX_SPECIES[suffix], name, seq) if mixed else (name, seq))
    built = {
        (a["species"], a["name"], a["sequence_gapped"].replace(".", ""))
        if mixed
        else (a["name"], a["sequence_gapped"].replace(".", ""))
        for a in alleles(db)
        if a["segment"] == seg
    }
    assert built - legacy == set()
    missing = {x[0] if not mixed else x[1] for x in legacy - built}
    allowed = set(MACAQUE_DROPOUTS) if db == "macaque" and seg == "V" else set()
    assert missing == allowed


def test_mixed_ids_prefixed_names_plain():
    rows = alleles("human+c57bl6")
    assert all(a["id"].startswith(("hs|", "mm|")) for a in rows)
    assert not any("|" in a["name"] for a in rows)
    by_id = {a["id"]: a for a in rows}
    assert by_id["hs|IGHV1-2*02"]["name"] == "IGHV1-2*02"


def test_human_tr_j_anchor_warnings():
    warnings = _core.germline_summary(str(germline.info("human").path))["warnings"]
    messages = {w["allele"]: w["message"] for w in warnings}
    for allele in ["TRAJ16*01", "TRAJ16*02", "TRAJ35*01", "TRBJ2-7*02", "TRGJP1*01", "TRGJP2*01"]:
        assert "no J Trp/Phe118 motif" in messages[allele]


def test_human_v_anchor():
    by_id = {a["id"]: a for a in alleles("human")}
    assert by_id["IGHV3-23*01"]["anchor"] == "285"


@pytest.mark.parametrize("db", DATABASES)
def test_builtins_load(db):
    s = _core.germline_summary(str(germline.info(db).path))
    assert s["name"] == db
    assert s["n_alleles"] == sum(s["counts"].values()) > 0


def test_macaque_v_anchored():
    def has_orf(seq):
        seq = seq.replace(".", "").upper()
        return any(
            not any(seq[i : i + 3] in ("TAA", "TAG", "TGA") for i in range(f, len(seq) - 2, 3))
            for f in range(3)
        )

    v = [a for a in alleles("macaque") if a["segment"] == "V" and has_orf(a["sequence_gapped"])]
    assert v
    missing = [a["name"] for a in v if not a["anchor"]]
    assert 1 - len(missing) / len(v) >= 0.99, f"unanchored ORF alleles: {missing}"


def test_macaque_gap_template_checksums_current():
    prov = germline.info("macaque").provenance
    seen = 0
    for p in prov:
        for f, sha in zip(p["files"], p["sha256"]):
            if f.startswith("gap-template:"):
                seen += 1
                assert sha == germline.info(f.split(":", 1)[1]).alleles_sha256, f"{f}: template rebuilt, macaque not"
    assert seen == 3


def test_macaque_never_a_default_template():
    assert "macaque" not in _core.GAP_DEFAULT_TEMPLATES
