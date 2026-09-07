"""Transparent, unmutated human TCR constructions through the public API."""

import hashlib
import json
import re
from dataclasses import dataclass
from pathlib import Path

import abstar
import polars as pl
import pytest
from abutils import Sequence
from Bio.Seq import Seq

from .helpers import read_fasta_records


TCR_DATA_DIR = Path(__file__).parent / "data" / "tcr"
PACKAGE_DIR = Path(abstar.__file__).parent


@dataclass(frozen=True, slots=True)
class TCRCase:
    sequence_id: str
    sequence: str
    locus: str
    allowed_v_calls: tuple[str, ...]
    allowed_d_calls: tuple[str | None, ...]
    allowed_j_calls: tuple[str, ...]
    junction: str
    productive: bool

    def as_sequence(self):
        return Sequence(self.sequence, id=self.sequence_id)


def load_tcr_cases() -> tuple[TCRCase, ...]:
    records = read_fasta_records(TCR_DATA_DIR / "sequences.fasta")
    definitions = json.loads((TCR_DATA_DIR / "cases.json").read_text())
    return tuple(
        TCRCase(
            sequence_id=item["sequence_id"], sequence=records[item["sequence_id"]],
            locus=item["locus"],
            allowed_v_calls=tuple(item["allowed_v_calls"]),
            allowed_d_calls=tuple(item["allowed_d_calls"]),
            allowed_j_calls=tuple(item["allowed_j_calls"]),
            junction=item["junction"], productive=item["productive"],
        )
        for item in definitions["cases"]
    )


@pytest.fixture(autouse=True)
def isolated_user_database(monkeypatch, tmp_path):
    monkeypatch.setenv("HOME", str(tmp_path))


def test_tcr_fixture_provenance_and_exact_construction():
    """Detect stale/fabricated alleles, boundaries, anchors, hashes or goldens."""
    definitions = json.loads((TCR_DATA_DIR / "cases.json").read_text())
    records = read_fasta_records(TCR_DATA_DIR / "sequences.fasta")
    assert definitions["receptor"] == "tcr"
    assert definitions["germline_database"] == "human"
    assert [case["locus"] for case in definitions["cases"]] == ["TRA", "TRB", "TRD", "TRG"]
    assert set(records) == {case["sequence_id"] for case in definitions["cases"]}
    for path, digest in definitions["source_files"].items():
        assert hashlib.sha256((PACKAGE_DIR.parent / path).read_bytes()).hexdigest() == digest
    root = PACKAGE_DIR / "germline_dbs" / "tcr" / "human"
    references = {
        (kind, segment): read_fasta_records(root / kind / f"{segment}.fasta")
        for kind in ("ungapped", "imgt_gapped") for segment in ("v", "d", "j", "c")
    }
    for item in definitions["cases"]:
        sequence = records[item["sequence_id"]]
        assert hashlib.sha256(sequence.encode()).hexdigest() == item["sequence_sha256"]
        assert ("d" in item["source_alleles"]) == (item["locus"] in ("TRB", "TRD"))
        for segment, source in item["source_alleles"].items():
            allele = source["allele"]
            assert allele.startswith(item["locus"])
            for kind in ("ungapped", "imgt_gapped"):
                bases = references[kind, segment][allele]
                assert hashlib.sha256(bases.encode()).hexdigest() == source[kind + "_sha256"]
            assert references["imgt_gapped", segment][allele].replace(".", "") == references["ungapped", segment][allele]
        reconstructed = ""
        for piece in item["construction"]:
            assert piece["query_start"] == len(reconstructed)
            if "allele" in piece:
                assert piece["allele"] == item["source_alleles"][piece["segment"]]["allele"]
                reference = references["ungapped", piece["segment"]][piece["allele"]]
                start, end = piece["germline_start"], piece["germline_end"]
                assert 0 <= start < end <= len(reference)
                reconstructed += reference[start:end]
            else:
                assert piece["segment"] in ("n1", "n2")
                reconstructed += piece["payload"]
            assert piece["query_end"] == len(reconstructed)
        assert reconstructed == sequence
        assert set(sequence) <= set("ACGT")
        vref = references["imgt_gapped", "v"][item["source_alleles"]["v"]["allele"]]
        assert vref[309:312] in ("TGT", "TGC")
        assert len(vref[:309].replace(".", "")) == item["anchor_evidence"]["v_imgt104_offset"] == item["junction_start"]
        jpiece = next(p for p in item["construction"] if p["segment"] == "j")
        jref = references["ungapped", "j"][jpiece["allele"]]
        anchors = [i for i in range(len(jref)) if re.match(r"TT[TC]GG[ACGT][ACGT]{3}GG[ACGT]", jref[i:])]
        assert anchors == [item["anchor_evidence"]["j_anchor_germline_offset"]]
        assert jpiece["query_start"] + anchors[0] + 3 == item["junction_end"]
        assert sequence[item["junction_start"]:item["junction_end"]] == item["junction"]
        assert str(Seq(item["junction"]).translate()) == item["junction_aa"]
        assert item["junction_aa"].startswith("C") and item["junction_aa"].endswith("F")
        assert item["junction_start"] % 3 == len(item["junction"]) % 3 == len(sequence) % 3 == 0
        assert str(Seq(sequence).translate()) == item["coding_translation"]
        assert "*" not in item["coding_translation"]
        assert item["productive"] is True
        vref = references["ungapped", "v"][item["source_alleles"]["v"]["allele"]]
        v_end = next((i for i, (query, germline) in enumerate(zip(sequence, vref)) if query != germline), len(vref))
        j_start = jpiece["query_start"]
        j_end = j_start + next((i for i, (query, germline) in enumerate(zip(sequence[j_start:], jref)) if query != germline), len(jref))
        assert {key: item["alignment_evidence"][key] for key in (
            "v_sequence_start", "v_sequence_end", "j_sequence_start", "j_sequence_end",
        )} == {"v_sequence_start": 0, "v_sequence_end": v_end,
              "j_sequence_start": j_start, "j_sequence_end": j_end}
        # Exact retained slices determine ambiguity independently of annotation.
        for segment in ("v", "d", "j", "c"):
            pieces = [p for p in item["construction"] if p["segment"] == segment]
            if not pieces:
                assert item[f"allowed_{segment}_calls"] == [None]
                continue
            piece = pieces[0]
            payload = sequence[piece["query_start"]:piece["query_end"]]
            ties = sorted({name.split("__")[0] for name, bases in references["ungapped", segment].items()
                           if name.startswith(item["locus"]) and bases[piece["germline_start"]:piece["germline_end"]] == payload})
            assert item[f"allowed_{segment}_calls"] == [",".join(ties)]


@pytest.mark.e2e
@pytest.mark.parametrize("case", load_tcr_cases(), ids=lambda case: case.locus)
def test_tcr_goldens(case):
    """Catch lost receptor/locus, calls, anchors, frame and strand propagation."""
    definition = next(item for item in json.loads((TCR_DATA_DIR / "cases.json").read_text())["cases"]
                      if item["sequence_id"] == case.sequence_id)
    reverse = str(Seq(case.sequence).reverse_complement())
    rows = abstar.run(
        [case.as_sequence(), Sequence(reverse, id=case.sequence_id + "_rc")],
        receptor="tcr", germline_database="human", n_processes=1,
        mmseqs_threads=1,
    )
    assert isinstance(rows, list) and len(rows) == 2
    for row, rev_comp, sequence_id in zip(rows, (False, True), (case.sequence_id, case.sequence_id + "_rc")):
        assert row.id == row["sequence_id"] == sequence_id
        assert row["annotation_status"] == "annotated"
        assert row["failure_reason"] is None
        assert row["rev_comp"] is rev_comp
        assert row["germline_database"] == "human"
        assert row["species"] == "homo_sapiens"
        assert row["locus"] == case.locus
        assert row["v_call"] in case.allowed_v_calls
        assert row["d_call"] in case.allowed_d_calls
        assert row["j_call"] in case.allowed_j_calls
        assert row["c_call"] in definition["allowed_c_calls"]
        assert row["junction"] == case.junction
        assert row["junction_aa"] == definition["junction_aa"]
        assert row["sequence_oriented"][definition["junction_start"]:definition["junction_end"]] == row["junction"]
        assert row["cdr3"] == case.junction[3:-3]
        assert row["cdr3_aa"] == definition["junction_aa"][1:-1]
        assert row["productive"] is case.productive
        assert row["productivity_issues"] == ""
        assert row["vj_in_frame"] is True
        for field in ("v_sequence_start", "v_sequence_end", "j_sequence_start", "j_sequence_end"):
            assert row[field] == definition["alignment_evidence"][field]
        assert row["sequence"] == case.sequence[:definition["alignment_evidence"]["j_sequence_end"]]
        assert row["sequence_vdjc"] == case.sequence
        assert row["v_identity"] == row["j_identity"] == 1.0
        assert row["v_mutations"] in (None, "")
        d_pieces = [piece for piece in definition["construction"] if piece["segment"] == "d"]
        if d_pieces:
            piece = d_pieces[0]
            assert row["d_sequence"] == row["d_germline"] == case.sequence[piece["query_start"]:piece["query_end"]]
            assert row["d_identity"] == 1.0
        else:
            assert row["d_gene"] is row["d_sequence"] is row["d_germline"] is None
        assert row["c_sequence"] == row["c_germline"] == case.sequence[row["c_sequence_start"]:row["c_sequence_end"]]
        assert row["c_identity"] == 1.0
        assert "row_id" not in row.annotations
        for suffix in ("", "_aa"):
            assert row[f"sequence_alignment{suffix}"].replace("-", "") == row[f"sequence{suffix}"]
            assert len(row[f"sequence_alignment{suffix}"]) == len(row[f"germline_alignment{suffix}"])
            assert len(row[f"sequence_vdjc{suffix}"]) == len(row[f"germline_vdjc{suffix}"])
            assert len(row[f"c_sequence{suffix}"]) == len(row[f"c_germline{suffix}"])
        for segment in ("v", "d", "j", "c"):
            if row[f"{segment}_call"] is not None:
                calls = row[f"{segment}_call"].split(",")
                assert calls == sorted(set(calls))
                assert all(call.startswith(case.locus) for call in calls)
    biological_fields = (
        "sequence", "sequence_aa", "sequence_vdjc", "sequence_vdjc_aa", "locus",
        "v_call", "d_call", "j_call", "c_call", "junction", "junction_aa",
        "productive", "productivity_issues",
        "v_sequence_start", "v_sequence_end", "d_sequence", "d_germline",
        "j_sequence_start", "j_sequence_end", "c_sequence_start", "c_sequence_end",
    )
    assert {key: rows[0][key] for key in biological_fields} == {key: rows[1][key] for key in biological_fields}


@pytest.mark.e2e
def test_tcr_receptor_and_coordinates_survive_chunked_project_output(tmp_path):
    """Inspect real assignment rows because receptor/anchor fields are internal."""
    from abstar.annotation.antibody import Antibody
    from abstar.annotation.annotator import annotate_single_sequence

    cases = load_tcr_cases()
    records = [record for case in cases for record in (
        case.as_sequence(),
        Sequence(str(Seq(case.sequence).reverse_complement()), id=case.sequence_id + "_rc"),
    )]
    abstar.run(
        records, project_path=str(tmp_path), receptor="tcr", germline_database="human",
        output_format="parquet", n_processes=2, chunksize=3, mmseqs_threads=1, debug=True,
    )
    rows = pl.read_parquet(tmp_path / "parquet" / "sequences.parquet").to_dicts()
    assert [row["sequence_id"] for row in rows] == [record.id for record in records]
    assert all(row["annotation_status"] == "annotated" and row["failure_reason"] is None for row in rows)
    assert all("row_id" not in row for row in rows)
    definitions = {item["sequence_id"]: item for item in json.loads((TCR_DATA_DIR / "cases.json").read_text())["cases"]}
    assignments = [row for path in sorted((tmp_path / "tmp").glob("chunk_*.parquet"))
                   if re.fullmatch(r"chunk_\d+\.parquet", path.name)
                   for row in pl.read_parquet(path).to_dicts()]
    assert [row["sequence_id"] for row in assignments] == [record.id for record in records]
    for assignment in assignments:
        sequence_id = assignment["sequence_id"]
        definition = definitions[sequence_id.removesuffix("_rc")]
        receptor = assignment.pop("receptor_type")
        assert receptor == "tcr"
        assert assignment["germline_database"] == "human"
        ab = Antibody(**assignment)
        ab.receptor_type = receptor
        ab = annotate_single_sequence(ab, germline_database="human")
        assert ab.receptor_type == "tcr"
        assert ab.locus == definition["locus"]
        assert (ab.junction_start, ab.junction_end) == (definition["junction_start"], definition["junction_end"])
        assert ab.junction == definition["junction"]
        assert ab.rev_comp is sequence_id.endswith("_rc")
        assert ab.productive is True and ab.productivity_issues == ""
        d_pieces = [piece for piece in definition["construction"] if piece["segment"] == "d"]
        if d_pieces:
            piece = d_pieces[0]
            assert (ab.d_sequence_start, ab.d_sequence_end) == (piece["query_start"], piece["query_end"])
            assert (ab.d_germline_start, ab.d_germline_end) == (piece["germline_start"], piece["germline_end"])
