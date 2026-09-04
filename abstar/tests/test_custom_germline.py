from types import SimpleNamespace

import pytest
from abutils import Sequence

from ..annotation.germline import get_germline_database_path
from ..core import germline


@pytest.fixture
def custom_genes(tmp_path):
    path = tmp_path / "genes.fasta"
    path.write_text(
        ">IGHVTEST*01\nACGT.ACGTACGT\n"
        ">IGHJTEST*01\nACGTACGTACGT\n"
    )
    return path


@pytest.fixture
def fake_mmseqs(monkeypatch):
    calls = []

    def make_database(input_file, output_file, debug=False):
        calls.append((input_file, output_file))
        with open(output_file, "w") as output:
            output.write("index")

    monkeypatch.setattr(germline, "_make_mmseqs_db", make_database)
    return calls


@pytest.mark.integration
def test_build_custom_database_clean_room_and_discovery(
    tmp_path, monkeypatch, custom_genes, fake_mmseqs
):
    monkeypatch.setenv("HOME", str(tmp_path))

    germline.build_germline_database(
        "custom",
        fastas=str(custom_genes),
        include_species_in_name=False,
        verbose=False,
    )

    database = tmp_path / ".abstar" / "germline_dbs" / "bcr" / "custom"
    assert get_germline_database_path("custom", "bcr") == str(database)
    assert (database / "imgt_gapped" / "v.fasta").exists()
    assert (database / "imgt_gapped" / "j.fasta").exists()
    assert (database / "ungapped" / "v.fasta").exists()
    assert (database / "ungapped" / "j.fasta").exists()
    assert not (database / "ungapped" / "d.fasta").exists()
    assert [path.rsplit("/", 1)[-1] for _, path in fake_mmseqs] == ["v", "j"]
    assert not list((database.parent).glob(".custom.staging-*"))


@pytest.mark.integration
def test_overwrite_removes_stale_files(
    tmp_path, monkeypatch, custom_genes, fake_mmseqs
):
    location = tmp_path / "databases"
    kwargs = dict(
        name="custom",
        fastas=str(custom_genes),
        location=str(location),
        include_species_in_name=False,
        verbose=False,
    )
    germline.build_germline_database(**kwargs)
    database = location / "bcr" / "custom"
    (database / "stale.txt").write_text("old")
    monkeypatch.setattr("builtins.input", lambda _: "yes")

    germline.build_germline_database(**kwargs)

    assert not (database / "stale.txt").exists()


@pytest.mark.integration
def test_failed_overwrite_preserves_existing_database(
    tmp_path, monkeypatch, custom_genes, fake_mmseqs
):
    location = tmp_path / "databases"
    kwargs = dict(
        name="custom",
        fastas=str(custom_genes),
        location=str(location),
        include_species_in_name=False,
        verbose=False,
    )
    germline.build_germline_database(**kwargs)
    database = location / "bcr" / "custom"
    marker = database / "existing.txt"
    marker.write_text("preserve me")
    monkeypatch.setattr("builtins.input", lambda _: "yes")
    monkeypatch.setattr(
        germline,
        "_make_mmseqs_db",
        lambda *args, **kwargs: (_ for _ in ()).throw(RuntimeError("index failed")),
    )

    with pytest.raises(RuntimeError, match="index failed"):
        germline.build_germline_database(**kwargs)

    assert marker.read_text() == "preserve me"
    assert not list((database.parent).glob(".custom.staging-*"))


@pytest.mark.parametrize(
    "sequences,receptor,error",
    [
        (
            [Sequence("ACGT", id="IGHV1*01"), Sequence("ACGT", id="IGHV1*01")],
            "bcr",
            "duplicate",
        ),
        ([Sequence("ACGT", id="TRAV1*01")], "tcr", "requires"),
        (
            [Sequence("ACGT", id="TRAV1*01"), Sequence("ACGT", id="TRAJ1*01"), Sequence("ACGT", id="TRAD1*01")],
            "tcr",
            "does not contain D",
        ),
        (
            [Sequence("ACGT!", id="IGHV1*01"), Sequence("ACGT", id="IGHJ1*01")],
            "bcr",
            "invalid nucleotide",
        ),
    ],
)
def test_validate_germlines_rejects_unsafe_inputs(sequences, receptor, error):
    with pytest.raises(ValueError, match=error):
        germline.validate_germlines(sequences, [], receptor)


@pytest.mark.parametrize(
    "identifier,receptor,reference_prefix",
    [("IGKVTEST*01", "bcr", "IGKV"), ("TRAVTEST*01", "tcr", "TRAV")],
)
def test_add_imgt_gaps_handles_mixed_inputs_with_locus_reference(
    monkeypatch, identifier, receptor, reference_prefix
):
    already_gapped = Sequence("A.CGT", id=identifier.replace("TEST", "GAPPED"))
    ungapped = Sequence("ACGT", id=identifier)
    lookups = []

    def get_reference(name, **kwargs):
        lookups.append((name, kwargs["receptor"]))
        return [Sequence("ACGT", id=f"{reference_prefix}REF*01")]

    monkeypatch.setattr(germline, "get_germline", get_reference)
    monkeypatch.setattr(
        germline.abutils.tl,
        "semiglobal_alignment",
        lambda *args, **kwargs: [
            SimpleNamespace(aligned_query="ACGT", aligned_target="ACGT")
        ],
    )

    result = germline.add_imgt_gaps(
        [already_gapped, ungapped], reference="human", receptor=receptor
    )

    assert result[0].sequence == "A.CGT"
    assert result[1].sequence == "ACGT"
    assert lookups == [(reference_prefix, receptor)]


def test_make_mmseqs_db_uses_checked_argument_lists(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(germline.abutils.bin, "get_path", lambda _: "/bin/mmseqs")

    def run(command, **kwargs):
        calls.append((command, kwargs))
        return SimpleNamespace(stdout="", stderr="")

    monkeypatch.setattr(germline.sp, "run", run)
    germline._make_mmseqs_db(str(tmp_path / "v.fasta"), str(tmp_path / "v"))

    assert calls[0][0] == [
        "/bin/mmseqs",
        "createdb",
        str(tmp_path / "v.fasta"),
        str(tmp_path / "v"),
    ]
    assert calls[1][0][:3] == ["/bin/mmseqs", "createindex", str(tmp_path / "v")]
    assert all(kwargs["check"] is True for _, kwargs in calls)
    assert all("shell" not in kwargs for _, kwargs in calls)
