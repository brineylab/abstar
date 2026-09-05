import hashlib
import os
import subprocess
from pathlib import Path
from types import SimpleNamespace

import pytest
from abutils import Sequence

from ..annotation.germline import get_germline_database_path
from ..core import germline


MMSEQS_COMPONENT_SUFFIXES = ("", ".dbtype", ".index", "_h", "_h.dbtype", "_h.index")


def _tree_snapshot(root):
    root = Path(root)
    if not root.exists():
        return None
    root_stat = root.lstat()
    snapshot = [(".", True, root_stat.st_mode, root_stat.st_size, root_stat.st_mtime_ns, None)]
    for path in sorted(root.rglob("*")):
        stat = path.lstat()
        digest = hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None
        snapshot.append(
            (str(path.relative_to(root)), path.is_dir(), stat.st_mode, stat.st_size,
             stat.st_mtime_ns, digest)
        )
    return tuple(snapshot)


@pytest.fixture(autouse=True)
def isolated_user_germline_root(tmp_path, monkeypatch):
    real_root = Path.home() / ".abstar" / "germline_dbs"
    before = _tree_snapshot(real_root)
    monkeypatch.setenv("HOME", str(tmp_path))
    yield
    assert _tree_snapshot(real_root) == before


@pytest.fixture
def custom_genes(tmp_path):
    path = tmp_path / "genes.fasta"
    path.write_text(
        ">IGHVTEST*01\nACGT.ACGTACGT\n"
        ">IGHJTEST*01\nACGTACGTACGT\n"
    )
    return path


@pytest.fixture
def complete_inputs(tmp_path):
    genes = tmp_path / "complete-genes.fasta"
    genes.write_text(
        ">IGHVTEST*01\nACGT.ACGTACGT\n"
        ">IGHDTEST*01\nACGT\n"
        ">IGHJTEST*01\nACGTACGTACGT\n"
    )
    constants = tmp_path / "constants.fasta"
    constants.write_text(">IGHMTEST*01\nACGTACGT\n")
    manifest = tmp_path / "manifest.txt"
    manifest.write_text("source: task-19 transactional fixture\nlicense: MIT\n")
    return genes, constants, manifest


def _write_fake_index(output_file):
    for suffix in MMSEQS_COMPONENT_SUFFIXES:
        Path(f"{output_file}{suffix}").write_text(f"index:{Path(output_file).name}{suffix}\n")


@pytest.fixture
def fake_mmseqs(monkeypatch):
    calls = []

    def make_database(input_file, output_file, debug=False):
        calls.append((input_file, output_file))
        _write_fake_index(output_file)

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
    replace = os.replace
    moves = []

    def record_replace(source, destination):
        source = Path(source)
        destination = Path(destination)
        if source == database:
            assert database.is_dir()
            assert destination.name.startswith("custom.backup-")
        elif source.name.startswith(".custom.staging-"):
            assert not database.exists()
            assert (source / "ungapped" / "v.fasta").is_file()
            assert (source / "mmseqs" / "j.index").is_file()
            assert destination == database
        moves.append((source, destination))
        replace(source, destination)

    monkeypatch.setattr(germline.os, "replace", record_replace)

    germline.build_germline_database(**kwargs)

    assert not (database / "stale.txt").exists()
    assert len(moves) == 2
    assert moves[0][0] == database
    assert moves[1][0].name.startswith(".custom.staging-")
    assert moves[1][1] == database


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


@pytest.mark.integration
def test_complete_database_is_validated_in_staging_before_atomic_publication(
    tmp_path, monkeypatch, complete_inputs, fake_mmseqs
):
    genes, constants, manifest = complete_inputs
    database_root = tmp_path / ".abstar" / "germline_dbs" / "bcr"
    database = database_root / "complete"
    validations = []
    publications = []
    validate = germline.validate_staged_database
    replace = os.replace

    def record_validation(staging_dir, receptor, **kwargs):
        assert Path(staging_dir).parent == database_root
        assert Path(staging_dir).name.startswith(".complete.staging-")
        validate(staging_dir, receptor, **kwargs)
        validations.append(Path(staging_dir))

    def record_publication(source, destination):
        source = Path(source)
        destination = Path(destination)
        if destination == database:
            assert validations == [source]
            assert not database.exists()
            publications.append((source, destination))
        replace(source, destination)

    monkeypatch.setattr(germline, "validate_staged_database", record_validation)
    monkeypatch.setattr(germline.os, "replace", record_publication)

    germline.build_germline_database(
        "complete",
        fastas=str(genes),
        constants=str(constants),
        manifest=str(manifest),
        include_species_in_name=False,
        verbose=False,
    )

    assert publications and publications[0][1] == database
    assert (database / "manifest.txt").read_text() == manifest.read_text()
    for segment in "vdjc":
        gapped = {sequence.id: sequence.sequence for sequence in
                  germline.abutils.io.read_fasta(str(database / "imgt_gapped" / f"{segment}.fasta"))}
        ungapped = {sequence.id: sequence.sequence for sequence in
                    germline.abutils.io.read_fasta(str(database / "ungapped" / f"{segment}.fasta"))}
        assert list(gapped) == list(ungapped)
        assert {identifier: sequence.replace(".", "") for identifier, sequence in gapped.items()} == ungapped
        for suffix in MMSEQS_COMPONENT_SUFFIXES:
            assert (database / "mmseqs" / f"{segment}{suffix}").is_file()
    assert not list(database_root.glob(".complete.staging-*"))
    assert not list(database_root.glob("complete.backup-*"))


@pytest.mark.parametrize("boundary", ("gapping", "v", "d", "j", "c"))
@pytest.mark.parametrize("replacing", (False, True), ids=("new", "replacement"))
def test_failure_at_each_build_boundary_cleans_staging_and_preserves_old_database(
    tmp_path, monkeypatch, complete_inputs, boundary, replacing
):
    genes, constants, manifest = complete_inputs
    database_root = tmp_path / ".abstar" / "germline_dbs" / "bcr"
    database = database_root / "custom"
    if replacing:
        database.mkdir(parents=True)
        (database / "marker.bin").write_bytes(b"existing-database\x00\xff")
        (database / "nested").mkdir()
        (database / "nested" / "state.txt").write_text("unchanged\n")
        before = _tree_snapshot(database)
        monkeypatch.setattr("builtins.input", lambda _: "yes")

    def failure(label):
        return germline.GermlineBuildExternalToolError(
            f"{label} failed",
            command=["/opt/task19-tool", label, "--checked"],
            returncode=31,
            stdout=f"{label}-stdout",
            stderr=f"{label}-stderr",
        )

    if boundary == "gapping":
        def fail_gapping(*args, **kwargs):
            try:
                raise ValueError("TASK19-GAPPING-CAUSE")
            except ValueError as cause:
                raise RuntimeError("TASK19-GAPPING-FAILURE") from cause

        monkeypatch.setattr(germline, "add_imgt_gaps", fail_gapping)
    else:
        monkeypatch.setattr(germline.abutils.bin, "get_path", lambda _: "/opt/mmseqs")

        def run(command, **kwargs):
            assert kwargs == {"check": True, "capture_output": True, "text": True}
            if command[1] == "createdb":
                _write_fake_index(command[3])
            elif command[1] == "createindex" and Path(command[2]).name == boundary:
                raise subprocess.CalledProcessError(
                    31, command, output=f"{boundary}-stdout", stderr=f"{boundary}-stderr"
                )
            return SimpleNamespace(stdout="", stderr="")

        monkeypatch.setattr(germline.sp, "run", run)

    expected_error = RuntimeError if boundary == "gapping" else germline.GermlineBuildExternalToolError
    with pytest.raises(expected_error) as captured:
        germline.build_germline_database(
            "custom", fastas=str(genes), constants=str(constants), manifest=str(manifest),
            include_species_in_name=False, verbose=False,
        )

    error = captured.value
    if boundary == "gapping":
        assert str(error) == "TASK19-GAPPING-FAILURE"
        assert isinstance(error.__cause__, ValueError)
        assert str(error.__cause__) == "TASK19-GAPPING-CAUSE"
    else:
        assert error.command[:2] == ("/opt/mmseqs", "createindex")
        assert Path(error.command[2]).name == boundary
        assert error.returncode == 31
        assert error.stdout == f"{boundary}-stdout"
        assert error.stderr == f"{boundary}-stderr"
        assert f"{boundary}-stdout" in str(error)
        assert f"{boundary}-stderr" in str(error)
    if replacing:
        assert _tree_snapshot(database) == before
    else:
        assert not database.exists()
    assert not list(database_root.glob(".custom.staging-*"))
    assert not list(database_root.glob("custom.backup-*"))


def test_incomplete_staged_index_is_rejected_before_publication(
    tmp_path, monkeypatch, custom_genes
):
    database_root = tmp_path / ".abstar" / "germline_dbs" / "bcr"
    database = database_root / "invalid"

    def incomplete_index(input_file, output_file, debug=False):
        _write_fake_index(output_file)
        if Path(output_file).name == "j":
            Path(f"{output_file}.index").unlink()

    monkeypatch.setattr(germline, "_make_mmseqs_db", incomplete_index)

    with pytest.raises(ValueError, match=r"incomplete staged germline database.*mmseqs/j\.index"):
        germline.build_germline_database(
            "invalid", fastas=str(custom_genes), include_species_in_name=False, verbose=False
        )

    assert not database.exists()
    assert not list(database_root.glob(".invalid.staging-*"))


def test_publication_failure_rolls_back_existing_database(
    tmp_path, monkeypatch, custom_genes, fake_mmseqs
):
    database_root = tmp_path / ".abstar" / "germline_dbs" / "bcr"
    database = database_root / "custom"
    database.mkdir(parents=True)
    (database / "marker.bin").write_bytes(b"original\x00database")
    before = _tree_snapshot(database)
    monkeypatch.setattr("builtins.input", lambda _: "yes")
    replace = os.replace

    def fail_publication(source, destination):
        if Path(source).name.startswith(".custom.staging-") and Path(destination) == database:
            raise OSError("TASK19-PUBLISH-FAILURE")
        replace(source, destination)

    monkeypatch.setattr(germline.os, "replace", fail_publication)

    with pytest.raises(OSError, match="TASK19-PUBLISH-FAILURE"):
        germline.build_germline_database(
            "custom", fastas=str(custom_genes), include_species_in_name=False, verbose=False
        )

    assert _tree_snapshot(database) == before
    assert not list(database_root.glob(".custom.staging-*"))
    assert not list(database_root.glob("custom.backup-*"))


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


def test_mmseqs_build_failure_retains_checked_process_diagnostics(monkeypatch):
    command = ["/opt/mmseqs", "createindex", "/tmp/custom/v", "/tmp/index"]

    def fail(actual, **kwargs):
        assert actual == command
        assert kwargs == {"check": True, "capture_output": True, "text": True}
        raise subprocess.CalledProcessError(
            23, actual, output="TASK19-STDOUT", stderr="TASK19-STDERR"
        )

    monkeypatch.setattr(germline.sp, "run", fail)

    with pytest.raises(RuntimeError) as captured:
        germline._run_mmseqs_command(command)

    error = captured.value
    assert error.command == tuple(command)
    assert error.returncode == 23
    assert error.stdout == "TASK19-STDOUT"
    assert error.stderr == "TASK19-STDERR"
    assert "TASK19-STDOUT" in str(error)
    assert "TASK19-STDERR" in str(error)


def test_mmseqs_launch_failure_retains_command(monkeypatch):
    command = ["/missing/mmseqs", "createdb", "/tmp/v.fasta", "/tmp/v"]

    def fail(actual, **kwargs):
        raise FileNotFoundError(2, "No such file or directory", actual[0])

    monkeypatch.setattr(germline.sp, "run", fail)

    with pytest.raises(germline.GermlineBuildExternalToolError) as captured:
        germline._run_mmseqs_command(command)

    assert captured.value.command == tuple(command)
    assert captured.value.returncode is None
    assert captured.value.stdout == ""
    assert "No such file or directory" in captured.value.stderr
