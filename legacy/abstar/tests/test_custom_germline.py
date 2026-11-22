import hashlib
import os
import stat
import subprocess
import warnings
from pathlib import Path
from types import SimpleNamespace

import pytest
from abutils import Sequence

from ..annotation.germline import get_germline_database_path
from ..core import germline


MMSEQS_COMPONENT_SUFFIXES = ("", ".dbtype", ".index", "_h", "_h.dbtype", "_h.index")


def _tree_snapshot(root):
    root = Path(root)
    if not os.path.lexists(root):
        return None
    root_stat = root.lstat()
    if stat.S_ISLNK(root_stat.st_mode):
        return ((".", "symlink", root_stat.st_mode, root_stat.st_size,
                 root_stat.st_mtime_ns, os.readlink(root)),)
    if stat.S_ISREG(root_stat.st_mode):
        return ((".", "file", root_stat.st_mode, root_stat.st_size,
                 root_stat.st_mtime_ns,
                 hashlib.sha256(root.read_bytes()).hexdigest()),)
    snapshot = [(".", "directory", root_stat.st_mode, root_stat.st_size,
                 root_stat.st_mtime_ns, None)]
    for path in sorted(root.rglob("*")):
        metadata = path.lstat()
        if stat.S_ISLNK(metadata.st_mode):
            kind = "symlink"
            detail = os.readlink(path)
        elif stat.S_ISREG(metadata.st_mode):
            kind = "file"
            detail = hashlib.sha256(path.read_bytes()).hexdigest()
        else:
            kind = "directory" if stat.S_ISDIR(metadata.st_mode) else "other"
            detail = None
        snapshot.append(
            (str(path.relative_to(root)), kind, metadata.st_mode, metadata.st_size,
             metadata.st_mtime_ns, detail)
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


@pytest.mark.parametrize("replacing", (False, True), ids=("new", "replacement"))
def test_gapping_failure_cleans_staging_and_preserves_old_database(
    tmp_path, monkeypatch, complete_inputs, replacing
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

    def fail_gapping(*args, **kwargs):
        try:
            raise ValueError("TASK19-GAPPING-CAUSE")
        except ValueError as cause:
            raise RuntimeError("TASK19-GAPPING-FAILURE") from cause

    monkeypatch.setattr(germline, "add_imgt_gaps", fail_gapping)

    with pytest.raises(RuntimeError) as captured:
        germline.build_germline_database(
            "custom", fastas=str(genes), constants=str(constants), manifest=str(manifest),
            include_species_in_name=False, verbose=False,
        )

    error = captured.value
    assert str(error) == "TASK19-GAPPING-FAILURE"
    assert isinstance(error.__cause__, ValueError)
    assert str(error.__cause__) == "TASK19-GAPPING-CAUSE"
    if replacing:
        assert _tree_snapshot(database) == before
    else:
        assert not database.exists()
    assert not list(database_root.glob(".custom.staging-*"))
    assert not list(database_root.glob("custom.backup-*"))


@pytest.mark.parametrize("step", ("createdb", "createindex"))
@pytest.mark.parametrize("boundary", tuple("vdjc"))
@pytest.mark.parametrize("replacing", (False, True), ids=("new", "replacement"))
def test_mmseqs_failure_at_each_segment_and_command_preserves_destination(
    tmp_path, monkeypatch, complete_inputs, step, boundary, replacing
):
    genes, constants, manifest = complete_inputs
    database_root = tmp_path / ".abstar" / "germline_dbs" / "bcr"
    database = database_root / "custom"
    if replacing:
        database.mkdir(parents=True)
        (database / "marker.bin").write_bytes(b"existing-database\x00\xff")
        before = _tree_snapshot(database)
        monkeypatch.setattr("builtins.input", lambda _: "yes")
    monkeypatch.setattr(germline.abutils.bin, "get_path", lambda _: "/opt/mmseqs")

    def run(command, **kwargs):
        assert kwargs == {"check": True, "capture_output": True, "text": True}
        output = command[3] if command[1] == "createdb" else command[2]
        segment = Path(output).name
        if command[1] == step and segment == boundary:
            raise subprocess.CalledProcessError(
                31, command, output=f"{step}-{boundary}-stdout",
                stderr=f"{step}-{boundary}-stderr",
            )
        if command[1] == "createdb":
            _write_fake_index(command[3])
        return SimpleNamespace(stdout="", stderr="")

    monkeypatch.setattr(germline.sp, "run", run)

    with pytest.raises(germline.GermlineBuildExternalToolError) as captured:
        germline.build_germline_database(
            "custom", fastas=str(genes), constants=str(constants), manifest=str(manifest),
            include_species_in_name=False, verbose=False,
        )

    error = captured.value
    assert error.command[:2] == ("/opt/mmseqs", step)
    output_position = 3 if step == "createdb" else 2
    assert Path(error.command[output_position]).name == boundary
    assert error.returncode == 31
    assert error.stdout == f"{step}-{boundary}-stdout"
    assert error.stderr == f"{step}-{boundary}-stderr"
    assert error.stdout in str(error)
    assert error.stderr in str(error)
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


@pytest.mark.parametrize("kind", ("file", "symlink", "dangling-symlink"))
def test_unsafe_destination_is_rejected_before_staging(
    tmp_path, custom_genes, fake_mmseqs, kind
):
    database_root = tmp_path / ".abstar" / "germline_dbs" / "bcr"
    database_root.mkdir(parents=True)
    database = database_root / "custom"
    if kind == "file":
        database.write_text("do not replace\n")
    else:
        target = tmp_path / ("target" if kind == "symlink" else "missing-target")
        if kind == "symlink":
            target.mkdir()
        database.symlink_to(target, target_is_directory=True)
    before = _tree_snapshot(database)

    with pytest.raises(ValueError, match="directory|symlink"):
        germline.build_germline_database(
            "custom", fastas=str(custom_genes), include_species_in_name=False,
            verbose=False,
        )

    assert _tree_snapshot(database) == before
    assert fake_mmseqs == []
    assert not list(database_root.glob(".custom.staging-*"))


def test_absent_destination_that_appears_during_build_is_not_touched(
    tmp_path, monkeypatch, custom_genes, fake_mmseqs
):
    database_root = tmp_path / ".abstar" / "germline_dbs" / "bcr"
    database = database_root / "custom"
    validate = germline.validate_staged_database

    def appear(staging_dir, receptor, **kwargs):
        validate(staging_dir, receptor, **kwargs)
        database.mkdir()
        (database / "owner.txt").write_text("other builder\n")

    monkeypatch.setattr(germline, "validate_staged_database", appear)

    with pytest.raises(RuntimeError, match="changed during build"):
        germline.build_germline_database(
            "custom", fastas=str(custom_genes), include_species_in_name=False,
            verbose=False,
        )

    assert (database / "owner.txt").read_text() == "other builder\n"
    assert not list(database_root.glob(".custom.staging-*"))


def test_approved_destination_identity_change_is_not_touched(
    tmp_path, monkeypatch, custom_genes, fake_mmseqs
):
    database_root = tmp_path / ".abstar" / "germline_dbs" / "bcr"
    database = database_root / "custom"
    database.mkdir(parents=True)
    (database / "approved.txt").write_text("approved original\n")
    displaced = database_root / "displaced"
    validate = germline.validate_staged_database
    monkeypatch.setattr("builtins.input", lambda _: "yes")

    def replace_identity(staging_dir, receptor, **kwargs):
        validate(staging_dir, receptor, **kwargs)
        os.replace(database, displaced)
        database.mkdir()
        (database / "owner.txt").write_text("new owner\n")

    monkeypatch.setattr(germline, "validate_staged_database", replace_identity)

    with pytest.raises(RuntimeError, match="changed during build"):
        germline.build_germline_database(
            "custom", fastas=str(custom_genes), include_species_in_name=False,
            verbose=False,
        )

    assert (database / "owner.txt").read_text() == "new owner\n"
    assert (displaced / "approved.txt").read_text() == "approved original\n"
    assert not list(database_root.glob(".custom.staging-*"))


@pytest.mark.parametrize(
    "component", ("ungapped/j.fasta", "mmseqs/j.index", "manifest.txt")
)
@pytest.mark.parametrize("defect", ("symlink", "empty"))
def test_staged_components_must_be_owned_regular_nonempty_files(
    tmp_path, monkeypatch, complete_inputs, fake_mmseqs, component, defect
):
    genes, constants, manifest = complete_inputs
    database_root = tmp_path / ".abstar" / "germline_dbs" / "bcr"
    database = database_root / "custom"
    validate = germline.validate_staged_database

    def corrupt(staging_dir, receptor, **kwargs):
        path = Path(staging_dir) / component
        path.unlink()
        if defect == "symlink":
            external = tmp_path / "external-component"
            external.write_text("external bytes\n")
            path.symlink_to(external)
        else:
            path.touch()
        validate(staging_dir, receptor, **kwargs)

    monkeypatch.setattr(germline, "validate_staged_database", corrupt)

    with pytest.raises(ValueError, match="regular non-symlink|empty|IDs differ"):
        germline.build_germline_database(
            "custom", fastas=str(genes), constants=str(constants), manifest=str(manifest),
            include_species_in_name=False, verbose=False,
        )

    assert not database.exists()
    assert not list(database_root.glob(".custom.staging-*"))


def test_first_replacement_rename_failure_preserves_primary_error_and_database(
    tmp_path, monkeypatch, custom_genes, fake_mmseqs
):
    database_root = tmp_path / ".abstar" / "germline_dbs" / "bcr"
    database = database_root / "custom"
    database.mkdir(parents=True)
    (database / "marker").write_text("original\n")
    before = _tree_snapshot(database)
    monkeypatch.setattr("builtins.input", lambda _: "yes")
    replace = os.replace

    def fail_first(source, destination):
        if Path(source) == database:
            raise OSError("TASK19-FIRST-RENAME")
        replace(source, destination)

    monkeypatch.setattr(germline.os, "replace", fail_first)

    with pytest.raises(OSError, match="TASK19-FIRST-RENAME"):
        germline.build_germline_database(
            "custom", fastas=str(custom_genes), include_species_in_name=False,
            verbose=False,
        )

    assert _tree_snapshot(database) == before
    assert not list(database_root.glob("custom.backup-*"))


def test_publish_and_rollback_failure_retains_both_errors_and_recovery_path(
    tmp_path, monkeypatch, custom_genes, fake_mmseqs
):
    database_root = tmp_path / ".abstar" / "germline_dbs" / "bcr"
    database = database_root / "custom"
    database.mkdir(parents=True)
    (database / "marker").write_text("recover me\n")
    monkeypatch.setattr("builtins.input", lambda _: "yes")
    replace = os.replace

    def fail_publish_and_rollback(source, destination):
        source = Path(source)
        if source.name.startswith(".custom.staging-"):
            raise OSError("TASK19-PUBLISH")
        if source.name.startswith("custom.backup-"):
            raise OSError("TASK19-ROLLBACK")
        replace(source, destination)

    monkeypatch.setattr(germline.os, "replace", fail_publish_and_rollback)

    with pytest.raises(germline.GermlinePublicationError) as captured:
        germline.build_germline_database(
            "custom", fastas=str(custom_genes), include_species_in_name=False,
            verbose=False,
        )

    error = captured.value
    assert "TASK19-PUBLISH" in str(error)
    assert "TASK19-ROLLBACK" in str(error)
    assert isinstance(error.primary_error, OSError)
    assert isinstance(error.rollback_error, OSError)
    assert Path(error.backup_path).is_dir()
    assert (Path(error.backup_path) / "marker").read_text() == "recover me\n"


def test_backup_cleanup_failure_logs_success_and_retained_backup_under_strict_warnings(
    tmp_path, monkeypatch, custom_genes, fake_mmseqs, caplog
):
    database_root = tmp_path / ".abstar" / "germline_dbs" / "bcr"
    database = database_root / "custom"
    database.mkdir(parents=True)
    (database / "marker").write_text("old database\n")
    monkeypatch.setattr("builtins.input", lambda _: "yes")
    rmtree = germline.shutil.rmtree

    def fail_backup_cleanup(path, *args, **kwargs):
        if Path(path).name.startswith("custom.backup-"):
            raise OSError("TASK19-CLEANUP")
        return rmtree(path, *args, **kwargs)

    monkeypatch.setattr(germline.shutil, "rmtree", fail_backup_cleanup)

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        with caplog.at_level("WARNING", logger=germline.__name__):
            result = germline.build_germline_database(
                "custom", fastas=str(custom_genes), include_species_in_name=False,
                verbose=False,
            )

    assert result is None
    assert database.is_dir()
    assert not (database / "marker").exists()
    backups = list(database_root.glob("custom.backup-*"))
    assert len(backups) == 1
    assert (backups[0] / "marker").read_text() == "old database\n"
    assert "TASK19-CLEANUP" in caplog.text
    assert str(backups[0]) in caplog.text


def test_logging_failure_cannot_turn_successful_publication_into_failure(
    tmp_path, monkeypatch, custom_genes, fake_mmseqs
):
    database_root = tmp_path / ".abstar" / "germline_dbs" / "bcr"
    database = database_root / "custom"
    database.mkdir(parents=True)
    (database / "marker").write_text("old database\n")
    monkeypatch.setattr("builtins.input", lambda _: "yes")
    rmtree = germline.shutil.rmtree

    def fail_backup_cleanup(path, *args, **kwargs):
        if Path(path).name.startswith("custom.backup-"):
            raise OSError("TASK19-CLEANUP")
        return rmtree(path, *args, **kwargs)

    monkeypatch.setattr(germline.shutil, "rmtree", fail_backup_cleanup)
    monkeypatch.setattr(
        germline.logger, "warning",
        lambda *args, **kwargs: (_ for _ in ()).throw(RuntimeError("broken handler")),
    )

    assert germline.build_germline_database(
        "custom", fastas=str(custom_genes), include_species_in_name=False,
        verbose=False,
    ) is None
    assert database.is_dir()
    assert not (database / "marker").exists()


def test_windows_lock_backend_is_balanced_and_closes_descriptor(
    tmp_path, monkeypatch
):
    calls = []

    class FakeMSVCRT:
        LK_LOCK = 1
        LK_UNLCK = 2

        @staticmethod
        def locking(descriptor, operation, size):
            calls.append((descriptor, operation, size))

    monkeypatch.setattr(germline, "_fcntl", None)
    monkeypatch.setattr(germline, "_msvcrt", FakeMSVCRT)

    with germline.database_build_lock(str(tmp_path), "custom"):
        descriptor = calls[0][0]
        os.fstat(descriptor)

    assert calls == [(descriptor, FakeMSVCRT.LK_LOCK, 1),
                     (descriptor, FakeMSVCRT.LK_UNLCK, 1)]
    with pytest.raises(OSError):
        os.fstat(descriptor)


def test_lock_symlink_is_rejected_before_windows_backend_open(
    tmp_path, monkeypatch
):
    target = tmp_path / "external-lock"
    target.write_text("external\n")
    (tmp_path / ".custom.lock").symlink_to(target)
    calls = []
    fake = SimpleNamespace(
        LK_LOCK=1, LK_UNLCK=2,
        locking=lambda *args: calls.append(args),
    )
    monkeypatch.setattr(germline, "_fcntl", None)
    monkeypatch.setattr(germline, "_msvcrt", fake)

    with pytest.raises(ValueError, match="lock must not be a symlink"):
        with germline.database_build_lock(str(tmp_path), "custom"):
            pass

    assert calls == []
    assert target.read_text() == "external\n"


def test_unlock_failure_after_publication_is_logged_and_descriptor_is_closed(
    tmp_path, monkeypatch, custom_genes, fake_mmseqs, caplog
):
    database = tmp_path / ".abstar" / "germline_dbs" / "bcr" / "custom"
    closed = []
    close = os.close

    def record_close(descriptor):
        closed.append(descriptor)
        close(descriptor)

    monkeypatch.setattr(
        germline, "_unlock_database_descriptor",
        lambda descriptor: (_ for _ in ()).throw(OSError("TASK19-UNLOCK")),
    )
    monkeypatch.setattr(germline.os, "close", record_close)

    with caplog.at_level("WARNING", logger=germline.__name__):
        result = germline.build_germline_database(
            "custom", fastas=str(custom_genes), include_species_in_name=False,
            verbose=False,
        )

    assert result is None
    assert database.is_dir()
    assert (database / "mmseqs" / "j.index").is_file()
    assert len(closed) == 1
    with pytest.raises(OSError):
        os.fstat(closed[0])
    assert "TASK19-UNLOCK" in caplog.text
    assert ".custom.lock" in caplog.text


def test_unlock_failure_during_body_error_preserves_body_and_closes_descriptor(
    tmp_path, monkeypatch, caplog
):
    closed = []
    close = os.close

    def record_close(descriptor):
        closed.append(descriptor)
        close(descriptor)

    monkeypatch.setattr(
        germline, "_unlock_database_descriptor",
        lambda descriptor: (_ for _ in ()).throw(OSError("TASK19-UNLOCK")),
    )
    monkeypatch.setattr(germline.os, "close", record_close)
    cause = ValueError("TASK19-BODY-CAUSE")
    body_error = RuntimeError("TASK19-BODY")
    body_error.__cause__ = cause

    with caplog.at_level("WARNING", logger=germline.__name__):
        with pytest.raises(RuntimeError) as captured:
            with germline.database_build_lock(str(tmp_path), "custom"):
                raise body_error

    assert captured.value is body_error
    assert captured.value.__cause__ is cause
    assert len(closed) == 1
    with pytest.raises(OSError):
        os.fstat(closed[0])
    assert "TASK19-UNLOCK" in caplog.text


@pytest.mark.parametrize("body_fails", (False, True), ids=("success", "body-error"))
def test_close_failure_is_logged_without_changing_body_outcome(
    tmp_path, monkeypatch, caplog, body_fails
):
    close = os.close
    closed = []

    def fail_after_close(descriptor):
        close(descriptor)
        closed.append(descriptor)
        raise OSError("TASK19-CLOSE")

    monkeypatch.setattr(germline.os, "close", fail_after_close)
    body_error = RuntimeError("TASK19-BODY")

    with caplog.at_level("WARNING", logger=germline.__name__):
        if body_fails:
            with pytest.raises(RuntimeError) as captured:
                with germline.database_build_lock(str(tmp_path), "custom"):
                    raise body_error
            assert captured.value is body_error
        else:
            with germline.database_build_lock(str(tmp_path), "custom"):
                pass

    assert len(closed) == 1
    assert "TASK19-CLOSE" in caplog.text
    assert ".custom.lock" in caplog.text


def test_header_only_staged_fasta_record_is_rejected(
    tmp_path, monkeypatch, custom_genes, fake_mmseqs
):
    database_root = tmp_path / ".abstar" / "germline_dbs" / "bcr"
    database = database_root / "custom"
    validate = germline.validate_staged_database

    def remove_sequence_content(staging_dir, receptor, **kwargs):
        for directory in ("imgt_gapped", "ungapped"):
            (Path(staging_dir) / directory / "j.fasta").write_text(">IGHJTEST*01\n")
        validate(staging_dir, receptor, **kwargs)

    monkeypatch.setattr(germline, "validate_staged_database", remove_sequence_content)

    with pytest.raises(ValueError, match="empty nucleotide"):
        germline.build_germline_database(
            "custom", fastas=str(custom_genes), include_species_in_name=False,
            verbose=False,
        )

    assert not database.exists()
    assert not list(database_root.glob(".custom.staging-*"))


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
    scratch = Path(calls[1][0][3])
    assert scratch.parent == tmp_path
    assert scratch.name.startswith(".v.index-")
    assert not scratch.exists()
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
