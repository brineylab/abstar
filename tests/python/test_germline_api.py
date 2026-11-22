import shutil
from pathlib import Path

import pytest

from abstar import errors, germline
from abstar.errors import AbstarError

DATA = Path(__file__).parent / "data"


def recipe(tmp_path, name, extra=""):
    """The fixture recipe and its data, copied to tmp_path under a new name."""
    d = tmp_path / f"recipe-{name}"
    shutil.copytree(DATA, d, dirs_exist_ok=True)
    text = (DATA / "ogrdb.toml").read_text().replace('name = "ogrdb-test"', f'name = "{name}"{extra}')
    (d / "ogrdb.toml").write_text(text)
    return d / "ogrdb.toml"


@pytest.fixture
def builtin_dir(tmp_path, monkeypatch):
    """A stand-in built-in directory holding a database named `human`."""
    d = tmp_path / "builtins"
    germline.build(recipe(tmp_path, "human"), output_dir=d / "human")
    monkeypatch.setattr(germline, "_BUILTIN_DIR", d)
    return d


def test_build_to_user_home(tmp_path, abstar_home):
    report = germline.build(recipe(tmp_path, "mydb"))
    path = tmp_path / "home" / "germline" / "mydb"
    assert report.path == path
    assert (path / "manifest.toml").is_file()
    assert report.name == "mydb"
    assert report.n_alleles == sum(report.counts.values()) > 0
    info = germline.info("mydb")
    assert info.kind == "user"
    assert info.path == path
    assert info.counts == report.counts
    assert info.species == ["human"]
    assert info.loci == ["IGH"]
    assert [(d.name, d.kind) for d in germline.list() if d.kind == "user"] == [("mydb", "user")]


def test_builtin_name_conflict(tmp_path, abstar_home, builtin_dir):
    assert germline.builtin_names() == ["human"]
    with pytest.raises(AbstarError) as e:
        germline.build(recipe(tmp_path, "human"))
    assert e.value.code == "BUILTIN_NAME_CONFLICT"
    assert e.value.exit_code == 2
    assert not abstar_home.exists()


def test_builtin_names_missing_dir(tmp_path, monkeypatch):
    monkeypatch.setattr(germline, "_BUILTIN_DIR", tmp_path / "nope")
    assert germline.builtin_names() == []


def test_exists_then_overwrite(tmp_path):
    r = recipe(tmp_path, "mydb")
    germline.build(r)
    with pytest.raises(AbstarError) as e:
        germline.build(r)
    assert e.value.code == "DATABASE_EXISTS"
    assert germline.build(r, overwrite=True).name == "mydb"


def test_resolve_order_and_not_found(tmp_path, abstar_home, builtin_dir):
    germline.build(recipe(tmp_path, "mydb"))
    # an explicit output_dir may reuse a built-in name
    germline.build(recipe(tmp_path, "human"), output_dir=abstar_home / "germline" / "human")
    # 1. an existing database directory wins, by path
    assert germline.resolve(builtin_dir / "human") == (builtin_dir / "human", "path")
    # 2. built-in before 3. user
    assert germline.resolve("human") == (builtin_dir / "human", "builtin")
    assert germline.resolve("mydb") == (abstar_home / "germline" / "mydb", "user")
    assert [(d.name, d.kind) for d in germline.list()] == [
        ("human", "builtin"), ("human", "user"), ("mydb", "user"),
    ]
    with pytest.raises(AbstarError) as e:
        germline.resolve("nope")
    assert e.value.code == "DATABASE_NOT_FOUND"
    assert "mydb" in e.value.message and "human" in e.value.message


def test_error_table():
    assert {k: v[0] for k, v in errors.ERROR_CODES.items()} == {
        "INPUT_NOT_FOUND": 2, "INPUT_UNREADABLE": 2, "INVALID_ARGUMENT": 2,
        "OUTPUT_NAME_CONFLICT": 2, "DATABASE_NOT_FOUND": 2, "DATABASE_INVALID": 2,
        "RECIPE_INVALID": 2, "DATABASE_EXISTS": 2, "BUILTIN_NAME_CONFLICT": 2,
        "MODEL_NOT_FOUND": 2, "MODEL_INVALID": 2,
        "NO_ALLELES_GAPPED": 2, "IO_ERROR": 1, "INTERNAL": 1, "STRICT_RECORD_FAILURE": 3,
    }
    assert all(isinstance(v[1], str) and v[1] for v in errors.ERROR_CODES.values())
    assert set(errors.EXIT_CODES) == {0, 1, 2, 3}


def test_database_source_by_name(tmp_path):
    germline.build(recipe(tmp_path, "base"))
    d = tmp_path / "donor"
    d.mkdir()
    (d / "genotype.txt").write_text("IGHJ6*03\n")
    (d / "donor.toml").write_text('name = "donor"\nkeep = "genotype.txt"\n[[source]]\ndatabase = "base"\n')
    report = germline.build(d / "donor.toml")
    assert report.restricted == ["IGH/J"]
    assert "IGH/J" not in report.passed_through and "IGH/V" in report.passed_through
    assert report.counts["IGH/J/human"] == 1


def test_core_error_maps_to_abstar_error(tmp_path):
    from abstar import _core
    with pytest.raises(_core.CoreError) as e:
        _core.germline_summary(str(tmp_path / "missing"))
    err = errors.from_core(e.value)
    assert (err.code, err.exit_code) == ("DATABASE_NOT_FOUND", 2)
    assert str(tmp_path / "missing") in err.message


def test_run_in_memory_is_arrow_stream(tmp_path):
    import polars as pl
    from abstar import _core
    germline.build(recipe(tmp_path, "mydb"))
    db = str(germline.resolve("mydb")[0])
    df = pl.DataFrame(_core.run_in_memory([], [("a", "ACGT"), ("b", "TTTT")], db, 1, False, None))
    assert df["sequence_id"].to_list() == ["a", "b"]
    empty = pl.DataFrame(_core.run_in_memory([], [], db, None, False, None))
    assert empty.height == 0 and "sequence_id" in empty.columns


def test_run_to_files_and_stdout(tmp_path, capfd):
    from abstar import _core
    germline.build(recipe(tmp_path, "mydb"))
    db = str(germline.resolve("mydb")[0])
    fasta = tmp_path / "reads.fasta"
    fasta.write_text(">r1\nACGT\n>r2\nGGGG\n")
    out = _core.run_to_files([str(fasta)], [], str(tmp_path / "out"), db, ["airr"], None, False, ["IGH"])
    assert sum(out["counts"].values()) == 2
    assert Path(out["inputs"][0]["outputs"]["airr"]).is_file()
    counts = _core.run_to_stdout([str(fasta)], [], db, None, False, None)
    assert sum(counts.values()) == 2
    assert capfd.readouterr().out.splitlines()[0].startswith("sequence_id\t")
    with pytest.raises(_core.CoreError) as e:
        _core.run_to_files([str(tmp_path / "nope.fa")], [], str(tmp_path / "o2"), db, ["airr"], None, False, None)
    assert e.value.args[0] == "INPUT_NOT_FOUND"


def test_run_to_stdout_broken_pipe(tmp_path):
    import ast, subprocess, sys
    germline.build(recipe(tmp_path, "mydb"))
    db = str(germline.resolve("mydb")[0])
    fasta = tmp_path / "many.fasta"
    fasta.write_text("".join(f">r{i}\nACGTACGTACGT\n" for i in range(50_000)))
    code = f"from abstar import _core; print(_core.run_to_stdout([{str(fasta)!r}], [], {db!r}, 1, False, None), file=__import__('sys').stderr)"
    p = subprocess.Popen([sys.executable, "-c", code], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    assert p.stdout.readline().startswith(b"sequence_id\t")
    p.stdout.close()
    err = p.stderr.read().decode()
    assert p.wait() == 0, err
    assert "Traceback" not in err
    assert sum(ast.literal_eval(err).values()) < 50_000  # stopped at the closed pipe


def test_recipe_name_traversal_rejected(tmp_path, abstar_home):
    victim = tmp_path / "x"
    victim.mkdir()
    (victim / "keep.txt").write_text("precious")
    r = recipe(tmp_path, "mydb")
    r.write_text(r.read_text().replace('name = "mydb"', 'name = "../x"'))
    before = sorted(tmp_path.rglob("*"))
    for overwrite in (False, True):
        with pytest.raises(AbstarError) as e:
            germline.build(r, overwrite=overwrite)
        assert e.value.code == "RECIPE_INVALID"
        assert "../x" in e.value.message
    assert sorted(tmp_path.rglob("*")) == before
    assert not abstar_home.exists()
    assert (victim / "keep.txt").read_text() == "precious"


def test_build_gap_v_default_templates(tmp_path):
    human = germline.resolve("human")[0] / "alleles.tsv"
    # strip a few shipped human IGHV alleles into a FASTA
    rows = [l.split("\t") for l in human.read_text().splitlines()]
    cols = rows[0]
    fa = ""
    for r in rows[1:]:
        d = dict(zip(cols, r))
        if d["id"] in ("IGHV1-2*02", "IGHV3-23*01"):
            fa += f">{d['id']}\n{d['sequence_gapped'].replace('.', '')}\n"
    assert fa.count(">") == 2
    (tmp_path / "v.fasta").write_text(fa)
    (tmp_path / "r.toml").write_text(
        'name = "gapped"\n[[source]]\nspecies = "human"\ngap_v = true\nv_fasta = ["v.fasta"]\n'
    )
    germline.build(tmp_path / "r.toml", output_dir=tmp_path / "out")
    files = germline.info(tmp_path / "out").provenance[0]["files"]
    for t in ("human", "c57bl6", "balbc"):
        assert f"gap-template:{t}" in files


def test_local_directory_does_not_shadow_builtin(tmp_path, monkeypatch, abstar_home):
    # #141 / 0010: a ./human database directory is used only when named as a path.
    work = tmp_path / "work"
    germline.build(recipe(tmp_path, "local"), output_dir=work / "human")
    monkeypatch.chdir(work)
    assert germline.resolve("human")[1] == "builtin"
    assert germline.resolve("./human") == (Path("human"), "path")
    assert germline.resolve(Path("human")) == (Path("human"), "path")
    # a bare name that is no built-in or user database still finds a local directory
    (work / "human").rename(work / "mine")
    assert germline.resolve("mine") == (Path("mine"), "path")
