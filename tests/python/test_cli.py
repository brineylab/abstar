import json
import subprocess
import sys
from pathlib import Path

from test_germline_api import recipe

ABSTAR = str(Path(sys.executable).parent / "abstar")
FA = ">r1\nACGTACGTACGTACGT\n>r2\nGGGGTTTTAAAACCCC\n>r3\nACACACACGTGTGTGT\n"


def abstar(*args, stdin=None):
    return subprocess.run([ABSTAR, *map(str, args)], input=stdin, capture_output=True, text=True)


def test_seq_to_stdout():
    p = abstar("run", "--seq", "ACGT")
    assert p.returncode == 0, p.stderr
    lines = p.stdout.splitlines()
    assert len(lines) == 2 and lines[0].startswith("sequence_id\t")
    assert lines[1].split("\t")[0] == "1"
    assert "annotated" in p.stderr  # summary on stderr


def test_quiet_suppresses_summary():
    p = abstar("run", "--seq", "ACGT", "-q")
    assert p.returncode == 0 and p.stderr == ""


def test_stdin():
    p = abstar("run", "-", stdin=FA)
    assert p.returncode == 0, p.stderr
    assert len(p.stdout.splitlines()) == 4


def test_broken_pipe(tmp_path):
    fa = tmp_path / "many.fasta"
    fa.write_text("".join(f">r{i}\nACGTACGTACGT\n" for i in range(10_000)))
    p = subprocess.Popen([ABSTAR, "run", str(fa), "--threads", "1"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    assert p.stdout.readline().startswith(b"sequence_id\t")
    p.stdout.close()
    err = p.stderr.read().decode()
    assert p.wait() == 0, err
    assert "Traceback" not in err


def test_output_dir_json_summary(tmp_path):
    fa = tmp_path / "in.fasta"
    fa.write_text(FA)
    p = abstar("run", fa, tmp_path / "out", "--json", "--format", "airr", "--format", "parquet")
    assert p.returncode == 0, p.stderr
    d = json.loads(p.stdout)
    assert sum(d["counts"].values()) == 3
    for path in d["inputs"][0]["outputs"].values():
        assert Path(path).is_file()


def test_output_dir_human_summary(tmp_path):
    fa = tmp_path / "in.fasta"
    fa.write_text(FA)
    p = abstar("run", fa, tmp_path / "out")
    assert p.returncode == 0 and p.stdout.strip() and p.stderr == ""


def test_parquet_without_output_dir():
    p = abstar("run", "--seq", "ACGT", "--format", "parquet")
    assert p.returncode == 2
    assert "error[INVALID_ARGUMENT]" in p.stderr


def test_json_error_shape(tmp_path):
    p = abstar("run", tmp_path / "missing.fasta", tmp_path / "out", "--json")
    assert p.returncode == 2
    assert json.loads(p.stderr)["error"]["code"] == "INPUT_NOT_FOUND"
    assert p.stdout == ""


def test_text_error(tmp_path):
    p = abstar("run", tmp_path / "missing.fasta")
    assert p.returncode == 2 and p.stderr.startswith("error[INPUT_NOT_FOUND]: ")
    assert "Traceback" not in p.stderr


def test_germline_commands(tmp_path):
    out = json.loads(abstar("germline", "list", "--json").stdout)
    assert any(d["name"] == "human" and d["kind"] == "builtin" for d in out)
    info = json.loads(abstar("germline", "info", "human", "--json").stdout)
    assert info["counts"]["IGH/V/human"] > 0
    rep = json.loads(abstar("germline", "build", recipe(tmp_path, "mydb"), "--json").stdout)
    assert rep["name"] == "mydb" and rep["n_alleles"] > 0
    assert abstar("germline", "list").returncode == 0
    assert abstar("germline", "info", "nope").returncode == 2


def test_no_seq_and_no_input():
    assert abstar("run").returncode == 2


def test_seq_with_output_dir(tmp_path):
    p = abstar("run", "--seq", "ACGT", tmp_path / "out")
    assert p.returncode == 0, p.stderr
    assert (tmp_path / "out" / "records.tsv").is_file()


def test_bad_threads():
    p = abstar("run", "--seq", "ACGT", "--threads", "-1")
    assert p.returncode == 2 and "INVALID_ARGUMENT" in p.stderr


def test_verbose_and_version():
    p = abstar("run", "--seq", "ACGT", "-v")
    assert "human" in p.stderr and p.stdout.count("\n") == 2
    assert abstar("--version").stdout.strip()
    assert abstar("manual").returncode == 0


def test_help_is_light():
    code = "import sys; from abstar.cli import app; from typer.testing import CliRunner; " \
           "r = CliRunner().invoke(app, ['--help']); assert r.exit_code == 0; " \
           "bad = [m for m in ('polars', 'abstar._core') if m in sys.modules]; assert not bad, bad"
    p = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert p.returncode == 0, p.stderr


def test_seq_ambiguous_positionals_rejected(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    fa = tmp_path / "reads.fasta"
    fa.write_text(FA)
    for args in (["run", fa, "--seq", "ACGT"], ["run", "a", "b", "--seq", "ACGT"]):
        p = abstar(*args)
        assert p.returncode == 2 and "error[INVALID_ARGUMENT]" in p.stderr and "mutually exclusive" in p.stderr
    assert not (tmp_path / "a").exists() and not (tmp_path / "b").exists() and fa.is_file()


def test_api_run_dash_reads_stdin():
    p = subprocess.run([sys.executable, "-c", "import abstar; print(abstar.run('-').height)"],
                       input=FA, capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    assert p.stdout.strip() == "3"


def test_build_overwrite_refuses_non_database_dir(tmp_path):
    victim = tmp_path / "victim"
    victim.mkdir()
    (victim / "thesis.docx").write_text("precious")
    p = abstar("germline", "build", recipe(tmp_path, "mydb"), "--output", victim, "--overwrite")
    assert p.returncode == 2 and "error[DATABASE_EXISTS]" in p.stderr and "not a germline database" in p.stderr
    assert (victim / "thesis.docx").read_text() == "precious"


def test_positional_is_always_a_path(tmp_path):
    for args in (["run", "ACGT"], ["run", "sample1"], ["run", "reads.txt", tmp_path / "out1"]):
        p = abstar(*args)
        assert p.returncode == 2 and p.stderr.startswith("error[INPUT_NOT_FOUND]: "), (args, p.stderr)
        assert p.stdout == ""
    assert not (tmp_path / "out1").exists()
