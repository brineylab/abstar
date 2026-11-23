import json
import subprocess
import sys
from pathlib import Path

from test_germline_api import recipe
from test_run_api import IGH_READ

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


def test_seq_real_read_is_annotated():
    p = abstar("run", "--seq", IGH_READ)
    assert p.returncode == 0, p.stderr
    head, row = (line.split("\t") for line in p.stdout.splitlines())
    rec = dict(zip(head, row))
    assert (rec["annotation_status"], rec["failure_reason"]) == ("annotated", "")
    assert rec["v_call"].startswith("IGHV1-2") and rec["j_call"] == "IGHJ4*02"
    assert (rec["rev_comp"], rec["productive"], rec["junction_aa"]) == ("F", "T", "CARDRITIFGVVIIPDYFDYW")


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


def _donor(tmp_path):
    abstar("germline", "build", recipe(tmp_path, "base"))
    d = tmp_path / "donor"
    d.mkdir()
    (d / "genotype.txt").write_text("IGHJ6*03\nIGHJ9*01\n")
    (d / "j.fasta").write_text(">IGHJ9*01\nACGTACGTACGTACGTACGT\n")
    (d / "donor.toml").write_text(
        'name = "donor"\nkeep = "genotype.txt"\n[[source]]\ndatabase = "base"\n'
        '[[source]]\nspecies = "human"\nj_fasta = ["j.fasta"]\n'
    )
    return d / "donor.toml"


def test_build_human_output(tmp_path):
    p = abstar("germline", "build", _donor(tmp_path))
    assert p.returncode == 0, p.stderr
    lines = p.stdout.splitlines()
    assert "restricted: IGH/J" in lines
    assert "passed through: IGH/V, IGH/C" in lines
    assert "warning: IGHJ9*01: no J Trp/Phe118 motif" in p.stderr
    assert "{" not in p.stderr


def test_list_survives_corrupt_user_database(abstar_home, tmp_path):
    abstar("germline", "build", recipe(tmp_path, "good"))
    bad = abstar_home / "germline" / "bad"
    bad.mkdir(parents=True)
    (bad / "manifest.toml").write_text("not = [toml")
    for leftover in ("good.staging-123", "good.old-456"):
        (abstar_home / "germline" / leftover).mkdir()
        (abstar_home / "germline" / leftover / "manifest.toml").write_text("")
    p = abstar("germline", "list", "--json")
    assert p.returncode == 0, p.stderr
    users = {d["name"]: d for d in json.loads(p.stdout) if d["kind"] == "user"}
    assert sorted(users) == ["bad", "good"]
    assert users["good"]["error"] is None and users["good"]["n_alleles"] > 0
    assert users["bad"]["n_alleles"] == 0 and "DATABASE_INVALID" in users["bad"]["error"]
    assert str(bad) in users["bad"]["error"]
    p = abstar("germline", "list")
    assert p.returncode == 0 and "bad" in p.stdout and "invalid" in p.stdout
    assert abstar("germline", "info", "good.staging-123").returncode == 2


def test_info_shows_manifest(tmp_path):
    abstar("germline", "build", _donor(tmp_path))
    d = json.loads(abstar("germline", "info", "donor", "--json").stdout)
    for k in ("built", "abstar_version", "alleles_sha256", "prefixes", "provenance", "warnings"):
        assert k in d, k
    assert len(d["alleles_sha256"]) == 64 and d["provenance"] and d["warnings"]
    p = abstar("germline", "info", "donor")
    assert p.returncode == 0
    for s in (d["built"], d["abstar_version"], d["alleles_sha256"], "provenance", "IGHJ9*01: no J Trp/Phe118 motif"):
        assert s in p.stdout, s


def _closed_stdout(*args):
    """Run abstar with a stdout pipe whose reader is already gone (`| head` that exited)."""
    import os

    r, w = os.pipe()
    os.close(r)
    try:
        return subprocess.run([ABSTAR, *map(str, args)], stdout=w, stderr=subprocess.PIPE, text=True)
    finally:
        os.close(w)


def test_closed_stdout_pipe(tmp_path):
    fa = tmp_path / "in.fasta"
    fa.write_text(FA)
    for args in (["manual"], ["germline", "list"], ["run", fa, tmp_path / "out"]):
        p = _closed_stdout(*args)
        assert p.returncode == 0, (args, p.returncode, p.stderr)
        assert "Traceback" not in p.stderr and "BrokenPipe" not in p.stderr, (args, p.stderr)


def test_closed_stdout_pipe_keeps_error_exit_code():
    # #142: buffered stdout plus an error; flushing into the closed pipe must not turn exit 2 into 0.
    import os

    code = (
        "import sys; from abstar import cli; print('partial output');"
        "sys.argv = ['abstar', 'germline', 'info', 'nope']; cli.main()"
    )
    r, w = os.pipe()
    os.close(r)
    try:
        p = subprocess.run([sys.executable, "-c", code], stdout=w, stderr=subprocess.PIPE, text=True)
    finally:
        os.close(w)
    assert p.returncode == 2, p.stderr
    assert p.stderr.startswith("error[DATABASE_NOT_FOUND]: ") and "Traceback" not in p.stderr

def test_unknown_locus_is_invalid_argument():
    p = abstar("run", "--seq", "ACGT", "--locus", "TRZ")
    assert p.returncode == 2 and p.stderr.startswith("error[INVALID_ARGUMENT]: ")
    assert "TRZ" in p.stderr and "IGH" in p.stderr and "germline database" not in p.stderr


def test_output_dir_is_a_file(tmp_path):
    fa = tmp_path / "in.fasta"
    fa.write_text(FA)
    p = abstar("run", fa, fa)
    assert p.returncode == 2 and p.stderr.startswith("error[INVALID_ARGUMENT]: "), p.stderr
    assert fa.read_text() == FA


def test_include_nonfunctional_flags(tmp_path):
    from test_run_api import kv229_read

    def v_call(*flags):
        p = abstar("run", "--seq", kv229_read(), *flags)
        assert p.returncode == 0, p.stderr
        head, row = (line.split("\t") for line in p.stdout.splitlines())
        return dict(zip(head, row))["v_call"]

    assert v_call() != "IGKV2-29*01" and v_call("--include-nonfunctional") == "IGKV2-29*01"
    p = abstar("run", "--seq", kv229_read(), tmp_path / "run", "--include-nonfunctional", "-q")
    assert p.returncode == 0, p.stderr
    assert json.loads((tmp_path / "run" / "run.json").read_text())["parameters"]["include_nonfunctional"] is True
    p = abstar("generate", tmp_path / "gen", "--n", "3", "--locus", "IGH", "--seed", "1", "--include-nonfunctional", "-q")
    assert p.returncode == 0, p.stderr
    assert json.loads((tmp_path / "gen" / "generate.json").read_text())["parameters"]["include_nonfunctional"] is True


def test_run_sequence_type_aa():
    from test_run_api import IGH_AA

    p = abstar("run", "--seq", IGH_AA, "--sequence-type", "aa")
    assert p.returncode == 0, p.stderr
    head, row = (line.split("\t") for line in p.stdout.splitlines())
    assert dict(zip(head, row))["v_gene"] == "IGHV1-2"
