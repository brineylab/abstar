import json
import subprocess
import tomllib

from test_cli import ABSTAR, _closed_stdout, abstar

from abstar import _core
from abstar.errors import ERROR_CODES
from abstar.manual import render


def test_fasta_to_stdout():
    p = abstar("generate", "--n", 3, "--locus", "IGH", "--seed", 1)
    assert p.returncode == 0, p.stderr
    lines = p.stdout.splitlines()
    assert len(lines) == 6
    assert lines[0] == ">g1" and lines[2] == ">g2" and set(lines[1]) <= set("ACGT")
    assert "n=3" in p.stderr and "seed=1" in p.stderr
    assert abstar("generate", "--n", 3, "--locus", "IGH", "--seed", 1).stdout == p.stdout
    q = abstar("generate", "--n", 3, "--locus", "IGH", "--seed", 1, "-q")
    assert q.returncode == 0 and q.stderr == "" and q.stdout == p.stdout


def test_closed_pipe():
    p = subprocess.Popen([ABSTAR, "generate", "--n", "100000", "--locus", "IGH"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    assert p.stdout.readline().startswith(b">g")
    p.stdout.close()
    err = p.stderr.read().decode()
    assert p.wait() == 0, err
    assert "Traceback" not in err
    p = _closed_stdout("generate", "--n", "1000", "--locus", "IGH")
    assert p.returncode == 0 and "Traceback" not in p.stderr, p.stderr


def test_output_dir_rerun_overwrites(tmp_path):
    out = tmp_path / "out"
    p = abstar("generate", out, "--n", 5, "--locus", "IGK", "--seed", 1, "--json")
    assert p.returncode == 0, p.stderr
    s = json.loads(p.stdout)
    assert s["n"] == 5 and s["seed"] == 1 and s["truth"] == {"airr": str(out / "generated.truth.tsv")}
    first = (out / "generated.fasta").read_text()
    p = abstar("generate", out, "--n", 2, "--locus", "IGK", "--seed", 2, "--format", "airr", "--format", "parquet")
    assert p.returncode == 0, p.stderr
    assert "generated.fasta" in p.stdout and p.stderr == ""
    assert (out / "generated.fasta").read_text().count(">") == 2 != first.count(">")
    assert (out / "generated.truth.parquet").is_file()
    assert len((out / "generated.truth.tsv").read_text().splitlines()) == 3
    assert json.loads((out / "generate.json").read_text())["seed"] == 2


def test_options_reach_the_generator(tmp_path):
    out = tmp_path / "out"
    p = abstar(
        "generate", out, "--n", 4, "--locus", "IGH", "--seed", 1, "--preset", "bnab", "--mutation", "0.05",
        "--long-cdr3", "0", "--indels", "0", "--truncate-5p", "uniform:0-5", "--truncate-3p", "2",
        "--constant", 30, "--productive-only", "--frameshift-indels", "--name", "x", "--threads", 2,
        "--model", "human-igh", "-g", "human", "-q",
    )
    assert p.returncode == 0, p.stderr
    assert p.stdout == "" and p.stderr == ""
    prm = json.loads((out / "generate.json").read_text())["parameters"]
    assert prm["preset"] == "bnab" and prm["mutation"] == "fixed:0.05" and prm["indel_rate"] == 0
    assert prm["truncate_5p"] == "uniform:0-5" and prm["truncate_3p"] == "fixed:2"
    assert prm["constant_len"] == 30 and prm["productive_only"] and prm["frameshift_indels"] and prm["threads"] == 2
    assert (out / "x.fasta").is_file()


def test_shm_table_option(tmp_path):
    p = abstar("generate", "--n", 2, "--locus", "IGH", "--seed", 1, "--shm-table", tmp_path / "missing")
    assert p.returncode == 2 and p.stderr.startswith("error[INVALID_ARGUMENT]: "), p.stderr


def test_errors():
    p = abstar("generate", "--n", 2, "--locus", "IGH", "--model", "nope", "--json")
    assert p.returncode == 2 and json.loads(p.stderr)["error"]["code"] == "MODEL_NOT_FOUND"
    p = abstar("generate", "--n", 2, "--locus", "IGH", "--format", "parquet")
    assert p.returncode == 2 and "error[INVALID_ARGUMENT]" in p.stderr
    p = abstar("generate", "--n", 2, "--locus", "IGH", "-g", "human+c57bl6")
    assert p.returncode == 2 and "error[INVALID_ARGUMENT]" in p.stderr


def test_model_commands():
    p = abstar("model", "list", "--json")
    assert p.returncode == 0, p.stderr
    models = json.loads(p.stdout)
    assert any(m["name"] == "human-igh" and m["kind"] == "builtin" for m in models)
    p = abstar("model", "list")
    assert p.returncode == 0 and "human-trb\tbuiltin\thuman\tTRB" in p.stdout
    p = abstar("model", "info", "human-igh")
    assert p.returncode == 0 and "IGH" in p.stdout and "Marcou" in p.stdout
    assert json.loads(abstar("model", "info", "human-igk", "--json").stdout)["locus"] == "IGK"
    p = abstar("model", "info", "nope")
    assert p.returncode == 2 and p.stderr.startswith("error[MODEL_NOT_FOUND]: ")


def test_manual_topics():
    text = render("truth-fields")
    for f in tomllib.loads(_core.truth_toml())["field"]:
        assert f"`{f['name']}`" in text
    assert "0025" in render("truth-fields", "v_call_valid")
    g = render("generate")
    for s in ("bnab", "`abstar generate`", "`abstar model list`", "--shm-table"):
        assert s in g, s
    for s in ("abstar generate", "abstar.generate(", "truth-fields"):
        assert s in render("examples"), s
    assert abstar("manual", "truth-fields").returncode == 0


def test_error_codes():
    assert ERROR_CODES["MODEL_NOT_FOUND"] == (2, "the recombination model does not exist")
    assert ERROR_CODES["MODEL_INVALID"] == (2, "the recombination model is malformed or does not fit the germline database")
