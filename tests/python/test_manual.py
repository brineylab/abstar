import re
import tomllib

import typer
from typer.testing import CliRunner

from abstar import _core, cli
from abstar.errors import ERROR_CODES, EXIT_CODES
from abstar.manual import render

runner = CliRunner()


def _walk(cmd):
    yield cmd
    for sub in getattr(cmd, "commands", {}).values():
        yield from _walk(sub)


def test_every_command_and_option_documented():
    text = render("commands")
    root = typer.main.get_command(cli.app)
    for cmd in _walk(root):
        if cmd is not root:
            assert cmd.name in text
        for p in cmd.params:
            if p.param_type_name == "argument":
                assert p.name.upper() in text  # an argument's `opts` is its internal name
            else:
                for o in p.opts + p.secondary_opts:
                    assert o in text, o


def test_every_field_documented():
    text = render("fields")
    for f in tomllib.loads(_core.fields_toml())["field"]:
        assert f["name"] in text


def test_field_with_decision():
    text = render("fields", "v_call_ties")
    assert "0008" in text
    assert "`v_call` is exactly one allele, always." in text
    assert "## Decision" in text


def test_errors_topic():
    text = render("errors")
    for code in ERROR_CODES:
        assert code in text
    for n in EXIT_CODES:
        assert re.search(rf"\| {n} \|", text)


def test_decisions_topic():
    assert "0008" in render("decisions")
    assert render("decisions", "0008").startswith("# 0008:")


def test_unknown_topic_exit(monkeypatch, capsys):
    import pytest

    monkeypatch.setattr("sys.argv", ["abstar", "manual", "bogus"])
    with pytest.raises(SystemExit) as e:
        cli.main()
    assert e.value.code == 2
    assert "fields" in capsys.readouterr().err


def test_unknown_subtopic():
    import pytest
    from abstar.errors import AbstarError

    with pytest.raises(AbstarError) as e:
        render("fields", "nope")
    assert e.value.code == "INVALID_ARGUMENT"


def test_cli_manual_stdout():
    r = runner.invoke(cli.app, ["manual"])
    assert r.exit_code == 0
    assert r.stdout.startswith("# abstar")
