import json
import shutil
import tomllib
from pathlib import Path

import pytest

import abstar
from abstar import _core
from abstar.errors import AbstarError
from abstar import GenerateSummary
from abstar.model import info as model_info, list as list_models, resolve as resolve_model

BUILTIN_MODELS = Path(abstar.__file__).parent / "models"


def _names():
    return [f["name"] for f in tomllib.loads(_core.truth_toml())["field"]]


def _code(fn):
    with pytest.raises(AbstarError) as e:
        fn()
    return e.value.code


def test_public_api_layout():
    import inspect
    import types

    assert isinstance(abstar.generate, types.FunctionType)
    params = inspect.signature(abstar.generate).parameters
    assert list(params)[:3] == ["n", "locus", "germline_db"] and "output_dir" in params
    assert abstar.generate.__doc__.startswith("Simulate `n` V(D)J reads")
    assert isinstance(abstar.model, types.ModuleType) and abstar.model.ModelInfo is type(model_info("human-igh"))
    assert abstar.GenerateSummary is GenerateSummary
    assert not hasattr(abstar.model, "generate")


def test_in_memory_columns_and_rows():
    df = abstar.generate(100, "IGH", seed=1)
    assert df.height == 100
    assert df.columns == _names()
    assert df["sequence_id"][0] == "g001" and set(df["locus"]) == {"IGH"}


def test_same_seed_same_frame():
    a = abstar.generate(50, "IGK", seed=7)
    assert a.equals(abstar.generate(50, "IGK", seed=7))
    assert not a.equals(abstar.generate(50, "IGK", seed=8))


def test_bnab_preset_mutation():
    m = abstar.generate(200, "IGH", seed=3, preset="bnab")["v_mutation_frequency"].mean()
    assert 0.15 <= m <= 0.35, m


def test_explicit_option_overrides_preset():
    df = abstar.generate(50, "IGH", seed=3, preset="bnab", mutation=0.0, indels=0.0)
    assert df["v_mutation_frequency"].max() == 0.0


def test_unknown_option_or_preset():
    assert _code(lambda: abstar.generate(5, "IGH", seed=1, bogus=1)) == "INVALID_ARGUMENT"
    assert _code(lambda: abstar.generate(5, "IGH", seed=1, preset="nope")) == "INVALID_ARGUMENT"
    assert _code(lambda: abstar.generate(5, "IGH", seed=1, mutation="weird:1")) == "INVALID_ARGUMENT"
    assert _code(lambda: abstar.generate(5, "IGH", seed=1, productive_only="yes")) == "INVALID_ARGUMENT"


def test_mutating_a_non_mutating_model_is_invalid():
    # human-tra has shm_model = "none" (G4: surfaced from Simulator::new as-is)
    assert _code(lambda: abstar.generate(5, "TRA", seed=1, mutation=0.1)) == "INVALID_ARGUMENT"
    df = abstar.generate(5, "TRA", seed=1, mutation=0.1, shm_model="motif")
    assert df["v_mutation_frequency"].max() > 0


def test_mixed_database_needs_model():
    assert _code(lambda: abstar.generate(5, "IGH", germline_db="human+c57bl6", seed=1)) == "INVALID_ARGUMENT"
    assert abstar.generate(5, "IGH", germline_db="human+c57bl6", model="human-igh", seed=1).height == 5


def test_model_not_found_and_locus_mismatch():
    assert _code(lambda: abstar.generate(5, "IGH", model="nope", seed=1)) == "MODEL_NOT_FOUND"
    assert _code(lambda: abstar.generate(5, "IGH", model="human-igk", seed=1)) == "INVALID_ARGUMENT"


def test_output_dir(tmp_path):
    s = abstar.generate(20, "IGL", seed=2, output_dir=tmp_path / "out", formats=["airr", "parquet"])
    assert isinstance(s, GenerateSummary)
    assert s.n == 20 and s.seed == 2 and s.productive <= s.attempts
    assert s.fasta == tmp_path / "out" / "generated.fasta"
    assert s.truth == {"airr": tmp_path / "out" / "generated.truth.tsv", "parquet": tmp_path / "out" / "generated.truth.parquet"}
    assert s.generate_json == tmp_path / "out" / "generate.json"
    for p in [s.fasta, s.generate_json, *s.truth.values()]:
        assert p.is_file()
    assert s.fasta.read_text().count(">") == 20
    doc = json.loads(s.generate_json.read_text())
    assert doc["model"]["name"] == "human-igl" and doc["database"]["name"] == "human"
    assert doc["parameters"]["preset"] is None and doc["seed"] == 2


def test_seed_drawn_and_recorded(tmp_path):
    s = abstar.generate(3, "IGH", output_dir=tmp_path, preset="memory")
    doc = json.loads(s.generate_json.read_text())
    assert doc["seed"] == s.seed and doc["parameters"]["preset"] == "memory"
    assert abstar.generate(3, "IGH", seed=s.seed, preset="memory")["sequence"].to_list() == [
        l for l in s.fasta.read_text().splitlines() if not l.startswith(">")
    ]


def test_bad_arguments():
    assert _code(lambda: abstar.generate(5, "IGH", formats=["csv"], output_dir="x")) == "INVALID_ARGUMENT"
    assert _code(lambda: abstar.generate(5, "IGH", threads=0)) == "INVALID_ARGUMENT"
    assert _code(lambda: abstar.generate(-1, "IGH")) == "INVALID_ARGUMENT"


def test_models_listed_and_resolved(abstar_home):
    names = [m.name for m in list_models()]
    assert names[:5] == ["human-igh", "human-igk", "human-igl", "human-tra", "human-trb"]
    m = model_info("human-igh")
    assert (m.kind, m.species, m.locus) == ("builtin", "human", "IGH") and m.sources
    assert m.path == BUILTIN_MODELS / "human-igh" / "model.toml"
    assert resolve_model(None, "human", "TRB") == (BUILTIN_MODELS / "human-trb" / "model.toml", "builtin")

    user = abstar_home / "models" / "mine"
    user.mkdir(parents=True)
    shutil.copy(BUILTIN_MODELS / "human-igk" / "model.toml", user / "model.toml")
    assert model_info("mine").kind == "user"
    assert model_info(user).kind == "path" and model_info(user / "model.toml").kind == "path"
    assert abstar.generate(3, "IGK", model="mine", seed=1).height == 3
    assert "mine" in [m.name for m in list_models()]


def test_user_model_may_not_shadow_builtin(abstar_home):
    d = abstar_home / "models" / "human-igh"
    d.mkdir(parents=True)
    shutil.copy(BUILTIN_MODELS / "human-igh" / "model.toml", d / "model.toml")
    assert _code(lambda: abstar.generate(3, "IGH", seed=1)) == "BUILTIN_NAME_CONFLICT"


def test_bnab_is_productive_only():
    df = abstar.generate(30, "IGH", seed=3, preset="bnab")
    assert df["productive"].all()
    assert not abstar.generate(30, "IGH", seed=3, preset="bnab", productive_only=False)["productive"].all()


def test_core_error_message_has_no_doubled_prefix():
    with pytest.raises(AbstarError) as e:
        abstar.generate(5, "IGH", seed=1, mutation="weird:1")
    assert e.value.message.startswith("bad distribution"), e.value.message


def test_invalid_model_file(tmp_path):
    bad = tmp_path / "model.toml"
    bad.write_text("not = [toml")
    assert _code(lambda: model_info(bad)) == "MODEL_INVALID"
    assert _code(lambda: abstar.generate(3, "IGH", model=bad, seed=1)) == "MODEL_INVALID"
