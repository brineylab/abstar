import pytest


@pytest.fixture(autouse=True)
def abstar_home(tmp_path, monkeypatch):
    """Redirect ABSTAR_HOME so no test touches the real ~/.abstar."""
    home = tmp_path / "home"
    monkeypatch.setenv("ABSTAR_HOME", str(home))
    return home
