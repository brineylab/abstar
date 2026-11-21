from importlib.metadata import version as _version
from typing import TYPE_CHECKING

__version__ = _version("abstar")

__all__ = ["run", "RunSummary", "InputSummary", "AbstarError", "germline", "generate", "GenerateSummary", "grade", "model"]

if TYPE_CHECKING:
    from . import germline, model
    from ._generate import GenerateSummary, generate
    from .grade import grade
    from .api import InputSummary, RunSummary, run
    from .errors import AbstarError


def __getattr__(name):
    # Lazy so `import abstar` loads no polars, typer or abstar._core.
    if name in ("run", "RunSummary", "InputSummary"):
        from . import api

        return getattr(api, name)
    if name in ("generate", "GenerateSummary"):
        from . import _generate

        return getattr(_generate, name)
    if name == "grade":
        from .grade import grade

        globals()["grade"] = grade  # the submodule import rebinds the name; keep the function
        return grade
    if name == "AbstarError":
        from .errors import AbstarError

        return AbstarError
    if name in ("germline", "model"):
        import importlib

        return importlib.import_module(f".{name}", __name__)
    raise AttributeError(f"module 'abstar' has no attribute {name!r}")
