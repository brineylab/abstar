from importlib.metadata import version as _version

__version__ = _version("abstar")

__all__ = ["run", "RunSummary", "InputSummary", "AbstarError", "germline"]


def __getattr__(name):
    # Lazy so `import abstar` loads no polars, typer or abstar._core.
    if name in ("run", "RunSummary", "InputSummary"):
        from . import api

        return getattr(api, name)
    if name == "AbstarError":
        from .errors import AbstarError

        return AbstarError
    if name == "germline":
        import importlib

        return importlib.import_module(".germline", __name__)
    raise AttributeError(f"module 'abstar' has no attribute {name!r}")
