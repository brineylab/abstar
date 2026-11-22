"""Sphinx config. The manual and decision records are generated from the real sources at build time."""
from pathlib import Path

import abstar
from abstar import manual

project = "abstar"
release = abstar.__version__
extensions = ["myst_parser"]
html_theme = "furo"
exclude_patterns = ["_build", "superpowers/**"]

_gen = Path(__file__).parent / "_generated"
_gen.mkdir(exist_ok=True)
(_gen / "manual.md").write_text(manual.render())
(_gen / "decisions.md").write_text(
    "\n\n".join(manual.render("decisions", i) for i in manual._decision_files())
)
