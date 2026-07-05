"""Sphinx configuration for bioseqkit."""

project = "bioseqkit"
author = "Jilai Cheng"
release = "0.1.0"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "myst_parser",
]

autodoc_typehints = "description"
templates_path = ["_templates"]
exclude_patterns = ["_build"]

html_theme = "alabaster"
source_suffix = {".rst": "restructuredtext", ".md": "markdown"}
