# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

from sphinx.ext import autodoc

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "tidesurf"
copyright = "2025, Jan Schleicher"
author = "Jan Schleicher"
release = "0.1.0"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx_copybutton",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
]

templates_path = ["_templates"]
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_book_theme"
html_theme_options = dict(
    use_repository_button=True,
    repository_url="https://github.com/janschleicher/tidesurf",
    repository_branch="main",
)
html_static_path = ["_static"]

# Automatic generation of API documentation
autosummary_generate = True
autodoc_default_options = {
    "member-order": "alphabetical",
    "show-inheritance": True,
}


def setup(app):
    def skip_member(app, what, name, obj, skip, options):
        # exclude attributes and methods added to enum by cython
        if name in [
            "real",
            "imag",
            "numerator",
            "denominator",
            "conjugate",
            "bit_length",
            "bit_count",
            "to_bytes",
            "from_bytes",
            "as_integer_ratio",
        ]:
            return True
        return None

    app.connect("autodoc-skip-member", skip_member)


class MockedClassDocumenter(autodoc.ClassDocumenter):
    def add_line(self, line: str, source: str, *lineno: int) -> None:
        if line == "   Bases: :py:class:`object`":
            return
        super().add_line(line, source, *lineno)


autodoc.ClassDocumenter = MockedClassDocumenter
