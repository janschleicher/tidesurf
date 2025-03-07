# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

from sphinx.ext import autodoc

import tidesurf

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "tidesurf"
copyright = "%Y, Jan Schleicher"
author = "Jan Schleicher"
release = ".".join(tidesurf.__version__.split(".")[:3])

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx_copybutton",
    "sphinx.ext.intersphinx",
    "sphinx.ext.napoleon",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "scanpydoc.definition_list_typed_field",
]

default_role = "literal"
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
html_title = "tidesurf"

# Automatic generation of API documentation
autosummary_generate = True
autodoc_default_options = {
    "member-order": "alphabetical",
    "show-inheritance": True,
}
autodoc_typehints = "both"
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_use_rtype = True
napoleon_use_param = True

intersphinx_mapping = dict(
    numpy=("https://numpy.org/doc/stable/", None),
    python=("https://docs.python.org/3", None),
    scipy=("https://docs.scipy.org/doc/scipy/", None),
)


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

    def link_replace(line):
        if line is not None:
            line = line.replace("np.", "~numpy.")
            line = line.replace("csr_matrix", "~scipy.sparse.csr_matrix")
            line = line.replace("List[", "~typing.List[")
            line = line.replace("Tuple[", "~typing.Tuple[")
            line = line.replace("Dict[", "~typing.Dict[")
            line = line.replace("Literal[", "~typing.Literal[")
            line = line.replace("Union[", "~typing.Union[")
            line = line.replace("Optional[", "~typing.Optional[")
        return line

    def autodoc_process_docstring(app, what, name, obj, options, lines):
        for i in range(len(lines)):
            lines[i] = link_replace(lines[i])

    def autodoc_process_signature(
        app, what, name, obj, options, signature, return_annotation
    ):
        return link_replace(signature), link_replace(return_annotation)

    app.connect("autodoc-skip-member", skip_member)
    app.connect("autodoc-process-docstring", autodoc_process_docstring)
    app.connect("autodoc-process-signature", autodoc_process_signature)


class MockedClassDocumenter(autodoc.ClassDocumenter):
    def add_line(self, line: str, source: str, *lineno: int) -> None:
        if line == "   Bases: :py:class:`object`":
            return
        super().add_line(line, source, *lineno)


autodoc.ClassDocumenter = MockedClassDocumenter
