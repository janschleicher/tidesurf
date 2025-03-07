import os

import pysam
from Cython.Build import cythonize
from setuptools import Extension, setup

# Set environment variable EXT_BUILD_DEBUG=1 to enable debug mode
IS_DEBUG = os.getenv("EXT_BUILD_DEBUG", False)

debug_macros = []
debug_cythonize_kw = dict(
    force=True,
    compiler_directives={
        "language_level": "3",
        "embedsignature": True,
        "annotation_typing": False,
    },
)
debug_include_path = []

if IS_DEBUG:
    print("Extension IS_DEBUG=True!")
    # Adding cython line trace for coverage report
    debug_macros += ("CYTHON_TRACE_NOGIL", 1), ("CYTHON_TRACE", 1)
    # Adding upper directory for supporting code coverage when running tests inside the cython package
    debug_include_path += [".."]
    # Some extra info for cython compilator
    debug_cythonize_kw.update(
        dict(
            gdb_debug=True,
            force=True,
            annotate=True,
        )
    )
    debug_cythonize_kw["compiler_directives"].update(
        {"linetrace": True, "profile": True, "binding": True}
    )

extensions = [
    Extension("tidesurf.enums", ["tidesurf/enums.pyx"]),
    Extension("tidesurf.transcript", ["tidesurf/transcript.py"]),
    Extension(
        "tidesurf.counter", ["tidesurf/counter.py"], include_dirs=pysam.get_include()
    ),
    Extension("tidesurf.main", ["tidesurf/main.py"], include_dirs=pysam.get_include()),
]

for ext in extensions:
    ext.define_macros = debug_macros

setup(
    name="tidesurf",
    ext_modules=cythonize(
        extensions,
        include_path=debug_include_path,
        **debug_cythonize_kw,
    ),
)
