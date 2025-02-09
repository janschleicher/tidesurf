# import numpy as np
import pysam
from Cython.Build import cythonize
from setuptools import Extension, setup

extensions = [
    Extension("tidesurf.enums", ["tidesurf/enums.pyx"]),
    Extension("tidesurf.transcript", ["tidesurf/transcript.py"]),
    Extension(
        "tidesurf.counter", ["tidesurf/counter.py"], include_dirs=pysam.get_include()
    ),
    Extension("tidesurf.main", ["tidesurf/main.py"], include_dirs=pysam.get_include()),
]

setup(
    name="tidesurf",
    ext_modules=cythonize(
        extensions,
        compiler_directives={
            "language_level": "3",
            "embedsignature": True,
            "annotation_typing": False,
        },
    ),
)
