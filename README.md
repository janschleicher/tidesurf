[![Python 3.10](https://img.shields.io/badge/python-3.10-blue.svg)](https://www.python.org/downloads/release/python-31015/)
[![Python 3.11](https://img.shields.io/badge/python-3.11-blue.svg)](https://www.python.org/downloads/release/python-31110/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Code style: ruff](https://img.shields.io/badge/code%20style-ruff-red)](https://github.com/astral-sh/ruff)
[![codecov](https://codecov.io/gh/janschleicher/tidesurf/branch/main/graph/badge.svg?token=dMenu3eZkX)](https://codecov.io/gh/janschleicher/tidesurf)
[![Python package](https://github.com/janschleicher/tidesurf/actions/workflows/python-package.yml/badge.svg?branch=main)](https://github.com/janschleicher/tidesurf/actions/workflows/python-package.yml)
![PyPI - Version](https://img.shields.io/pypi/v/tidesurf)

# tidesurf

This repository provides a Tool for IDentification and Enumeration of Spliced and Unspliced Read Fragments using Python.

## Installation

### From PyPI

Set up a virtual environment using Conda with Python version >=3.10 and activate it:

    conda create -n <envName> python=3.10
    conda activate <envName>

Install the package from PyPI:
    
    pip install tidesurf

### Latest version from GitHub

Clone the repository:

    git clone git@github.com:janschleicher/tidesurf.git

Change into the directory and install with pip:
    
    cd tidesurf
    pip install -e .

## Usage

```
```

## Contributing

For contributing, you should install `tidesurf` in development mode:

    pip install -e ".[dev]"

This will install the additional dependencies `ruff` and `pytest`, which are used for formatting and code style, and testing, respectively.
Please run these before commiting new code.