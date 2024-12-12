[![Python 3.10](https://img.shields.io/badge/python-3.10-blue.svg)](https://www.python.org/downloads/release/python-31015/)
[![Python 3.11](https://img.shields.io/badge/python-3.11-blue.svg)](https://www.python.org/downloads/release/python-31110/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Code style: ruff](https://img.shields.io/badge/code%20style-ruff-red)](https://github.com/astral-sh/ruff)
[![codecov](https://codecov.io/gh/janschleicher/tidesurf/branch/main/graph/badge.svg?token=dMenu3eZkX)](https://codecov.io/gh/janschleicher/tidesurf)
[![Python package](https://github.com/janschleicher/tidesurf/actions/workflows/python-package.yml/badge.svg?branch=main)](https://github.com/janschleicher/tidesurf/actions/workflows/python-package.yml)

# tidesurf

This repository provides a Tool for IDentification and Enumeration of Spliced and Unspliced Read Fragments using Python.

## Installation

Set up a virtual environment using Conda with Python version >=3.10

    conda create -n <envName> python=3.10

Clone the repository

    git clone git@github.com:janschleicher/tidesurf.git

Change into the directory and install with pip
    
    cd tidesurf
    pip install -e .

## Usage

```
usage: tidesurf [-h] [-v] [--orientation {sense,antisense}] [-o OUTPUT]
                [--filter_cells] [--whitelist WHITELIST | --num_umis NUM_UMIS]
                SAMPLE_DIR GTF_FILE

Program: tidesurf (Tool for IDentification and Enumeration of Spliced and Unspliced Read Fragments)
Version: 0.1.dev31+g30b513e.d20241211

positional arguments:
  SAMPLE_DIR            Sample directory containing Cell Ranger output.
  GTF_FILE              GTF file with transcript information.

options:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  --orientation {sense,antisense}
                        Orientation of reads with respect to transcripts. For
                        10x Genomics, use 'sense' for three prime and
                        'antisense' for five prime.
  -o OUTPUT, --output OUTPUT
                        Output directory.
  --filter_cells        Filter cells based on a whitelist.
  --whitelist WHITELIST
                        Whitelist for cell filtering. Set to 'cellranger' to
                        use barcodes in the sample directory. Alternatively,
                        provide a path to a whitelist.
  --num_umis NUM_UMIS   Minimum number of UMIs for filtering a cell.
```

## Contributing

For contributing, you should install `tidesurf` in development mode:

    pip install -e ".[dev]"

This will install the additional dependencies `ruff` and `pytest`, which are used for formatting and code style, and testing, respectively.
Please run these before commiting new code.