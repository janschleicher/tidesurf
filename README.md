[![Python 3.10](https://img.shields.io/badge/python-3.10-blue.svg)](https://www.python.org/downloads/release/python-31015/)
[![Python 3.11](https://img.shields.io/badge/python-3.11-blue.svg)](https://www.python.org/downloads/release/python-31110/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Code style: ruff](https://img.shields.io/badge/code%20style-ruff-red)](https://github.com/astral-sh/ruff)
[![codecov](https://codecov.io/gh/janschleicher/tidesurf/branch/main/graph/badge.svg?token=dMenu3eZkX)](https://codecov.io/gh/janschleicher/tidesurf)
[![Python package](https://github.com/janschleicher/tidesurf/actions/workflows/python-package.yml/badge.svg?branch=main)](https://github.com/janschleicher/tidesurf/actions/workflows/python-package.yml)
[![PyPI - Version](https://img.shields.io/pypi/v/tidesurf)](https://pypi.org/project/tidesurf/)
[![Documentation Status](https://readthedocs.org/projects/tidesurf/badge/?version=latest)](https://tidesurf.readthedocs.io/latest/?badge=latest)

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
usage: tidesurf [-h] [-v] [--orientation {sense,antisense}] [-o OUTPUT]
                [--no_filter_cells]
                [--whitelist WHITELIST | --num_umis NUM_UMIS]
                [--min_intron_overlap MIN_INTRON_OVERLAP]
                [--multi_mapped_reads]
                SAMPLE_DIR GTF_FILE

Program: tidesurf (Tool for IDentification and Enumeration of Spliced and Unspliced Read Fragments)
Version: 0.2.0

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
  --no_filter_cells     Do not filter cells.
  --whitelist WHITELIST
                        Whitelist for cell filtering. Set to 'cellranger' to
                        use barcodes in the sample directory. Alternatively,
                        provide a path to a whitelist.
  --num_umis NUM_UMIS   Minimum number of UMIs for filtering a cell.
  --min_intron_overlap MIN_INTRON_OVERLAP
                        Minimum number of bases that a read must overlap with
                        an intron to be considered intronic.
  --multi_mapped_reads  Take reads mapping to multiple genes into account
                        (default: reads mapping to more than one gene are
                        discarded).
```

## Contributing

For contributing, you should install `tidesurf` in development mode:

    pip install -e ".[dev]"

This will install the additional dependencies `ruff` and `pytest`, which are used for formatting and code style, and testing, respectively.
Please run these before commiting new code.

## Citation

If you use `tidesurf` in your research, please cite the following publication:

Schleicher, J.T., and Claassen, M. (2025).
Accurate quantification of spliced and unspliced transcripts for single-cell RNA sequencing with tidesurf.
_bioRxiv_ 2025.01.28.635274; DOI: [10.1101/2025.01.28.635274](https://doi.org/10.1101/2025.01.28.635274).