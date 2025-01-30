Installation
============

Default installation
--------------------

To use **tidesurf**, we recommend installing it in a fresh virtual environment with Python version :math:`\geq` 3.10.
You can set up the environment as follows:

.. code-block:: console

    conda create -n <env_name> python=3.10
    conda activate <env_name>

PyPI
~~~~

You can install the package from PyPI:

.. code-block:: console

    pip install tidesurf

GitHub
~~~~~~

To install the latest development version, clone the repository from GitHub:

.. code-block:: console

    git clone git@github.com:janschleicher/tidesurf.git

Change into the directory and install with pip:

.. code-block:: console

    cd tidesurf
    pip install -e .

Contributing
------------

To install the development version for contributing to the project, you need to install the optional dependencies ``ruff`` and ``pytest``.
They are used for code formatting and testing, respectively, and should be run before committing changes.

.. code-block:: console

    pip install -e ".[dev]"