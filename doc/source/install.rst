.. _installation:

Installation
============

:mod:`metaseq` relies on the standard `Scientific Python Stack
<http://www.scipy.org/stackspec.html>`_, specifically, `NumPy
<http://www.numpy.org/>`_, `SciPy <http://www.scipy.org/index.html>`_,
`matplotlib <http://matplotlib.org/>`_, and `pandas
<http://pandas.pydata.org/>`_. It also relies on standard genomics tools like
`BEDTools <http://bedtools.readthedocs.org/en/latest/>`_, `samtools
<http://samtools.sourceforge.net/>`_, and others.

Installing each of these prerequisites would be time-consuming and frustrating.
However `conda <http://conda.pydata.org/docs/intro.html>`_ combined with the
`bioconda project <https://github.com/bioconda/recipes>`_ makes it easy and
straightforward.


Easy installation method
------------------------

The following assumes you have `conda
<http://conda.pydata.org/docs/intro.html>`_ installed. `conda` is
a cross-platform package manager that started on Python but now includes other
languages and binaries as well.  There are two options for installing conda:
the `Anaconda Python distribution
<https://store.continuum.io/cshop/anaconda/>`_, which includes the full
scientific Python stack, or the smaller, barebones `Miniconda
<http://conda.pydata.org/miniconda.html>`_.


The commands below work with either Linux or Mac OSX. See `managing
environments <http://conda.pydata.org/docs/using/envs.html>`_ to decide if you
want to use the main `conda` environment or an isolated environment. If you
want to use an isolated environment, use `create -n my-env-name` instead of
`install` in the commands below.

**Option 1:** Install the `metaseq` package and Python dependencies

.. code-block:: bash

    conda install --channel bioconda metaseq

**Option 2:** Install the `metaseq` package and Python dependencies, plus BEDTools, samtools,
htslib, and the UCSC tools bigWigSummary, bedGraphToBigWig, bedToBigBed, and
bigBedToBed. This is a good option if you don't already maintain separate
installations of these tools

.. code-block:: bash

    conda install --channel bioconda metaseq-all


.. _windows:

Windows
~~~~~~~
It is difficult to do bioinformatics work on Windows. The most convenient
option will be to run Linux on your Windows machine (via dual-boot or virtual
machine), and follow the Linux instructions above to install metaseq and
requirements.
