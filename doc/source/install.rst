Installation
============

:mod:`metaseq` relies on the standard `Scientific Python Stack
<http://www.scipy.org/stackspec.html>`_.  Specifically, `NumPy
<http://www.numpy.org/>`_, `matplotlib
<http://matplotlib.org/>`_, and `pandas <http://pandas.pydata.org/>`_.

If you do not already have these packages installed, the easiest way to get
them is to download the free `Anaconda Python Distribution
<https://store.continuum.io/cshop/anaconda/>`_, which comes with these packages
(and more) already installed.

.. note::

    The `Anaconda Python Distribution
    <https://store.continuum.io/cshop/anaconda/>`_ allows you to install these
    packages even if you do not have admin or root permission on your machine.


Step 1. Install non-Python prerequisites
----------------------------------------
The following non-Python programs are needed:

* A C/C++ compiler
* `BEDTools`, `samtools`, and `Tabix`

The `installation page for pybedtools
<https://pythonhosted.org/pybedtools/main.html>`_ has installation instructions
for these prerequisites.  Alternatively, see the :ref:`from_scratch` section below for
installation scripts.


Step 2. Install Python packages
-------------------------------

Option 1 (recommended): Install from PyPI
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The most robust method for installing :mod:`metaseq` is to do a 2-stage
installation.  First, ensure the base prerequisistes are installed:

::

    pip install Cython numpy pycurl


Then install :mod:`metaseq`, which will install any remaining dependencies:

::

    pip install metaseq

If you are not using the Anaconda Python Distribution, you may need to
be root in order to successfully run the above commands.

Option 2: Install from source
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
::

    git clone https://github.com/daler/metaseq.git
    cd metaseq
    pip install -r requirements.txt
    python setup.py develop




.. _from_scratch:

Installing from scratch
-----------------------

Continuous integration tests are run on `Travis-CI
<https://travis-ci.org/daler/metaseq>`_.  If you have trouble installing
dependencies above, the following scripts allow you to re-create a full test
environment on Ubuntu 12.04 LTS.

The path names are relative to the top-level source directory.

`.travis.yml`
~~~~~~~~~~~~~
This is the configuration file for Travis-CI; when running on a normal machine
the commands in the "install" list should be run from the command line.  Then
the `.travis-test.sh` script (below) can be run.

You should set the `TRAVIS_BUILD_DIR` environment variable to wherever you want
the dependencies to be installed.

.. literalinclude:: ../../.travis.yml
    :language: yaml

`.travis-test.sh`
~~~~~~~~~~~~~~~~~
This script in turn calls installation scripts for `bedtools`, `samtools`, and
`tabix` (below).

.. literalinclude:: ../../.travis-test.sh
    :language: bash

`.install-bedtools2.sh`
~~~~~~~~~~~~~~~~~~~~~~~
.. literalinclude:: ../../.install-bedtools2.sh
    :language: bash

`.install-tabix.sh`
~~~~~~~~~~~~~~~~~~~
.. literalinclude:: ../../.install-tabix.sh
    :language: bash

`.install-samtools.sh`
~~~~~~~~~~~~~~~~~~~~~~
.. literalinclude:: ../../.install-samtools.sh
    :language: bash


Python dependencies
-------------------
Here is a list of dependencies, and why they are needed by :mod:`metaseq`:

:NumPy:
    fast arrays; required for pandas, scikit-learn, matplotlib, and much of
    metaseq.

:Cython:
    Dependency for pybedtools, bx-python, statsmodels, pandas

:SciPy:
    dependency


:matplotlib:
    interactive plotting

:pandas:
    fast tablular data

:bx-python:
    Accessing bigWig and bigBed files

:pysam:
    Accessing BAM files

:scikit-learn (optional):
    Clustering

:gffutils:
    GFF and GTF manipulation

:pybedtools:
    BED, BAM, GTF, GFF, VCF manipulation

:urlgrabber:
    For downloading example data

:PyYAML:
    Config files
