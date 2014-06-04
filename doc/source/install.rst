Installation
============

:mod:`metaseq` relies on the standard `Scientific Python Stack
<http://www.scipy.org/stackspec.html>`_.  Specifically, `NumPy
<http://www.numpy.org/>`_, `matplotlib <http://matplotlib.org/>`_, and `pandas
<http://pandas.pydata.org/>`_. It also relies on standard genomics tools like
`BEDTools <http://bedtools.readthedocs.org/en/latest/>`_, `samtools
<http://samtools.sourceforge.net/>`_, and others.

These prerequisites can make :mod:`metaseq` difficult to install, so we provide
an installation script that handles all of the complicated work.


Easy installation method
------------------------
If you are just trying out :mod:`metaseq`, the best way to do so is to download
the installation script.  This script works on Mac OSX and Linux, and will:

    - Download and install `BEDTools
      <http://bedtools.readthedocs.org/en/latest/>`_, `samtools
      <http://samtools.sourceforge.net/>`_, `tabix
      <http://samtools.sourceforge.net/tabix.shtml>`_, and the UCSC tools
      `bigWigSummary, bigWigToBedGraph, and bedGraphToBigWig
      <http://hgdownload.cse.ucsc.edu/admin/exe/>`_

    - Download and install `Miniconda
      <http://conda.pydata.org/miniconda.html>`_, which sets up an isolated
      Python environment that is separate from anything you might already have
      installed

    - Create an isolated Python environment with `conda
      <http://conda.pydata.org/docs/examples/create.html>`_

    - Download and install prerequisites for :mod:`metaseq` into the test
      environment

    - Download and install :mod:`metaseq` itself into the test environment


The script will tell you what it's doing, and at the end will prompt you if you
want to add the installation locations to your PATH variable (if you're not
sure what this is, then you should say "yes").  It will print a README.txt file
with the results and some additional instructions to finalize the installation.


.. note::

    On Mac OSX, you will first need to install `Xcode
    <https://developer.apple.com/xcode/>`_, which provides C and C++ compilers.
    You can get Xcode for free directly from Apple, and the version to get
    depends on the version of OSX you are running.  You may have to register
    for a free developer account.

Mac OSX
~~~~~~~
Paste the following commands in a Terminal window to download the script and
perform the installation using default settings::

    curl -O https://raw.githubusercontent.com/daler/metaseq/master/create-metaseq-test-environment.sh
    bash create-metaseq-test-environment.sh

Linux
~~~~~
On Linux, `wget` is usually available by default instead of `curl`.  So paste
these commands into a terminal instead::

    wget https://raw.githubusercontent.com/daler/metaseq/master/create-metaseq-test-environment.sh
    bash create-metaseq-test-environment.sh

Customizing
~~~~~~~~~~~
If you want to customize the installation locations, specify versions, or only
install a subset of the prerequisites, you can view the help with::

    bash create-metaseq-test-environment.sh -h

Uninstalling
~~~~~~~~~~~~

Uninstalling is straightforward -- as long as you used the default locations,
then simply remove the `miniconda` directory and the `tools` directory.


Detailed installation
---------------------
If you do not want to use the installation script described above, you can
install the components needed by :mod:`metaseq` manually.  Note however that
the installation script offers lots of flexibility, and allows you to
mix-and-match components.  For example, the script can install just the
genomics tools but not Miniconda or metaseq itself. This can help speed up the
installation process. See the help for that script for details.

Step 1: Non-python programs
~~~~~~~~~~~~~~~~~~~~~~~~~~~
The following non-Python programs are needed:

* A C/C++ compiler
* `BEDTools`, `samtools`, and `Tabix`
* bigWigSummary, bigWigToBedGraph, bedGraphToBigWig

The `installation page for pybedtools
<https://pythonhosted.org/pybedtools/main.html>`_ has installation instructions
for these prerequisites.  Alternatively, see the :ref:`from_scratch` section below for
installation scripts.


Step 2. Install Python packages
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Option 1 (recommended): Install from PyPI
+++++++++++++++++++++++++++++++++++++++++
The most robust method for installing :mod:`metaseq` is to do a 2-stage
installation.  First, ensure the base prerequisites are installed.  If any of
these are installed, a message will be printed on the screen indicating so.  Note that the
Anaconda Python Distribution comes with these packages, so you don't
necessarily need to run this:

::

    pip install Cython numpy pycurl


Then install :mod:`metaseq`, which will install any remaining dependencies:

::

    pip install metaseq

If you are not using the Anaconda Python Distribution, you may need to
be root in order to successfully run the above commands.

Option 2: Install from source
+++++++++++++++++++++++++++++
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


Python dependencies
-------------------
Here is a list of dependencies, and why they are needed by :mod:`metaseq`:

:NumPy:
    fast arrays; required for pandas, scikit-learn, matplotlib, and much of
    metaseq.

:Cython:
    dependency for pybedtools, bx-python, statsmodels, pandas

:matplotlib:
    interactive plotting

:pandas:
    fast tablular data

:bx-python:
    accessing bigWig and bigBed files

:pysam:
    accessing BAM files

:scikit-learn (optional):
    clustering

:gffutils:
    GFF and GTF manipulation

:pybedtools:
    BED, BAM, GTF, GFF, VCF manipulation

:urlgrabber:
    downloading example data

:PyYAML:
    config files
