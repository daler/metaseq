Installation
============

:mod:`metaseq` relies on the standard `Scientific Python Stack
<http://www.scipy.org/stackspec.html>`_, specifically, `NumPy
<http://www.numpy.org/>`_, `SciPy <http://www.scipy.org/index.html>`_,
`matplotlib <http://matplotlib.org/>`_, and `pandas
<http://pandas.pydata.org/>`_. It also relies on standard genomics tools like
`BEDTools <http://bedtools.readthedocs.org/en/latest/>`_, `samtools
<http://samtools.sourceforge.net/>`_, and others.

These prerequisites can make :mod:`metaseq` difficult to install, so we provide
an installation script that handles all of the complicated work.


Easy installation method
------------------------
If you are just trying out :mod:`metaseq`, the best way to do so is to download
the installation script (see the appropriate section below for your operating
system, :ref:`mac` or :ref:`linux`).  This script works on Mac OSX and Linux, and will:

    - Download and install `BEDTools
      <http://bedtools.readthedocs.org/en/latest/>`_, `samtools
      <http://samtools.sourceforge.net/>`_, `tabix
      <http://samtools.sourceforge.net/tabix.shtml>`_, and the UCSC tools
      `bigWigSummary, bigWigToBedGraph, and bedGraphToBigWig
      <http://hgdownload.cse.ucsc.edu/admin/exe/>`_

    - Download and install `Miniconda
      <http://conda.pydata.org/miniconda.html>`_ (a slimmed-down version of the
      `Anaconda Python distribution
      <https://store.continuum.io/cshop/anaconda/>`_) which sets up an isolated
      Python environment that is separate from anything you might already have
      installed (see :ref:`whyconda` for more details on this).

    - Create an isolated Python environment with `conda
      <http://conda.pydata.org/docs/examples/create.html>`_

    - Download and install prerequisites for :mod:`metaseq` into the test
      environment

    - Download and install :mod:`metaseq` itself into the test environment


The script will tell you what it's doing, and at the end will prompt you if you
want to add the installation locations to your PATH variable (if you're not
sure what this is, then you should say "yes").  It will print a README.txt file
with the results and some additional instructions to finalize the installation.



.. warning::

    The installation script depends on several external servers (UCSC, github,
    PyPI) beyond our immediate control.  If the script seems to hang for more
    than a couple of minutes, or you get unxplained error messages, please use
    Ctrl-C to abort and try again later.

    If you run into difficulties that are not solved by re-running the script,
    please `open an issue on github <https://github.com/daler/metaseq/issues>`_
    describing the details of the problem.

In the end, you will have a complete scientific Python installation, along with
some commonly-used genomics tools.

.. _mac:

Mac OSX
~~~~~~~

.. note::

    On Mac OSX, you will first need to install `Xcode
    <https://developer.apple.com/xcode/>`_, which provides C and C++ compilers.
    You can get Xcode for free directly from Apple, and the version to get
    depends on the version of OSX you are running.  Note that you may have to
    register for a free developer account.

To download the script and perform the installation using default settings,
paste the following two commands in a Terminal window.  The first command
downloads the script, and the second command runs it::

    curl -O https://raw.githubusercontent.com/daler/metaseq/master/create-metaseq-test-environment.sh

::

    bash create-metaseq-test-environment.sh

.. _linux:

Linux
~~~~~
.. note::

    On Linux, you will need a C and C++ compiler as well as the zlib
    development libraries, which don't come installed by default.  In Ubuntu,
    the following command should install these for you::

        sudo apt-get install build-essential zlib1g-dev

On Linux, `wget` is usually available by default instead of `curl`.  So paste
these two commands into a terminal instead to perform the installation using
default settings.  The first command downloads the script, and the second
command runs it::

    wget https://raw.githubusercontent.com/daler/metaseq/master/create-metaseq-test-environment.sh

::

    bash create-metaseq-test-environment.sh

Customizing
~~~~~~~~~~~
If you want to customize the installation locations, specify versions, or only
install a subset of the prerequisites, you can view the help with::

    bash create-metaseq-test-environment.sh -h

Uninstalling
~~~~~~~~~~~~

Uninstalling is straightforward.  **Assuming you used the default locations:**

* Delete ``~/miniconda/envs/metaseq-test`` to uninstall just the test
  environment.
* Delete ``~/miniconda`` to uninstall the test environment and all of miniconda.
* Delete ``~/tools`` to uninstall the genomics tools.  Specifically, the
  installation script creates the following directories and files within
  `~/tools`:

    * ``bedtools<VERSION>/``  (where BEDTools is installed)
    * ``samtools<VERSION>/``  (where samtools is installed)
    * ``tabix<VERSION>/``  (where tabix is installed)
    * ``ucsc/`` (where bigWigSummary and other UCSC programs are installed)
    * ``logs/``  (any logs from the installation process)
    * ``README.txt`` (post-installation instructions)
    * ``miniconda-paths`` (describes where miniconda was installed)
    * ``paths`` (describes where genomics tools were installed)

* Optionally, if you added anything to your PATH, you can delete the relevant
  lines in your `~/.bashrc` or `~/.bash_profile` file, but this is not strictly
  necessary if these directories are deleted.



Custom installation
-------------------
Even if you do not want to use the default full installation script described
above, it can still be useful to install the individual components.  See the
help for that script for the full details, but useful flags are:

* `-M` disables the miniconda installation
* `-i` controls which genomics tools are installed
* `-g` controls which :mod:`metaseq` version to install (specified as tags or commits from
  github).  The special tag "disable" will disable installation of metaseq.

Some example use-cases:

* Only install BEDTools::

    bash create-metaseq-test-environment.sh -M -i "bedtools" -g disable

* Install just the latest commit of metaseq into your system-wide Python
  installation (note: you will need to run the script with sudo priviliges,
  since it uses `pip install`)::

    bash create-metaseq-test-environment.sh -M -i "" -g master

* Same thing, but install it into the test environment::

    bash create-metaseq-test-environment.sh -i "" -g master


Manual installation
-------------------

Step 1: Non-python programs
~~~~~~~~~~~~~~~~~~~~~~~~~~~
The following non-Python programs are needed:

* A C and C++ compiler
* `BEDTools`, `samtools`, and `Tabix`
* bigWigSummary, bigWigToBedGraph, bedGraphToBigWig

If you don't already have them installed, the installation script described
above is the easiest way to get these.


Step 2. Install Python packages
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Option 1: Install from PyPI
+++++++++++++++++++++++++++
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



Footnotes
---------

.. _whyconda:

Miniconda instead of virtualenv?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the past, the standard way of creating isolated environments was to use
`virtualenv <http://virtualenv.readthedocs.org/en/latest/>`_.  The standard
procedure is to create a blank environment, and `pip install` all necessary
requirements, essentially installing everything from scratch.  However, for
packages like :mod:`metaseq` with many dependencies, installing and compiling
from scratch can take a lot of time.

Recently, the Anaconda Python distribution has provided another way of creating
isolated environments.  It has made it much easier to install the scienfific
Python stack because it provides pre-compiled versions of numpy, scipy,
matplotlib, and other hard-to-install packages.  This drastically reduces the
amount of time it takes to set up an isolated environment.

We decided to use Miniconda (a slimmed-down version of Anaconda) for the
:mod:`meteaseq` installation script because it provides the user with an
isolated environment in a fraction of the time of a full virtualenv
installation, and does not require a FORTRAN compiler for installing scipy.

Tests
~~~~~
After every change to :mod:`metaseq`, tests are run by the Travis-CI continuous
integration service.  You can always check the status by visiting
https://travis-ci.org/daler/metaseq/.  These tests are run by setting up the
test environment in Ubuntu 12.04 using the `create-metaseq-test-environment.sh`
script described above.
