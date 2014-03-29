Installation
============

From PyPI
---------
::

    pip install metaseq

From source
-----------
::

    git clone https://github.com/daler/metaseq.git
    cd metaseq
    sudo pip install -r requirements.txt
    sudo python setup.py develop

From scratch
------------

Continuous integration tests are run on Travis-CI.  If you have trouble
installing dependencies above, the following scripts allow you to re-create a full test
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


