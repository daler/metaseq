#!/bin/bash

# Tests the installation of requirements via bioconda channel and via pip for
# both Python 2 and 3.
#
# Runs main tests and doctests for Python 2 and 3.

set -e
set -x

HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"


# clone into a separate directory just for this python version
src=/tmp/metaseq
mkdir $src

# travis-ci pulls the repo using --depth=50 which creates a shallow clone.
# Getting rid of the `shallow` file lets us clone it elsewhere, avoiding
# the "fatal: attempt to fetch/clone from a shallow repository" error.
#
# The reason we're cloning in the first place is to avoid root doing
# anything in the existing directory -- especially creating the sdist
rm -f $HERE/../.git/shallow
git clone $HERE/.. $src
cd $src

# extract version
VERSION=$(python -c 'exec(open("metaseq/version.py").read());print(__version__)')

# ------------------------------------------------------------------------
# Install prerequisites with conda -- otherwise it'll take a long time
# (but install metaseq from sdist)
conda create -y -n metaseq-env -c bioconda python=2 --file requirements.txt --file bioconda-requirements.txt
set +x; source activate metaseq-env; set -x
python setup.py clean sdist
pip install dist/metaseq-${VERSION}.tar.gz

conda install -y nose
export METASEQ_PROCESSES=1
nosetests metaseq/test/test.py

# Install tools and run doctests
conda install -y --file docs-requirements.txt
(cd doc && make clean && make html)
set +x; source deactivate; set -x

