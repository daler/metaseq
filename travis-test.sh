#!/bin/bash

set -e
set -x

conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda

conda create -y -n metaseq-env python=2 --file requirements.txt --file bioconda-requirements.txt
set +x; source activate metaseq-env; set -x

VERSION=$(python -c 'exec(open("metaseq/version.py").read());print(__version__)')
python setup.py clean sdist
pip install dist/metaseq-${VERSION}.tar.gz

conda install -y nose
export METASEQ_PROCESSES=1
nosetests metaseq/test/test.py

# Install tools and run doctests
conda install -y --file docs-requirements.txt
(cd doc && make clean && make html)
set +x; source deactivate; set -x
