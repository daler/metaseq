#!/bin/bash

#
# This script ensures that an isolated installation of metaseq and genomics
# tools can be successfully installed and passes all tests.
#
# It is intended to be used just before releasing a new version compiled by
# setup.py sdist (as opposed to using the current state of the repository, as
# is done in the travis-ci tests).
#
# Takes ~5 mins to run.
#
set -e

# Get rid of any source distributions and other version of miniconda or tools
# installed by the create-metaseq-test-environment script.
rm -rf "$(pwd)/tools" "$(pwd)/miniconda" "$(pwd)/dist"

# Extract version info
VERSION=$(cut -f 3 -d " " $(pwd)/metaseq/version.py | sed 's/"//g')

# Create sdist
python setup.py sdist

# Install tools and miniconda into this directory.
# `-g disable` disables installing metaseq itself, which we'll do later
./create-metaseq-test-environment.sh -d "$(pwd)/tools" -t -m "$(pwd)/miniconda" -e metaseq-test -v -g disable

# Source the paths
source "$(pwd)/tools/miniconda-paths"
source "$(pwd)/tools/paths"

# activate the new environment
source activate metaseq-test
echo "$(pwd)"

# install the sdist version
pip install dist/metaseq-${VERSION}.tar.gz

# finally run the tests...
nosetests -v metaseq/test/test.py
