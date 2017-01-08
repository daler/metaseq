#!/bin/bash

# Sometimes when running python setup.py sdist upload, MAINIFEST.in can catch
# extras in the development source dir that haven't been tested. This ensures
# that the only things making it into the source distribution has been commited
# to the repo.
HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
rm -rf /tmp/metaseq
git clone $HERE/.. /tmp/metaseq
cd /tmp/metaseq
python setup.py sdist upload
