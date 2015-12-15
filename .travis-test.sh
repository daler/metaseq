#!/bin/bash

set -e

# This limits the number of CPUs used in test.py.
export METASEQ_PROCESSES=1

nosetests -v metaseq/test/test.py

(cd doc && make html)
