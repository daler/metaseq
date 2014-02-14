#!/bin/bash

set -e

# Copied from pybedtools travis setup scripts.
TABIX_VERSION=0.2.6
SAMTOOLS_VERSION=0.1.19
./.install-tabix.sh $TABIX_VERSION ${TRAVIS_BUILD_DIR}
./.install-bedtools2.sh ${TRAVIS_BUILD_DIR}
./.install-samtools.sh $SAMTOOLS_VERSION ${TRAVIS_BUILD_DIR}
export PATH=${TRAVIS_BUILD_DIR}/tabix-${TABIX_VERSION}:$PATH
export PATH=${TRAVIS_BUILD_DIR}/samtools-${SAMTOOLS_VERSION}:$PATH
export PATH=${TRAVIS_BUILD_DIR}/bedtools2/bin:$PATH
echo $PATH


# This limits the number of CPUs used in test.py.
export METASEQ_PROCESSES=1

# Prepare UCSC's bigWigSummary
mkdir -p ${TRAVIS_BUILD_DIR}/ucsc
wget \
    http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigSummary \
    -O ${TRAVIS_BUILD_DIR}/ucsc/bigWigSummary
chmod +x ${TRAVIS_BUILD_DIR}/ucsc/bigWigSummary
export PATH=${TRAVIS_BUILD_DIR}/ucsc:$PATH


nosetests -v metaseq/test/test.py
