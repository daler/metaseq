Overview
========
This directory contains the materials needed to recreate the IPython notebook
used for the SciPy 2013 talk, "metaseq: a Python framework for integrating
high-throughput sequencing analyses"

Getting the data
================
The `data/get-data.py` script will download publicly available data from the
ENCODE project, from mouse mm9 assembly.  To save on disk space and execution
time in the IPython notebook, only data from chr7 are kept.

This script operates on the raw data provided by the ENCODE project.  There are
lots of commented-out data sets, which you can un-comment at the expense of
download and processing time, but by default only the following data are
retrieved:


* reads from bam files of mapped reads for Gata1 and input ChIP-seq in G1E-ER cells
* Download full BED files of called peaks for Gata1 in G1E-ER cells
* Annotations from Ensembl for mm9
* counts tables from PMID:23390196 (Ldb1 knockout vs WT RNA-seq)


The data are then processed: bam files are indexed, converted to bedGraph, and
then to bigWig; BED files are filtered to only contain chr7 peaks, then
converted to bigBed; the GTF file is filtered to include only chr7 annotations,
then converted to a gfftutils database; the counts tables are run through DESeq
(in R) to create final results tables.

Requirements
------------
Non-Python requirements
~~~~~~~~~~~~~~~~~~~~~~~
* samtools
* BEDTools
* R
* DESeq package for R
* bedToBigBed (from UCSC)
* bedGraphToBigWig (from UCSC)
* unix tools wget, cut, gunzip

Python requirements
-------------------
These are the same requirements for `metaseq`:

* gffutils
* pybedtools
* bx-python
* pysam
* matplotlib
* scipy
* scikit-learn
* numpy

Run the script
--------------
Change to the `data` directory, and run::

    $ python get-data.py


IPython Notebook
================
The IPython notebook, `metaseq-demo.ipynb`, can be used directly.

Slides
======

In order to create slides, the dev version of IPython was used::

    git clone https://github.com/ipython/ipython.git \
    && cd ipython \
    && sudo python setup.py develop

I've included the version of nbconvert that I was able to get custom CSS
working with.

To make the slides::

    make slides

And then open `nbconvert/metaseq-demo_slides.html`.
