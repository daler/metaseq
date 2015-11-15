.. metachip documentation master file, created by
   sphinx-quickstart on Sat Oct  8 14:09:02 2011.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. include:: README.rst

Some use-cases for :mod:`metaseq`:

* Create NumPy arrays (which can be plotted as heatmaps) where each row
  represents a feature (say, TSS +/- 1kb) at fairly high speed.  The actual
  speed is highly dependent on your data, but a rule of thumb is that it
  takes about 10s for every 10,000 2-kb features.  Once you have the heatmap,
  you can sort, cluster, zoom, and pan to explore interesting groups.  See
  the :ref:`CTCF example` for more on this as well as the
  :mod:`metaseq.genomic_signal` module.

* "Mini-browsers".  Instead of tweaking something, then laboriously converting
  to bedGraph, then bigWig, uploading your tracks to the genome browser, then
  searching for the region of interest, you can can just connect a mini-browser
  to your raw data and spawn figures locally.  See the :ref:`example_session` for
  more on this as well as the :mod:`metaseq.minibrowser` API docs.

* Cluster genes based on the spatial distribution of ChIP-seq peaks around
  their TSSs.

* Scatter plot of DESeq results (basemeana vs basemeanb) where points are
  colored according to the number of ChIP peaks in the gene.  This, too, can be
  attached to mini-browsers, enabling you to click on a point to see the
  genomic signal.  See the :mod:`metaseq.results_table` module for more on
  this.

* Pie charts of where peaks fall within annotated genes -- TSS, poly-A
  site, intron, exon, etc.  See the :mod:`metaseq.integration` module for more
  on this.

Where possible, the inputs are standard formats -- BED, GFF, GTF, BAM, SAM,
DESeq results as saved from R, or even arbitrary tab-delimited data files that
have a header.  If you take the time to convert to bigWig or bigBed,
performance will be improved.



Contents:

.. toctree::
    :maxdepth: 2

    install
    running-the-examples
    example_session
    example_session_2
    autodocs
    changes
