Metaseq
=======
.. note:: 

    The main documentation for `metachip` can be found at
    http://packages.python.org/metaseq/

The goal of :mod:`metaseq` is to tie together lots of existing software into
a framework for exploring genomic data.  It focuses on flexibility and
interactive exploration of data.

:mod:`metaseq` is sort of like a glue that brings together lot of
great software commonly used for genomic data analysis.  It gives you the tools you
need to do analyses that are difficult to do in any other way that I know of.
It ties together other packages like:

* :mod:`Pysam` [1]_ for access to BAM and SAM files
* :mod:`bx-python` [2]_ for access to bigWig and bigBed files
* :mod:`pybedtools` [3]_, which integrates BEDTools into Python, for fast "genome
  algebra"
* :mod:`matplotlib` [4]_ for rich, interactive, highly customizable plotting
* :mod:`gffutils` [5]_ for hierarchical access to gene models
* :mod:`Cython` [6]_ for speed
* :mod:`scikits.learn` [7]_ for fast clustering
* :mod:`NumPy` [8]_ for large, fast arrays




:mod:`metaseq` also:

* Supplies a fast "rebinning" routine that is compiled in C for speed
* Accepts bigBed, bigWig, BAM and SAM files as input for random access of
  continuous data (and with Tabix, accepts BED, GFF, GTF, VCF as well)
* Support for tunable (i.e., tradeoffs between CPU, RAM, and I/O) parallel
  processing of data accessed from the file formats mentioned above
* Provides a framework for "mini-browsers", zoomable and pannable Python-only
  figures that show genomic signal and gene models and are spawned by clicking
  on features of interest
* Integrates RNA-seq and ChIP-seq data by using gffutils to convert between
  gene IDs and coordinates as well as handle full gene models and plotting.


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
  to your raw data and spawn figures locally.  See the :ref:`CTCF example` for
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



.. [1] http://code.google.com/p/pysam/
.. [2] https://bitbucket.org/james_taylor/bx-python/wiki/Home
.. [3] http://packages.python.org/pybedtools
.. [4] http://matplotlib.sourceforge.net/
.. [5] http://packages.python.org/gffutils
.. [6] http://cython.org/
.. [7] http://scikit-learn.org/stable/
.. [8] http://www.scipy.org/
