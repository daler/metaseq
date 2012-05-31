Metaseq
=======
.. note:: 

    The main documentation for `metaseq` can be found at
    http://packages.python.org/metaseq/

The goal of `metaseq` is to tie together lots of existing software into
a framework for exploring genomic data.  It focuses on flexibility and
interactive exploration of data.

`metaseq` is sort of like a glue that brings together lot of great software
commonly used for genomic data analysis.  It gives you the tools you need to do
analyses that are difficult to do in any other way that I know of.  It ties
together other packages like:

* `Pysam` [1]_ for access to BAM and SAM files
* `bx-python` [2]_ for access to bigWig and bigBed files
* `pybedtools` [3]_, which integrates BEDTools into Python, for fast "genome
  algebra"
* `matplotlib` [4]_ for rich, interactive, highly customizable plotting
* `gffutils` [5]_ for hierarchical access to gene models
* `Cython` [6]_ for speed
* `scikits.learn` [7]_ for fast clustering
* `NumPy` [8]_ for large, fast arrays


`metaseq` also:

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


.. [1] http://code.google.com/p/pysam/
.. [2] https://bitbucket.org/james_taylor/bx-python/wiki/Home
.. [3] http://packages.python.org/pybedtools
.. [4] http://matplotlib.sourceforge.net/
.. [5] http://packages.python.org/gffutils
.. [6] http://cython.org/
.. [7] http://scikit-learn.org/stable/
.. [8] http://www.scipy.org/
