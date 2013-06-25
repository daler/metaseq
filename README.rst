`Metaseq`
=========
.. note:: 

    The main documentation for `metaseq` can be found at
    http://packages.python.org/metaseq/

    See http://packages.python.org/metaseq/example_session.html for an example
    session.


The goal of `metaseq` is to tie together lots of existing software into
a framework for exploring genomic data.  It focuses on flexibility and
interactive exploration of data.

There are three main components:

1. Tools for accessing genomic data, based on interval, and returning NumPy
   arrays (see the `metaseq.genomic_signal` module, and
   `metaseq.integration.chipseq` for more specific uses).  An example use is to
   extract histone modification ChIP-seq signal from called peaks of
   a transcription factor ChIP-seq experiment, and plot them as heatmaps or as
   summarized line plots.

2. Tools for working with tables of processed results, keyed by gene ID (see
   `metaseq.results_table`).  An example use is to plot interactive
   scatterplots of RNA-seq results, with differentially
   expressed genes highlighted along with genes that contain a ChIP-seq peak
   within 1kb of their TSS.

3. A customizable, persistent, indexed mapping from gene ID -> genomic
   coordinates (implemented with `gffutils`) for combining the interval and
   ID-based components.

Other parts of `metaseq` combine these three components in other useful ways.
`metaseq` also:

* Provides parallel access for getting the genomic signal over many intervals
  simulataneously

* Offers a range of options for binning data or otherwise modifying the
  returned genomic signal

* Accepts bigBed, bigWig, BAM and SAM files as input for random access of
  continuous data (and with Tabix, accepts BED, GFF, GTF, VCF as well)

* Integrates RNA-seq and ChIP-seq data by using `gffutils` to convert between
  gene IDs and coordinates as well as handle full gene models and plotting.

* With the `metaseq.minibrowser` module, provides pure Python genome
  browser-like functionality in `matplotlib` figures.  Connecting these as
  callbacks to scatterplots will open a figure window showing the genomic
  region around the gene whose point was clicked.
