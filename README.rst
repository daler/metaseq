Metaseq
=======

The main documentation for `metaseq` can be found at http://packages.python.org/metaseq/

Briefly, the goal of `metaseq` is to tie together lots of existing software into
a framework for exploring genomic data.  It focuses on flexibility and
interactive exploration and plotting of disparate genomic data sets.

The example session shows the the process of creating a plot like this,
starting with BAM and GTF files:

.. figure:: demo.png

    Top: Heatmap of ATF3 ChIP-seq signal over transcription start sites (TSS) on
    chr17 in human K562 cells.  Middle: average ChIP enrichment over all TSSs
    +/- 1kb, with 95% CI band.  Bottom: Integration with ATF3 knockdown RNA-seq
    results, showing differential enrichment over transcripts that went up,
    down, or were unchanged upon ATF3 knockdown.

`metaseq` offers:

* A format-agnostic API for accessing "genomic signal" that allows you to work
  with BAM, BED, VCF, GTF, GFF, bigBed, and bigWig using the same API.

* Parallel data access from the file formats mentioned above

* "Mini-browsers", zoomable and pannable Python-only  figures that show genomic
  signal and gene models and are spawned by clicking on features of interest

* A wrapper around pandas.DataFrames to simplify working with tabular results
  data, (think DESeq results tables)

* Integrates data keyed by genomic interval (think BAM or BED files) with data
  keyed by gene ID (e.g., Cufflinks or DESeq results tables)

