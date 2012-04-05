Metaseq
=======
`metaseq` is a framework that makes it easy to work with a combination of
ChIP-seq, RNA-seq, or anything-seq.

Example use-cases:

    * Matrix plots:

        Generate an MxN matrix (e.g., M = genes and N = bins, TSS +/- 1kb) showing the
        ChIP-seq signal of of up-regulated genes from an RNA-seq experiment
        compared to down-regulated genes

        The signal is calculated directly from the IP and input BAM files, and
        scaling to library size normalizes the data

        Since the result is a NumPy array, it's trivial to, say, sort by gene
        expression or to take column averages.

    * Cluster genes based on the spatial distribution of ChIP-seq peaks around
      their TSSs.

        A matrix can be generated using the same interface as for BAM files,
        using a BED or bigBed file instead, again resulting in a NumPy array
        for further analysis

    * Scatter plot of DESeq results (basemeana vs basemeanb) where points are
      colored according to the number of ChIP peaks in the gene

        Using callback functions in matplotlib, clicking on a point can print
        gene information, or can even open a window showing a mini
        genome-browser view of that gene

    * Pie charts of where peaks fall within annotated genes -- TSS, poly-A
      site, intron, exon, etc

Where possible, the inputs are standard formats -- BED, GFF, GTF, BAM, SAM,
DESeq results as saved from R, or even arbitrary tab-delimited data files that
have a header.  If you take the time to convert to bigWig or bigBed,
performance will be improved

`metaseq` stands on the shoulders of a large body of existing Python packages
to provide a flexible, intuitive, and powerful framework.  These packages
include:

    * `pysam` for parsing BAM files

    * `bx-python` for access to bigWig and bigBed files

    * `matplotlib` for fast, interactive, and extremely flexible plotting

    * `numpy` for fast arrays

    * `scikits.learn` for clustering (specifically, the minibatch k-means
      implementation)

    * `pybedtools`, an interface to the `BEDTools` suite, for flexible and
      powerful "genome algebra" and feature-level manipulation

    * `gffutils`, a lightweight database framework for navigating the hierarchy
      of gene annotations (exon -> transcript -> gene) in a portable format

    * `Cython` for implementing computationally intensive code in C

In addition to providing the "glue" between these various packages, `metaseq`
also supplies binning code written in Cython and helpers for parallelizing data
access.
