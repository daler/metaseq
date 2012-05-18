Inputs
======
:mod:`metachip` is useful for exploring relationships between intervals (e.g.,
annotated genes; ChIP-seq peaks) and numeric data (e.g., ChIP-seq signal;
RNA-seq expression data).

Imagine you did a ChIP-seq experiment on your favorite protein, and you have
binding sites.  Then you do an RNA-seq experiment where you knock down that
protein.  You'd like to know things like:

* which of the differentially expressed genes have a peak in near their promoter
* what the genome-wide average ChIP-seq signal at gene TSSs looks like compared to the
  average ChIP-seq signal at TTSs
* the same thing, but only on differentially expressed genes
* the same thing, but only on genes down-regulated upon knockdown
* different "promoter region" definitions

All of this can be done in other ways, or by tying together multiple scripts.
But :mod:`metachip` allows all of this to be done within one framework.


:mod:`metachip` :
* ties together several other packages for manipulating genomic features,
  :mod:`pybedtools` and :mod:`gffutils`.
* Supplies a fast "rebinning" routine that is compiled in C for speed
* Accepts BAM, SAM, BED, GTF, GFF, VCF files as input for intervals
* Accepts bigBed, bigWig, BAM files as input for continuous data (and with
  Tabix, accepts BED, GFF, GTF, VCF as well)
