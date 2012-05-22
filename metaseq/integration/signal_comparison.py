"""
Module for working with comparisons between two genomic_signal objects, e.g.,
genome-wide correlation
"""
import sys
import pybedtools
import itertools
import numpy as np


def compare(signal1, signal2, features, outfn, comparefunc=np.subtract,
        batchsize=5000, array_kwargs=None, verbose=False):
    """
    Compares two genomic signal objects and outputs results as a bedGraph file.
    Can be used for entire genome-wide comparisons due to its parallel nature.

    Typical usage would be to create genome-wide windows of equal size to
    provide as `features`::

        windowsize = 10000
        features = pybedtools.BedTool().window_maker(
           genome='hg19', w=windowsize)

    You will usually want to choose bins for the array based on the final
    resolution you would like. Say you would like 10-bp bins in the final
    bedGraph; using the example above you would use array_kwargs={'bins':
    windowsize/10}.  Or, for single-bp resolution (beware: file will be large),
    use {'bins': windowsize}.

    Here's how it works.  This function:

        * Takes `batchsize` features at a time from `features`

        * Constructs normalized (RPMMR) arrays in parallel for each input
          genomic signal object for those `batchsize` features

        * Applies `comparefunc` (np.subtract by default) to the arrays to get
          a "compared" (e.g., difference matrix by default) for the `batchsize`
          features.

        * For each row in this matrix, it outputs each nonzero column as
          a bedGraph format line in `outfn`

    `comparefunc` is a function with the signature::

        def f(x, y):
            return z

    where `x` and `y` will be arrays for `signal1` and `signal2` (normalized to
    RPMMR) and `z` is a new array.  By default this is np.subtract, but another
    common `comparefunc` might be a log2-fold-change function::

        def lfc(x, y):
            return np.log2(x / y)

    :param signal1: A genomic_signal object
    :param signal2: Another genomic_signal object
    :param features: An iterable of pybedtools.Interval objects. A list will be
        created for every `batchsize` features, so you need enough memory for
        this.
    :param comparefunc: Function to use to compare arrays (default is
        np.subtract)
    :param outfn: String filename to write bedGraph file
    :param batchsize: Number of features (each with length `windowsize` bp) to
        process at a time
    :param array_kwargs: Kwargs passed directly to genomic_signal.array.  Needs
        `processes` and `chunksize` if you want parallel processing
    :param verbose: Be noisy
    """
    fout = open(outfn, 'w')
    fout.write('track type=bedGraph\n')

    i = 0
    this_batch = []
    for feature in features:
        if i <= batchsize:
            this_batch.append(feature)
            i += 1
            continue

        if verbose:
            print 'working on batch of %s' % batchsize
            sys.stdout.flush()

        arr1 = signal1.array(this_batch, **array_kwargs).astype(float)
        arr2 = signal2.array(this_batch, **array_kwargs).astype(float)
        arr1 /= signal1.million_mapped_reads()
        arr2 /= signal2.million_mapped_reads()
        compared = comparefunc(arr1, arr2)

        for feature, row in itertools.izip(this_batch, compared):
            start = feature.start
            bins = len(row)
            binsize = len(feature) / len(row)

            # Quickly move on if nothing here.  speed increase prob best for
            # sparse data
            if sum(row) == 0:
                continue

            for j in range(0, len(row)):
                score = row[j]
                stop = start + binsize
                if score != 0:
                    fout.write('\t'.join([
                        feature.chrom,
                        str(start),
                        str(stop),
                        str(score)]) + '\n')
                start = start + binsize
        this_batch = []
        i = 0
    fout.close()


if __name__ == "__main__":
    import metaseq
    ip_bam = metaseq.genomic_signal(
            metaseq.example_filename(
                'wgEncodeUwTfbsK562CtcfStdAlnRep1.bam'), 'bam')
    control_bam = metaseq.genomic_signal(
            metaseq.example_filename(
                'wgEncodeUwTfbsK562InputStdAlnRep1.bam'), 'bam')

    BINSIZE = 10
    WINDOWSIZE = 10000
    BINS = WINDOWSIZE / BINSIZE
    features = pybedtools.BedTool()\
            .window_maker(genome='hg19', w=WINDOWSIZE)\
            .filter(lambda x: x.chrom == 'chr19')

    result = compare(
            signal1=ip_bam,
            signal2=control_bam,
            features=features,
            outfn='diffed.bedgraph',
            array_kwargs=dict(bins=BINS, processes=6, chunksize=50),
            verbose=True)
