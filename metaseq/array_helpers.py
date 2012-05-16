import pybedtools
import numpy as np
import multiprocessing
import itertools
import pysam
import genomic_signal
import sys
from rebin import rebin, float_rebin
from helpers import chunker
import filetype_adapters


def _local_count(reader, feature, stranded=False):
    """
    The count of genomic signal (typcially BED features) found within an
    interval.

    Usually this only makes sense for BED or BAM (not bigWig) files.

    :param feature: pybedtools.Interval object
    :param stranded: If `stranded=True`, then only counts signal on the same
        strand as `feature`.
    """
    chrom = feature.chrom
    start = feature.start
    stop = feature.stop
    if stranded:
        strand = feature.strand
    else:
        strand = '.'

    count = 0
    for al in reader[feature]:
        if stranded and al.strand != strand:
            continue
        count += 1
    return count


def _local_coverage(reader, features, read_strand=None,
        fragment_size=1, shift_width=0, bins=None, use_score=False):
    """
    Computes a 1D vector of coverage at the coordinates for each feature in
    `features`, extending each read by `fragmentsize` bp.

    Depending on the arguments provided, this method can return a vector
    containing values from a single feature or from concatenated features.

    An example of the latter case: `features` can be a 3-tuple of
    pybedtools.Intervals representing (TSS + 1kb upstream, gene, TTS + 1kb
    downstream) and `bins` can be [100, 1000, 100].  This will return a vector
    of length 1200 containing the three genomic intervals binned into 100,
    1000, and 100 bins respectively.  Note that is up to the caller to
    construct the right axes labels in the final plot!

    :param features:
        A pybedtools.Interval object or an iterable of pybedtools.Interval
        objects.

    :param bins:
        * If `bins` is None, then each value in the returned array will
          correspond to one bp in the genome.

        * If `features` is a single Interval, then `bins` is an integer or
          None.

        * If `features` is an iterable of Intervals, `bins` is an iterable of
          integers of the same length as `features`.

    :param fragment_size:
        Integer. Each item from the genomic signal (e.g., reads from a BAM
        file) will be extended this many bp in the 3' direction.  Higher
        fragment sizes will result in smoother signal

    :param shift_width:
        Integer. Each item from the genomic signal (e.g., reads from a BAM
        file) will be shifted this number of bp in the 3' direction.  This can
        be useful for reconstructing a ChIP-seq profile, using the shift width
        determined from the peak-caller (e.g., modeled `d` in MACS)

    :param read_strand:
        String. If `read_strand` is one of "+" or "-", then only reads on that
        strand will be considered and reads on the opposite strand ignored.
        Useful for plotting genomic signal for stranded libraries


    :param use_score:
        If True, then each bin will contain the sum of the *score* attribute of
        genomic features in that bin instead of the *number* of genomic
        features falling within each bin.

    :rtype: NumPy array

    If a feature has a "-" strand attribute, then the resulting profile will be
    *relative to a minus-strand feature*.  That is, the resulting profile will
    be reversed.

    Returns arrays `x` and `y`.  `x` is in genomic coordinates, and `y` is
    the coverage at each of those coordinates after extending fragments.

    The total number of reads is guaranteed to be the same no matter how
    it's binned.

    (with ideas from
    http://www-huber.embl.de/users/anders/HTSeq/doc/tss.html)
    """

    if not (isinstance(features, list) or isinstance(features, tuple)):
        if bins is not None:
            if not isinstance(bins, int):
                raise ValueError(
                        "bins must be an int, got %s" % type(bins))
        features = [features]
        bins = [bins]
    else:
        if not len(bins) == len(features):
            raise ValueError(
                    "bins must have same length as feature list")

    # To keep speed high and memory low, there are different rebin funcs
    # depending on data type you need
    if not use_score:
        score = 1
        dtype = "int"
        rebin_func = rebin
    else:
        dtype = "float"
        rebin_func = float_rebin

    profiles = []
    xs = []
    for feature, nbin in zip(features, bins):
        chrom = feature.chrom
        start = feature.start
        stop = feature.stop
        strand = feature.strand

        # Extend the window to catch reads that would extend into the
        # requested window
        _fs = fragment_size or 0
        pos = pybedtools.Interval(
                chrom,
                max(start - _fs - shift_width, 0),
                stop + _fs + shift_width,
                )
        feature_size = stop - start

        # start off with an array of zeros to represent the feature
        profile = np.zeros(feature_size, dtype=dtype)

        for al in reader[pos]:

            if read_strand:
                if al.strand != read_strand:
                    continue

            # Shift fragment by modeled distance, if specified.
            if al.strand == '-':
                al.start -= shift_width
                al.stop -= shift_width
            else:
                al.start += shift_width
                al.stop += shift_width

            # Extend 3' by fragment size
            if fragment_size is not None:
                if al.strand == '-':
                    al.start = al.stop - fragment_size
                else:
                    al.stop = al.start + fragment_size

            # Convert to 0-based coords that can be used as indices into
            # array, making sure not to overflow the window.
            start_ind = al.start - start
            stop_ind = al.stop - start
            start_ind = max(start_ind, 0)
            stop_ind = min(stop_ind, feature_size)
            if start_ind >= feature_size or stop_ind < 0:
                continue
            # Finally, increment profile
            if use_score:
                score = float(al.score)

            profile[start_ind:stop_ind] += score

        # If no bins, return genomic coords
        if nbin is None:
            x = np.arange(start, stop)

        # Otherwise do the downsampling; resulting x is stll in genomic
        # coords
        else:

            #profile, x = signal.resample(profile, bins, x)
            xi = np.linspace(
                    start, stop - (stop - start) / float(nbin), nbin)
            profile = rebin_func(profile, nbin)
            x = xi

        # Minus-strand profiles should be flipped left-to-right.
        if strand == '-':
            profile = profile[::-1]
        xs.append(x)
        profiles.append(profile)

    return np.hstack(xs), np.hstack(profiles)


def _array_parallel(fn, cls, genelist, chunksize=25, processes=1, **kwargs):
    """
    Returns an array of genes in `genelist`, using `bins` bins.

    `genelist` is a list of pybedtools.Interval objects

    Splits `genelist` into pieces of size `chunksize`, creating an array
    for each chunk and merging ret

    A chunksize of 25-100 seems to work well on 8 cores.
    """
    pool = multiprocessing.Pool(processes)
    chunks = list(chunker(genelist, chunksize))

    # pool.map can only pass a single argument to the mapped function, so you
    # need this trick for passing multiple arguments; idea from
    # http://stackoverflow.com/questions/5442910/
    #               python-multiprocessing-pool-map-for-multiple-arguments
    #
    arrays = pool.map(
                _array_star,
                itertools.izip(
                    itertools.repeat(fn),
                    itertools.repeat(cls),
                    chunks,
                    itertools.repeat(kwargs)))
    stacked_arrays = np.row_stack(arrays)
    del arrays
    return stacked_arrays


def _array_star(args):
    """
    Unpacks the tuple `args` and calls _array.  Needed to pass multiple args to
    a pool.map-ed function
    """
    fn, cls, genelist, kwargs = args
    return _array(fn, cls, genelist, **kwargs)


def _array(fn, cls, genelist, **kwargs):
    """
    Returns a "meta-feature" array, with len(genelist) rows and `bins`
    cols.  Each row contains the number of reads falling in each bin of
    that row's modified feature.
    """
    reader = cls(fn)
    biglist = []
    if 'bins' in kwargs:
        if isinstance(kwargs['bins'], int):
            kwargs['bins'] = [kwargs['bins']]

    for gene in genelist:
        if not isinstance(gene, (list, tuple)):
            gene = [gene]
        coverage_x, coverage_y = _local_coverage(reader, gene,
                                                     **kwargs)
        biglist.append(coverage_y)
    return np.array(biglist)


def _local_coverage_bigwig(bigwig, features, bins=None):
    """
    Returns matrix of coverage of `features` using `bins` -- see
    :func:`metaseq.array_helpers._local_coverage` for more info.
    """

    if not (isinstance(features, list) or isinstance(features, tuple)):
        if bins is not None:
            if not isinstance(bins, int):
                raise ValueError(
                        "bins must be an int, got %s" % type(bins))
        features = [features]
        bins = [bins]
    else:
        if not len(bins) == len(features):
            raise ValueError(
                    "bins must have same length as feature list")

    profiles = []
    xs = []
    for feature, nbin in zip(features, bins):
        chrom = feature.chrom
        start = feature.start
        stop = feature.stop

        s = bigwig.summarize_from_full(chrom, start, stop, nbin)
        x = np.linspace(s.start, s.end, s.size)
        s.sum_data = np.ma.masked_where(s.sum_data == 0,
                                        s.sum_data,
                                        copy=False)

        profile = s.sum_data / s.valid_count

        if feature.strand == '-':
            profile = profile[::-1]

        xs.append(x)
        profiles.append(profile)

    return np.hstack(xs), np.hstack(profiles)
