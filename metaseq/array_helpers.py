import pybedtools
import numpy as np
import multiprocessing
import itertools
import pysam
import genomic_signal
import sys
from rebin import rebin
from helpers import chunker
import helpers
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
    if isinstance(feature, basestring):
        # assume it's in chrom:start-stop format
        chrom, coords = feature.split(':')
        start, stop = coords.split('-')
        feature = pybedtools.create_interval_from_list([chrom, start, stop])
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


def _local_coverage(reader, features, read_strand=None, fragment_size=None,
                    shift_width=0, bins=None, use_score=False, accumulate=True,
                    preserve_total=False, verbose=False):
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

    Parameters
    ----------
    features : various
        Can be a single interval, an iterable yielding intervals, or an
        iterable-of-iterables.

        Intervals must have chrom, start, and stop attributes.

        If a single interval, then return a 1-D array for that interval.

        If an iterable of single intervals, then return an array, one row for
        each interval.

        If an iterable of iterables, then the number of nested intervals must
        match the number of bins provided.

    bins : None, int, list
        * If `bins` is None, then each value in the returned array will
          correspond to one bp in the genome.

        * If `features` is a single Interval, then `bins` is an integer or
          None.

        * If `features` is an iterable of Intervals, `bins` is an iterable of
          integers of the same length as `features`.

    fragment_size : None or int
        If not None, then each item from the genomic signal (e.g., reads from
        a BAM file) will be extended `fragment_size` bp in the 3' direction.
        Higher fragment sizes will result in smoother signal.

    shift_width : int
        Each item from the genomic signal (e.g., reads from a BAM
        file) will be shifted `shift_width` bp in the 3' direction.  This can
        be useful for reconstructing a ChIP-seq profile, using the shift width
        determined from the peak-caller (e.g., modeled `d` in MACS)

    read_strand : None or str
        If `read_strand` is one of "+" or "-", then only items from the genomic
        signal (e.g., reads from a BAM file) on that strand will be considered
        and reads on the opposite strand ignored.  Useful for plotting genomic
        signal for stranded libraries

    use_score : bool
        If True, then each bin will contain the sum of the *score* attribute of
        genomic features in that bin instead of the *number* of genomic
        features falling within each bin.

    accumulate : bool
        If False, then only record *that* there was something there, rather
        than acumulating reads.  This is useful for making matrices with called
        peaks.

    preserve_total : bool
        If True, re-scales the returned value so that each binned row's total
        is equal to the sum of the original, un-binned data.  The units of the
        returned array will be in "total per bin".  This is useful for, e.g.,
        counting reads in features.  If `preserve_total` is False, then the
        returned array will have units of "density"; this is more generally
        useful and is the default behavior.

    :rtype: NumPy array

    If a feature has a "-" strand attribute, then the resulting profile will be
    *relative to a minus-strand feature*.  That is, the resulting profile will
    be reversed.

    Returns arrays `x` and `y`.  `x` is in genomic coordinates, and `y` is
    the coverage at each of those coordinates after extending fragments.

    The total number of reads is guaranteed to be the same no matter how it's
    binned.

    (with ideas from
    http://www-huber.embl.de/users/anders/HTSeq/doc/tss.html)

    """
    if isinstance(features, basestring):
        features = helpers.tointerval(features)

    if not ((isinstance(features, list) or isinstance(features, tuple))):
        if bins is not None:
            if not isinstance(bins, int):
                raise ValueError(
                    "bins must be an int, got %s" % type(bins))
        features = [features]
        bins = [bins]
    else:
        if bins is None:
            bins = [None for i in features]
        if not len(bins) == len(features):
            raise ValueError(
                "bins must have same length as feature list")
    # nomenclature:
    #   "window" is region we're getting data for
    #   "alignment" is one item in that region
    #
    profiles = []
    xs = []
    for window, nbin in zip(features, bins):
        window = helpers.tointerval(window)
        chrom = window.chrom
        start = window.start
        stop = window.stop
        strand = window.strand

        # Extend the window to catch reads that would extend into the
        # requested window
        _fs = fragment_size or 0
        padded_window = pybedtools.Interval(
            chrom,
            max(start - _fs - shift_width, 0),
            stop + _fs + shift_width,
        )
        window_size = stop - start

        # start off with an array of zeros to represent the window
        profile = np.zeros(window_size, dtype=float)

        for interval in reader[padded_window]:

            if read_strand:
                if interval.strand != read_strand:
                    continue

            # Shift interval by modeled distance, if specified.
            if shift_width:
                if interval.strand == '-':
                    interval.start -= shift_width
                    interval.stop -= shift_width
                else:
                    interval.start += shift_width
                    interval.stop += shift_width

            # Extend fragment size from 3'
            if fragment_size:
                if interval.strand == '-':
                    interval.start = interval.stop - fragment_size
                else:
                    interval.stop = interval.start + fragment_size

            # Convert to 0-based coords that can be used as indices into
            # array
            start_ind = interval.start - start

            # If the feature goes out of the window, then only include the part
            # that's inside the window
            start_ind = max(start_ind, 0)

            # Same thing for stop
            stop_ind = interval.stop - start
            stop_ind = min(stop_ind, window_size)

            # Skip if the feature is shifted outside the window. This can
            # happen with large values of `shift_width`.
            if start_ind >= window_size or stop_ind < 0:
                continue

            # Finally, increment profile
            if use_score:
                score = float(interval.score)
            else:
                score = 1

            if accumulate:
                if verbose:
                    print '%s-%s += %s' % (start_ind, stop_ind, score)
                profile[start_ind:stop_ind] += score
            else:
                profile[start_ind:stop_ind] = score

        # If no bins, return genomic coords
        if nbin is None:
            x = np.arange(start, stop)

        # Otherwise do the downsampling; resulting x is stll in genomic
        # coords
        else:

            #xi = np.linspace(
            #        start, stop - (stop - start) / float(nbin), nbin)
            xi, profile = rebin(x=np.arange(start, stop), y=profile, nbin=nbin)
            if not accumulate:
                nonzero = profile != 0
                profile[profile != 0] = 1
            x = xi

        # Minus-strand profiles should be flipped left-to-right.
        if strand == '-':
            profile = profile[::-1]
        xs.append(x)
        if preserve_total:
            if nbin is not None:
                scale = window_size / float(nbin)
                profile *= scale
        profiles.append(profile)

    return np.hstack(xs), np.hstack(profiles)


def _local_coverage_bigwig(bigwig, features, bins=None, accumulate=True,
                           preserve_total=False):
    """
    Returns matrix of coverage of `features` using `bins` -- see
    :func:`metaseq.array_helpers._local_coverage` for more info.
    """
    if isinstance(features, basestring):
        features = helpers.tointerval(features)
    if not (isinstance(features, list) or isinstance(features, tuple)):
        if bins is not None:
            if not isinstance(bins, int):
                raise ValueError(
                    "bins must be an int, got %s" % type(bins))
        features = [features]
        bins = [bins]
    else:
        if bins is None:
            bins = [None for i in features]
        if not len(bins) == len(features):
            raise ValueError(
                "bins must have same length as feature list")

    profiles = []
    xs = []
    for window, nbin in zip(features, bins):
        window = helpers.tointerval(window)
        chrom = window.chrom
        start = window.start
        stop = window.stop
        strand = window.strand
        profile = bigwig.summarize(window, bins=(nbin or len(window)))


        # If no bins, return genomic coords
        if nbin is None:
            x = np.arange(start, stop)

        # Otherwise do the downsampling; resulting x is stll in genomic
        # coords
        else:

            #xi = np.linspace(
            #        start, stop - (stop - start) / float(nbin), nbin)
            xi, profile = rebin(x=np.arange(start, stop), y=profile, nbin=nbin)
            x = xi
        if not accumulate:
            nonzero = profile != 0
            profile[profile != 0] = 1

        # Minus-strand profiles should be flipped left-to-right.
        if strand == '-':
            profile = profile[::-1]
        xs.append(x)
        if preserve_total:
            scale = window_size / float(nbin)
            profile *= scale
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
    return pool.map(
        _array_star,
        itertools.izip(
            itertools.repeat(fn),
            itertools.repeat(cls),
            chunks,
            itertools.repeat(kwargs)))


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
    _local_coverage_func = cls.local_coverage
    biglist = []
    if 'bins' in kwargs:
        if isinstance(kwargs['bins'], int):
            kwargs['bins'] = [kwargs['bins']]

    for gene in genelist:
        if not isinstance(gene, (list, tuple)):
            gene = [gene]
        coverage_x, coverage_y = _local_coverage_func(
            reader, gene, **kwargs)

        biglist.append(coverage_y)
    return biglist
