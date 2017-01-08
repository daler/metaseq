import pybedtools
import numpy as np
import multiprocessing
import itertools
import pysam
import sys
from helpers import chunker
import helpers
from helpers import rebin
import filetype_adapters


class ArgumentError(Exception):
    pass


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
        feature = helpers.tointerval(feature)
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
                    preserve_total=False, method=None, function="mean",
                    zero_inf=True, zero_nan=True, processes=None,
                    stranded=True, verbose=False):
    """
    Returns a binned vector of coverage.

    Computes a 1D vector of coverage at the coordinates for each feature in
    `features`, extending each read by `fragmentsize` bp.

    Some arguments cannot be used for bigWig files due to the structure of
    these files.  The parameters docstring below indicates whether or not an
    argument can be used with bigWig files.

    Depending on the arguments provided, this method can return a vector
    containing values from a single feature or from concatenated features.

    An example of the flexibility afforded by the latter case:

        `features` can be a 3-tuple of pybedtools.Intervals representing (TSS
        + 1kb upstream, gene, TTS + 1kb downstream) and `bins` can be [100,
        1000, 100].  This will return a vector of length 1200 containing the
        three genomic intervals binned into 100, 1000, and 100 bins
        respectively.  Note that is up to the caller to construct the right
        axes labels in the final plot!

    Parameters
    ----------
    features : str, interval-like object, or list

        Can be a single interval or an iterable yielding intervals.

        Interval-like objects must have chrom, start, and stop attributes, and
        optionally a strand attribute.  One exception to this that if
        `features` is a single string, it can be of the form "chrom:start-stop"
        or "chrom:start-stop[strand]".

        If `features` is a single interval, then return a 1-D array for that
        interval.

        If `features` is an iterable of intervals, then return a 1-D
        array that is a concatenation of signal for these intervals.

        Available for bigWig.

    bins : None, int, list
        If `bins` is None, then each value in the returned array will
        correspond to one bp in the genome.

        If `features` is a single Interval, then `bins` is an integer or None.

        If `features` is an iterable of Intervals, `bins` is an iterable of
        integers of the same length as `features`.

        Available for bigWig.

    fragment_size : None or int
        If not None, then each item from the genomic signal (e.g., reads from
        a BAM file) will be extended `fragment_size` bp in the 3' direction.
        Higher fragment sizes will result in smoother signal.  Not available
        for bigWig.

    shift_width : int
        Each item from the genomic signal (e.g., reads from a BAM
        file) will be shifted `shift_width` bp in the 3' direction.  This can
        be useful for reconstructing a ChIP-seq profile, using the shift width
        determined from the peak-caller (e.g., modeled `d` in MACS). Not
        available for bigWig.

    read_strand : None or str
        If `read_strand` is one of "+" or "-", then only items from the genomic
        signal (e.g., reads from a BAM file) on that strand will be considered
        and reads on the opposite strand ignored.  Useful for plotting genomic
        signal for stranded libraries. Not available for bigWig.

    stranded : bool
        If True, then the profile will be reversed for features whose strand
        attribute is "-".

    use_score : bool
        If True, then each bin will contain the sum of the *score* attribute of
        genomic features in that bin instead of the *number* of genomic
        features falling within each bin. Not available for bigWig.

    accumulate : bool
        If False, then only record *that* there was something there, rather
        than acumulating reads.  This is useful for making matrices with called
        peaks. Available for bigWig.

    preserve_total : bool
        If True, re-scales the returned value so that each binned row's total
        is equal to the sum of the original, un-binned data.  The units of the
        returned array will be in "total per bin".  This is useful for, e.g.,
        counting reads in features.  If `preserve_total` is False, then the
        returned array will have units of "density"; this is more generally
        useful and is the default behavior.  Available for bigWig, but not when
        using method="ucsc_summarize".

    method : str; one of [ "summarize" | "get_as_array" | "ucsc_summarize" ]
        Only used for bigWig.  The method specifies how data are extracted from
        the bigWig file.  "summarize" is the default.  It's quite fast, but may
        yield slightly different results when compared to running this same
        function on the BAM file from which the bigWig was created.

        "summarize" uses bx-python.  The values returned will not be exactly
        the same as the values returned when local_coverage is called on a BAM,
        BED, or bigBed file, but they will be close.  This method is quite
        fast, and is the default when bins is not None.

        "get_as_array" uses bx-python, but does a separate binning step.  This
        can be slower than the other two methods, but the results are exactly
        the same as those from a BAM, BED, or bigBed file.  This method is
        always used if bins=None.

        "ucsc_summarize" is an alternative version of "summarize".  It uses the
        UCSC program `bigWigSummary`, which must already installed and on your
        path.

    function : str; one of ['sum' | 'mean' | 'min' | 'max' | 'std']
        Determine the nature of the values returned. Only valid if `method` is
        "summarize" or "ucsc_summarize", which also implies bigWig. Default is
        "mean". If `method="ucsc_summarize", then there is an additional option
        for function, "coverage", which returns the percent of region that is
        covered.

    zero_inf, zero_nan : bool
        Only used for bigWig. If either are True, sets any missing or inf
        values to zero before returning.

        If `method="ucsc_summarize"`, missinv values are always reported as
        zero. If `method="get_as_array"`, missing values always reported as
        nan.

        Values can be -inf, inf, or nan for missing values when
        `method="summarize"` according to the following table:

        ========== ========================
        `function` missing values appear as
        ========== ========================
        "sum"      0
        "mean"     nan
        "min"      inf
        "max"      -inf
        "std"      nan
        ========== ========================

    processes : int or None
        The feature can be split across multiple processes.

    Returns
    -------

    1-d NumPy array


    Notes
    -----
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
    # bigWig files are handled differently, so we need to know if we're working
    # with one; raise exeception if a kwarg was supplied that's not supported.
    if isinstance(reader, filetype_adapters.BigWigAdapter):
        is_bigwig = True
        defaults = (
            ('read_strand', read_strand, None),
            ('fragment_size', fragment_size, None),
            ('shift_width', shift_width, 0),
            ('use_score', use_score, False),
            ('preserve_total', preserve_total, False),
        )
        for name, check, default in defaults:
            if (
                ((default is None) and (check is not default))
                or
                (check != default)
            ):
                raise ArgumentError(
                    "Argument '%s' not supported for bigWig" % name)

        if method == 'ucsc_summarize':
            if preserve_total:
                raise ArgumentError(
                    "preserve_total=True not supported when using "
                    "method='ucsc_summarize'")
    else:
        is_bigwig = False

    if isinstance(reader, filetype_adapters.BamAdapter):
        if use_score:
            raise ArgumentError("Argument 'use_score' not supported for "
                                "bam")

    # e.g., features = "chr1:1-1000"
    if isinstance(features, basestring):
        features = helpers.tointerval(features)

    if not ((isinstance(features, list) or isinstance(features, tuple))):
        if bins is not None:
            if not isinstance(bins, int):
                raise ArgumentError(
                    "bins must be an int, got %s" % type(bins))
        features = [features]
        bins = [bins]
    else:
        if bins is None:
            bins = [None for i in features]
        if not len(bins) == len(features):
            raise ArgumentError(
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

        if not is_bigwig:
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

                # If the feature goes out of the window, then only include the
                # part that's inside the window
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
                    if preserve_total:
                        profile[start_ind:stop_ind] += (
                            score / float((stop_ind - start_ind)))
                    else:
                        profile[start_ind:stop_ind] += score

                else:
                    profile[start_ind:stop_ind] = score

        else:  # it's a bigWig
            profile = reader.summarize(
                window,
                method=method,
                function=function,
                bins=(nbin or len(window)),
                zero_inf=zero_inf,
                zero_nan=zero_nan,
                )

        # If no bins, return genomic coords
        if (nbin is None):
            x = np.arange(start, stop)

        # Otherwise do the downsampling; resulting x is stll in genomic
        # coords
        else:
            if preserve_total:
                total = float(profile.sum())
            if not is_bigwig or method == 'get_as_array':
                xi, profile = rebin(
                    x=np.arange(start, stop), y=profile, nbin=nbin)
                if not accumulate:
                    nonzero = profile != 0
                    profile[profile != 0] = 1
                x = xi

            else:
                x = np.linspace(start, stop - 1, nbin)

        # Minus-strand profiles should be flipped left-to-right.
        if stranded and strand == '-':
            profile = profile[::-1]
        xs.append(x)
        if preserve_total:
            scale = profile.sum() / total
            profile /= scale
        profiles.append(profile)

    stacked_xs = np.hstack(xs)
    stacked_profiles = np.hstack(profiles)
    del xs
    del profiles
    return stacked_xs, stacked_profiles


def _array_parallel(fn, cls, genelist, chunksize=250, processes=1, **kwargs):
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
    results = pool.map(
        func=_array_star,
        iterable=itertools.izip(
            itertools.repeat(fn),
            itertools.repeat(cls),
            chunks,
            itertools.repeat(kwargs)))
    pool.close()
    pool.join()
    return results

def _count_array_parallel(fn, cls, genelist, chunksize=250, processes=1, **kwargs):
    pool = multiprocessing.Pool(processes)
    chunks = list(chunker(genelist, chunksize))
    results = pool.map(
        func=_count_array_star,
        iterable=itertools.izip(
            itertools.repeat(fn),
            itertools.repeat(cls),
            chunks,
            itertools.repeat(kwargs)))
    pool.close()
    pool.join()
    return results


def _count_array_star(args):
    fn, cls, genelist, kwargs = args
    return _count_array(fn, cls, genelist, **kwargs)

def _count_array(fn, cls, genelist, **kwargs):
    reader = cls(fn)
    _local_count_func = cls.local_count
    biglist = []
    for gene in genelist:
        c = _local_count_func(
            reader, gene, **kwargs)
        biglist.append(c)
    return biglist

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
