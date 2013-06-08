from copy import copy
import os
import sys
import time
import scipy
import matplotlib
import pybedtools
import subprocess
import re

coord_re = re.compile(
    r"""
    (?P<chrom>.+):
    (?P<start>\d+)-
    (?P<stop>\d+)
    (?:\[(?P<strand>.)\])?""", re.VERBOSE)

def gfffeature_to_interval(feature):
    return pybedtools.create_interval_from_list(feature.tostring().split('\t'))


def chunker(f, n):
    """
    Utility function to split iterable `f` into `n` chunks
    """
    f = iter(f)
    x = []
    while 1:
        if len(x) < n:
            try:
                x.append(f.next())
            except StopIteration:
                if len(x) > 0:
                    yield tuple(x)
                break
        else:
            yield tuple(x)
            x = []


def nice_colormap(z, invalid=None):
    """
    Dynamically scales the midpoint to the median of the positive values.

    Returns a colormap ready for imshow or pcolor or whatever.
    """
    norm = matplotlib.colors.Normalize()
    z = z.copy()
    norm(z)

    # Set max to 99th percentile
    norm.vmax = scipy.stats.scoreatpercentile(z.ravel(), 99)
    zeropoint = norm(0)

    # split from zero to max(z) into chunks
    dcolor = (1 - zeropoint) / 3

    # midpoint of color change is median of positive values.
    medpoint = norm(scipy.stats.scoreatpercentile(z[z > 0].ravel(), 50))

    # construct that bitch
    cdict = {
               'red':  ((0.0, 0.0, 0.0),
                       (zeropoint, 1.0, 1.0),
                       (medpoint, 1., 1.0),
                       (1.0, 1.0, 1.0)),

             'green': ((0.0, 0.0, 0.0),
                       (zeropoint, 1.0, 1.0),
                       (medpoint, .9, .9),
                       (1.0, 0.0, 0.0)),

             'blue':  ((0.0, 0.0, 1.0),
                       (zeropoint, 1.0, 1.0),
                       (medpoint, .0, .0),
                       (1.0, 0.0, 0.0))
            }
    cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap', cdict, 256)

    # NaNs and stuff will be medium gray.
    if invalid is not None:
        cmap.set_bad(invalid, 1.0)

    return norm, cmap


def bam2bigwig(bam, bigwig, genome, scale=1e6, verbose=False):
    """
    Uses BEDTools to go from BAM to bedgraph, then bedGraphToBigWig to get the
    final bigwig.
    """
    if scale is not None:
        cmds = ['samtools', 'view', '-F', '0x4', '-c', bam]
        p = subprocess.Popen(cmds, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
        stdout, stderr = p.communicate()
        total_reads = float(stdout)
        reads_per_scale = total_reads / scale
        if verbose:
            sys.stderr.write('%s total reads\n' % total_reads)
            sys.stderr.flush()

    chromsizes = pybedtools.chromsizes_to_file(pybedtools.chromsizes(genome))

    t0 = time.time()
    bedgraph = pybedtools.BedTool(bam)\
            .genome_coverage(bg=True, g=chromsizes, scale=scale)\
            .moveto('bedgraph.bedgraph')
    print bedgraph.fn
    if verbose:
        sys.stderr.write('Completed bedGraph in %.1fs\n' % (time.time() - t0))
        sys.stderr.flush()

    cmds = ['bedGraphToBigWig', bedgraph.fn, chromsizes, bigwig]
    p = subprocess.Popen(cmds, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()

    if verbose:
        sys.stderr.write('Completed bigWig %s\n' % bigwig)
        sys.stderr.flush()


def bedgraph2bigwig(bedgraph, bigwig, genome, verbose=False):
    """
    Create a bigWig from `bedgraph`.

    :param bedgraph: Input filename of bedgraph
    :param bigwig: Output filename of bigWig to create
    :param genome: String assembly name of genome
    :param verbose: Print messages to stderr
    """
    chromsizes = pybedtools.chromsizes_to_file(pybedtools.chromsizes(genome))
    cmds = ['bedGraphToBigWig', bedgraph, chromsizes, bigwig]
    p = subprocess.Popen(cmds, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    if verbose:
        sys.stderr.write('Completed bigWig %s\n' % bigwig)
        sys.stderr.flush()

def data_dir():
    """
    Returns the data directory that contains example files for tests and
    documentation.
    """
    return os.path.join(os.path.dirname(__file__), 'test', 'data')


def example_filename(fn):
    """
    Return a bed file from the pybedtools examples directory.  Use
    :func:`list_example_files` to see a list of files that are included.
    """
    fn = os.path.join(data_dir(), fn)
    if not os.path.exists(fn):
        raise ValueError("%s does not exist" % fn)
    return fn


def split_feature(f, n):
    """
    Split an interval into `n` roughly equal portions
    """
    if not isinstance(n, int):
        raise ValueError('n must be an integer')
    orig_feature = copy(f)
    step = (f.stop - f.start) / n
    for i in range(f.start, f.stop, step):
        f = copy(orig_feature)
        start = i
        stop = min(i + step, orig_feature.stop)
        f.start = start
        f.stop = stop
        yield f
        if stop == orig_feature.stop:
            break

def tointerval(s):
    """
    If string, then convert to an interval; otherwise just return the input
    """
    if isinstance(s, basestring):
        m = coord_re.search(s)
        if m.group('strand'):
            return pybedtools.create_interval_from_list([
                m.group('chrom'),
                m.group('start'),
                m.group('stop'),
                '.',
                '0',
                m.group('strand')])
        else:
            return pybedtools.create_interval_from_list([
                m.group('chrom'),
                m.group('start'),
                m.group('stop'),
            ])
    return s
