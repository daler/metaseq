from copy import copy
import os
import sys
import glob
import time
import matplotlib
import pybedtools
import subprocess
import re
import numpy as np


def rebin(x, y, nbin):
    xi = np.linspace(x.min(), x.max(), nbin)
    return xi, np.interp(xi, x, y)


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


def list_example_files(pattern=None, full_path=False):
    if pattern is None:
        files = os.listdir(data_dir())
    else:
        files = glob.glob(os.path.join(data_dir(), pattern))
    if full_path:
        return files
    else:
        return [os.path.basename(i) for i in files]


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

coord_re = re.compile(
    r"""
    (?P<chrom>.+):
    (?P<start>\d+)-
    (?P<stop>\d+)
    (?:\[(?P<strand>.)\])?""", re.VERBOSE)


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
