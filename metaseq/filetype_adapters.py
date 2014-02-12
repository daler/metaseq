"""
This module provides classes that make a file format conform to a uniform API.
These are not generally needed by end-users, rather, they are used internally
by higher-level code like :mod:`metaseq.genomic_signal`.

File-type adapters accept a filename of the appropriate format (which is not
checked) as the only argument to their constructor.

Subclasses must define __getitem__ to accept a pybedtools.Interval and return
an iterator of pybedtools.Intervals

Subclasses must define make_fileobj(), which returns an object to be iterated
over in __getitem__
"""
from bx.bbi.bigbed_file import BigBedFile
from bx.bbi.bigwig_file import BigWigFile
import numpy as np
import subprocess
import pysam
import pybedtools
import os
import sys

strand_lookup = {16: '-', 0: '+'}


class BaseAdapter(object):
    """
    Base class for filetype adapters
    """
    def __init__(self, fn):
        self.fn = fn
        self.fileobj = None
        self.fileobj = self.make_fileobj()

    def __getitem__(self, key):
        raise ValueError('Subclasses must define __getitem__')

    def make_fileobj(self):
        raise ValueError('Subclasses must define make_fileobj')


class BamAdapter(BaseAdapter):
    """
    Adapter that provides random access to BAM objects using Pysam
    """
    def __init__(self, fn):
        super(BamAdapter, self).__init__(fn)

    def make_fileobj(self):
        return pysam.Samfile(self.fn, 'rb')

    def __getitem__(self, key):
        iterator = self.fileobj.fetch(
            key.chrom,
            key.start,
            key.stop)
        for r in iterator:
            start = r.pos
            curr_end = r.pos
            for op, bp in r.cigar:
                start = curr_end
                curr_end += bp
                if op == 0:
                    interval = pybedtools.Interval(
                        self.fileobj.references[r.rname],
                        start,
                        curr_end,
                        strand=strand_lookup[r.flag & 0x0010])
                    interval.file_type = 'bed'
                    yield interval


class BedAdapter(BaseAdapter):
    """
    Adapter that provides random access to BED files via Tabix
    """
    def __init__(self, fn):
        super(BedAdapter, self).__init__(fn)

    def make_fileobj(self):
        obj = pybedtools.BedTool(self.fn)
        if not obj._tabixed():
            obj = obj.sort().tabix(in_place=False, force=False, is_sorted=True)
            self.fn = obj.fn
        return obj

    def __getitem__(self, key):
        bt = self.fileobj.tabix_intervals(
            '%s:%s-%s' % (key.chrom, key.start, key.stop))
        for i in bt:
            yield i
        del bt


class BigBedAdapter(BaseAdapter):
    """
    Adapter that provides random access to bigBed files via bx-python
    """
    def __init__(self, fn):
        super(BigBedAdapter, self).__init__(fn)

    def make_fileobj(self):
        return BigBedFile(open(self.fn))

    def __getitem__(self, key):
        chrom = key.chrom
        start = key.start
        stop = key.end
        bx_intervals = self.fileobj.get(chrom, start, stop)
        if bx_intervals is None:
            raise StopIteration
        for i in bx_intervals:
            interval = pybedtools.create_interval_from_list(i.fields)
            interval.file_type = 'bed'
            yield interval


class BigWigAdapter(BaseAdapter):
    """
    Adapter that provides random access to bigWig files bia bx-python
    """
    def __init__(self, fn):
        super(BigWigAdapter, self).__init__(fn)

    def make_fileobj(self):
        return self.fn

    def __getitem__(self, key):
        raise NotImplementedError(
            "__getitem__ not implemented for %s" % self.__class__.__name__)

    def summarize(self, interval, bins=None, method='summarize', function='mean'):

        # We may be dividing by zero in some cases, which raises a warning in
        # NumPy based on the IEEE 754 standard (see
        # http://docs.scipy.org/doc/numpy/reference/generated/
        #       numpy.seterr.html)
        #
        # That's OK -- we're expecting that to happen sometimes. So temporarily
        # disable this error reporting for the duration of this method.
        orig = np.geterr()['invalid']
        np.seterr(invalid='ignore')

        if (bins is None) or (method == 'get_as_array'):
            bw = BigWigFile(open(self.fn))
            s = bw.get_as_array(
                interval.chrom,
                interval.start,
                interval.stop,)
            if s is None:
                s = np.zeros((interval.stop - interval.start,))
            else:
                s[np.isnan(s)] = 0

        elif method == 'ucsc_summarize':
            if function in ['mean', 'min', 'max', 'std', 'coverage']:
                return self.ucsc_summarize(interval, bins, function=function)
            else:
                raise ValueError('function "%s" not supported by UCSC\'s bigWigSummary')

        else:
            bw = BigWigFile(open(self.fn))
            s = bw.summarize(
                interval.chrom,
                interval.start,
                interval.stop, bins)
            if s is None:
                s = np.zeros((bins,))
            else:
                if function == 'sum':
                    s = s.sum_data
                if function == 'mean':
                    s = s.sum_data / s.valid_count
                    s[np.isnan(s)] = 0
                if function == 'min':
                    s = s.min_val
                    s[np.isinf(s)] = 0
                if function == 'max':
                    s = s.max_val
                    s[np.isinf(s)] = 0
                if function == 'std':
                    s = (s.sum_squares / s.valid_count)
                    s[np.isnan(s)] = 0

        # Reset NumPy error reporting
        np.seterr(divide=orig)
        return s

    def ucsc_summarize(self, interval, bins=None, function='mean'):
        if bins is None:
            bins = len(interval)
        y = np.zeros(bins)

        cmds = [
            'bigWigSummary',
            self.fn,
            interval.chrom,
            str(interval.start),
            str(interval.stop),
            str(bins),
            '-type=%s' % function]
        p = subprocess.Popen(
            cmds,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

        def gen():
            try:
                for line in p.stdout:
                    yield line
            finally:
                if p.poll() is None:
                    return
                else:
                    p.wait()
                    err = p.stderr.read().strip()
                    if p.returncode not in (0, None):
                        if err.startswith('no data'):
                            return
                        raise ValueError(
                            "cmds: %s: %s" %
                            (' '.join(cmds), p.stderr.read()))
                    if len(err) != 0:
                        sys.stderr.write(err)

        for line in gen():
            for i, x in enumerate(line.split('\t')):
                try:
                    y[i] = float(x)
                except ValueError:
                    pass
        return np.array(y)
