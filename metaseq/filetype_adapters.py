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
from bx.intervals.io import StrandFormatError
import numpy as np
import subprocess
import pysam
import pybedtools
import os
import sys
from textwrap import dedent

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
            str(key.chrom),
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
        return pybedtools.IntervalFile(self.fn)

    def __getitem__(self, key):
        for i in self.fileobj.all_hits(key):
            yield i


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
        try:
            bx_intervals = self.fileobj.get(chrom, start, stop)
        except StrandFormatError:
            raise NotImplementedError(dedent(
                """
                It appears you have a version of bx-python where bigBed files
                are temporarily unsupported due to recent changes in the
                bx-python dependency. In the meantime, please convert bigBed to
                BAM like this:

                    bigBedToBed {0} tmp.bed
                    bedtools bedtobam -i tmp.bed > {0}.bam

                and create a genomic signal object using this {0}.bam file.
                """.format(self.fn)))
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

    def summarize(self, interval, bins=None, method='summarize',
                  function='mean', zero_inf=True, zero_nan=True):
        """
        Parameters
        ----------

        interval : object
            Object with chrom (str), start (int) and stop (int) attributes.

        bins : int or None
            Number of bins; if None, bins will be the length of the interval

        method : summarize | ucsc_summarize | get_as_array
            "summarize" and "get_as_array" use bx-python; "ucsc_summarize" uses
            bigWigSummarize. See other notes in docstring for
            metaseq.array_helpers._local_coverage. If None, defaults to
            "summarize".

        function : mean | min | max | std | coverage
            Determines the nature of the summarized values. Ignored if
            `method="get_as_array"`; "coverage" is only valid if method is
            "ucsc_summarize".

        zero_inf, zero_nan : bool
            If `zero_inf` is True, set any inf or -inf to zero before
            returning. If `zero_nan` is True, set any nan values to zero before
            returning.
        """

        if method is None:
            method = 'summarize'

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
                if zero_nan:
                    s[np.isnan(s)] = 0
                if zero_inf:
                    s[np.isinf(s)] = 0

        elif method == 'ucsc_summarize':
            if function in ['mean', 'min', 'max', 'std', 'coverage']:
                return self.ucsc_summarize(interval, bins, function=function)
            else:
                raise ValueError('function "%s" not supported by UCSC\'s'
                                 'bigWigSummary')

        elif method == 'summarize':
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
                elif function == 'mean':
                    s = s.sum_data / s.valid_count
                    if zero_nan:
                        s[np.isnan(s)] = 0
                elif function == 'min':
                    s = s.min_val
                    if zero_inf:
                        s[np.isinf(s)] = 0
                elif function == 'max':
                    s = s.max_val
                    if zero_inf:
                        s[np.isinf(s)] = 0
                elif function == 'std':
                    s = (s.sum_squares / s.valid_count)
                    if zero_nan:
                        s[np.isnan(s)] = 0
                else:
                    raise ValueError(
                            'function "%s" not supported by bx-python'
                            % function
                    )
        else:
            raise ValueError("method '%s' not in [summarize, ucsc_summarize, get_as_array]" % method)

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
