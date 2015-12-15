"""
The classes in this module enable random access to a variety of file formats
(BAM, bigWig, bigBed, BED) using a uniform syntax, and allow you to compute
coverage across many features in parallel or just a single feature.

Using classes in the :mod:`metaseq.integration` and :mod:`metaseq.minibrowser`
modules, you can connect these objects to matplotlib figures that show a window
into the data, making exploration easy and interactive.

Generally, the :func:`genomic_signal` function is all you need -- just provide
a filename and the format and it will take care of the rest, returning
a genomic signal of the proper type.

Adding support for a new format is straightforward:

    * Write a new adapter for the format in :mod:`metaseq.filetype_adapters`
    * Subclass one of the existing classes below, setting the `adapter`
      attribute to be an instance of this new adapter
    * Add the new class to the `_registry` dictionary to enable support for the
      file format.

Note that to support parallel processing and to avoid repeating code, these
classes delegate their local_coverage methods to the
:func:`metaseq.array_helpers._local_coverage` function.
"""

import os
import sys
import subprocess

import numpy as np
from bx.bbi.bigwig_file import BigWigFile

import pybedtools

from array_helpers import _array, _array_parallel, _local_coverage, \
    _local_count, _count_array, _count_array_parallel
import filetype_adapters
import helpers
from helpers import rebin


def supported_formats():
    """
    Returns list of formats supported by metaseq's genomic signal objects.
    """
    return _registry.keys()


def genomic_signal(fn, kind):
    """
    Factory function that makes the right class for the file format.

    Typically you'll only need this function to create a new genomic signal
    object.

    :param fn: Filename
    :param kind:
        String.  Format of the file; see
        metaseq.genomic_signal._registry.keys()
    """
    try:
        klass = _registry[kind.lower()]
    except KeyError:
        raise ValueError(
            'No support for %s format, choices are %s'
            % (kind, _registry.keys()))
    m = klass(fn)
    m.kind = kind
    return m


class BaseSignal(object):
    """
    Base class to represent objects from which genomic signal can be
    calculated/extracted.

    `__getitem__` uses the underlying adapter the instance was created with
    (e.g., :class:`metaseq.filetype_adapters.BamAdapter` for
    a :class:`BamSignal` object).
    """
    def __init__(self, fn):
        self.fn = fn

    def array(self, features, processes=None, chunksize=1, ragged=False,
              **kwargs):
        """
        Creates an MxN NumPy array of genomic signal for the region defined by
        each feature in `features`, where M=len(features) and N=(bins or
        feature length)

        Parameters
        ----------
        features : iterable of interval-like objects
            An iterable of interval-like objects; see docstring for
            `local_coverage` method for more details.

        processes : int or None
            If not None, then create the array in parallel, giving each process
            chunks of length `chunksize` to work on.

        chunksize : int
            `features` will be split into `chunksize` pieces, and each piece
            will be given to a different process. The optimum value is
            dependent on the size of the features and the underlying data set,
            but `chunksize=100` is a good place to start.

        ragged : bool
            If False (default), then return a 2-D NumPy array.  This requires
            all rows to have the same number of columns, which you get when
            supplying `bins` or if all features are of uniform length.  If
            True, then return a list of 1-D NumPy arrays

        Notes
        -----
        Additional keyword args are passed to local_coverage() which performs
        the work for each feature; see that method for more details.
        """
        if processes is not None:
            arrays = _array_parallel(
                self.fn, self.__class__, features, processes=processes,
                chunksize=chunksize, **kwargs)
        else:
            arrays = _array(self.fn, self.__class__, features, **kwargs)
        if not ragged:
            stacked_arrays = np.row_stack(arrays)
            del arrays
            return stacked_arrays
        else:
            return arrays

    def local_coverage(self, features, *args, **kwargs):
        processes = kwargs.pop('processes', None)
        if not processes:
            return _local_coverage(self.adapter, features, *args, **kwargs)

        if isinstance(features, (list, tuple)):
            raise ValueError(
                "only single features are supported for parallel "
                "local_coverage")

        # we don't want to have self.array do the binning
        bins = kwargs.pop('bins', None)

        # since if we got here processes is not None, then this will trigger
        # a parallel array creation
        features = helpers.tointerval(features)
        x = np.arange(features.start, features.stop)
        features = list(helpers.split_feature(features, processes))
        ys = self.array(
            features, *args, bins=None, processes=processes, ragged=True,
            **kwargs)
        # now we ravel() and re-bin
        y = np.column_stack(ys).ravel()
        if bins:
            xi, yi = rebin(x, y, bins)
            del x, y
            return xi, yi
        return x, y

    local_coverage.__doc__ = _local_coverage.__doc__


class BigWigSignal(BaseSignal):
    def __init__(self, fn):
        """
        Class for operating on bigWig files
        """
        super(BigWigSignal, self).__init__(fn)
        self.adapter = filetype_adapters.BigWigAdapter(fn)


class IntervalSignal(BaseSignal):
    def __init__(self, fn):
        """
        Abstract class for bed, BAM and bigBed files.
        """
        BaseSignal.__init__(self, fn)

    def local_count(self, *args, **kwargs):
        return _local_count(self.adapter, *args, **kwargs)

    local_count.__doc__ = _local_count.__doc__

    def count_array(self, features, processes=None, chunksize=1,  **kwargs):
        if processes is not None:
            arrays = _count_array_parallel(
                self.fn, self.__class__, features, processes=processes,
                chunksize=chunksize, **kwargs)
        else:
            arrays = _count_array(self.fn, self.__class__, features, **kwargs)
        return np.concatenate(arrays)


class BamSignal(IntervalSignal):
    def __init__(self, fn):
        """
        Class for operating on BAM files.
        """
        BaseSignal.__init__(self, fn)
        self._readcount = None
        self.adapter = filetype_adapters.BamAdapter(self.fn)

    def genome(self):
        """
        "genome" dictionary ready for pybedtools, based on the BAM header.
        """
        # This gets the underlying pysam Samfile object
        f = self.adapter.fileobj
        d = {}
        for ref, length in zip(f.references, f.lengths):
            d[ref] = (0, length)
        return d

    def mapped_read_count(self, force=False):
        """
        Counts total reads in a BAM file.

        If a file self.bam + '.scale' exists, then just read the first line of
        that file that doesn't start with a "#".  If such a file doesn't exist,
        then it will be created with the number of reads as the first and only
        line in the file.

        The result is also stored in self._readcount so that the time-consuming
        part only runs once; use force=True to force re-count.

        Parameters
        ----------
        force : bool
            If True, then force a re-count; otherwise use cached data if
            available.
        """
        # Already run?
        if self._readcount and not force:
            return self._readcount

        if os.path.exists(self.fn + '.mmr') and not force:
            for line in open(self.fn + '.mmr'):
                if line.startswith('#'):
                    continue
                self._readcount = float(line.strip())
                return self._readcount

        cmds = ['samtools',
                'view',
                '-c',
                '-F', '0x4',
                self.fn]
        p = subprocess.Popen(
            cmds, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = p.communicate()
        if stderr:
            sys.stderr.write('samtools says: %s' % stderr)
            return None
        mapped_reads = int(stdout)

        # write to file so the next time you need the lib size you can access
        # it quickly
        if not os.path.exists(self.fn + '.mmr'):
            fout = open(self.fn + '.mmr', 'w')
            fout.write(str(mapped_reads) + '\n')
            fout.close()

        self._readcount = mapped_reads
        return self._readcount


class BigBedSignal(IntervalSignal):
    def __init__(self, fn):
        """
        Class for operating on bigBed files.
        """
        IntervalSignal.__init__(self, fn)
        self.adapter = filetype_adapters.BigBedAdapter(fn)


class BedSignal(IntervalSignal):
    def __init__(self, fn):
        """
        Class for operating on BED files.
        """
        IntervalSignal.__init__(self, fn)
        self.adapter = filetype_adapters.BedAdapter(fn)


_registry = {
    'bam': BamSignal,
    'bed': BedSignal,
    'gff': BedSignal,
    'gtf': BedSignal,
    'vcf': BedSignal,
    'bigwig': BigWigSignal,
    'bigbed': BigBedSignal,
}
