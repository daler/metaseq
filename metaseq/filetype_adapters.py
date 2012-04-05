"""
File-type adapters accept a filename of the appropriate format (not checked) as
the only argument to their constructor.

Subclasses must define __getitem__ to accept a pybedtools.Interval and return
an iterator of pybedtools.Intervals

Subclasses must define make_fileobj(), which returns an object to be iterated
over in __getitem__
"""
from bx.bbi.bigbed_file import BigBedFile
import pysam
import pybedtools

strand_lookup = {16: '-', 0: '+'}


class BaseAdapter(object):
    def __init__(self, fn):
        self.fn = fn
        self.fileobj = None
        self.fileobj = self.make_fileobj()

    def __getitem__(self, key):
        raise ValueError('Subclasses must define __getitem__')

    def make_fileobj(self):
        raise ValueError('Subclasses must define make_fileobj')


class BamAdapter(BaseAdapter):
    def __init__(self, fn):
        super(BamAdapter, self).__init__(fn)

    def make_fileobj(self):
        return pysam.Samfile(self.fn)

    def __getitem__(self, key):
        iterator = self.fileobj.fetch(
                key.chrom,
                key.start,
                key.stop)
        for r in iterator:
            interval = pybedtools.Interval(
                self.fileobj.references[r.rname],
                r.pos,
                r.pos + r.qend,
                strand=strand_lookup[r.flag & 0x0010])
            interval.file_type = 'bed'
            yield interval


class BedAdapter(BaseAdapter):
    def __init__(self, fn):
        super(BedAdapter, self).__init__(fn)

    def make_fileobj(self):
        
        obj = pybedtools.BedTool(self.fn)
        if not obj._tabixed():
            obj = obj.sort().tabix(in_place=True, force=True)
            self.fn = obj.fn
        return obj

    def __getitem__(self, key):
        bt = self.fileobj.tabix_intervals(
                '%s:%s-%s' % (key.chrom, key.start, key.stop))
        for i in bt:
            yield i
        del bt


class BigBedAdapter(BaseAdapter):
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
            interval = pybedtools.Interval(
                    i.chrom, i.start, i.end, strand=i.strand)
            interval.file_type = 'bed'
            yield interval
