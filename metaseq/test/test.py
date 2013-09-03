"""
Many of these tests use the minimal test/data/gdc.bed file which has just
enough complexity to be useful in testing corner cases.  When reading through
the tests, it's useful to have that file open to understand what's happening.
"""
import metaseq
import multiprocessing
from metaseq.array_helpers import ArgumentError
import numpy as np
from nose.tools import assert_raises
gs = {}
for kind in ['bed', 'bam', 'bigbed', 'bigwig']:
    gs[kind] = metaseq.genomic_signal(metaseq.example_filename('gdc.%s' % kind), kind)


def test_tointerval():
    assert metaseq.helpers.tointerval("chr2L:1-10[-]").strand == '-'
    assert metaseq.helpers.tointerval("chr2L:1-10[+]").strand == '+'
    assert metaseq.helpers.tointerval("chr2L:1-10").strand == '.'


def test_local_count():

    def check(kind, coord, expected, stranded):
        result = gs[kind].local_count(coord, stranded=stranded)
        assert result == expected, (kind, coord, result)

    for kind in ['bam', 'bigbed', 'bed']:
        for coord, expected, stranded in (
            ('chr2L:1-80', 3, False),       #  easy case
            ('chr2L:1000-3000', 0, False),  #  above upper boundary
            ('chr2L:1-9', 0, False),        #  below lower boundary
            ('chr2L:71-73[-]', 2, False),   #  unstranded = 2
            ('chr2L:71-73[-]', 1, True),    #  stranded = 1
            ('chr2L:70-71', 2, False),      #  pathological corner case
            ('chr2L:75-76', 0, False),      #  pathological corner case
        ):
            yield check, kind, coord, expected, stranded


def test_local_coverage_stranded():
    def check(kind, coord, expected):
        result = gs[kind].local_coverage(coord)
        assert np.all(result[0] == expected[0]) and np.all(result[1] == expected[1]), (kind, coord, result)

    for kind in ['bam', 'bigbed', 'bed', 'bigwig']:
        for coord, expected in (
            ('chr2L:1-20[-]',
             (
                 np.array([1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15, 16, 17, 18, 19]),
                 np.array([0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 1., 1., 1., 1., 0., 0., 0., 0., 0.])[::-1],
                 # note reverse------------------------------------------------------------------------^^^^^^
             ),
            ),
            ('chr2L:68-76[-]',
             (
                 np.array([68, 69, 70, 71, 72, 73, 74, 75]),
                 np.array([0., 0., 2., 2., 2., 2., 2., 0.])[::-1],
                 # note reverse----------------------------^^^^^^

             ),
            ),
        ):
            yield check, kind, coord, expected


def test_local_coverage_shifted():
    def check(kind, coord, shift_width, expected):
        result = gs[kind].local_coverage(coord, shift_width=shift_width)
        assert np.all(result[0] == expected[0]) and np.all(result[1] == expected[1]), (kind, coord, result)

    for kind in ['bam', 'bigbed', 'bed']:
        for coord, shift_width, expected in (
            ('chr2L:1-20', -2,
             (
                 np.array([1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15, 16, 17, 18, 19]),
                 np.array([0., 0., 0., 0., 0., 0., 0., 1., 1., 1., 1., 1., 0., 0., 0., 0., 0., 0., 0.]),
             ),
            ),
            # this one is complex, because the minus-strand read shifts left,
            # and the plus-strand shifts right.
            ('chr2L:68-76', 1,
             (
                 np.array([68, 69, 70, 71, 72, 73, 74, 75]),
                 np.array([0., 1., 1., 2., 2., 2., 1., 1.]),
             ),
            ),

            # shift the reads all the way out of the window...
            ('chr2L:68-76', 10,
             (
                 np.array([68, 69, 70, 71, 72, 73, 74, 75]),
                 np.array([0., 0., 0., 0., 0., 0., 0., 0.]),
             ),
            ),
        ):
            yield check, kind, coord, shift_width, expected


def test_local_coverage_read_strand():
    """
    checks stranded full binning

    excludes bigwig since strand doesn't make sense for that format.
    """
    def check(kind, coord, read_strand, expected):
        result = gs[kind].local_coverage(coord, read_strand=read_strand)
        assert np.all(result[0] == expected[0]) and np.all(result[1] == expected[1]), (kind, coord, result)

    for kind in ['bam', 'bigbed', 'bed']:
        for coord, read_strand, expected in (
            ('chr2L:1-20', '+',
             (
                 np.array([1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15, 16, 17, 18, 19]),
                 np.array([0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 1., 1., 1., 1., 0., 0., 0., 0., 0.]),
             ),
            ),
            ('chr2L:1-20', '-',
             (
                 np.array([1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15, 16, 17, 18, 19]),
                 np.array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]),
             ),
            ),
        ):
            yield check, kind, coord, read_strand, expected


def test_local_coverage_fragment_size():
    def check(kind, coord, fragment_size, expected):
        result = gs[kind].local_coverage(coord, fragment_size=fragment_size)
        assert np.all(result[0] == expected[0]) and np.all(result[1] == expected[1]), (kind, coord, result)

    for kind in ['bam', 'bigbed', 'bed']:
        for coord, fragment_size, expected in (
            ('chr2L:1-20', 7,
             (
                 np.array([1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15, 16, 17, 18, 19]),
                 np.array([0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 1., 1., 1., 1., 1., 1., 0., 0., 0.]),
             ),
            ),

            ('chr2L:68-76', 6,
             (
                 np.array([68, 69, 70, 71, 72, 73, 74, 75]),
                 np.array([0., 1., 2., 2., 2., 2., 2., 1.]),
             ),
            ),

            ('chr2L:68-76', 1,
             (
                 np.array([68, 69, 70, 71, 72, 73, 74, 75]),
                 np.array([0., 0., 1., 0., 0., 0., 1., 0.]),
             ),
            ),
        ):
            yield check, kind, coord, fragment_size, expected


def test_local_coverage_score():
    def check(kind, coord, expected):
        result = gs[kind].local_coverage(coord, use_score=True)
        assert np.all(result[0] == expected[0]) and np.all(result[1] == expected[1]), (kind, coord, result)

    for kind in ['bigbed', 'bed']:
        for coord, expected in (
            ('chr2L:1-20',
             (
                 np.array([1,  2,  3,  4,  5,  6,  7,  8,  9,  10,   11,   12,   13,   14,   15, 16, 17, 18, 19]),
                 np.array([0., 0., 0., 0., 0., 0., 0., 0., 0., 255., 255., 255., 255., 255., 0., 0., 0., 0., 0.]),
             ),
            ),
            ('chr2L:68-76',
             (
                 np.array([68, 69, 70,   71,   72,   73,   74,   75]),
                 np.array([0., 0., 510., 510., 510., 510., 510., 0.]),
             ),
            ),
        ):
            yield check, kind, coord, expected


def test_local_coverage_full():
    """generator of tests for local coverage

    ensures that all formats are consistent in their results when retrieving
    the full un-binned data.
    """
    def check(kind, coord, processes, expected):
        result = gs[kind].local_coverage(coord, processes=processes)
        assert np.all(result[0] == expected[0]) and np.all(result[1] == expected[1]), (kind, coord, result)

    for kind in ['bam', 'bigbed', 'bed', 'bigwig']:
        for coord, expected in (
            ('chr2L:1-20',
             (
                 np.array([1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15, 16, 17, 18, 19]),
                 np.array([0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 1., 1., 1., 1., 0., 0., 0., 0., 0.]),
             ),
            ),
            ('chr2L:68-76',
             (
                 np.array([68, 69, 70, 71, 72, 73, 74, 75]),
                 np.array([0., 0., 2., 2., 2., 2., 2., 0.]),
             ),
            ),
            ('chr2L:568-576',
             (
                 np.array([568, 569, 570, 571, 572, 573, 574, 575]),
                 np.array([0., 0.,   0.,  0.,  0.,  0.,  0.,  0.]),
             ),
            ),
        ):
            for processes in [None, multiprocessing.cpu_count()]:

                yield check, kind, coord, processes, expected


def test_local_coverage_binned():
    """generator of tests for local coverage

    ensures that all formats are consistent in their results when retrieving
    binned data.
    """
    def check(kind, coord, processes, expected):
        if kind == 'bigwig':
            result = gs[kind].local_coverage(coord, bins=8, method='get_as_array', processes=processes)
        else:
            result = gs[kind].local_coverage(coord, bins=8, processes=processes)
        try:
            assert np.allclose(result[0], expected[0]) and np.allclose(result[1], expected[1])
        except:
            print (kind, coord, result, expected)
            raise

    for kind in ['bam', 'bigbed', 'bed', 'bigwig']:
        for coord, expected in (
            ('chr2L:1-20',
             (
                 np.array([ 1., 3.57142857, 6.14285714, 8.71428571, 11.28571429, 13.85714286, 16.42857143, 19.]),
                 np.array([ 0., 0.,         0.,         0.,         1.,          1.,          0.,          0. ]),
             ),
            ),
            ('chr2L:68-76',
             (
                 np.array([68, 69, 70, 71, 72, 73, 74, 75]),
                 np.array([0., 0., 2., 2., 2., 2., 2., 0.]),
             ),
            ),
        ):
            for processes in [None, multiprocessing.cpu_count()]:
                yield check, kind, coord, processes, expected


def test_array_binned():
    def check(kind, coord, processes, expected):
        if kind == 'bigwig':
            result = gs[kind].array(coord, bins=8, method='get_as_array', processes=processes)
        else:
            result = gs[kind].array(coord, bins=8, processes=processes)
        try:
            assert np.allclose(result, expected)
        except:
            print (kind, coord, result, expected)
            raise

    for kind in ['bam', 'bigbed', 'bed', 'bigwig']:
        for coord, expected in (
            (['chr2L:1-20'],
             np.array([[0., 0., 0., 0., 1., 1., 0., 0. ]]),
            ),
            (['chr2L:1-20', 'chr2L:1-20[-]'],
             np.array([[0., 0., 0., 0., 1., 1., 0., 0. ],
                       [0., 0., 1., 1., 0., 0., 0., 0. ]]),
            ),
            (['chr2L:68-76'],
             np.array([[0., 0., 2., 2., 2., 2., 2., 0.]]),
             ),
        ):
            for processes in [None, multiprocessing.cpu_count()]:
                yield check, kind, coord, processes, expected


def test_array_ragged():
    def check(kind, coord, processes, expected):
        if kind == 'bigwig':
            result = gs[kind].array(coord, method='get_as_array', ragged=True, processes=processes)
        else:
            result = gs[kind].array(coord, processes=processes, ragged=True)
        try:
            if isinstance(result, list):
                for i, j in zip(result, expected):
                    assert np.allclose(i, j)
        except:
            print (kind, coord, result, expected)
            raise

    for kind in ['bam', 'bigbed', 'bed', 'bigwig']:
        for coord, expected in (
            (['chr2L:1-20'],
             np.array([[0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  1.,  1.,  1.,  1.,   1.,  0.,  0.,  0.,  0.,  0.]])
            ),
            (
                [['chr2L:1-20'], ['chr2L:1-19']],
                [
                    np.array([ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  1.,  1.,  1.,  1.,   1.,  0.,  0.,  0.,  0.,  0.]),
                    np.array([ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  1.,  1.,  1.,  1.,   1.,  0.,  0.,  0.,  0., ])
                ]
            ),

            ([['chr2L:68-76', 'chr2L:68-76'], 'chr2L:68-76'],
             [
                 np.array([0., 0., 2., 2., 2., 2., 2., 0., 0., 0., 2., 2., 2., 2., 2., 0.]),
                 np.array([0., 0., 2., 2., 2., 2., 2., 0.])
             ],
            ),
        ):
            for processes in [None, multiprocessing.cpu_count()]:
                yield check, kind, coord, processes, expected


def test_bigwig_methods():
    """
    using bx-python's `summarize` gives different results than UCSC's summarize
    and summarizing from the original BAM.

    possible bx-python bug, or just a resolution thing?  Haven't worked that
    out yet, so test results here.
    """
    # y0 is identical to the above test, and should be the same as
    # local_coverage for bam, bigbed, and bed as well.
    #
    # y1 may not always be identical . . . I think it's due to resolution of
    # the summarized file, especially when working with such short genomic
    # intervals?
    x0 = np.array([1., 3.57142857, 6.14285714, 8.71428571, 11.28571429, 13.85714286, 16.42857143, 19.])
    y0 = np.array([ 0., 0.,         0.,         0.,         1.,          1.,          0.,          0. ]),
    y1 = np.array([ 0., 0.,         0.,         0.,         1.,          1.,          1.,          0. ]),
    #                                                         note difference here ---^
    x, y = gs['bigwig'].local_coverage('chr2L:1-20', bins=8, method='get_as_array')
    assert np.allclose(x0, x)
    assert np.allclose(y0, y)

    x, y = gs['bigwig'].local_coverage('chr2L:1-20', bins=8, method='summarize')
    assert np.allclose(x0, x)
    assert np.allclose(y1, y)

    x, y = gs['bigwig'].local_coverage('chr2L:1-20', bins=8, method='ucsc_summarize')
    assert np.allclose(x0, x)
    assert np.allclose(y0, y)


def test_nonbigwig_kwargs():
    """
    test non-bigwig kwargs with a bigwig file

    these kwargs are not valid for bigwigs because they don't make sense for
    the format.
    """
    nonwigs = {
        'read_strand': '+',
        'fragment_size': 200,
        'shift_width': 80,
        'use_score': True,
    }
    for nw in nonwigs.items():
        kwargs = dict((nw,))
        try:
            assert_raises(
                ArgumentError, gs['bigwig'].local_coverage, 'chr2L:1-20',
                **kwargs)
        except AssertionError:
            print kwargs
            raise

    assert_raises(ArgumentError, gs['bigwig'].local_coverage, 'chr2L:1-20',
                  method='ucsc_summarize', preserve_total=True)


def test_local_coverage_kwarg_errors():
    """
    ensures the argument-checking is working
    """
    def check(kind, kwargs):
        assert_raises(
            ArgumentError,
            gs[kind].local_coverage, **kwargs)

    items = [
        ('bam', dict(bins='a', features='chr2L:1-20')),
        ('bam', dict(use_score=True, features='chr2L:1-20')),
        ('bam', dict(bins=[8, 10], features=['chr2L:1-20'])),

    ]

    for kind, kwargs in items:
        yield check, kind, kwargs


def test_local_coverage_full_multifeature():
    def check(kind, coord, expected):
        result = gs[kind].local_coverage(coord)
        assert np.all(result[0] == expected[0]) and np.all(result[1] == expected[1]), (kind, coord, result)

    for kind in ['bam', 'bigbed', 'bed', 'bigwig']:
        for coord, expected in (
            (['chr2L:1-20', 'chr2L:68-76'],
             (
                 np.array([1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 68, 69, 70, 71, 72, 73, 74, 75]),
                 np.array([0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 1., 1., 1., 1., 0., 0., 0., 0., 0., 0., 0., 2., 2., 2., 2., 2., 0.]),
             ),
            ),
        ):
            yield check, kind, coord, expected


def test_errors():
    "these things should raise errors"
    def check(error, callable_obj, args, kwargs):
        assert_raises(error, callable_obj, *args, **kwargs)

    class X(metaseq.filetype_adapters.BaseAdapter):
        def make_fileobj(self):
            return None

    items = [
        (ValueError, metaseq.filetype_adapters.BaseAdapter, (metaseq.example_filename('gdc.bed'),), {}),
        (NotImplementedError, metaseq.filetype_adapters.BigWigAdapter(metaseq.example_filename('gdc.bigwig')).__getitem__, (0,), {}),
        (ValueError, X("").__getitem__, (0,), {}),
        (ValueError, gs['bam'].local_coverage, ['chr2L:1-5', 'chr2L:1-5'], dict(processes=multiprocessing.cpu_count())),
    ]
    for error, callable_obj, args, kwargs in items:
        yield check, error, callable_obj, args, kwargs

def test_supported_formats():
    assert set(metaseq._genomic_signal.supported_formats()) \
        == set(['bam', 'bigwig', 'bed', 'gff', 'gtf', 'vcf', 'bigbed'])
    assert_raises(ValueError, metaseq.genomic_signal, "", 'unsupported')

def test_bam_genome():
    assert gs['bam'].genome() == {'chr2L': (0, 23011544L)}

def test_bam_mmr():
    import os
    fn = metaseq.example_filename('gdc.bam.scale')
    if os.path.exists(fn):
        os.unlink(fn)
    assert gs['bam'].million_mapped_reads(force=True) == 8e-6
    gs['bam']._readcount = None
    assert gs['bam'].million_mapped_reads(force=False) == 8e-6

