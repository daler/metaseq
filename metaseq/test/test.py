import metaseq
import numpy as np
from nose.tools import assert_raises
gs = {}
for kind in ['bed', 'bam', 'bigbed', 'bigwig']:
    gs[kind] = metaseq.genomic_signal(metaseq.example_filename('gdc.%s' % kind), kind)


def test_local_count():

    def check_local_count(kind, coord, expected, stranded):
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
            yield check_local_count, kind, coord, expected, stranded




def test_local_coverage_full():
    """generator of tests for local coverage

    ensures that all formats are consistent in their results when retrieving
    the full un-binned data.
    """
    def check_local_coverage(kind, coord, expected):
        result = gs[kind].local_coverage(coord)
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
        ):
            yield check_local_coverage, kind, coord, expected

def test_local_coverage_binned():
    """generator of tests for local coverage

    ensures that all formats are consistent in their results when retrieving
    the full un-binned data.
    """
    def check_local_coverage(kind, coord, expected):
        if kind == 'bigwig':
            result = gs[kind].local_coverage(coord, bins=8, method='get_as_array')
        else:
            result = gs[kind].local_coverage(coord, bins=8)
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
            yield check_local_coverage, kind, coord, expected


def test_nonbigwig_kwargs():
    "exceptions raised when using non-bigwig kwargs with a bigwig file?"
    nonwigs = {
        'read_strand': '+',
        'fragment_size': 200,
        'shift_width': 80,
        'use_score': True,
    }
    for nw in nonwigs.items():
        kwargs = dict((nw,))
        try:
            assert_raises(ValueError, gs['bigwig'].local_coverage, 'chr2L:1-20', **kwargs)
        except AssertionError:
            print kwargs
            raise
