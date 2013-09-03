"""
low-level tests
"""
import os
import multiprocessing
import logging


import sys
import datetime
from pprint import pprint as pp
import numpy as np
import pybedtools
import metaseq
from nose.tools import assert_raises
PROCESSES = 8
def test_example_data_exists():
    assert os.path.exists(metaseq.example_filename('x.bam'))
    assert os.path.exists(metaseq.example_filename('gdc.bam'))


class ResultsTable_features_test(object):
    def setup(self):
        deseq_fn = metaseq.example_filename('ex.deseq')
        db_fn = metaseq.example_filename('dmel-all-r5.33-cleaned.gff.db')
        self.d = metaseq.ResultsTable(deseq_fn, db_fn)

    def test_integer_indexing(self):
        # don't go out of range
        d = self.d[0]
        assert_raises(IndexError, d.__getitem__, 1)
        assert len(self.d[0:5]) == 5

    def test_mutability_of_indexed(self):
        # indexes into self should be copies
        d = self.d[0]
        assert self.d.id[0] == d.id[0]
        d.id[0] = 'asdf'
        print d.id
        assert d.id[0] == 'asdf'
        assert d.id[0] != self.d.id[0]

    def test_enrichment_inds(self):
        # inds should be as long as self
        assert len(self.d.enriched()) \
                == len(self.d.disenriched()) \
                == len(self.d) \
                == len(self.d.data)

    def test_enrichment(self):
        d = self.d[:]
        assert not np.all(d.padj < 0.05)
        d.padj[0] = 0.05
        ind = d.enriched()
        enriched = d[ind]

        print enriched.padj
        assert np.all(enriched.padj < 0.05)
        #assert sum(enriched.padj == 0.05) == 1
        assert not np.any(enriched.padj > 0.05)

        assert np.all(enriched.log2foldchange > 0)
        enriched = self.d[self.d.enriched(pval=0.1)]
        assert np.all(enriched.padj < 0.1)

    def test_gene_ind(self):
        assert_raises(KeyError, self.d.gene_data, 'asdf')


class GenomicSignal_local_count_test(object):
    # BAM file looks like this:
    """
    None	0	chr2L	11	255	5M	*	0	0	CGACA	IIIII	NM:i:0	NH:i:1
    None	16	chr2L	71	255	5M	*	0	0	TTCTC	IIIII	NM:i:0	NH:i:1
    None	0	chr2L	71	255	5M	*	0	0	GAGAA	IIIII	NM:i:0	NH:i:1
    None	16	chr2L	141	255	5M	*	0	0	CACCA	IIIII	NM:i:0	NH:i:1
    None	0	chr2L	141	255	5M	*	0	0	TGGTG	IIIII	NM:i:0	NH:i:1
    None	16	chr2L	151	255	5M	*	0	0	GTTCA	IIIII	NM:i:0	NH:i:1
    None	0	chr2L	161	255	5M	*	0	0	GATAA	IIIII	NM:i:0	NH:i:1
    None	0	chr2L	211	255	5M	*	0	0	AAATA	IIIII	NM:i:0	NH:i:1
    """
    def setup(self):
        self.m = metaseq.genomic_signal(
                metaseq.example_filename('gdc.bam'), kind='bam')
        line = '[%s] %s\n' % (datetime.datetime.now(), self.__class__.__name__)
        print line
        sys.stdout.flush()

    def teardown(self):
        pass

    def test_count(self):
        interval = pybedtools.create_interval_from_list(
                ['chr2L', '0', '100', '.', '.', '+'])
        assert self.m.local_count(interval) == 3
        assert self.m.local_count(interval, stranded=True) == 2

class GenomicSignal_local_coverage_test(object):
    # BAM file looks like this:
    """
    None	0	chr2L	11	255	5M	*	0	0	CGACA	IIIII	NM:i:0	NH:i:1
    None	16	chr2L	71	255	5M	*	0	0	TTCTC	IIIII	NM:i:0	NH:i:1
    None	0	chr2L	71	255	5M	*	0	0	GAGAA	IIIII	NM:i:0	NH:i:1
    None	16	chr2L	141	255	5M	*	0	0	CACCA	IIIII	NM:i:0	NH:i:1
    None	0	chr2L	141	255	5M	*	0	0	TGGTG	IIIII	NM:i:0	NH:i:1
    None	16	chr2L	151	255	5M	*	0	0	GTTCA	IIIII	NM:i:0	NH:i:1
    None	0	chr2L	161	255	5M	*	0	0	GATAA	IIIII	NM:i:0	NH:i:1
    None	0	chr2L	211	255	5M	*	0	0	AAATA	IIIII	NM:i:0	NH:i:1
    """

    # ------------------------------------------------------------------------------
    # Note that SAM files like this are 1-based, but metaseq, like BAM, is 0-based.
    # ------------------------------------------------------------------------------

    def setup(self):
        self.m = metaseq.genomic_signal(
                metaseq.example_filename('gdc.bam'), kind='bam')
        line = '[%s] %s\n' % (datetime.datetime.now(), self.__class__.__name__)
        print line
        sys.stdout.flush()
        pass

    def teardown(self):
        pass

    def test_bintotal_invariance(self):
        # read count should always be the same, no matter how you bin things.
        feature = pybedtools.create_interval_from_list(['chr2L', '0', '20', '.', '.', '+'])
        x0, y0 = self.m.local_coverage(feature, fragment_size=10, bins=40)
        x1, y1 = self.m.local_coverage(feature, fragment_size=10, bins=20)
        x2, y2 = self.m.local_coverage(feature, fragment_size=10, bins=10)
        x3, y3 = self.m.local_coverage(feature, fragment_size=10, bins=5)
        x4, y4 = self.m.local_coverage(feature, fragment_size=10, bins=7)
        assert y0.sum() == y1.sum() == y2.sum() == y3.sum() == y4.sum()

    def test_identity_bins(self):
        # Same number of bins as bp in the feature
        feature = pybedtools.create_interval_from_list(['chr2L', '0', '20', '.', '.', '+'])
        x, y = self.m.local_coverage(feature, fragment_size=10, bins=20)
        pp(zip(x, y))
        assert np.allclose(
                np.array(zip(x, y)),
                np.array(
            [(0, 0),
             (1, 0),
             (2, 0),
             (3, 0),
             (4, 0),
             (5, 0),
             (6, 0),
             (7, 0),
             (8, 0),
             (9, 0),
             (10, 1),
             (11, 1),
             (12, 1),
             (13, 1),
             (14, 1),
             (15, 1),
             (16, 1),
             (17, 1),
             (18, 1),
             (19, 1)]))

    def test_doubled_bins(self):
        #Feature is 20 bp, but ask for 40 bins.  Should space out bins with zeros
        feature = pybedtools.create_interval_from_list(['chr2L', '0', '20', '.', '.', '+'])
        x, y = self.m.local_coverage(feature, fragment_size=10, bins=40)
        pp(zip(x, y))
        assert np.allclose(
                np.array(zip(x, y)),
                np.array(

                [(0.0, 0),
                 (0.5, 0),
                 (1.0, 0),
                 (1.5, 0),
                 (2.0, 0),
                 (2.5, 0),
                 (3.0, 0),
                 (3.5, 0),
                 (4.0, 0),
                 (4.5, 0),
                 (5.0, 0),
                 (5.5, 0),
                 (6.0, 0),
                 (6.5, 0),
                 (7.0, 0),
                 (7.5, 0),
                 (8.0, 0),
                 (8.5, 0),
                 (9.0, 0),
                 (9.5, 0),
                 (10.0, 1),
                 (10.5, 0),
                 (11.0, 1),
                 (11.5, 0),
                 (12.0, 1),
                 (12.5, 0),
                 (13.0, 1),
                 (13.5, 0),
                 (14.0, 1),
                 (14.5, 0),
                 (15.0, 1),
                 (15.5, 0),
                 (16.0, 1),
                 (16.5, 0),
                 (17.0, 1),
                 (17.5, 0),
                 (18.0, 1),
                 (18.5, 0),
                 (19.0, 1),
                 (19.5, 0)]),
                )

    def test_halfed_bins(self):
        feature = pybedtools.create_interval_from_list(['chr2L', '0', '20', '.', '.', '+'])
        x, y = self.m.local_coverage(feature, fragment_size=10, bins=10)
        pp(zip(x, y))
        assert np.allclose(
                np.array(zip(x, y)),
                np.array(
                [(0, 0),
                 (2, 0),
                 (4, 0),
                 (6, 0),
                 (8, 0),
                 (10, 2),
                 (12, 2),
                 (14, 2),
                 (16, 2),
                 (18, 2)]),
                )

    def test_strange_bins(self):
        feature = pybedtools.create_interval_from_list(['chr2L', '0', '20', '.', '.', '+'])
        x, y = self.m.local_coverage(feature, fragment_size=10, bins=7)
        pp(zip(x, y))
        assert np.allclose(
                np.array(zip(x, y)),
                np.array(
                [(0.0, 0),
                 (2.8571428571428572, 0),
                 (5.7142857142857144, 0),
                 (8.5714285714285712, 2),
                 (11.428571428571429, 3),
                 (14.285714285714286, 3),
                 (17.142857142857142, 2)]))

    def test_shift_of_0(self):
        feature = pybedtools.create_interval_from_list(['chr2L', '60', '80', '.', '.', '+'])
        x, y = self.m.local_coverage(feature, fragment_size=5, shift_width=0)
        pp(zip(x,y))
        assert zip(x, y) == \
            [(60, 0),
             (61, 0),
             (62, 0),
             (63, 0),
             (64, 0),
             (65, 0),
             (66, 0),
             (67, 0),
             (68, 0),
             (69, 0),
             (70, 2),
             (71, 2),
             (72, 2),
             (73, 2),
             (74, 2),
             (75, 0),
             (76, 0),
             (77, 0),
             (78, 0),
             (79, 0)]

    def test_shiftwith_of_1(self):
        #
        # Reads on opposite strands shift oppositely...
        #
        #      |||||   original
        #
        #     |||||    minus strand read, leftshift 1
        #       |||||  plus strand read, rightshift 1
        #     1122211
        #
        feature = pybedtools.create_interval_from_list(['chr2L', '60', '80', '.', '.', '+'])
        x, y = self.m.local_coverage(feature, fragment_size=5, shift_width=1)
        pp(zip(x, y))
        assert zip(x, y) == \
            [(60, 0),
             (61, 0),
             (62, 0),
             (63, 0),
             (64, 0),
             (65, 0),
             (66, 0),
             (67, 0),
             (68, 0),
             (69, 1),
             (70, 1),
             (71, 2),
             (72, 2),
             (73, 2),
             (74, 1),
             (75, 1),
             (76, 0),
             (77, 0),
             (78, 0),
             (79, 0)]

    def test_shiftwidth_and_fragmentsize(self):
        #
        # Reads on opposite strands shift oppositely...
        #
        #      |||||   original
        #
        #    +|||||    minus strand read, leftshift 1 and additional to left (3')
        #       |||||+ plus strand read, rightshift 1 and additional to right (3')
        #    111222111
        feature = pybedtools.create_interval_from_list(['chr2L', '60', '80', '.', '.', '+'])
        x, y = self.m.local_coverage(feature, fragment_size=6, shift_width=1)
        pp(zip(x, y))
        assert zip(x, y) == \
            [(60, 0),
             (61, 0),
             (62, 0),
             (63, 0),
             (64, 0),
             (65, 0),
             (66, 0),
             (67, 0),
             (68, 1),
             (69, 1),
             (70, 1),
             (71, 2),
             (72, 2),
             (73, 2),
             (74, 1),
             (75, 1),
             (76, 1),
             (77, 0),
             (78, 0),
             (79, 0)]

    def test_shiftwidth_of_1_plus_only(self):
        # The plus-strand read should shift right by 1
        feature = pybedtools.create_interval_from_list(['chr2L', '0', '20', '.', '.', '+'])
        x, y = self.m.local_coverage(feature, fragment_size=5, shift_width=1)
        pp(zip(x, y))
        assert zip(x, y) == \
            [(0, 0),
             (1, 0),
             (2, 0),
             (3, 0),
             (4, 0),
             (5, 0),
             (6, 0),
             (7, 0),
             (8, 0),
             (9, 0),
             (10, 0),
             (11, 1),
             (12, 1),
             (13, 1),
             (14, 1),
             (15, 1),
             (16, 0),
             (17, 0),
             (18, 0),
             (19, 0)]

    def test_plus_feature(self):
        # first make one where fragments are exactly as long as reads

        feature = pybedtools.create_interval_from_list(['chr2L', '0', '20', '.', '.', '+'])
        x, y = self.m.local_coverage(feature, fragment_size=5)
        pp(zip(x,y))
        assert zip(x, y) == \
                [(0, 0),
                 (1, 0),
                 (2, 0),
                 (3, 0),
                 (4, 0),
                 (5, 0),
                 (6, 0),
                 (7, 0),
                 (8, 0),
                 (9, 0),
                 (10, 1),
                 (11, 1),
                 (12, 1),
                 (13, 1),
                 (14, 1),
                 (15, 0),
                 (16, 0),
                 (17, 0),
                 (18, 0),
                 (19, 0)]

    def test_minus_feature(self):
        # minus strand flips the profile
        feature = pybedtools.create_interval_from_list(['chr2L', '0', '20', '.', '.', '-'])
        xm, ym = self.m.local_coverage(feature, fragment_size=5)
        pp(zip(xm,ym))
        assert zip(xm, ym) == \
                [(0, 0),
                 (1, 0),
                 (2, 0),
                 (3, 0),
                 (4, 0),
                 (5, 1),
                 (6, 1),
                 (7, 1),
                 (8, 1),
                 (9, 1),
                 (10, 0),
                 (11, 0),
                 (12, 0),
                 (13, 0),
                 (14, 0),
                 (15, 0),
                 (16, 0),
                 (17, 0),
                 (18, 0),
                 (19, 0)]

    def test_minus_feature_is_reversed_profile(self):
        feature1 = pybedtools.create_interval_from_list(['chr2L', '0', '20', '.', '.', '+'])
        feature2 = pybedtools.create_interval_from_list(['chr2L', '0', '20', '.', '.', '-'])
        x, y = self.m.local_coverage(feature1, fragment_size=5)
        xm, ym = self.m.local_coverage(feature2, fragment_size=5)
        pp(list(enumerate(zip(ym, y))))
        assert list(ym) == list(y[::-1])
        assert list(xm) == list(x)

    def test_plus_feature_minus_reads(self):
        # the plus-strand read from 10-14 should not appear -- so all zeros
        feature = pybedtools.create_interval_from_list(['chr2L', '0', '20', '.', '.', '+'])
        x, y = self.m.local_coverage(feature, read_strand='-', fragment_size=5)
        pp(zip(x, y))
        assert zip(x, y) == \
            [(0, 0),
             (1, 0),
             (2, 0),
             (3, 0),
             (4, 0),
             (5, 0),
             (6, 0),
             (7, 0),
             (8, 0),
             (9, 0),
             (10, 0),
             (11, 0),
             (12, 0),
             (13, 0),
             (14, 0),
             (15, 0),
             (16, 0),
             (17, 0),
             (18, 0),
             (19, 0)]

    def test_truncated(self):
        # what happens if we stop halfway through a read?
        feature = pybedtools.create_interval_from_list(['chr2L', '0', '12', '.', '.', '+'])
        x, y = self.m.local_coverage(feature, fragment_size=5)
        pp(zip(x, y))
        assert zip(x, y) == \
            [(0, 0),
             (1, 0),
             (2, 0),
             (3, 0),
             (4, 0),
             (5, 0),
             (6, 0),
             (7, 0),
             (8, 0),
             (9, 0),
             (10, 1),
             (11, 1)]

    def test_fragmentsize(self):
        feature = pybedtools.create_interval_from_list(['chr2L', '0', '25', '.', '.', '+'])
        x, y = self.m.local_coverage(feature, fragment_size=10)
        pp(zip(x, y))
        assert zip(x, y) == \
            [(0, 0),
             (1, 0),
             (2, 0),
             (3, 0),
             (4, 0),
             (5, 0),
             (6, 0),
             (7, 0),
             (8, 0),
             (9, 0),
             (10, 1),
             (11, 1),
             (12, 1),
             (13, 1),
             (14, 1),
             (15, 1),
             (16, 1),
             (17, 1),
             (18, 1),
             (19, 1),
             (20, 0),
             (21, 0),
             (22, 0),
             (23, 0),
             (24, 0)]


class GenomicSignal_array_test(object):
    # BAM file looks like this:
    """
    None	0	chr2L	11	255	5M	*	0	0	CGACA	IIIII	NM:i:0	NH:i:1
    None	16	chr2L	71	255	5M	*	0	0	TTCTC	IIIII	NM:i:0	NH:i:1
    None	0	chr2L	71	255	5M	*	0	0	GAGAA	IIIII	NM:i:0	NH:i:1
    None	16	chr2L	141	255	5M	*	0	0	CACCA	IIIII	NM:i:0	NH:i:1
    None	0	chr2L	141	255	5M	*	0	0	TGGTG	IIIII	NM:i:0	NH:i:1
    None	16	chr2L	151	255	5M	*	0	0	GTTCA	IIIII	NM:i:0	NH:i:1
    None	0	chr2L	161	255	5M	*	0	0	GATAA	IIIII	NM:i:0	NH:i:1
    None	0	chr2L	211	255	5M	*	0	0	AAATA	IIIII	NM:i:0	NH:i:1
    """

    # ------------------------------------------------------------------------------
    # Note that SAM files like this are 1-based, but metaseq, like BAM, is 0-based.
    # ------------------------------------------------------------------------------

    def setup(self):
        self.m = metaseq.genomic_signal(
                metaseq.example_filename('gdc.bam'), kind='bam')
        pass

    def teardown(self):
        pass

    def simple_array_test(self):
        features = [
                (0, 20),
                (60, 80),
                (200, 220)]
        features = [pybedtools.Interval('chr2L', *i) for i in features]
        arr = self.m.array(features, bins=20, fragment_size=5)
        assert np.all(arr == np.array(
                [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0],
                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0],
                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0]]))

    def shift_1_simple_array_test(self):
        #Shifting the feature to the right by one should shift the values to the left by one.
        features = [
                (0, 20),
                (61, 81),
                (200, 220)]
        features = [pybedtools.Interval('chr2L', *i) for i in features]
        arr = self.m.array(features, bins=20, fragment_size=5)
        print arr
        assert np.all(arr == np.array(
                [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0],
                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0],
                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0]]))


    def parallel_array_test(self):
        features = [
                (0, 20),
                (61, 81),
                (200, 220)]
        features = [pybedtools.Interval('chr2L', *i) for i in features]
        arr0 = self.m.array(features, bins=20, fragment_size=5)
        assert np.all(arr0 == np.array(
                [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0],
                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0],
                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0]]))

        arr1 = self.m.array(
                features, bins=20, fragment_size=5, chunksize=1, processes=PROCESSES)
        assert np.all(arr0 == arr1)

        # mix up the chunksize and processes
        arr2 = self.m.array(
                features, bins=20, fragment_size=5, chunksize=3, processes=PROCESSES)
        assert np.all(arr0 == arr2)

        # use more features and test for identity again
        features *= 1000
        print len(features)
        arr0 = self.m.array(features, bins=20, fragment_size=5)
        arr1 = self.m.array(
                features, bins=20, fragment_size=5, chunksize=5, processes=PROCESSES)
        print arr0.shape
        print arr1.shape
        print (arr0 != arr1)
        assert np.all(arr0 == arr1)


class GenomicSignal_bigBed_local_coverage_test(GenomicSignal_local_coverage_test):
    """
    chr2L   10      15      None    255     +
    chr2L   70      75      None    255     -
    chr2L   70      75      None    255     +
    chr2L   140     145     None    255     -
    chr2L   140     145     None    255     +
    chr2L   150     155     None    255     -
    chr2L   160     165     None    255     +
    chr2L   210     215     None    255     +
    """
    def setup(self):
        self.m = metaseq.genomic_signal(
                metaseq.example_filename('gdc.bigbed'), kind='bigbed')


class GenomicSignal_bigBed_array_test(GenomicSignal_array_test):
    def setup(self):
        self.m = metaseq.genomic_signal(
                metaseq.example_filename('gdc.bigbed'), kind='bigbed')

