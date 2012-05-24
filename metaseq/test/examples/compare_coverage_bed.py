"""
Benchmark script to compare metaseq to coverageBed on a very specific task
(counting reads in windows).

Currently metaseq is ~18x faster (coverageBed = 66s, metaseq = 3.5s), but there
are discrepancies where metaseq is not counting as many reads as bedtools for
some windows.  Across the 5913 windows on chr19,

 * 904 (15%) of windows are off by one
 * 19 are off by 5 reads
 * 20 are off by anywhere from 9 to 41 reads

Diagnostic plots are generated at the end of the script.

TODO: figure out what's causing the discrepancies (open vs closed intervals?
Binning artifact? CIGAR operations?)
"""
import os
import sys
import time
import numpy as np
import metaseq
import pybedtools

from matplotlib import pyplot as plt

bam_fn = metaseq.example_filename('wgEncodeUwTfbsK562CtcfStdAlnRep1.bam')

if not os.path.exists(bam_fn):
    raise ValueError(
            'Please run download_data.py in test/data dir to retrieve ENCODE '
            'data used for examples')


# Construct 10kb windows, but subset to only use chr19 (to speed up the test)
print 'creating windows...'
sys.stdout.flush()
windows = pybedtools.BedTool()\
        .window_maker(genome='hg19', w=10000)\
        .filter(lambda x: x.chrom == 'chr19')\
        .saveas()


def run_bedtools():

    # set up a BAM-based BedTool
    bt = pybedtools.BedTool(bam_fn)

    print 'pybedtools coverageBed starting...',
    sys.stdout.flush()
    t0 = time.time()

    # Need to sort to maintain consistency with metaseq
    bt_result = bt.coverage(windows, counts=True).sort()
    bt_array = np.array([i[-1] for i in bt_result], dtype=int)

    t1 = time.time()
    print 'completed in %.2fs' % (t1 - t0)
    sys.stdout.flush()
    return bt_array


def run_metaseq():
    # set up a BamSignal object
    m = metaseq.genomic_signal(
            metaseq.example_filename('wgEncodeUwTfbsK562CtcfStdAlnRep1.bam'),
            kind='bam')

    print 'metaseq starting...',
    sys.stdout.flush()
    t0 = time.time()

    # Tweak processes and chunksize as needed to balance CPUs and I/O.
    PROCESSES = 6
    CHUNKSIZE = 100

    # the trick is to use a single bin...
    ms_array = m.array(
            windows, processes=PROCESSES, chunksize=CHUNKSIZE, bins=1)

    t1 = time.time()
    print 'completed in %.2fs' % (t1 - t0)
    sys.stdout.flush()
    return ms_array.ravel()


if __name__ == "__main__":
    bt_array = run_bedtools()
    ms_array = run_metaseq()

    # Diagnostic plots
    #
    # Scatter plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(np.log1p(bt_array), np.log1p(ms_array), 'k.', alpha=0.2)
    ax.set_xlabel('log(coverageBed counts + 1)')
    ax.set_ylabel('log(metaseq counts + 1)')
    ax.set_title('Read counts across 10kb-windows in chr19')

    # Histogram of differences
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    diffed = ms_array - bt_array
    diffed_nonzero = diffed[diffed != 0]
    counts, bins, patches = ax2.hist(diffed_nonzero)
    ax2.set_xlabel(
            'Nonzero difference in read counts\n '
            'for each window\n'
            '(metaseq minus bedtools)')
    ax2.set_ylabel('Number of windows')
    fig2.subplots_adjust(bottom=0.2)

    # Table of differences
    print 'bin\tcount'
    for cnt, bn in zip(counts[::-1], bins[::-1]):
        print '%s\t%s' % (bn, cnt)

    plt.show()
