"""
Practical testing grounds to see what sorts of features are needed.  Heavily
commented to serve as interim documentation.

Use the download_data.py script in the test/data dir to get ENCODE CTCF
ChIP-seq data.

Different modes -- TSS, intron, peaks.

Each one generates features of interest, and then grabs the raw data from the
BAM files to construct the signal from that region.  Doing that for thousands
of features results in a matrix.

Cluster the matrix, sort the clusters, adjust heatmap . . . and then add
a strip of dots along the left side that, when zoomed in, can be clicked to
spawn a mini-browser that shows the IP and input signal for that row as well as
all nearby genes.
"""

import os

import numpy as np
from matplotlib import pyplot as plt

import pybedtools
import metaseq
from metaseq.integration import chipseq
from metaseq import colormap_adjust

# Edit the settings files to make tweaks
import atf3_peaks_settings as settings
import atf3_peaks_helpers as helpers

# global list that will store spawned figs, so you can use close_figs() to
# close them all
FIGS = []


def close_figs():
    """
    Convenience function to close all mini-browser figures
    """
    for fig in FIGS:
        plt.close(fig)

# Choices for RUN_TYPE are:
# * 'intron': all introns of all genes on the selected chromosomes
# * 'TSS'   : gene-level TSSs, +/- upstream and downstream bp
# * 'peaks' : peaks from ENCODE; acts as a positive control on the numbers

RUN_TYPE = 'TSS'

try:
    chip = chipseq.Chipseq(
            ip_bam=metaseq.example_filename(
                'wgEncodeHaibTfbsK562Atf3V0416101AlnRep1.bam'
                ),
            control_bam=metaseq.example_filename(
                'wgEncodeHaibTfbsK562RxlchV0416101AlnRep1.bam'
                ),
            dbfn=metaseq.example_filename(
                'Homo_sapiens.GRCh37.66.cleaned.gtf.db')
            )
except ValueError:
    raise ValueError("please use the download_data.py script in the "
                     "data directory")


if RUN_TYPE == "TSS":
    # Gets all genes on selected chroms, then applies the TSS modifier and
    # saves the results
    tss_fn = 'example_tsses.gtf'
    if not os.path.exists(tss_fn):
        features = pybedtools.BedTool(helpers.gene_generator())\
                .filter(helpers.chromfilter)\
                .each(helpers.TSS, upstream=settings.UPSTREAM,
                        downstream=settings.DOWNSTREAM)\
                .saveas(tss_fn)
    else:
        features = pybedtools.BedTool(tss_fn)

elif RUN_TYPE == "intron":
    # Gets all genes and exons on selected chroms, then subtracts exons from
    # genes.
    intron_fn = 'example_introns.gtf'
    if not os.path.exists(intron_fn):
        features = pybedtools.BedTool(helpers.intron_generator())\
                .filter(helpers.chromfilter)\
                .saveas(intron_fn)
    else:
        features = pybedtools.BedTool(intron_fn)

elif RUN_TYPE == 'peaks':
    # Extends the ENCODE peaks and filters out ones with pvals higher than 1e-5
    features = pybedtools.BedTool(peaks)\
            .filter(helpers.chromfilter)\
            .filter(helpers.peak_filter)\
            .each(helpers.peak_extender)\
            .saveas()

# This does most of the work -- given the list of `features`, we send chunks of
# 50 features to each of 8 processes, binning reads in the BAM file in to
# settings.BINS bins and extending the reads 3'-war by settings.FRAGMENT_SIZE
chip.diff_array(
        features=features,
        array_kwargs=dict(processes=8, chunksize=50, bins=settings.BINS,
            fragment_size=settings.FRAGMENT_SIZE),
        )

# Nice colormap centered on zero that doesn't get too saturated on the
# negatives if they're not as extreme as the positives
cmap = colormap_adjust.smart_colormap(
        chip.diffed_array.min(), chip.diffed_array.max())


# Calculate the TIP scores for all features (see Cheng et al. 2001,
# Bioinformatics 27(23):3221-3227)
row_order = np.argsort(metaseq.plotutils.tip_zscores(chip.diffed_array))

# Indices to use if we want mini-batch k-means clustering
# row_order, breaks = metaseq.plotutils.clustered_sortind(d, k=10)

# x-axis for plots, also used for the extent of the matrix
x = np.linspace(-settings.UPSTREAM, settings.DOWNSTREAM, settings.BINS)

# Make and show the 4-panel fig.
chip.plot(x, row_order=row_order, imshow_kwargs=dict(cmap=cmap))
plt.show()
