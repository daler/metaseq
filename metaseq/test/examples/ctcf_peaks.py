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
import matplotlib
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from scikits.statsmodels.sandbox.stats.multicomp import fdrcorrection0
import pybedtools
from pybedtools.contrib import plotting
import gffutils
from gffutils.helpers import asinterval
from gffutils.contrib.plotting import Gene
import metaseq
from metaseq import colormap_adjust
import ctcf_peaks_settings as settings
import ctcf_peaks_helpers as helpers

# global list that will store spawned figs, so you can use close_figs() to
# close them all
FIGS = []


def close_figs():
    for fig in FIGS:
        plt.close(fig)

# Choices for RUN_TYPE are:
# * 'intron': all introns of all genes on the selected chromosomes
# * 'TSS'   : gene-level TSSs, +/- upstream and downstream bp
# * 'peaks' : peaks from ENCODE; acts as a positive control on the numbers

RUN_TYPE = 'TSS'

# Genomic signal objects provide convenient random and parallel access to BAM
# files.
try:
    ip = metaseq.genomic_signal(
            metaseq.example_filename('wgEncodeUwTfbsK562CtcfStdAlnRep1.bam'),
            'bam')
    inp = metaseq.genomic_signal(
            metaseq.example_filename('wgEncodeUwTfbsK562InputStdAlnRep1.bam'),
            'bam')
    peaks = metaseq.example_filename(
            'wgEncodeUwTfbsK562CtcfStdPkRep1.narrowPeak.gz')
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

# Create a NumPy array with len(features) rows and `bins` columns, and use
# `processes` cores to do so.
kwargs = dict(features=features, processes=8, chunksize=50, bins=settings.BINS,
        fragment_size=settings.FRAGMENT_SIZE)
ip_arr = ip.array(**kwargs)
inp_arr = inp.array(**kwargs)

# Scale arrays to library size.
# Note that for speed, ip_arr and inp_arr are both integer arrays so they need
# to be converted to float first before scaling.
ip_arr = ip_arr.astype(float) / ip.million_mapped_reads()
inp_arr = inp_arr.astype(float) / inp.million_mapped_reads()


# ip_arr and inp_arr now have units of "reads per million mapped reads", so we
# can subtract input from IP to get the number enriched reads per million
# mapped reads for each bin for each gene
diffed = ip_arr - inp_arr

# Create a nice red/blue colormap
cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
'test', ['#cde3ef', '#FFFFFF', '#b11902'], N=256)

# Convert `diffed` to log scale, but modified in such a way that values that
# were initially negative are still negative.
d = metaseq.plotutils.nice_log(diffed)

# Center red-blue colorbar on zero
vmin = d.min()
vmax = d.max()
cmap = colormap_adjust.cmap_center_point_adjust(cmap, [vmin, vmax], 0)

# Calculate the TIP scores and FDR for all features (see Cheng et al. 2001,
# Bioinformatics 27(23):3221-3227)
fdr = metaseq.plotutils.tip_fdr(d)


if RUN_TYPE == 'TSS':
    # Note that since `d` has retained the sort order of the original features
    # used to create it, we can now iterate through the features and TIP FDR
    # values to determine which have significant regulatory scores
    THRESH = 0.05
    sig_genes = []
    for feature, feature_fdr in zip(features, fdr):
        if feature_fdr < THRESH:
            sig_genes.append(feature['gene_id'])


# TODO: "minibrowser" really ought to be built-in to metaseq . . . use mixins
# for flexibility?
#
# see minibrowser.py for an implementation and details
def minibrowser(feature):
    """
    Given a feature, spawns a new figure showing the local coverage of that
    feature in IP, input, and diffed.

    Look up in the FeatureDB to see what genes are nearby; plot those in
    a second set of axes along with the singal for that entire region.

    Since the goal is to double-check the matrix, we show exactly the data used
    for the feature area, and possibly-differently-binned data for the rest of
    the locus.

    Vertical dotted lines indicate the extent of the original feature.
    """
    fig = plt.figure(figsize=(8, 4))
    FIGS.append(fig)
    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])

    ax = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1], sharex=ax)

    # Create add gene models first so that we know how far out to reach for BAM
    # data.
    extent = [feature.start, feature.stop]
    nearby_genes = settings.G.overlapping_features(
            feature.chrom, feature.start, feature.stop, featuretype='gene')
    ybase = 0
    ngenes = 0
    for nearby_gene in nearby_genes:
        ngenes += 1
        extent.extend([nearby_gene.start, nearby_gene.stop])
        gene_collection = Gene(
                settings.G,
                nearby_gene,
                transcripts=['mRNA'],
                cds=['CDS'],
                utrs=['exon'],
                ybase=ybase,
                color="0.5", picker=5)
        gene_collection.name = nearby_gene.id
        gene_collection.add_to_ax(ax2)
        ybase += gene_collection.max_y

    def gene_callback(event):
        """
        prints the gene's name when it's clicked on
        """
        artist = event.artist
        print artist.name

    fig.canvas.mpl_connect('pick_event', gene_callback)

    xmin = min(extent)
    xmax = max(extent)
    ymax = ngenes

    # 1% padding seems to work well
    padding = (xmax - xmin) * 0.01
    ax2.axis('tight')

    # Make a new feature to represent the region plus surrounding genes
    interval = pybedtools.create_interval_from_list(feature.fields)
    interval.start = xmin - padding
    interval.stop = xmax + padding

    # Important -- otherwise, the entire signal will be flipped and won't be
    # consistent with the un-flipped gene models.
    interval.strand = "."

    # Use the bp-per-bins resolution of the original feature
    bp_per_bin = float(len(feature)) / settings.BINS
    bins = int(len(interval) / bp_per_bin)

    # Get the signal
    x, ipsig = ip.local_coverage(interval, bins=bins, fragment_size=200)
    x,  inputsig = inp.local_coverage(interval, bins=bins, fragment_size=200)

    # Scale
    ipsig = ipsig.astype(float) / ip.million_mapped_reads()
    inputsig = inputsig.astype(float) / inp.million_mapped_reads()

    # Normalize
    diffed = ipsig - inputsig

    # And plot
    ax.plot(x, ipsig, 'r-', label='IP')
    ax.plot(x, inputsig, color='.5', linestyle='--', label='input')
    ax.plot(x, diffed, 'b-', label='enrichment')
    ax.set_ylabel('reads per million mapped reads')
    ax.axis('tight')
    ax.legend(loc='best', prop=dict(size=10), borderpad=0.5, labelspacing=0.2,
            frameon=False)

    # Delimits data from which the original matrix row came from
    ax.axvline(feature.start, color='k', linestyle=":")
    ax.axvline(feature.stop, color='k', linestyle=":")

    ax.xaxis.set_visible(False)
    ax2.yaxis.set_visible(False)
    ax2.set_xlabel('Genomic coords (%s)' % feature.chrom)


def callback(event):
    """
    Callback for the left-hand-side "strip axes" so that clicking on a point
    spawns a mini-browser figure
    """
    artist = event.artist
    ind = artist.ind
    limit = 5
    browser = True
    if len(event.ind) > limit:
        print "More than %s genes selected; not spawning browsers" % limit
        browser = False
    for i in event.ind:
        feature = artist.features[ind[::-1][i]]
        print feature,
        if browser:
            minibrowser(feature)

# Indices to use if we want to sort by TIP fdr
ind = np.argsort(fdr)

# Indices to use if we want mini-batch k-means clustering
ind, breaks = metaseq.plotutils.clustered_sortind(d, k=10)

# Convert features to a list we can access later
feature_list = list(features)

# This function sets up a figure with convenient axes, see the docstring for
# more
axes_info = metaseq.plotutils.matrix_and_line_shell(strip=True)
fig, matrix_ax, line_ax, strip_ax, cbar_ax = axes_info

# x-axis for plots
x = np.linspace(-settings.UPSTREAM, settings.DOWNSTREAM, settings.BINS)

# We'll need this to get the right axes on the heatmap
extent = (-settings.UPSTREAM, settings.DOWNSTREAM, 0, diffed.shape[0])

# Plot the heatmap
mappable = matrix_ax.imshow(
        d[ind],                  # re-order by `ind`
        aspect='auto',           # fill in existing axes
        cmap=cmap,               # the cmap we made
        extent=extent,           # to get the correct axes
        interpolation='nearest'  # prevent blurring when zoomed in
        )

# Colorbar on the provided axes
plt.colorbar(mappable, cbar_ax)

# Plot column (axis=0) averages in the lower panel
line_ax.plot(x, diffed.mean(axis=0))

# Plot dummy points in the "strip" axes; these will be assigned a callback
line, = strip_ax.plot(
        np.zeros((d.shape[0],)),
        np.arange(d.shape[0]) + 0.5,
        color='.5',
        markeredgewidth=0,
        marker='o',
        linestyle='None',
        picker=5)

# Attach info to the artist so that the callback can access it.
line.features = feature_list
line.ind = ind

# Assign the callback
fig.canvas.mpl_connect('pick_event', callback)

# Axes clean-up
matrix_ax.axis('tight')
strip_ax.xaxis.set_visible(False)
matrix_ax.yaxis.set_visible(False)
matrix_ax.xaxis.set_visible(False)

# Uncomment for testing -- zooms in on the top 50 so you don't have to manually
# zoom for testing the callback
#matrix_ax.axis(ymax= 50, ymin=0)

plt.show()
