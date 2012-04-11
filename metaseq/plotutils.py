from itertools import groupby
import matplotlib
from matplotlib import pyplot as plt
import numpy as np
from scipy import stats

def prepare_logged(x, y):
    """
    scales x and y so that they have zero values at the global min
    """
    xi = np.log2(x)
    yi = np.log2(y)

    xv = np.isfinite(xi)
    yv = np.isfinite(yi)

    global_min = min(xi[xv].min(), yi[yv].min())
    global_max = max(xi[xv].max(), yi[yv].max())

    xi[~xv] = global_min
    yi[~yv] = global_min

    return xi, yi


def matrix_and_line_shell(figsize=(5, 12), strip=False):
    """
    Helper function to construct an empty figure that has space for a matrix,
    a summary line plot directly below it, a colorbar axis, and an optional
    "strip" axis that parallels the matrix (and shares its y-axis) where data
    can be added to create callbacks.

    Returns a tuple of (fig, matrix_ax, line_ax, strip_ax, colorbar_ax).

    If `strip` is False, then the returned `strip_ax` will be None.
    """
    fig = plt.figure(figsize=figsize)

    # Constants to keep track
    if strip:
        STRIP_COLS = 1
    else:
        STRIP_COLS = 0
    ROWS = 4
    COLS = 8 + STRIP_COLS
    MAT_COLS = 7
    MAT_ROWS = 3
    LINE_ROWS = ROWS - MAT_ROWS

    mat_ax = plt.subplot2grid(
            shape=(ROWS, COLS),
            loc=(0, STRIP_COLS),
            rowspan=MAT_ROWS,
            colspan=MAT_COLS,
            )

    line_ax = plt.subplot2grid(
            shape=(ROWS, COLS),
            loc=(MAT_ROWS, STRIP_COLS),
            rowspan=LINE_ROWS,
            colspan=MAT_COLS,
            sharex=mat_ax)

    if strip:
        strip_ax = plt.subplot2grid(
                shape=(ROWS, COLS),
                loc=(0, 0),
                rowspan=MAT_ROWS,
                colspan=STRIP_COLS,
                sharey=mat_ax,
                )
    else:
        strip_ax = None

    cax = plt.subplot2grid(
            shape=(ROWS, COLS),
            loc=(ROWS - MAT_ROWS, MAT_COLS + STRIP_COLS),
            rowspan=1,
            colspan=1,
            )

    fig.subplots_adjust(hspace=0.1, wspace=0.2, right=0.88, left=0.23)
    return fig, mat_ax, line_ax, strip_ax, cax


def clustered_sortind(x, k=10, scorefunc=None):
    """
    Uses MiniBatch k-means clustering to cluster matrix `x` into `k` groups of
    rows.  The clusters are then sorted by `scorefunc` -- by default, the max
    peak height when all rows in a cluster are averaged, or
    cluster.mean(axis=0).max().

    Returns the index that will sort the rows of `x` and a list of "breaks".
    `breaks` is essentially a cumulative row count for each cluster boundary.
    In other words, after plotting the array you can use axhline on each
    "break" to plot the cluster boundary.

    If `k` is a list or tuple, iteratively try each one and select the best
    with the lowest mean distance from cluster centers.
    """
    try:
        from sklearn.cluster import MiniBatchKMeans
    except ImportError:
        raise ImportError('please install scikits.learn for '
                          'clustering.')

    # If integer, do it once and we're done
    if isinstance(k, int):
        best_k = k

    else:
        mean_dists = {}
        for _k in k:
            mbk = MiniBatchKMeans(init='k-means++', k=_k, chunk_size=45)
            mbk.fit(x)
            mean_dists[_k] = mbk.transform(x).mean()
        best_k = sorted(mean_dists.items(), key=lambda x: x[1])[-1][0]

    mbk = MiniBatchKMeans(init='k-means++', k=best_k, chunk_size=45)
    mbk.fit(x)
    k = best_k
    labels = mbk.labels_
    scores = np.zeros(labels.shape, dtype=float)

    if not scorefunc:
        def scorefunc(x):
            return x.mean(axis=0).max()

    for label in range(k):
        ind = labels == label
        score = scorefunc(x[ind, :])
        scores[ind] = score

    pos = 0
    breaks = []
    ind = np.argsort(scores)
    for k, g in groupby(labels[ind]):
        pos += len(list(g))
        breaks.append(pos)

    return ind, breaks


def input_ip_plots(iparr, inputarr, diffed, x, sort_ind,
                   prefix=None, limits1=(None, None), limits2=(None, None),
                   hlines=None, vlines=None):

    """
    All-in-one plotting function to make a 5-panel (IP, input, and diffed;
    2 line plots showing averages) figure out of pre-computed data.

    `iparr` and `inputarr` are arrays constructed from MetaChip.array().

    `diffed` is the difference of the two, possibly re-scaled.

    `x` is the extent to use -- for TSSs, maybe something like
    np.linspace(-1000, 1000, bins), or just bin IDs like np.arange(bins).

    `sort_ind` is the row order for each of the 3 panels -- usually interesting
    to use clustered_sortind.

    `prefix` is used to prefix plot titles with '%(prefix)s IP", etc

    `limits1` is a tuple passed to the Normalize function for IP and input.

    `limits2` is a tuple passed tot he Normalize function for the diffed array

    `hlines` is a list of (position, kwarg) tuples for plotting horizontal
    lines.  Kwargs are passed directly to axhline. Useful for delimiting
    clusters, if you used clustered_sortind.

    `vlines` is a list of (position, kwargs) tuples.  A vertical line will be
    plotted at each position using kwargs.
    """

    # global min and max
    gmin = min(iparr.min(), inputarr.min())
    gmax = max(iparr.max(), inputarr.max())

    fig = plt.figure(figsize=(10, 10))

    # 3 arrays, 2 line plots, a gene strip, and 2 colorbars.  Plots share the
    # axes that make sense
    #
    # 3 arrays
    ax1 = plt.subplot2grid((9, 9), (0, 0),
            colspan=3, rowspan=6)
    ax2 = plt.subplot2grid((9, 9), (0, 3),
            colspan=3, rowspan=6, sharex=ax1, sharey=ax1)
    ax3 = plt.subplot2grid((9, 9), (0, 6),
            colspan=3, rowspan=6, sharex=ax1, sharey=ax1)

    # 2 line plots
    ax4 = plt.subplot2grid((9, 9), (6, 3), colspan=3, rowspan=3, sharex=ax1)
    ax5 = plt.subplot2grid((9, 9), (6, 6), colspan=3, rowspan=3, sharex=ax1)

    # 2 colorbars
    cax1 = plt.Axes(fig, rect=(0.05, 0.25, 0.25, 0.025))
    cax2 = plt.Axes(fig, rect=(0.05, 0.15, 0.25, 0.025))

    # For nice imshow axes
    extent = (min(x), max(x), 0, diffed.shape[0])

    cm = matplotlib.cm.gist_gray
    cm.set_bad('k')
    cm.set_over('r')
    cm.set_under('b')

    limits1 = list(limits1)
    limits2 = list(limits2)

    all_base = np.column_stack((iparr.ravel(), inputarr.ravel())).ravel()

    if limits1[0] is None:
        limits1[0] = stats.scoreatpercentile(
                all_base, 1. / all_base.size)
    if limits1[1] is None:
        limits1[1] = stats.scoreatpercentile(
                all_base, 100 - 1. / all_base.size)
    if limits2[0] is None:
        limits2[0] = stats.scoreatpercentile(
                diffed.ravel(), 1. / all_base.size)
    if limits2[1] is None:
        limits2[1] = stats.scoreatpercentile(
                diffed.ravel(), 100 - 1. / all_base.size)

    del all_base

    imshow_kwargs = dict(
                        interpolation='nearest',
                        aspect='auto',
                        cmap=cm,
                        norm=matplotlib.colors.Normalize(*limits1),
                        extent=extent,
                        origin='lower')

    # modify kwargs for diffed (by changing the normalization)
    diffed_kwargs = imshow_kwargs.copy()
    diffed_kwargs['norm'] = matplotlib.colors.Normalize(*limits2)

    # IP
    mappable1 = ax1.imshow(iparr[sort_ind, :], **imshow_kwargs)

    # input
    mappable2 = ax2.imshow(inputarr[sort_ind, :], **imshow_kwargs)

    # diffed
    mappable3 = ax3.imshow((diffed)[sort_ind, :], **diffed_kwargs)

    # IP and input line plot with vertical line
    ax4.plot(x, inputarr.mean(axis=0), color='k', linestyle='--',
            label='input')
    ax4.plot(x, iparr.mean(axis=0), color='k', label='ip')
    ax4.axvline(0, color='k', linestyle=':')

    # Diffed line plot with vertical line
    ax5.plot(x, diffed.mean(axis=0), 'k', label='enrichment')
    ax5.axvline(0, color='k', linestyle=':')

    # Colorbars
    cbar1 = fig.colorbar(mappable1, cax1, orientation='horizontal')
    cbar2 = fig.colorbar(mappable3, cax2, orientation='horizontal')
    fig.add_axes(cax1)
    fig.add_axes(cax2)

    # labeling...
    ax1.set_ylabel('features')
    plt.setp(ax2.get_yticklabels(), visible=False)
    plt.setp(ax3.get_yticklabels(), visible=False)
    ax4.set_xlabel('bp')
    ax4.set_ylabel('mean reads per million mapped reads')
    ax5.set_xlabel('bp')
    cax1.set_xlabel('Reads per million mapped reads')
    cax2.set_xlabel('Enrichment (RPMMR)')

    if prefix is None:
        prefix = ""
    ax1.set_title('%s IP' % prefix)
    ax2.set_title('%s input' % prefix)
    ax3.set_title('Difference')

    # diffed line plot should have y ax on right
    ax5.yaxis.set_ticks_position('right')
    ax5.yaxis.set_label_position('right')
    ax5.set_ylabel('enriched reads per million mapped reads')

    # Legends
    ax4.legend(loc='best', frameon=False)
    ax5.legend(loc='best', frameon=False)

    # Make sure everybody snaps to xmin/xmax
    for ax in [ax1, ax2, ax3, ax4, ax5]:
        ax.axis(xmin=extent[0], xmax=extent[1])

    if not hlines:
        hlines = []
    if not vlines:
        vlines = []

    for ax in [ax1, ax2, ax3]:
        for pos, kwargs in hlines:
            ax.axhline(pos, **kwargs)
        for pos, kwargs in vlines:
            ax.axvline(pos, **kwargs)

    fig.subplots_adjust(bottom=0.05, top=0.95, hspace=0.75, wspace=0.9)

    return fig
