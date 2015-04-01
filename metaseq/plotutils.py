"""
Module with handy utilities for plotting genomic signal
"""
import itertools
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import mlab
import numpy as np
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import gridspec
import colormap_adjust
from scipy import stats


def ci_plot(x, arr, conf=0.95, ax=None, line_kwargs=None, fill_kwargs=None):
    """
    Plots the mean and 95% ci for the given array on the given axes

    Parameters
    ----------
    x : 1-D array-like
        x values for the plot

    arr : 2-D array-like
        The array to calculate mean and std for

    conf : float [.5 - 1]
        Confidence interval to use

    ax : matplotlib.Axes
        The axes object on which to plot

    line_kwargs : dict
        Additional kwargs passed to Axes.plot

    fill_kwargs : dict
        Additiona kwargs passed to Axes.fill_between
    """
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    line_kwargs = line_kwargs or {}
    fill_kwargs = fill_kwargs or {}

    m, lo, hi = ci(arr, conf)
    ax.plot(x, m, **line_kwargs)
    ax.fill_between(x, lo, hi, **fill_kwargs)
    return ax


def imshow(arr, x=None, ax=None, vmin=None, vmax=None, percentile=True,
           strip=False, features=None, conf=0.95, sort_by=None,
           line_kwargs=None, fill_kwargs=None, imshow_kwargs=None, figsize=(5, 12),
           width_ratios=(4, 1), height_ratios=(4, 1),
           subplot_params=dict(wspace=0.1, hspace=0.1),
           subset_by=None, subset_order=None,):
    """
    Do-it-all function to help with plotting heatmaps

    Parameters
    ----------
    arr : array-like

    x : 1D array
        X values to use.  If None, use range(arr.shape[1])

    ax : matplotlib.Axes
        If not None, then only plot the array on the provided axes.  This will
        ignore any additional arguments provided that apply to figure-level
        configuration or to the average line plot.  For example, `figsize`,
        `width_ratios`, `height_ratios`, `subplot_params`, `line_kwargs`, and
        `fill_kwargs` will all be ignored.

    vmin, vmax : float

    percentile : bool
        If True, then treat values for `vmin` and `vmax` as percentiles rather
        than absolute values.

    strip : bool
        Include a strip plot alongside the array

    features : pybedtools.BedTool or string filename
        Features used to construct the array

    conf : float
        Confidence interval to use in line plot.

    sort_by : array-like
        Use the provided array to sort the array (e.g., an array of expression
        values).  This array will be argsorted to get the proper order.

    line_kwargs, fill_kwargs : dict
        Passed directly to `ci_plot`.

    figsize : tuple
        (Width, height) of the figure to create.

    imshow_kwargs : dict
        Passed directly to matplotlib.pyplot.imshow.  By default, arguments
        used are `origin='lower'`, `aspect="auto"` and a colormap from
        colormap_adjust.smart_colormap generated using the provided `vmin` and
        `vmax`.

    width_ratios, height_ratios: tuple
        These tuples are passed to the `new_shell` function.  The default
        values set up a 2x2 configuration of panels for heatmap, line plot,
        colorbar axes, and optional strip plot.  However modifying
        `width_ratios` or `height_ratios` can be used to create more or fewer panels.

    subplot_params : dict
        Passed to Figure.subplots_adjust

    subset_by : array
        An array of any type (but usually int or str) that contains a class
        label for each row in the heatmap array.  For example, to subset by
        expression, an array the values of "up", "down", or "unchanged" at each
        of the positions could be provided.

        Note that the heatmap array is first sorted by `sort_by` and then split
        into groups according to `subset_by`, so each subset remains sorted by
        `sort_by`.

    subset_order : list-like
        This provides the order in which the subsets are plotted.  Since the
        default imshow arguments contain `origin="lower"`, these will be
        plotted in order starting at the bottom of the heatmap.

    """
    if ax is None:
        fig = new_shell(
            figsize=figsize,
            strip=strip,
            subplot_params=subplot_params,
            width_ratios=width_ratios,
            height_ratios=height_ratios)

    if x is None:
        x = np.arange(arr.shape[1] + 1)

    if percentile:
        if vmin is None:
            vmin = arr.min()
        else:
            vmin = mlab.prctile(arr.ravel(), vmin)
        if vmax is None:
            vmax = arr.max()
        else:
            vmax = mlab.prctile(arr.ravel(), vmax)
    else:
        if vmin is None:
            vmin = arr.min()
        if vmax is None:
            vmax = arr.max()

    cmap = colormap_adjust.smart_colormap(vmin, vmax)
    _imshow_kwargs = dict(origin='lower', cmap=cmap, vmin=vmin, vmax=vmax,
                          aspect='auto')
    if imshow_kwargs is not None:
        _imshow_kwargs.update(imshow_kwargs)

    # previously we did an argsort first; with subsetting we don't want to do
    # that yet....
    #if sort_by is not None:
    #    ind = np.argsort(sort_by)
    #else:
    #    ind = np.arange(arr.shape[0])

    if sort_by is None:
        sort_by = np.arange(arr.shape[0])

    if ax is None:
        array_ax = fig.array_axes
    else:
        array_ax = ax

    # If not provided, assume all in the same subset.
    if subset_by is None:
        subset_by = np.zeros(arr.shape[0])

    # Ensure always array, since we're doing indexing tricks
    if not isinstance(subset_by, np.ndarray):
        subset_by = np.array(subset_by)

    # If not provided, use sorted order
    if subset_order is None:
        subset_order = sorted(np.unique(subset_by))

    inds = []
    for cls in subset_order:
        subset_ind = np.nonzero(subset_by == cls)[0]
        subset_sort_by = sort_by[subset_ind]
        subset_argsort_by = np.argsort(subset_sort_by)
        inds.append(subset_ind[subset_argsort_by])
    ind = np.concatenate(inds)

    mappable = array_ax.imshow(
        arr[ind, :],
        extent=(x.min(), x.max(), 0, arr.shape[0]),
        **_imshow_kwargs
    )

    if line_kwargs is None:
        line_kwargs = {}
    if fill_kwargs is None:
        fill_kwargs = {}

    if isinstance(line_kwargs, dict):
        line_kwargs = [line_kwargs]
    if isinstance(fill_kwargs, dict):
        fill_kwargs = [fill_kwargs]

    _line_kwargs = itertools.cycle(line_kwargs)
    _fill_kwargs = itertools.cycle(fill_kwargs)

    if ax is None:
        plt.colorbar(mappable, fig.cax)
        for subset_ind, label, _lkw, _fkw in zip(inds, subset_order, _line_kwargs, _fill_kwargs):
            ci_plot(
                x,
                arr[subset_ind],
                ax=fig.line_axes,
                line_kwargs=_lkw,
                fill_kwargs=_fkw,
            )
        return fig
    else:
        return ax.figure


def add_labels_to_subsets(ax, subset_by, subset_order, text_kwargs=None,
                          add_hlines=True, hline_kwargs=None):
    """
    Helper function for adding labels to subsets within a heatmap.

    Assumes that imshow() was called with `subsets` and `subset_order`.

    Parameters
    ----------
    ax : matplotlib.Axes
        The axes to label.  Generally you can use `fig.array_axes` attribute of
        the Figure object returned by `metaseq.plotutils.imshow`.

    subset_by, subset_order : array, list
        See `metaseq.plotutils.imshow()` docstring; these should be the same
        `subsets` and `subset_order` that were provided to that function.
    """

    _text_kwargs = dict(transform=ax.get_yaxis_transform())
    if text_kwargs:
        _text_kwargs.update(text_kwargs)

    _hline_kwargs = dict(color='k')
    if hline_kwargs:
        _hline_kwargs.update(hline_kwargs)
    pos = 0
    for label in subset_order:
        ind = subset_by == label
        last_pos = pos
        pos += sum(ind)
        if add_hlines:
            ax.axhline(pos, **_hline_kwargs)
        ax.text(
            1.1,
            last_pos + (pos - last_pos)/2.0,
            label,
            **_text_kwargs)


def calculate_limits(array_dict, method='global', percentiles=None, limit=()):
    """
    Calculate limits for a group of arrays in a flexible manner.

    Returns a dictionary of calculated (vmin, vmax), with the same keys as
    `array_dict`.

    Useful for plotting heatmaps of multiple datasets, and the vmin/vmax values
    of the colormaps need to be matched across all (or a subset) of heatmaps.

    Parameters
    ----------
    array_dict : dict of np.arrays

    method : {'global', 'independent', callable}
        If method="global", then use the global min/max values across all
        arrays in array_dict.  If method="independent", then each array will
        have its own min/max calcuated.  If a callable, then it will be used to
        group the keys of `array_dict`, and each group will have its own
        group-wise min/max calculated.


    percentiles : None or list
        If not None, a list of (lower, upper) percentiles in the range [0,100].
    """
    if percentiles is not None:
        for percentile in percentiles:
            if not 0 <= percentile <= 100:
                raise ValueError("percentile (%s) not between [0, 100]")

    if method == 'global':
        all_arrays = np.concatenate(
            [i.ravel() for i in array_dict.itervalues()]
        )
        if percentiles:
            vmin = mlab.prctile(
                all_arrays, percentiles[0])
            vmax = mlab.prctile(
                all_arrays, percentiles[1])

        else:
            vmin = all_arrays.min()
            vmax = all_arrays.max()
        d = dict([(i, (vmin, vmax)) for i in array_dict.keys()])

    elif method == 'independent':
        d = {}
        for k, v in array_dict.iteritems():
            d[k] = (v.min(), v.max())

    elif hasattr(method, '__call__'):
        d = {}
        sorted_keys = sorted(array_dict.keys(), key=method)
        for group, keys in itertools.groupby(sorted_keys, method):
            keys = list(keys)
            all_arrays = np.concatenate([array_dict[i] for i in keys])
            if percentiles:
                vmin = mlab.prctile(
                    all_arrays, percentiles[0])
                vmax = mlab.prctile(
                    all_arrays, percentiles[1])
            else:
                vmin = all_arrays.min()
                vmax = all_arrays.max()
            for key in keys:
                d[key] = (vmin, vmax)
    return d


def ci(arr, conf=0.95):
    """
    Column-wise confidence interval.

    Parameters
    ----------
    arr : array-like

    conf : float
        Confidence interval

    Returns
    -------
    m : array
        column-wise mean
    lower : array
        lower column-wise confidence bound
    upper : array
        upper column-wise confidence bound
    """
    m = arr.mean(axis=0)
    n = len(arr)
    se = arr.std(axis=0) / np.sqrt(n)
    h = se * stats.t._ppf((1 + conf) / 2., n - 1)
    return m, m - h, m + h


def nice_log(x):
    """
    Uses a log scale but with negative numbers.

    :param x: NumPy array
    """
    neg = x < 0
    xi = np.log2(np.abs(x) + 1)
    xi[neg] = -xi[neg]
    return xi


def tip_zscores(a):
    """
    Calculates the "target identification from profiles" (TIP) zscores
    from Cheng et al. 2001, Bioinformatics 27(23):3221-3227.

    :param a: NumPy array, where each row is the signal for a feature.
    """
    weighted = a * a.mean(axis=0)
    scores = weighted.sum(axis=1)
    zscores = (scores - scores.mean()) / scores.std()
    return zscores


def fdrcorrection(pvals, alpha=0.05, method='indep'):
    '''
    NOTE: This function was copied from
    statsmodels.sandbox.stats.multicomp.fdrcorrection0, from statsmodels
    version 0.5.0.

    This is to avoid requiring all of statsmodels to be a dependency for
    metaseq, just for this function.




    pvalue correction for false discovery rate

    This covers Benjamini/Hochberg for independent or positively correlated and
    Benjamini/Yekutieli for general or negatively correlated tests. Both are
    available in the function multipletests, as method=`fdr_bh`, resp.
    `fdr_by`.

    Parameters
    ----------
    pvals : array_like
        set of p-values of the individual tests.
    alpha : float
        error rate
    method : {'indep', 'negcorr')

    Returns
    -------
    rejected : array, bool
        True if a hypothesis is rejected, False if not
    pvalue-corrected : array
        pvalues adjusted for multiple hypothesis testing to limit FDR

    Notes
    -----

    If there is prior information on the fraction of true hypothesis, then
    alpha should be set to alpha * m/m_0 where m is the number of tests, given
    by the p-values, and m_0 is an estimate of the true hypothesis.  (see
    Benjamini, Krieger and Yekuteli)

    The two-step method of Benjamini, Krieger and Yekutiel that estimates the
    number of false hypotheses will be available (soon).

    Method names can be abbreviated to first letter, 'i' or 'p' for fdr_bh and
    'n' for fdr_by.

    '''
    pvals = np.asarray(pvals)

    pvals_sortind = np.argsort(pvals)
    pvals_sorted = pvals[pvals_sortind]
    sortrevind = pvals_sortind.argsort()

    if method in ['i', 'indep', 'p', 'poscorr']:
        ecdffactor = _ecdf(pvals_sorted)
    elif method in ['n', 'negcorr']:
        cm = np.sum(1./np.arange(1, len(pvals_sorted)+1))  # corrected this
        ecdffactor = _ecdf(pvals_sorted) / cm
#    elif method in ['n', 'negcorr']:
#        cm = np.sum(np.arange(len(pvals)))
#        ecdffactor = ecdf(pvals_sorted)/cm
    else:
        raise ValueError('only indep and necorr implemented')
    reject = pvals_sorted <= ecdffactor*alpha
    if reject.any():
        rejectmax = max(np.nonzero(reject)[0])
        reject[:rejectmax] = True

    pvals_corrected_raw = pvals_sorted / ecdffactor
    pvals_corrected = np.minimum.accumulate(pvals_corrected_raw[::-1])[::-1]
    pvals_corrected[pvals_corrected > 1] = 1
    return reject[sortrevind], pvals_corrected[sortrevind]


def tip_fdr(a, alpha=0.05):
    """
    Returns adjusted TIP p-values for a particular `alpha`.

    (see :func:`tip_zscores` for more info)

    :param a: NumPy array, where each row is the signal for a feature
    :param alpha: False discovery rate

    """
    zscores = tip_zscores(a)
    pvals = stats.norm.pdf(zscores)
    rejected, fdrs = fdrcorrection(pvals)
    return fdrs


def prepare_logged(x, y):
    """
    Transform `x` and `y` to a log scale while dealing with zeros.

    This function scales `x` and `y` such that the points that are zero in one
    array are set to the min of the other array.

    When plotting expression data, frequently one sample will have reads in
    a particular feature but the other sample will not.  Expression data also
    tends to look better on a log scale, but log(0) is undefined and therefore
    cannot be shown on a plot.  This function allows these points to be shown,
    piled up along one side of the plot.

    :param x,y: NumPy arrays
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


def new_shell(figsize=(5, 12), strip=False, height_ratios=(4, 1),
              width_ratios=(4, 1), subplot_params=None):
    if subplot_params is None:
        subplot_params = {}
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(
        len(height_ratios),
        len(width_ratios),
        height_ratios=height_ratios,
        width_ratios=width_ratios,
        **subplot_params)
    fig.array_axes = plt.subplot(gs[0, 0])
    if strip:
        fig.strip_axes = plt.subplot(gs[0, 1], sharey=fig.array_axes)
    fig.line_axes = plt.subplot(gs[1, 0], sharex=fig.array_axes)
    fig.cax = plt.subplot(gs[1, 1])
    fig.gs = gs
    return fig


def matrix_and_line_shell(figsize=(5, 12), strip=False):
    """
    Helper function to construct an empty figure that has space for a matrix,
    a summary line plot directly below it, a colorbar axis, and an optional
    "strip" axis that parallels the matrix (and shares its y-axis) where data
    can be added to create callbacks.

    Returns a tuple of (fig, matrix_ax, line_ax, strip_ax, colorbar_ax) that
    can then be used to plot upon.

    :param figsize: Tuple of (width, height), in inches, for the figure to be
    :param strip: If `strip` is False, then the returned `strip_ax` will be
        None and no strip axes will be created.
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
    Uses MiniBatch k-means clustering to cluster matrix into groups.

    Each cluster of rows is then sorted by `scorefunc` -- by default, the max
    peak height when all rows in a cluster are averaged, or
    cluster.mean(axis=0).max().

    Returns the index that will sort the rows of `x` and a list of "breaks".
    `breaks` is essentially a cumulative row count for each cluster boundary.
    In other words, after plotting the array you can use axhline on each
    "break" to plot the cluster boundary.

    If `k` is a list or tuple, iteratively try each one and select the best
    with the lowest mean distance from cluster centers.

    :param x: Matrix whose rows are to be clustered
    :param k: Number of clusters to create or a list of potential clusters; the
        optimum will be chosen from the list
    :param scorefunc: Optional function for sorting rows within clusters.  Must
        accept a single argument of a NumPy array.
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
            mbk = MiniBatchKMeans(init='k-means++', n_clusters=_k)
            mbk.fit(x)
            mean_dists[_k] = mbk.transform(x).mean()
        best_k = sorted(mean_dists.items(), key=lambda x: x[1])[-1][0]

    mbk = MiniBatchKMeans(init='k-means++', n_clusters=best_k)
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
    for k, g in itertools.groupby(labels[ind]):
        pos += len(list(g))
        breaks.append(pos)

    return ind, breaks


def new_clustered_sortind(x, k=10, row_key=None, cluster_key=None):
    """
    Uses MiniBatch k-means clustering to cluster matrix into groups.

    Each cluster of rows is then sorted by `scorefunc` -- by default, the max
    peak height when all rows in a cluster are averaged, or
    cluster.mean(axis=0).max().

    Returns the index that will sort the rows of `x` and a list of "breaks".
    `breaks` is essentially a cumulative row count for each cluster boundary.
    In other words, after plotting the array you can use axhline on each
    "break" to plot the cluster boundary.

    If `k` is a list or tuple, iteratively try each one and select the best
    with the lowest mean distance from cluster centers.

    :param x: Matrix whose rows are to be clustered
    :param k: Number of clusters to create or a list of potential clusters; the
        optimum will be chosen from the list
    :param row_key:
        Optional function to act as a sort key for sorting rows within
        clusters.  Signature should be `scorefunc(a)` where `a` is a 1-D NumPy
        array.
    :param cluster_key:
        Optional function for sorting clusters.  Signature is `clusterfunc(a)`
        where `a` is a NumPy array containing all rows of `x` for cluster `i`.
        It must return a single value.
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
            mbk = MiniBatchKMeans(init='k-means++', n_clusters=_k)
            mbk.fit(x)
            mean_dists[_k] = mbk.transform(x).mean()
        best_k = sorted(mean_dists.items(), key=lambda x: x[1])[-1][0]

    mbk = MiniBatchKMeans(init='k-means++', n_clusters=best_k)
    mbk.fit(x)
    k = best_k
    labels = mbk.labels_
    scores = np.zeros(labels.shape, dtype=float)

    if cluster_key:
        # It's easier for calling code to provide something that operates on
        # a cluster level, but here it's converted to work on a label level
        # that looks in to the array `x`.
        def _cluster_key(i):
            return cluster_key(x[labels == i, :])
        sorted_labels = sorted(range(k), key=_cluster_key)
    else:
        # Otherwise just use them as-is.
        sorted_labels = range(k)

    if row_key:
        # Again, easier to provide a function to operate on a row.  But here we
        # need it to accept an index
        def _row_key(i):
            return row_key(x[i, :])

    final_ind = []
    breaks = []
    pos = 0
    for label in sorted_labels:
        # which rows in `x` have this label
        label_inds = np.nonzero(labels == label)[0]
        if row_key:
            label_sort_ind = sorted(label_inds, key=_row_key)
        else:
            label_sort_ind = label_inds
        for li in label_sort_ind:
            final_ind.append(li)
        pos += len(label_inds)
        breaks.append(pos)

    return np.array(final_ind), np.array(breaks)


def input_ip_plots(iparr, inputarr, diffed, x, sort_ind,
                   prefix=None, limits1=(None, None), limits2=(None, None),
                   hlines=None, vlines=None):

    """
    All-in-one plotting function to make a 5-panel figure.

    Panels are IP, input, and diffed; plus 2 line plots showing averages.

    :param iparr, inputarr: NumPy arrays constructed by a genomic_signal object
    :param diffed: Difference of `iparr` and `inputarr`, but can be some other
                   transformation.
    :param x: Extent to use -- for TSSs, maybe something like
        np.linspace(-1000, 1000, bins), or for just bin IDs, something like
        `np.arange(bins)`.

    :param sort_ind: row order for each of the 3 panels -- usually interesting
        to use `clustered_sortind` or `tip_zscores`

    :param prefix: Used to prefix plot titles with '%(prefix)s IP", etc
    :param limits1: Tuple passed to the Normalize function for IP and input.
    :param limits2: Tuple passed tot he Normalize function for the diffed array
    :param hlines: List of (position, kwarg) tuples for plotting horizontal
        lines.  Kwargs are passed directly to axhline. Useful for delimiting
        clusters, if you used `clustered_sortind` and have both `row_order` and
        `breaks`.
    :param vlines: List of (position, kwargs) tuples.  A vertical line will be
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
    ax1 = plt.subplot2grid(
        (9, 9), (0, 0), colspan=3, rowspan=6)
    ax2 = plt.subplot2grid(
        (9, 9), (0, 3), colspan=3, rowspan=6, sharex=ax1, sharey=ax1)
    ax3 = plt.subplot2grid(
        (9, 9), (0, 6), colspan=3, rowspan=6, sharex=ax1, sharey=ax1)

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
        limits1[0] = mlab.prctile(
            all_base, 1. / all_base.size)
    if limits1[1] is None:
        limits1[1] = mlab.prctile(
            all_base, 100 - 1. / all_base.size)
    if limits2[0] is None:
        limits2[0] = mlab.prctile(
            diffed.ravel(), 1. / all_base.size)
    if limits2[1] is None:
        limits2[1] = mlab.prctile(
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


def _updatecopy(orig, update_with, keys=None, override=False):
    """
    Update a copy of dest with source.  If `keys` is a list, then only update
    with those keys.
    """
    d = orig.copy()
    if keys is None:
        keys = update_with.keys()
    for k in keys:
        if k in update_with:
            if k in d and not override:
                continue
            d[k] = update_with[k]
    return d


def _clean(z):
    """
    Return a version of z that only has finite values
    """
    return z[np.isfinite(z)]


class MarginalHistScatter(object):
    def __init__(self, ax, hist_size=0.6, pad=0.05):
        """
        Class to enable incremental appending of scatterplots, each of which
        generate additional marginal histograms.
        """
        self.scatter_ax = ax
        self.fig = ax.figure

        self.divider = make_axes_locatable(self.scatter_ax)

        self.top_hists = []
        self.right_hists = []
        self.hist_size = hist_size
        self.pad = pad
        self.xfirst_ax = None
        self.yfirst_ax = None

        # will hold histogram data
        self.hxs = []
        self.hys = []

    @property
    def xmax(self):
        return self.scatter_ax.dataLim.xmax

    @property
    def ymax(self):
        return self.scatter_ax.dataLim.ymax

    @property
    def xmin(self):
        return self.scatter_ax.dataLim.xmin

    @property
    def ymin(self):
        return self.scatter_ax.dataLim.ymin

    @property
    def limits(self):
        return (self.xmin, self.xmax, self.ymin, self.ymax)

    def append(self, x, y, scatter_kwargs, hist_kwargs=None, xhist_kwargs=None,
               yhist_kwargs=None, num_ticks=3, labels=None, hist_share=False,
               marginal_histograms=True):
        """
        Adds a new scatter to self.scatter_ax as well as marginal histograms
        for the same data, borrowing addtional room from the axes.

        Parameters
        ----------

        x, y : array-like
            Data to be plotted

        scatter_kwargs : dict
            Keyword arguments that are passed directly to scatter().

        hist_kwargs : dict
            Keyword arguments that are passed directly to hist(), for both the
            top and side histograms.

        xhist_kwargs, yhist_kwargs : dict
            Additional, margin-specific kwargs for the x or y histograms
            respectively.  These are used to update `hist_kwargs`

        num_ticks : int
            How many tick marks to use in each histogram's y-axis

        labels : array-like
            Optional NumPy array of labels that will be set on the collection
            so that they can be accessed by a callback function.

        hist_share : bool
            If True, then all histograms will share the same frequency axes.
            Useful for showing relative heights if you don't want to use the
            hist_kwarg `normed=True`

        marginal_histograms : bool
            Set to False in order to disable marginal histograms and just use
            as a normal scatterplot.
        """
        scatter_kwargs = scatter_kwargs or {}
        hist_kwargs = hist_kwargs or {}
        xhist_kwargs = xhist_kwargs or {}
        yhist_kwargs = yhist_kwargs or {}
        yhist_kwargs.update(dict(orientation='horizontal'))

        # Plot the scatter
        coll = self.scatter_ax.scatter(x, y, **scatter_kwargs)
        coll.labels = labels

        if not marginal_histograms:
            return

        xhk = _updatecopy(hist_kwargs, xhist_kwargs)
        yhk = _updatecopy(hist_kwargs, yhist_kwargs)

        axhistx = self.divider.append_axes(
            'top', size=self.hist_size,
            pad=self.pad, sharex=self.scatter_ax, sharey=self.xfirst_ax)

        axhisty = self.divider.append_axes(
            'right', size=self.hist_size,
            pad=self.pad, sharey=self.scatter_ax, sharex=self.yfirst_ax)

        axhistx.yaxis.set_major_locator(
            MaxNLocator(nbins=num_ticks, prune='both'))

        axhisty.xaxis.set_major_locator(
            MaxNLocator(nbins=num_ticks, prune='both'))

        if not self.xfirst_ax and hist_share:
            self.xfirst_ax = axhistx

        if not self.yfirst_ax and hist_share:
            self.yfirst_ax = axhisty

        # Keep track of which axes are which, because looking into fig.axes
        # list will get awkward....
        self.top_hists.append(axhistx)
        self.right_hists.append(axhisty)

        # Scatter will deal with NaN, but hist will not.  So clean the data
        # here.
        hx = _clean(x)
        hy = _clean(y)

        self.hxs.append(hx)
        self.hys.append(hy)

        # Only plot hists if there's valid data
        if len(hx) > 0:
            if len(hx) == 1:
                _xhk = _updatecopy(orig=xhk, update_with=dict(bins=[hx[0], hx[0]]), keys=['bins'])
                axhistx.hist(hx, **_xhk)
            else:
                axhistx.hist(hx, **xhk)
        if len(hy) > 0:
            if len(hy) == 1:
                _yhk = _updatecopy(orig=yhk, update_with=dict(bins=[hy[0], hy[0]]), keys=['bins'])
                axhisty.hist(hy, **_yhk)
            else:
                axhisty.hist(hy, **yhk)

        # Turn off unnecessary labels -- for these, use the scatter's axes
        # labels
        for txt in axhisty.get_yticklabels() + axhistx.get_xticklabels():
            txt.set_visible(False)

        for txt in axhisty.get_xticklabels():
            txt.set_rotation(-90)

    def add_legends(self, xhists=True, yhists=False, scatter=True, **kwargs):
        """
        Add legends to axes.
        """
        axs = []
        if xhists:
            axs.extend(self.hxs)
        if yhists:
            axs.extend(self.hys)
        if scatter:
            axs.extend(self.ax)

        for ax in axs:
            ax.legend(**kwargs)
