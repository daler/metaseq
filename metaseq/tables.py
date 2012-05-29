"""
Module for working with tables of data, typically with at least one column
containing gene accessions.  This module makes heavy use of `pandas`.

In order to maximize compatibility with pandas.DataFrame and pandas.Panel
classes, this module uses a functional approach rather than a class-based
approach (see also the discussion cautioning subclassing pandas.DataFrames at
https://groups.google.com/forum/#!topic/pystatsmodels/hIg9oEpHpWI)

"""
import sys
import pandas
from pandas import read_table
import numpy as np
from matplotlib import pyplot as plt

import metaseq
import gffutils
from gffutils.helpers import asinterval
import pybedtools


def skiprows(fn):
    """
    Figure out the number of comment-lines to skip at the start of file `fn`.
    """
    for i, line in enumerate(open(fn)):
        if line.startswith('#'):
            continue
        else:
            break
    return i


def clean_unknown(df, db, index_col='id'):
    """
    Removes rows for genes not found in the database (e.g., "ambiguous"
    fields from DESeq results)
    """
    ind = []
    for i, gene_id in enumerate(df[index_col]):
        try:
            db[gene_id]
            ind.append(i)
        except gffutils.FeatureNotFoundError:
            continue
    return df.ix[ind, :]


def enriched(df, pval=0.05, sig_column='padj', en_column='log2FoldChange',
        ind=True, exclude_nans=False):
    """
    Identify enriched genes.

    :param pval: Alpha to use
    :param sig_column: Which column to use for significance
    :param en_column: Which column to use for enrichment
    :param ind: Boolean; if True will return just the index but if False,
        will return a copy
    """
    idx = (df[sig_column] < pval) & (df[en_column] > 0)

    if exclude_nans:
        idx = idx & df[en_column].notnull()

    if ind:
        return idx
    else:
        return df[idx]


def scatter(x, y, xfunc=None, yfunc=None, xscale=None,
        yscale=None, xlab=None, ylab=None, genes_to_highlight=None,
        label_genes=False, general_kwargs=dict(color="k", alpha=0.2,
            linewidths=0), marginal=True, offset_kwargs={}, label_kwargs=None,
        ax=None, one_to_one=None, callback=None, xlab_prefix=None, ylab_prefix=None):
    """
    Do-it-all method for making annotated scatterplots.

    :param x, y:
        Variables to plot -- say, "df.baseMeanA" and "df.baseMeanB"

    :param xfunc, yfunc:
        Functions to apply to `xvar` and `yvar` respectively. Default is
        log2; set to None to have no transformation.

    :param xlab, ylab:
        Labels for x and y axes; default is to use function names for
        `xfunc` and `yfunc` and variable names `xvar` and `yvar`, e.g.,
        "log2(baseMeanA)"

    :param ax:
        If `ax=None`, then makes a new fig and returns the Axes object,
        otherwise, plots onto `ax`

    :param general_kwargs:
        Kwargs for matplotlib.scatter; specifies how all points look

    :param genes_to_highlight:
        Provides lots of control to colors.  It is a list of (`ind`,
        `kwargs`) tuples, where each `ind` specifies genes to plot with
        `kwargs`.  Each dictionary updates a copy of `general_kwargs`. If
        `genes_to_highlight` has a "name" kwarg, this must be a list that't
        the same length as `ind`.  It will be used to label the genes in
        `ind` using `label_kwargs`.

    :param marginal:
        Boolean, toggles display of non-finite cases where `x` or `y` is
        nonzero but the other one is zero (and `xfunc` or `yfunc` are log)

    :param callback:
        Function to call upon clicking a point. Default is to print the
        gene name, but an example of another useful callback would be
        a mini-browser connected to a genomic_signal object from which the
        expression data were calculated.

    :param one_to_one:
        If not None, a dictionary of matplotlib.plot kwargs that will be
        used to plot a 1:1 line.

    :param label_kwargs:
        Kwargs for labeled genes (e.g., dict=(style='italic')).  Will only
        be used if an entry in `genes_to_highlight` has a `name` key.

    :param offset_kwargs:
        Kwargs to be passed to matplotlib.transforms.offset_copy, used for
        adjusting the positioning of gene labels in relation to the actual
        point.

    :param xlab_prefix, ylab_prefix:
        Optional label prefix that will be added to the beginning of `xlab`
        and/or `ylab`.
    """
    # Construct defaults---------------------------------------------------
    def identity(x):
        return x.copy()

    if xlab_prefix is None:
        xlab_prefix = ""

    if ylab_prefix is None:
        ylab_prefix = ""

    if xlab is None:
        try:
            xname = x.name
        except AttributeError:
            xname = 'x'

        if xfunc is not None:
            xlab = xlab_prefix + "%s(%s)" % (xfunc.__name__, xname)
        else:
            xlab = xlab_prefix + "%s" % xname

    if ylab is None:
        try:
            yname = y.name
        except AttributeError:
            yname = 'y'
        if yfunc is not None:
            ylab = ylab_prefix + "%s(%s)" % (yfunc.__name__, yname)
        else:
            ylab = ylab_prefix + "%s" % yname

    if xfunc is None:
        xfunc = identity

    if yfunc is None:
        yfunc = identity

    if general_kwargs is None:
        general_kwargs = {}

    if genes_to_highlight is None:
        genes_to_highlight = []

    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    if label_kwargs is None:
        label_kwargs = dict(horizontalalignment='right',
                verticalalignment='center', style='italic',
                bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))

    # ---------------------------------------------------------------------

    xi = xfunc(x)
    yi = yfunc(y)

    xv = np.isfinite(xi.astype(float))
    yv = np.isfinite(yi.astype(float))

    global_min = min(xi[xv].min(), yi[yv].min())
    global_max = max(xi[xv].max(), yi[yv].max())

    if marginal:
        xi[~xv] = global_min
        yi[~yv] = global_min

    # Plot everybody
    ax.scatter(xi, yi, picker=5, **general_kwargs)

    # one-to-one line, if kwargs were specified
    if one_to_one:
        ax.plot([global_min, global_max],
                [global_min, global_max],
                **one_to_one)

    # plot any specially-highlighted genes, and label if specified
    for ind, kwargs in genes_to_highlight:
        names = kwargs.pop('names', None)
        updated_kwargs = general_kwargs.copy()
        updated_kwargs.update(kwargs)
        ax.scatter(xi[ind], yi[ind], **updated_kwargs)

        if names:
            transOffset = matplotlib.transforms.offset_copy(ax.transData,
                    fig=ax.figure, **offset_kwargs)

            for xii, yii, name in zip(xi[ind], yi[ind], names):
                ax.text(xii,
                        yii,
                        name,
                        transform=transOffset,
                        **label_kwargs)

    # register callback
    if callback is not None:
        ax.figure.canvas.mpl_connect('pick_event', callback)

    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)

    return ax


def dataframe_features(df, db):
    """
    Returns a pybedtools.BedTool object representing the indexed genes.

    This function does not check to ensure a gene exists in the database;
    please use :func:`clean_unknown` on `df` first.

    :param df: A pandas.DataFrame object, indexed on gene IDs
    :param db: A gffutils.FeatureDB object
    """
    def generator():
        for gene_id in df.index:
            yield asinterval(db[gene_id])

    return pybedtools.BedTool(generator())


def deseq_dataframe(fn, cols=None, db=None, index_col='id'):
    """
    Returns a "cleaned" pandas.DataFrame from DESeq results.

    :param fn: Filename
    :param cols: If not None, an iterable of columns to include
    :param db: If not None, a gffutils.FeatureDB object which will be used
        to strip from the results any "genes" not included in the database.
        This will remove things like the DESeq "ambiguous" and "no_feature"
        rows.
    :param index_col: Which column to use as the index column (default is
        "id")

    """
    data = read_table(fn, sep='\t', skiprows=skiprows(fn))
    if cols:
        data = data.ix[:, cols]
    if db:
        data = clean_unknown(data, db, index_col=index_col)
    data = data.set_index(index_col)
    return data


def panel_splom(p, val='baseMeanA', hist_kwargs=None, **kwargs):
    """
    Scatterplot matrix of datasets in a panel.

    :param p: A pandas.Panel object
    :param val: Which value to plot for each item in the Panel
    :param hist_kwargs: If not None, a histogram will be plotted along the
        diagonal.
    """
    if hist_kwargs is None:
        hist_kwargs = {}

    try:
        xfunc = kwargs['xfunc']
    except KeyError:
        xfunc = lambda x: x

    try:
        yfunc = kwargs['yfunc']
    except KeyError:
        yfunc = lambda y: y

    fig = plt.figure()
    nrows = len(p.items)
    ncols = len(p.items)
    axind = 1
    for i in p.items:
        for j in p.items:
            ax = fig.add_subplot(nrows, ncols, axind)
            if (i == j) and hist_kwargs:
                ax.hist(xfunc(p.ix[i][val]), **hist_kwargs)
            else:
                scatter(p.ix[i][val], p.ix[j][val], ax=ax, xlab_prefix=i + " ", ylab_prefix=j + " ", **kwargs)
            axind += 1


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import numpy as np

    dbfn = metaseq.example_filename('Homo_sapiens.GRCh37.66.cleaned.gtf.db')
    db = gffutils.FeatureDB(dbfn)

    p = pandas.Panel(
            {
                'uninduced_1': deseq_dataframe(
                    metaseq.example_filename('GSM847565_SL2585.table'),
                    index_col='id', db=db),

                'induced_1': deseq_dataframe(
                    metaseq.example_filename('GSM847566_SL2592.table'),
                    index_col='id', db=db),
                'uninduced_2': deseq_dataframe(
                    metaseq.example_filename('GSM847567_SL4337.table'),
                    index_col='id', db=db),

                'induced_2': deseq_dataframe(
                    metaseq.example_filename('GSM847568_SL4326.table'),
                    index_col='id', db=db),
                })

    from pybedtools.featurefuncs import TSS

    features = dataframe_features(p['uninduced_1'], db)
    strong_peaks = pybedtools.BedTool(
            metaseq.example_filename('wgEncodeHaibTfbsK562Atf3V0416101PkRep1.broadPeak.gz'))\
                    .filter(lambda x: int(x.score) > 800)

    features_with_peaks = features\
            .each(TSS, upstream=1, downstream=1)\
            .intersect(strong_peaks, c=True)

    peak_ind = np.array([int(i[-1]) > 0 for i in features_with_peaks])


    def callback(event):
        ind = event.ind
        for i in event.ind:
            print p.ix[:, i]
            print db[p.major_axis[i]].attributes
            sys.stdout.flush()

    scatter(
            x=p['uninduced_1'].fpkm,
            y=p['induced_1'].fpkm,
            xfunc=np.log2,
            yfunc=np.log2,
            callback=callback,
            genes_to_highlight=[
                (peak_ind, dict(color='r', alpha=0.5)),
                ]
            )


    """
    panel_splom(p,
            val='fpkm',
            xfunc=np.log2, yfunc=np.log2, one_to_one=dict(color='r', linestyle='--'),
            callback=callback,
            genes_to_highlight=[
                (peak_ind, dict(color='r')),
                ]
            )
    """

    plt.show()
