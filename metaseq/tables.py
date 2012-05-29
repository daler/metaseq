"""
Module for working with tables of data, typically with at least one column
containing gene accessions.

In order to maximize compatibility with pandas.DataFrame and pandas.Panel
classes, this module uses a functional approach rather than a class-based
approach.
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
    Figure out the number of comment-lines to skip in file `fn`.
    """
    for i, line in enumerate(open(fn)):
        if line.startswith('#'):
            continue
        else:
            break
    return i


def dataframe(fn, cols=None):
    """
    Reads in a tab-delimited file.

    Lines at the beginning of the file starting with a "#" are skipped.

    :param fn: Filename to read
    :param cols: Optional list of columns to include.  Any columns not in this
        list will not be included.
    """
    data = read_table(fn, sep='\t', skiprows=skiprows(fn))
    if cols:
        data = data.ix[:, cols]
    data.index = data['id']
    return data


def monkeypatch_class(name, bases, namespace):
    """
    From GvR,
    http://mail.python.org/pipermail/python-dev/2008-January/076194.html
    """
    assert len(bases) == 1, "Exactly one base class required"
    base = bases[0]
    for name, value in namespace.iteritems():
        if name != "__metaclass__":
            setattr(base, name, value)
    return base


class DESeqResults(pandas.DataFrame):
    __metaclass__ = monkeypatch_class


    def clean_unknown(self, db):
        """
        Removes rows for genes not found in the database (e.g., "ambiguous"
        fields from DESeq results)
        """
        ind = []
        for i, gene_id in enumerate(self.data.id):
            try:
                db[gene_id]
                ind.append(i)
            except gffutils.FeatureNotFoundError:
                continue
        return self.data.ix[ind, :]

    def enriched(self, pval=0.05, sig_column='padj', en_column='log2FoldChange', ind=True):
        """
        Identify enriched genes.

        :param pval: Alpha to use
        :param sig_column: Which column to use for significance
        :param en_column: Which column to use for enrichment
        :param ind: Boolean; if True will return just the index but if False,
            will return a copy
        """
        idx = (self[sig_column] < pval) & (self[en_column] > 0)
        if ind:
            return idx
        else:
            return self[idx]

    def scatter(self, xvar, yvar, xfunc=None, yfunc=None, xscale=None,
            yscale=None, xlab=None, ylab=None, genes_to_highlight=None,
            label_genes=False, pvals=None, general_kwargs=dict(color="k",
                alpha=0.2, linewidths=0), marginal=True, offset_kwargs={},
            label_kwargs=None, featuretypes=None, GFFDB=None, ax=None,
            one_to_one=None, callback=None):
        """
        Do-it-all method for making scatterplots.

        :df: A pandas.DataFrame object

        :param xvar, yvar:
            Variables to plot -- say, "baseMeanA" and "baseMeanB"

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
        """
        # Construct defaults---------------------------------------------------
        def identity(x):
            return x.copy()

        if xlab is None:
            if xfunc is not None:
                xlab = "%s(%s)" % (xfunc.__name__, xvar)
            else:
                xlab = "%s" % xvar

        if ylab is None:
            if yfunc is not None:
                ylab = "%s(%s)" % (yfunc.__name__, yvar)
            else:
                ylab = "%s" % yvar


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

        if featuretypes is None:
            featuretypes = {}

        # ---------------------------------------------------------------------

        x = self[xvar]
        y = self[yvar]

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

    def features_from_genes(self, func=None, *args, **kwargs):
        """
        Retrieves each gene from the database and stores it in the `features`
        column.

        Since pybedtools.Intervals are themselves iterable, we need to do some
        extra work to make sure an array is created with one entire interval
        per index.
        """
        if func is None:
            func = lambda x: x


        features = []
        dummy = pybedtools.Interval('None', 1, 1)
        for gene_id in self['id']:
            try:
                gene = func(asinterval(self.db[gene_id]), *args, **kwargs)
                features.append(gene)
            except gffutils.FeatureNotFoundError:
                features.append([dummy])
                continue

        x = np.empty(len(features), dtype=object)
        for i, feature in enumerate(features):
            x[i] = feature
        return x


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import numpy as np
    cols = ['id', 'baseMean', 'baseMeanA', 'baseMeanB', 'foldChange', 'log2FoldChange', 'pval', 'padj']

    fn = metaseq.example_filename('rrp6-s2-polyA.final.summary')
    fn2 = metaseq.example_filename('rrp40-s2-polyA.final.summary')

    dbfn = metaseq.example_filename('dmel-all-r5.33-cleaned.gff.db')


    #z = DESeqResults(fn, dbfn=dbfn)
    #x = DESeqResults(fn2, dbfn=dbfn)

    x = DESeqResults(dataframe(fn))
    y = DESeqResults(dataframe(fn2))

    p = pandas.Panel({'x': x, 'y': y})

    plt.show()
