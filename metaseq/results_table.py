"""
Module for working with tables of data, typically with at least one columns
containing gene accessions.

Results tables saved from DESeq are a perfect example of this, but any data can
be used.
"""

import sys
import tempfile
from matplotlib.mlab import csv2rec
import pybedtools
import gffutils
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
from scipy import stats
from rpy2.robjects import r
import rpy2.robjects as RO
import random
from gffutils.helpers import asinterval


class ResultsTable(object):
    def __init__(self, data, dbfn=None, id_column='id', csv2rec_kwargs=None):
        """
        Generic class for handling tables of data.

        :param data: If a string, assume it's a filename and load that using
            `csv2rec_kwargs`. Otherwise, assume it's a record array.
        :param dbfn: Filename for a `gffutils.FeatureDB`. Optional, but really
            handy.
        :param id_column: Which column contains gene accessions that can be
            looked up in the `gffutils.FeatureDB`.
        :param csv2rec_kwargs: Kwargs passed to `matplotlib.mlab.csv2rec`.
            Default is dict(delimiter="\\t", missing="NA").
        """
        if csv2rec_kwargs is None:
            csv2rec_kwargs = dict(delimiter='\t', missing='NA')

        if isinstance(data, basestring):
            data = csv2rec(data, **csv2rec_kwargs)

        self.id_column = id_column
        self.data = data
        self.dbfn = dbfn
        self.gffdb = gffutils.FeatureDB(dbfn)
        self._cached_lookup = None

    @property
    def colnames(self):
        return self.data.dtype.names

    def __repr__(self):
        return "<%s instance with %s items>" % (self.__class__.__name__,
                                                len(self.data))

    def __str__(self):
        fields = ['chrom', 'source', 'featuretype', 'start', 'end', 'score',
                'strand', 'frame', 'attributes']
        s = []
        for i, item in enumerate(self.data):
            d = dict(zip(self.colnames, item))
            d['_index'] = i
            try:
                feature = self.gffdb[item.id]
                d.update(zip(fields, feature.tostring().strip().split('\t')))
            except gffutils.FeatureNotFoundError:
                d.update({'attributes': 'Feature not found'})
            for key, val in d.items():
                s.append('%s: %s' % (key, val))
        return '\n'.join(s)

    def __len__(self):
        return len(self.data)

    def __getitem__(self, ind):
        orig_kwargs = dict(dbfn=self.dbfn, id_column=self.id_column)
        if isinstance(ind, int):
            if ind > len(self) - 1:
                raise IndexError
            new_instance = self.__class__(
                    self.data[ind:ind + 1].copy(), **orig_kwargs)
        else:
            new_instance = self.__class__(self.data[ind].copy(), **orig_kwargs)
        return new_instance

    def __getattr__(self, attr):
        if attr in self.__dict__:
            return getattr(self, attr)
        return getattr(self.data, attr)



    def strip_unknown_features(self):
        """
        Remove features not found in the `gffutils.FeatureDB`.  This will
        typically include 'ambiguous', 'no_feature', etc, but can also be
        useful if the database was created from a different one than was used
        to create the table.
        """
        ind = []
        for i, gene_id in enumerate(self.id):
            try:
                self.gffdb[gene_id]
                ind.append(i)
            except gffutils.FeatureNotFoundError:
                pass
        ind = np.array(ind)
        return self[ind]

    def random_subset(self, n, idx=True):
        """
        Random subset of all rows

        :param n: Number of rows to return
        :param idx: If True, return the index; if False, returns a subsetted
            version.
        """
        ind = random.sample(xrange(len(self.data)), n)
        if idx:
            return ind
        return self[ind]

    def sorted_by(self, attr, absolute=False, reverse=False):
        """
        Re-sort by an attribute and return a copy.

        :param attr: Attribute to sort by; must be a column in `self.colnames`
        :param absolute: If True, then ignore sign when sorting
        :param reverse: If True, highest values are first
        """
        vals = getattr(self, attr)
        if absolute:
            vals = abs(vals)
        ind = np.argsort(vals)
        if reverse:
            ind = ind[::-1]
        return self[ind]

    def histogram_of_hits(self, bed, field='log2foldchange',
            log=False, labels=None):
        """
        Plots a histogram of data values indicated by `field` for all genes
        with and without peaks in `bed`.
        """
        have_peaks = self.genes_with_peak(bed, as_ind=True)
        hits = getattr(self, field)[have_peaks]
        misses = getattr(self, field)[~have_peaks]
        if log:
            hits = np.log2(hits)
            misses = np.log2(misses)

        hits = hits[np.isfinite(hits)]
        misses = misses[np.isfinite(misses)]

        fig = plt.figure()
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)
        mmax = max(hits.max(), misses.max())
        bins = np.linspace(0, mmax, 50)
        kwargs = dict(bins=bins, color=(.5, .5, .5))

        ax1.hist(hits, **kwargs)
        ax2.hist(misses, **kwargs)

        ax1.set_title('genes with peaks in promoter')
        ax2.set_title('genes without peaks in promoter')

        ax1.set_xlabel(field)
        ax2.set_xlabel(field)

        ax1.set_ylabel('Number of genes\n(total=%s)' % len(hits),
                       va='center')
        ax2.set_ylabel('Number of genes\n(total=%s)' % len(misses),
                       va='center')

        fig.subplots_adjust(hspace=0.5)

        # Null hypothesis is that the two samples come from populations having
        # the same location (Sokal & Rolf p.427).
        #
        # Result is one-tailed pval; multiply by 2 to get two-tailed pval.
        s = stats.mannwhitneyu(hits, misses)
        results = {'U': s[0],
                   'pval': s[1]}
        ax2.text(x=.7,
                 y=.7,
                 s='Mann-Whitney U\np=%.3f' % s[1],
                 transform=ax2.transAxes)
        return results, fig

    def features(self, ignore_unknown=False):
        """
        Generator of currently-selected features.

        Looks up each feature in the attached `gffutils.FeatureDB` and converts
        it into a `pybedtools.Interval` object for use with `pybedtools`.
        Raises a warning if you haven't yet attached a `gffutils.FeatureDB` to
        this instance.

        :param ignore_unknown: If `ignore_unknown=False` then an exception will
            be raised if a feature cannot be found; if `ignore_unknown=True`
            then silently ignore these cases. Consider using the
            `strip_unknown_features()` method to handle these cases up front.
        """
        if not self.gffdb:
            raise ValueError('Please attach a GFF database created by '
                             'gffutils by setting the .gffdb attribute to the '
                             'database\'s path.')

        for i in self.data[self.id_column]:
            try:
                yield asinterval(self.gffdb[i])
            except gffutils.FeatureNotFoundError:
                if ignore_unknown:
                    continue
                else:
                    raise gffutils.FeatureNotFoundError('%s not found' % i.id)

    def align_with(self, other):
        """
        Ensure identical sorting of this object's data with another.

        Returns `self`, sorted the same way as `other`.

        :param other: Another instance of a ResultsTable or ResultsTable
            subclass.
        """
        ind = self.gene_ind(other[other.id_column])
        return self[ind]

    def scatter(self, x, y, xfunc=None, yfunc=None, xscale=None,
            yscale=None, xlab=None, ylab=None, genes_to_highlight=None,
            label_genes=False, general_kwargs=dict(color="k", alpha=0.2,
                linewidths=0), marginal=True, offset_kwargs={},
            label_kwargs=None, ax=None, one_to_one=None, callback=None,
            xlab_prefix=None, ylab_prefix=None):
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

    def genes_in_common(self, other):
        """
        Return a list of shared IDs.

        :param other: List of gene IDs, or another similar object
        """
        these = set(self[self.id_column].tolist())
        if isinstance(other, ResultsTable):
            those = set(other[other.id_column].tolist())
        else:
            those = set(other)
        common = these.intersection(those)
        return list(common)

    def gene_ind(self, genes, idx=True):
        """
        Returns an array of indices for `genes`.

        Useful for matching up two ResultsTable instances that are not
        guaranteed to have the same gene order (though they should have the
        same total gene set)

        :param genes: An iterable of feature accessions that are in the
            accession column.
        """
        # make a dictionary mapping current gene to index.
        if not self._cached_lookup:
            self._cached_lookup = dict(zip(self.data.id, np.arange(len(self.data))))
        ind = []
        for gene in genes:
            ind.append(self._cached_lookup[gene])
        ind = np.array(ind)
        if idx:
            return ind
        return self[ind]


class DESeqResults(ResultsTable):
    def __init__(self, *args, **kwargs):
        """
        Subclass of :class:`ResultsTable` specifically for working with DESeq
        results.
        """
        super(DESeqResults, self).__init__(*args, **kwargs)

    def strip_deseq_nongenes(self):
        """
        DESeq adds "no_feature", "not_aligned", etc. features.  This method
        removes them for better plotting.
        """
        to_remove = [
                'no_feature',
                'not_aligned',
                'alignment_not_unique',
                'ambiguous',
                'too_low_aQual']
        remove_ind = self.gene_ind(to_remove)
        keep_ind = []
        for i in range(len(self)):
            if i not in remove_ind:
                keep_ind.append(i)
        return self[keep_ind]

    def enriched(self, pval=0.05, column='padj', idx=True):
        """
        Enriched genes at `pval` significance.

        :param pval: Alpha to use as a cutoff
        :param column: Column to apply cutoff to
        :param idx: If True, return the index; if False, returns a subsetted
            version.
        """
        ind1 = self.data[column] <= pval
        ind2 = self.log2foldchange > 0
        ind = ind1 & ind2
        if idx:
            return ind
        return self[ind]

    def disenriched(self, pval=0.05, column='padj', idx=True):
        """
        Disenriched genes at `pval` significance.

        :param pval: Alpha to use as a cutoff
        :param column: Column to apply cutoff to
        :param idx: If True, return the index; if False, returns a subsetted
            version.
        """
        ind1 = self.data[column] <= pval
        ind2 = self.log2foldchange < 0
        ind = ind1 & ind2
        if idx:
            return ind
        return self[ind]

    def nonsig(self, pval=0.05, column='padj', idx=True):
        """
        Non-significant genes (that were still expressed at some level)

        :param pval: Alpha to use as a cutoff
        :param column: Column to apply cutoff to
        :param idx: If True, return the index; if False, returns a subsetted
            version.
        """
        ind = (self.data[column] > pval) & (self.basemean > 0)
        if idx:
            return ind
        return self[ind]

    def random_nonsig(self, n, pval=0.05, column='padj'):
        """
        Random subset of nonsignifcant genes at `pval` significance
        """
        # get inds of nonsig as bool
        ind1 = self.nonsig(pval=pval, column=column)

        # indices as integers
        ind2 = np.nonzero(ind1)[0]

        # indices to keep
        keepers = random.sample(ind2, n)

        # array of all False
        final_ind = np.ones_like(ind1) == 0

        # Only set the keepers to True
        final_ind[keepers] = True

        if idx:
            return final_ind
        return self[final_ind]


def hypergeom(m, n, n1, n2):
    """
    m = overlapping genes
    n = total genes that could be sampled
    n1 = number of genes in set 1
    n2 = number of genes in set 2
    """
    if m == 0:
        return 1.0
    return r['phyper'](min(n1, n2), n1, n - n1, n2)[0] \
            - r['phyper'](m - 1, n1, n - n1, n2)[0]


def hypergeom_scipy(m, n, n1, n2, p=False):
    """
    Given gene counts `n1` and `n2`, each drawn from `n` total genes, return
    the probability that `m` genes would be shared due random chance alone.

    e.g.,

    n1 = 100  # significantly enriched genes from sample 1
    n2 = 50   # significantly enriched genes from sample 2
    n = 15000 # total number of genes that could be sampled
    m = 10  # number of genes that overlap in the two lists

    See: http://www.nslij-genetics.org/wli/pub/ieee-embs06.pdf

    Thanks to Brent Pedersen (https://github.com/brentp/bio-playground) for
    implementation.
    >>> hypergeom(1, 1000, 1000, 1000) # has to be shared.
    1.0

    >>> all(hypergeom(i, 1000, 1000, 1000) == 1.0 for i in range(100))
    True

    >>> hypergeom(1, 30000, 20, 20)
    0.013253396616299651

    >>> hypergeom(2, 30000, 20, 20)
    7.9649366037104485e-05

    >>> hypergeom(11, 30000, 20, 20)
    4.516176321800458e-11

    >>> hypergeom(10, 30000, 20, 20) # very low prob.
    4.516176321800458e-11

    >>> hypergeom(20, 30000, 20, 20) # very low chance that all are shared.
    4.516176321800458e-11

    """
    if m <= 0:
        return 1.0
    mmin = m - 1
    mmax = min(n1, n2)
    return stats.hypergeom.cdf(mmax, n, n1, n2) \
            - stats.hypergeom.cdf(mmin, n, n1, n2)
