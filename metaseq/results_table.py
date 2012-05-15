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
from helpers import gfffeature_to_interval


class ResultsTable(object):
    def __init__(self, data, dbfn=None):
        """
        If `data` is a string, assume it's a filename and load that.
        Otherwise, assume it's a record array.
        """
        if isinstance(data, basestring):
            data = csv2rec(data, delimiter='\t', missing='NA')
        self.data = data
        self.dbfn = dbfn
        self.gffdb = gffutils.FeatureDB(dbfn)
        self._gene_ind_cache = None

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
        """
        Return a copy of self, attaching a GFFDB etc
        """
        if isinstance(ind, int):
            if ind > len(self) - 1:
                raise IndexError
            new_instance = ResultsTable(self.data.copy()[ind:ind + 1],
                    dbfn=self.dbfn)
        else:
            new_instance = ResultsTable(self.data.copy()[ind], dbfn=self.dbfn)
        return new_instance

    def __getattr__(self, attr):
        try:
            return getattr(self.data, attr)
        except AttributeError:
            raise AttributeError('ResultsTable has no attribute "%s"' % attr)

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

    def strip_unknown_features(self):
        """
        Remove features not found in the GFFDB.  This will typically include
        'ambiguous', 'no_feature', etc, but can also be useful if the GFFDB was
        created from a different one than was used to count reads
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

    def gene_data(self, gene_id):
        if self._gene_ind_cache is None:
            self._gene_ind_cache = {}
            for i, name in enumerate(self.id):
                self._gene_ind_cache[name] = i
        try:
            return self[self._gene_ind_cache[gene_id]]
        except KeyError:
            raise KeyError('No gene %s found' % gene_id)

    def enriched(self, pval=0.05, column='padj'):
        """
        Index of enriched genes at `pval` significance
        """
        ind1 = self.data[column] <= pval
        ind2 = self.log2foldchange > 0
        return ind1 & ind2

    def disenriched(self, pval=0.05, column='padj'):
        """
        Index of disenriched genes at `pval` significance
        """
        ind1 = self.data[column] <= pval
        ind2 = self.log2foldchange < 0
        return ind1 & ind2

    def nonsig(self, pval=0.05, column='padj'):
        """
        Index of non-significant genes at `pval` significance
        """
        return (self.data[column] > pval) & (self.basemean > 0)

    def random_subset(self, n):
        """
        Random subset of all available genes
        """
        return random.sample(xrange(len(self.data)), n)

    def random_nonsig(self, pval, n):
        """
        Random subset of nonsignifcant genes at `pval` significance
        """
        # get inds of nonsig as bool
        ind1 = self.nonsig(pval)

        # indices as integers
        ind2 = np.nonzero(ind1)[0]

        # indices to keep
        keepers = random.sample(ind2, n)

        # array of all False
        final_ind = np.ones_like(ind1) == 0

        # Only set the keepers to True
        final_ind[keepers] = True

        return final_ind

    @property
    def colnames(self):
        return self.data.dtype.names

    def sorted_by(self, attr, absolute=False, reverse=False):
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

    def features_to_gff(self, fn=None):
        """
        Save features as a GFF file on disk as `fn`. If `fn` is None, then save
        to a temp file.  Always returns the filename.
        """
        if fn is None:
            fn = tempfile.mktemp()
        fout = open(fn, 'w')
        for gene_id in self.data.id:
            try:
                gene = self.gffdb[gene_id]
                fout.write(str(gene) + '\n')
            except gffutils.FeatureNotFoundError:
                pass
        fout.close()
        return fn

    def features(self):
        """
        Generator of currently-selected features from the GFF database.  Raises
        a warning if you haven't yet attached a GFFDB to this instance.
        """
        if not self.gffdb:
            raise ValueError('Please attach a GFF database created by '
                             'gffutils by setting the .gffdb attribute to the '
                             'database\'s path.')

        for i in self.data.id:
            try:
                yield gfffeature_to_interval(self.gffdb[i])
            except gffutils.FeatureNotFoundError:
                yield None
                #sys.stderr.write('no data for "%s", skipping...\n' % i)
                #sys.stderr.flush()

    def genes_with_peak(self, bed, as_ind=False):
        """
        Given a `bed` file (say, of peaks), search for peaks in promoters of
        the currently-selected-and-sorted genes.

        Returns a list of ALL genes, with their score attribute set to the
        number of peaks.

        If `as_ind` is True, then don't return a list of genes -- just return
        True/False if a peak was found.
        """
        hits = []

        x = pybedtools.IntervalFile(bed)
        for feature in self.features():
            c = 0
            if feature is not None:
                if isinstance(feature, gffutils.Feature):
                    interval = gfffeature_to_Interval(feature)
                else:
                    interval = feature
                c += x.count_hits(interval)
            if as_ind:
                hits.append(c > 0)
            else:
                hits.append(feature)
        if as_ind:
            return np.array(hits)
        return hits

    def genes_in_common(self, other):
        """
        Return a list of genes shared between this ResultsTable instance and
        another one.  If `other` is a ResultsTable object, then it will get the
        id field; if it's a list-like object, it will be converted to a set.
        """
        these = set(self.data.id.tolist())
        if isinstance(other, ResultsTable):
            those = set(other.data.id.tolist())
        else:
            those = set(other)
        common = these.intersection(those)
        return list(common)

    def gene_ind(self, genes):
        """
        Returns an array of indices for `genes`.  Useful for matching up two
        ResultsTable instances that are not guaranteed to have the same gene
        order.

        If a gene in the list cannot be found, then the index will be NaN so
        that the returned index and the length of `genes` will always be the
        same.
        """
        # make a dictionary mapping current gene to index.
        d = dict(zip(self.data.id, np.arange(len(self.data))))
        ind = []
        for gene in genes:
            try:
                ind.append(d[gene])
            except KeyError:
                ind.append(np.NaN)
        return np.array(ind)

    def chisq(self, bed, obs_ind, exp_ind, as_string=False,
            title=""):
        """
        Chi-square test of gene with peaks in BED file `bed`.

        `obs_ind` will typically be self.enriched(); `exp_ind` will typically
        be self.nonsig().

        If `as_string` is True, then return a nicely formatted string instead
        of a dictionary.  This will include the optional `title`, too.
        """
        # construct temp files of the genes of interest...

        obs_w_peak = sum(self[obs_ind].genes_with_peak(bed,
            as_ind=True))
        exp_w_peak = sum(self[exp_ind].genes_with_peak(bed,
            as_ind=True))
        total_obs = float(sum(obs_ind))
        total_exp = float(sum(exp_ind))
        x = RO.FloatVector((obs_w_peak, total_obs - obs_w_peak))
        p = RO.FloatVector((exp_w_peak, total_exp - exp_w_peak))
        res = r['chisq.test'](x=x, p=p, rescale_p=True)
        pval = res.rx('p.value')[0][0]
        d = {
               'total_obs': total_obs,
               'total_exp': total_exp,
              'obs_w_peak': obs_w_peak,
              'exp_w_peak': exp_w_peak,
               'obs_ratio': obs_w_peak / total_obs,
               'exp_ratio': exp_w_peak / total_exp,
                    'pval': pval}
        if not as_string:
            return d

        d['title'] = title
        d['underline'] = '-' * (len(title) + 1)
        return """
%(title)s:
%(underline)s
      \tobs\texp
total \t%(total_obs)i\t%(total_exp)i
w/peak\t%(obs_w_peak)s\t%(exp_w_peak)s
ratio \t%(obs_ratio).4f\t%(exp_ratio).4f
%(underline)s
pval: %(pval)s
%(underline)s
""" % d

    def hgpval(self, other, totalgenes=None):
        """
        Returns the hypergeometric pval of this ResultsTable object with
        another.  `totalgenes` is the number of possible genes, or the total
        number of "gene" features in the GFF database if None.

        Low pvals means the sets of genes share significantly more genes than
        would be expected by chance alone.
        """
        if totalgenes is None:
            totalgenes = self.gffdb.count_features_of_type('gene')
        n1 = len(self.data)

        if isinstance(other, ResultsTable):
            n2 = len(other.data)
            fb = frozenset(other.data.id)
        else:
            n2 = len(other)
            fb = frozenset(other)

        fa = frozenset(self.data.id)
        m = len(fa.intersection(fb))
        pval = hypergeom(m, totalgenes, n1, n2)
        results = {
                           'pval': pval,
                    'total genes': totalgenes,
                           'set1': n1,
                           'set2': n2,
                        'overlap': m
                  }
        return results

    def scatter(self, xvar, yvar, xfunc=np.log2, yfunc=np.log2, xscale=None,
                yscale=None, xlab=None, ylab=None, genes_to_highlight=None,
                label_genes=False, pvals=None, outfn=None, general_kwargs=None,
                force_gene_pval_color=True, color_by_variance=False,
                marginal=True, offset_kwargs={}, label_kwargs=None,
                featuretypes=None, GFFDB=None, ax=None, one_to_one=None,
                callback=None):
        """
        Do-it-all method for making scatterplots.

        `xvar` and `yvar` are the variables to plot -- say, "basemeana" and
        "basemeanb"

        `xfunc` and `yfunc` are the functions to apply to x and y respectively.

        `xscale` and `yscale`, if not None, will be multiplied before applying
        xfunc and yfunc.

        If ax=None, then makes a new fig and returns the Axes object,
        otherwise, plots onto `ax`

        `general_kwargs` specifies how all points look

        `genes_to_highlight` provides lots of control to colors.  It is a list
        of (ind, kwargs) tuples, where each ind specifies genes to plot with
        `kwargs`.  Each dictionary updates a copy of `general_kwargs`.

        If `genes_to_highlight` has a "name" kwarg, this must be a list that't
        the same length as `ind`.  It will be used to label the genes in `ind`

        `marginal` toggles display of non-finite cases where x or y is nonzero
        but the other one is zero (and `xfunc` or `yfunc` are log)

        `callback` is a function to call upon clicking a point.  Default is to
        print the GFFFeature of the clicked point.

        """
        # Construct defaults---------------------------------------------------
        if general_kwargs is None:
            general_kwargs = {}

        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)

        if label_kwargs is None:
            label_kwargs = dict(horizontalalignment='right',
                    verticalalignment='center', style='italic',
                    bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))

        if featuretypes is None:
            featuretypes = {}

        if xlab is None:
            xlab = "%s(%s)" % (xfunc.__name__, xvar)

        if ylab is None:
            ylab = "%s(%s)" % (yfunc.__name__, yvar)

        if callback is None:
            def callback(event):
                for ind in event.ind:
                    sys.stderr.write(str(self[ind]))
                    sys.stderr.flush()

        # ---------------------------------------------------------------------

        x = self.data[xvar].copy()
        y = self.data[yvar].copy()
        if xscale is not None:
            x *= xscale

        if yscale is not None:
            y *= scale

        xi = xfunc(x)
        yi = yfunc(y)

        xv = np.isfinite(xi)
        yv = np.isfinite(yi)

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
        ax.figure.canvas.mpl_connect('pick_event', callback)

        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)

        return ax

    def filter_features(self, features):
        """
        Given a list of pybedtools.Interval features (say, from the output of
        metaseq.MetaChip.array_parallel), return the subset of those features
        that are also in this ResultsTable instance.

        For example,

        enriched = self[self.enriched()]
        enriched_tss = enriched.filter_gene_list(all_tss_features)

        Note that this matches **names**, so your accessions need to match
        up....
        """
        these_features = set(self.id)

        def in_this_instance(feature):
            return feature.name in these_features
        return filter(in_this_instance, features)


def rank_plot(genelists, height=0.5, spacing=0.2):
    """
    strip plot.  `genelists` is a list.  Each item in the list contains the
    following dictionary:

        genelist: list of features with .npeaks attribute
        colorfunc: function to determine box color given score.
        kwargs: Dictionary of kwargs that are passed to BrokenBarHCollection.
                Default is to pass only the generated colors.
        label: label to use as ytick label

    Example usage::

        >>> def colorfunc1(x):
        ...     if x.score > 0:
        ...        return 'r'
        ...    else:
        ...        return 'w'
        >>> def colorfunc2(x):
        ...     if x.score > 0:
        ...        return 'b'
        ...    else:
        ...        return 'w'
        >>> genelists = []
        >>> genelists.append(dict(genelist=hits1, colorfunc=colorfunc1,
        ...                  kwargs=dict(linewidths=0.3)))
        >>> genelists.append(dict(genelist=hits2, colorfunc=colorfunc2,
        ...                  kwargs=dict(linewidths=0.3)))
        >>> rank_plot(genelists)

    """
    fig = plt.figure(figsize=(12, 1.8), facecolor='w')
    ax = fig.add_subplot(111)
    ymin = 0
    yticks = []
    yticklabels = []
    for i, data in enumerate(genelists[::-1]):
        genelist = data['genelist']
        try:
            colorfunc = data['colorfunc']
        except KeyError:
            def colorfunc(x):
                if x is not None:
                    if x.score > 0:
                        return 'r'
                return 'w'
        try:
            kwargs = data['kwargs'].copy()
        except KeyError:
            kwargs = {}

        colors = [colorfunc(k) for k in genelist]
        if 'facecolors' not in kwargs:
            kwargs['facecolors'] = colors

        xranges = [(j, 1) for j in range(len(genelist))]
        ywidth = height
        ax.add_collection(
                matplotlib.collections.BrokenBarHCollection(xranges,
                                                            (ymin, ywidth),
                                                            **kwargs))
        yticks.append(ymin + height / 2.)
        try:
            yticklabels.append(data['label'])
        except KeyError:
            yticklabels.append('')
        ymin += height + spacing
    ax.axis('tight')
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklabels)
    ax.set_xlabel('ranked genes')
    fig.subplots_adjust(bottom=0.5, top=0.8)
    return fig


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
