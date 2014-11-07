import copy
from textwrap import dedent
import numpy as np
import pandas
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as patches
import matplotlib
import plotutils
from matplotlib.transforms import blended_transform_factory
from matplotlib.collections import EventCollection
import gffutils
import pybedtools
from pybedtools import featurefuncs

_base_doc = """%s
The underlying pandas.DataFrame is always available with the `data`
attribute.

Any attributes not explicitly in this class will be looked for in the
underlying pandas.DataFrame.

Parameters
----------
data : string or pandas.DataFrame
    If string, assumes it's a filename and calls
    pandas.read_table(data, **import_kwargs).

db : string or gffutils.FeatureDB
    Optional database that can be used to generate features

import_kwargs : dict
    These arguments will be passed to pandas.read_table() if `data` is
    a filename.
"""


class ResultsTable(object):
    __doc__ = _base_doc % dedent(
        """
        Wrapper around a pandas.DataFrame that adds additional functionality.
        """)
    def __init__(self, data, db=None, import_kwargs=None):
        if isinstance(data, basestring):
            import_kwargs = import_kwargs or {}
            data = pandas.read_table(data, **import_kwargs)
        if not isinstance(data, pandas.DataFrame):
            raise ValueError("`data` is not a pandas.DataFrame")
        self.data = data

        self._kwargs = dict(db=db, import_kwargs=import_kwargs)
        self.attach_db(db)
        self._cached_features = None

    def __getattr__(self, attr):
        if attr in self.__dict__:
            return getattr(self, attr)
        return getattr(self.data, attr)

    def __getitem__(self, attr):
        if isinstance(attr, basestring):
            return self.data.__getitem__(attr)
        else:
            return self.__class__(self.data.__getitem__(attr), **self._kwargs)

    def update(self, dataframe):
        """
        Updates the current data with a new dataframe.

        This extra step is required to get around the fancy pandas.DataFrame
        indexing (like .ix, .iloc, etc).
        """
        return self.__class__(dataframe, **self._kwargs)

    def copy(self):
        data = self.data.copy(deep=True)
        return self.__class__(data, db=self.db, import_kwargs=self._kwargs)

    def __repr__(self):
        s = []
        s.append("<%s instance, wrapping the following:"
                 % self.__class__.__name__)
        s.append('')
        s.extend('\t' + i for i in repr(self.data).splitlines(False))
        s.append('>')
        return '\n'.join(s)

    def attach_db(self, db):
        """
        Attach a gffutils.FeatureDB for access to features.

        Useful if you want to attach a db after this instance has already been
        created.

        Parameters
        ----------
        db : gffutils.FeatureDB
        """
        if db is not None:
            if isinstance(db, basestring):
                db = gffutils.FeatureDB(db)
            if not isinstance(db, gffutils.FeatureDB):
                raise ValueError(
                    "`db` must be a filename or a gffutils.FeatureDB")
        self._kwargs['db'] = db
        self.db = db

    def features(self, ignore_unknown=False):
        """
        Generator of features.

        If a gffutils.FeatureDB is attached, returns a pybedtools.Interval for
        every feature in the dataframe's index.

        Parameters
        ----------
        ignore_unknown : bool
            If True, silently ignores features that are not found in the db.
        """
        if not self.db:
            raise ValueError("Please attach a gffutils.FeatureDB")
        for i in self.data.index:
            try:
                yield gffutils.helpers.asinterval(self.db[i])
            except gffutils.FeatureNotFoundError:
                if ignore_unknown:
                    continue
                else:
                    raise gffutils.FeatureNotFoundError('%s not found' % i)

    def reindex_to(self, x, attribute="Name"):
        """
        Returns a copy that only has rows corresponding to feature names in x.

        Parameters
        ----------
        x : str or pybedtools.BedTool
            BED, GFF, GTF, or VCF where the "Name" field (that is, the value
            returned by feature['Name']) or any arbitrary attribute

        attribute : str
            Attribute containing the name of the feature to use as the index.
        """
        names = [i[attribute] for i in x]
        new = self.copy()
        new.data = new.data.reindex(names)
        return new

    def five_prime(self, upstream=1, downstream=0):
        """
        Creates a BED/GFF file of the 5' end of each feature represented in the
        table and returns the resulting pybedtools.BedTool object. Needs an
        attached database.

        Parameters
        ----------
        upstream, downstream : int
            Number of basepairs up and downstream to include
        """
        return pybedtools.BedTool(self.features())\
            .each(featurefuncs.five_prime, upstream, downstream)\
            .saveas()

    def three_prime(self, upstream=0, downstream=1):
        """
        Creates a BED/GFF file of the 3' end of each feature represented in the
        table and returns the resulting pybedtools.BedTool object. Needs an
        attached database.

        Parameters
        ----------
        upstream, downstream : int
            Number of basepairs up and downstream to include
        """
        return pybedtools.BedTool(self.features())\
            .each(featurefuncs.three_prime, upstream, downstream)\
            .saveas()

    TSS = five_prime
    TTS = three_prime

    def align_with(self, other):
        """
        Align the dataframe's index with another.
        """
        return self.__class__(self.data.reindex_like(other), **self._kwargs)

    def genes_in_common(self, other):
        """
        Convenience method for getting the genes found in both dataframes.
        """
        return self.index & other.index

    def __and__(self, other):
        return self.index & other.index

    def __or__(self, other):
        return self.index | other.index

    def __sub__(self, other):
        return self.index - other.index

    def __len__(self):
        return len(self.data)

    def scatter(self, x, y, xfunc=None, yfunc=None, xscale=None, yscale=None,
                xlab=None, ylab=None, genes_to_highlight=None,
                label_genes=False, marginal_histograms=False,
                general_kwargs=dict(color="k", alpha=0.2, picker=True),
                general_hist_kwargs=None, offset_kwargs={}, label_kwargs=None,
                ax=None, one_to_one=None, callback=None, xlab_prefix=None,
                ylab_prefix=None, sizefunc=None, hist_size=0.3, hist_pad=0.0,
                nan_offset=0.015, pos_offset=0.99, linelength=0.01,
                neg_offset=0.005, figure_kwargs=None):
        """
        Do-it-all method for making annotated scatterplots.

        Parameters
        ----------

        x, y : array-like
            Variables to plot.  Must be names in self.data's DataFrame.  For
            example, "baseMeanA" and "baseMeanB"

        xfunc, yfunc : callable
            Functions to apply to `xvar` and `yvar` respectively. Default is
            log2; set to None to have no transformation.

        xlab, ylab : string
            Labels for x and y axes; default is to use function names for
            `xfunc` and `yfunc` and variable names `xvar` and `yvar`, e.g.,
            "log2(baseMeanA)"

        ax : None or Axes object
            If `ax=None`, then makes a new fig and returns the Axes object,
            otherwise, plots onto `ax`

        general_kwargs : dict
            Kwargs for matplotlib.scatter; specifies how all points look

        genes_to_highlight : list of (index, dict) tuples
            Provides lots of control to colors.  It is a list of (`ind`,
            `kwargs`) tuples, where each `ind` specifies genes to plot with
            `kwargs`.  Each dictionary updates a copy of `general_kwargs`. If
            `genes_to_highlight` has a "name" kwarg, this must be a list that't
            the same length as `ind`.  It will be used to label the genes in
            `ind` using `label_kwargs`.

        callback : callable
            Function to call upon clicking a point. Must accept a single
            argument which is the gene ID. Default is to print the gene name,
            but an example of another useful callback would be a mini-browser
            connected to a genomic_signal object from which the expression data
            were calculated.

        one_to_one : None or dict
            If not None, a dictionary of matplotlib.plot kwargs that will be
            used to plot a 1:1 line.

        label_kwargs : dict
            Kwargs for labeled genes (e.g., dict=(style='italic')).  Will only
            be used if an entry in `genes_to_highlight` has a `name` key.

        offset_kwargs : dict
            Kwargs to be passed to matplotlib.transforms.offset_copy, used for
            adjusting the positioning of gene labels in relation to the actual
            point.

        xlab_prefix, ylab_prefix : str
            Optional label prefix that will be added to the beginning of `xlab`
            and/or `ylab`.

        hist_size : float
            Size of marginal histograms

        hist_pad : float
            Spacing between marginal histograms

        nan_offset, pos_offset, neg_offset : float
            Offset, in units of "fraction of axes" for the NaN, +inf, and -inf
            "rug plots"

        linelength : float
            Line length for the rug plots

        """
        _x = self.data[x]
        _y = self.data[y]

        # Construct defaults---------------------------------------------------
        def identity(x):
            return x.copy()

        # Axis label setup
        if xlab_prefix is None:
            xlab_prefix = ""

        if ylab_prefix is None:
            ylab_prefix = ""

        if xlab is None:
            xlab = x

            if xfunc is not None:
                xlab = xlab_prefix + "%s(%s)" % (xfunc.__name__, str(x))
            else:
                xlab = xlab_prefix + "%s" % (str(x))

        if ylab is None:
            ylab = y
            if yfunc is not None:
                ylab = ylab_prefix + "%s(%s)" % (yfunc.__name__, str(y))
            else:
                ylab = ylab_prefix + "%s" % (str(y))

        if xfunc is None:
            xfunc = identity

        if yfunc is None:
            yfunc = identity

        if general_kwargs is None:
            general_kwargs = {}

        if general_hist_kwargs is None:
            general_hist_kwargs = {}

        if genes_to_highlight is None:
            genes_to_highlight = []

        if ax is None:
            if figure_kwargs is None:
                figure_kwargs = {}
            fig = plt.figure(**figure_kwargs)
            ax = fig.add_subplot(111)

        if label_kwargs is None:
            label_kwargs = dict(
                horizontalalignment='right',
                verticalalignment='center',
                style='italic',
                bbox=dict(facecolor='w', edgecolor='None', alpha=0.5)
            )

        # Clean data ---------------------------------------------------------
        xi = xfunc(_x)
        yi = yfunc(_y)

        # handle inf, -inf, and NaN
        x_is_pos_inf = np.isinf(xi) & (xi > 0)
        x_is_neg_inf = np.isinf(xi) & (xi < 0)
        x_is_nan = np.isnan(xi)
        y_is_pos_inf = np.isinf(yi) & (yi > 0)
        y_is_neg_inf = np.isinf(yi) & (yi < 0)
        y_is_nan = np.isnan(yi)

        # Indexes for valid values
        x_valid = ~(x_is_pos_inf | x_is_neg_inf | x_is_nan)
        y_valid = ~(y_is_pos_inf | y_is_neg_inf | y_is_nan)

        # global min/max
        gmin = max(xi[x_valid].min(), yi[y_valid].min())
        gmax = min(xi[x_valid].max(), yi[y_valid].max())

        # Convert any integer indexes into boolean, and create a new list of
        # genes to highlight.  This handles optional hist kwargs.

        # We'll compile a new list of genes to highlight.
        _genes_to_highlight = []

        for block in genes_to_highlight:
            ind = block[0]

            # Convert to boolean
            if ind.dtype != 'bool':
                new_ind = (np.zeros_like(xi) == 0)
                new_ind[ind] = True
                _genes_to_highlight.append(
                    tuple([new_ind] + list(block[1:]))
                )

            # If it's a DataFrame, we only want the boolean values;
            else:
                if hasattr(ind, 'values'):
                    ind = ind.values
                _genes_to_highlight.append(
                    tuple([ind] + list(block[1:]))
                )

        # Now we remove any genes from in allind (which will be plotted using
        # `general_kwargs`) that will be plotted by genes_to_highlight.  This
        # avoids double-plotting.
        allind = np.zeros_like(xi) == 0
        for block in _genes_to_highlight:
            ind = block[0]
            allind[ind] = False

        # Copy over the color and alpha if they're not specified
        general_hist_kwargs = plotutils._updatecopy(
            orig=general_hist_kwargs, update_with=general_kwargs,
            keys=['color', 'alpha'])

        # Put the non-highlighted genes at the beginning of _genes_to_highlight
        # list so we can just iterate over one list
        _genes_to_highlight.insert(
            0,
            (allind, general_kwargs, general_hist_kwargs)
        )

        # Set up the object that will handle the marginal histograms
        self.marginal = plotutils.MarginalHistScatter(
            ax, hist_size=hist_size, pad=hist_pad)

        # Set up kwargs for x and y rug plots
        rug_x_kwargs = dict(
            linelength=linelength,
            transform=blended_transform_factory(ax.transData, ax.transAxes),
        )
        rug_y_kwargs = dict(
            linelength=linelength,
            transform=blended_transform_factory(ax.transAxes, ax.transData),
            orientation='vertical',
        )

        # EventCollection objects need a color as a 3-tuple, so set up
        # a converter here.
        color_converter = matplotlib.colors.ColorConverter().to_rgb

        # Plot the one-to-one line, if kwargs were specified
        if one_to_one:
            ax.plot([gmin, gmax],
                    [gmin, gmax],
                    **one_to_one)

        # Plot 'em all, and label if specified

        # In order to avoid calling the callback function multiple times when
        # we have overlapping genes to highlight (e.g., a gene that is both
        # upregulated AND has a peak), keep track of everything that's been
        # added so far.
        self._seen = np.ones_like(xi) == 0
        for block in _genes_to_highlight:
            ind = block[0]
            kwargs = block[1]

            if len(block) == 3:
                hist_kwargs = block[2]
            else:
                hist_kwargs = {}

            names = kwargs.pop('names', None)
            _marginal_histograms = (
                kwargs.pop('marginal_histograms', False) or
                marginal_histograms)

            updated_kwargs = plotutils._updatecopy(
                orig=kwargs, update_with=general_kwargs)

            updated_hist_kwargs = plotutils._updatecopy(
                orig=hist_kwargs, update_with=general_hist_kwargs)
            updated_hist_kwargs = plotutils._updatecopy(
                orig=updated_hist_kwargs, update_with=kwargs,
                keys=['color', 'alpha'], override=True)

            xhist_kwargs = updated_kwargs.pop('xhist_kwargs', None)
            yhist_kwargs = updated_kwargs.pop('yhist_kwargs', None)

            self.marginal.append(
                xi[ind & x_valid & y_valid],
                yi[ind & x_valid & y_valid],
                scatter_kwargs=dict(**updated_kwargs),
                hist_kwargs=updated_hist_kwargs,
                xhist_kwargs=xhist_kwargs,
                yhist_kwargs=yhist_kwargs,
                marginal_histograms=_marginal_histograms,
            )

            # This is important for callbacks: here we grab the last-created
            # collection,
            coll = self.marginal.scatter_ax.collections[-1]
            coll.df = self.data
            coll.ind = ind & x_valid & y_valid

            color = color_converter(updated_kwargs['color'])
            rug_x_kwargs['color'] = color
            rug_y_kwargs['color'] = color

            # Note: if both x and y are not valid, then they will not be on the
            # plot.
            items = [
                # top rug, y is +inf and x is valid
                (xi, ind & x_valid & y_is_pos_inf, pos_offset, rug_x_kwargs),

                # one of the bottom rugs, where y is NaN
                (xi, ind & x_valid & y_is_nan, nan_offset, rug_x_kwargs),

                # bottom rug, y is -inf
                (xi, ind & x_valid & y_is_neg_inf, neg_offset, rug_x_kwargs),

                # right rug, x is +inf
                (yi, ind & y_valid & x_is_pos_inf, pos_offset, rug_y_kwargs),

                # one of the left rugs; x is NaN
                (yi, ind & y_valid & x_is_nan, nan_offset, rug_y_kwargs),

                # left rug, x is -inf
                (yi, ind & y_valid & x_is_neg_inf, neg_offset, rug_y_kwargs),
            ]
            for values, index, offset, kwargs in items:
                coll = EventCollection(
                    values[index], lineoffset=offset, **kwargs)
                coll.df = self.data
                coll.ind = index
                ax.add_collection(coll)

            if names:
                transOffset = matplotlib.transforms.offset_copy(
                    ax.transData, fig=ax.figure, **offset_kwargs)

                for xii, yii, name in zip(xi[ind], yi[ind], names):
                    ax.text(xii,
                            yii,
                            name,
                            transform=transOffset,
                            **label_kwargs)

        # register callback
        if callback is None:
            callback = self._default_callback

        def wrapped_callback(event):
            for _id in self._id_callback(event):
                callback(_id)

        ax.figure.canvas.mpl_connect('pick_event', wrapped_callback)

        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        ax.axis('tight')

        return ax

    def radviz(self, column_names, transforms=dict(), **kwargs):
        """
        Radviz plot.

        Useful for exploratory visualization, a radviz plot can show
        multivariate data in 2D.  Conceptually, the variables (here, specified
        in `column_names`) are distributed evenly around the unit circle.  Then
        each point (here, each row in the dataframe) is attached to each
        variable by a spring, where the stiffness of the spring is proportional
        to the value of corresponding variable.  The final position of a point
        represents the equilibrium position with all springs pulling on it.

        In practice, each variable is normalized to 0-1 (by subtracting the
        mean and dividing by the range).

        This is a very exploratory plot.  The order of `column_names` will
        affect the results, so it's best to try a couple different orderings.
        For other caveats, see [1].

        Additional kwargs are passed to self.scatter, so subsetting, callbacks,
        and other configuration can be performed using options for that method
        (e.g., `genes_to_highlight` is particularly useful).

        Parameters
        ----------
        column_names : list
            Which columns of the dataframe to consider.  The columns provided
            should only include numeric data, and they should not contain any
            NaN, inf, or -inf values.

        transforms : dict
            Dictionary mapping column names to transformations that will be
            applied just for the radviz plot.  For example, np.log1p is
            a useful function. If a column name is not in this dictionary, it
            will be used as-is.

        ax : matplotlib.Axes
            If not None, then plot the radviz on this axes.  If None, then
            a new figure will be created.

        kwargs : dict
            Additional arguments are passed to self.scatter.  Note that not all
            possible kwargs for self.scatter are necessarily useful for
            a radviz plot (for example, margninal histograms would not be
            meaningful).

        Notes
        -----
        This method adds two new variables to self.data: "radviz_x" and
        "radviz_y".  It then calls the self.scatter method, using these new
        variables.

        The data transformation was adapted from the
        pandas.tools.plotting.radviz function.

        References
        ----------
        [1]  Hoffman,P.E. et al. (1997) DNA visual and analytic data mining. In
             the Proceedings of the IEEE Visualization. Phoenix, AZ, pp.
             437-441.
        [2] http://www.agocg.ac.uk/reports/visual/casestud/brunsdon/radviz.htm
        [3] http://pandas.pydata.org/pandas-docs/stable/visualization.html\
                #radviz
        """
        # make a copy of data
        x = self.data[column_names].copy()

        for k, v in transforms.items():
            x[k] = v(x[k])

        def normalize(series):
            mn = min(series)
            mx = max(series)
            return (series - mn) / (mx - mn)

        df = x.apply(normalize)

        to_plot = []

        n = len(column_names)

        s = np.array([(np.cos(t), np.sin(t))
                      for t in [2.0 * np.pi * (i / float(n))
                                for i in range(n)]])

        for i in range(len(x)):
            row = df.irow(i).values
            row_ = np.repeat(np.expand_dims(row, axis=1), 2, axis=1)
            to_plot.append((s * row_).sum(axis=0) / row.sum())

        x_, y_ = zip(*to_plot)
        self.data['radviz_x'] = x_
        self.data['radviz_y'] = y_

        ax = self.scatter('radviz_x', 'radviz_y', **kwargs)

        ax.add_patch(patches.Circle((0.0, 0.0), radius=1.0, facecolor='none'))
        for xy, name in zip(s, column_names):
            ax.add_patch(patches.Circle(xy, radius=0.025, facecolor='gray'))
            if xy[0] < 0.0 and xy[1] < 0.0:
                ax.text(xy[0] - 0.025, xy[1] - 0.025, name,
                        ha='right', va='top', size='small')
            elif xy[0] < 0.0 and xy[1] >= 0.0:
                ax.text(xy[0] - 0.025, xy[1] + 0.025, name,
                        ha='right', va='bottom', size='small')
            elif xy[0] >= 0.0 and xy[1] < 0.0:
                ax.text(xy[0] + 0.025, xy[1] - 0.025, name,
                        ha='left', va='top', size='small')
            elif xy[0] >= 0.0 and xy[1] >= 0.0:
                ax.text(xy[0] + 0.025, xy[1] + 0.025, name,
                        ha='left', va='bottom', size='small')

        ax.axis('equal')
        return ax

    def _id_callback(self, event):
        # event.ind is the index into event's x and y data.
        #
        # event.artist.ind is the index of the entire artist into the original
        # dataframe.
        subset_df = event.artist.df.ix[event.artist.ind]
        for i in event.ind:
            _id = subset_df.index[i]
            yield _id

    def _default_callback(self, i):
        print self.data.ix[i]

    def strip_unknown_features(self):
        """
        Remove features not found in the `gffutils.FeatureDB`.  This will
        typically include 'ambiguous', 'no_feature', etc, but can also be
        useful if the database was created from a different one than was used
        to create the table.
        """
        if not self.db:
            return self
        ind = []
        for i, gene_id in enumerate(self.data.index):
            try:
                self.db[gene_id]
                ind.append(i)
            except gffutils.FeatureNotFoundError:
                pass
        ind = np.array(ind)
        return self.__class__(self.data.ix[ind], **self._kwargs)

    def genes_with_peak(self, peaks, transform_func=None, split=False,
                        intersect_kwargs=None, id_attribute='ID', *args,
                        **kwargs):
        """
        Returns a boolean index of genes that have a peak nearby.

        Parameters
        ----------
        peaks : string or pybedtools.BedTool
            If string, then assume it's a filename to a BED/GFF/GTF file of
            intervals; otherwise use the pybedtools.BedTool object directly.

        transform_func : callable
            This function will be applied to each gene object returned by
            self.features().  Additional args and kwargs are passed to
            `transform_func`. For example, if you're looking for peaks within
            1kb upstream of TSSs, then pybedtools.featurefuncs.TSS would be
            a useful `transform_func`, and you could supply additional kwargs
            of `upstream=1000` and `downstream=0`.

            This function can return iterables of features, too. For example,
            you might want to look for peaks falling within the exons of
            a gene.  In this case, `transform_func` should return an iterable
            of pybedtools.Interval objects.  The only requirement is that the
            `name` field of any feature matches the index of the dataframe.

        intersect_kwargs : dict
            kwargs passed to pybedtools.BedTool.intersect.

        id_attribute : str
            The attribute in the GTF or GFF file that contains the id of the
            gene. For meaningful results to be returned, a gene's ID be also
            found in the index of the dataframe.

            For GFF files, typically you'd use `id_attribute="ID"`.  For GTF
            files, you'd typically use `id_attribute="gene_id"`.
        """
        def _transform_func(x):
            """
            In order to support transform funcs that return a single feature or
            an iterable of features, we need to wrap it
            """
            result = transform_func(x)
            if isinstance(result, pybedtools.Interval):
                result = [result]
            for i in result:
                if i:
                    yield result

        intersect_kwargs = intersect_kwargs or {}
        if not self._cached_features:
            self._cached_features = pybedtools\
                .BedTool(self.features())\
                .saveas()

        if transform_func:
            if split:
                features = self._cached_features\
                    .split(_transform_func, *args, **kwargs)
            else:
                features = self._cached_features\
                    .each(transform_func, *args, **kwargs)

        else:
            features = self._cached_features

        hits = list(set([i[id_attribute] for i in features.intersect(
            peaks, **intersect_kwargs)]))
        return self.data.index.isin(hits)


class DifferentialExpressionResults(ResultsTable):

    __doc__ = _base_doc % dedent("""
    A ResultsTable subclass for working with differential expression results.

    Adds methods for up/down regulation, ma_plot, and sets class variables for
    which columns should be considered for pval, log fold change, and mean
    values. This class acts as a parent for subclasses like DESeqResults,
    EdgeRResults, and others/
    """)
    pval_column = 'padj'
    lfc_column = 'log2FoldChange'
    mean_column = 'baseMean'

    def __init__(self, data, db=None, header_check=True, **kwargs):
        import_kwargs = kwargs.pop('import_kwargs', {})
        if header_check and isinstance(data, basestring):
            comment_char = import_kwargs.get('comment', '#')
            for i, line in enumerate(open(data)):
                if line[0] != comment_char:
                    break
            import_kwargs['skiprows'] = i
        import_kwargs['na_values'] = ['nan']

        import_kwargs['index_col'] = import_kwargs.pop('index_col', 0)
        super(DifferentialExpressionResults, self).__init__(
            data=data, db=db, import_kwargs=import_kwargs, **kwargs)

    def changed(self, thresh=0.05, idx=True):
        """
        Changed features.

        {threshdoc}
        """
        ind = self.data[self.pval_column] <= thresh
        if idx:
            return ind
        return self[ind]

    def unchanged(self, thresh=0.05, idx=True):
        """
        Changed features.

        {threshdoc}
        """
        ind = (
            (self.data[self.pval_column] > thresh)
            | np.isnan(self.data[self.pval_column])
        )
        if idx:
            return ind
        return self[ind]

    def enriched(self, thresh=0.05, idx=True):
        """
        Enriched features.

        {threshdoc}
        """
        return self.upregulated(thresh=thresh, idx=idx)

    def upregulated(self, thresh=0.05, idx=True):
        """
        Upregulated features.

        {threshdoc}
        """
        ind = (
            (self.data[self.pval_column] <= thresh)
            & (self.data[self.lfc_column] > 0)
        )
        if idx:
            return ind
        return self[ind]

    def downregulated(self, thresh=0.05, idx=True):
        """
        Downregulated features.

        {threshdoc}
        """
        ind = (
            (self.data[self.pval_column] <= thresh)
            & (self.data[self.lfc_column] < 0)
        )
        if idx:
            return ind
        return self[ind]

    def disenriched(self, thresh=0.05, idx=True):
        """
        Disenriched features.

        {threshdoc}
        """
        return self.downregulated(thresh=thresh, idx=idx)

    def ma_plot(self, thresh, up_kwargs=None, dn_kwargs=None,
                zero_line=None, **kwargs):
        """
        MA plot

        Plots the average read count across treatments (x-axis) vs the log2
        fold change (y-axis).

        Additional kwargs are passed to self.scatter (useful ones might include
        `genes_to_highlight`)

        Parameters
        ----------
        thresh : float
            Features with values <= `thresh` will be highlighted in the plot.

        up_kwargs, dn_kwargs : None or dict
            Kwargs passed to matplotlib's scatter(), used for styling up/down
            regulated features (defined by `thresh` and `col`)

        zero_line : None or dict
            Kwargs passed to matplotlib.axhline(0).

        """
        genes_to_highlight = kwargs.pop('genes_to_highlight', [])
        genes_to_highlight.append(
            (self.upregulated(thresh),
             up_kwargs or dict(color='r')))
        genes_to_highlight.append(
            (self.downregulated(thresh),
             dn_kwargs or dict(color='b')))
        if zero_line is None:
            zero_line = {}
        x = self.mean_column
        y = self.lfc_column

        if 'xfunc' not in kwargs:
            kwargs['xfunc'] = np.log
        ax = self.scatter(
            x=x,
            y=y,
            genes_to_highlight=genes_to_highlight,
            **kwargs)
        if zero_line:
            ax.axhline(0, **zero_line)
        return ax

    threshdoc = """
    Parameters
    ----------
    thresh : float
        Only features with <= `thresh` will be returned

    idx : bool
        If True, a boolean index will be returned.  If False, a new object will
        be returned that has been subsetted.
    """
    enriched.__doc__ = enriched.__doc__.format(threshdoc=threshdoc)
    disenriched.__doc__ = disenriched.__doc__.format(threshdoc=threshdoc)
    upregulated.__doc__ = upregulated.__doc__.format(threshdoc=threshdoc)
    downregulated.__doc__ = downregulated.__doc__.format(threshdoc=threshdoc)


class EdgeRResults(DifferentialExpressionResults):
    __doc__ = _base_doc % dedent(
        """
        Class for working with results from edgeR.

        Just like a DifferentialExpressionResults object, but sets the
        pval_column, lfc_column, and mean_column to the names used in edgeR's
        output.
        """)
    pval_column = 'FDR'
    lfc_column = 'logFC'
    mean_column = 'logCPM'


class DESeqResults(DifferentialExpressionResults):
    __doc__ = _base_doc % dedent(
        """
        Class for working with results from DESeq.

        Just like a DifferentialExpressionResults object, but sets the
        pval_column, lfc_column, and mean_column to the names used in edgeR's
        output.
        """)
    def colormapped_bedfile(self, genome, cmap=None):
        """
        Create a BED file with padj encoded as color

        Features will be colored according to adjusted pval (phred
        transformed).  Downregulated features have the sign flipped.

        Parameters
        ----------
        cmap : matplotlib colormap
            Default is matplotlib.cm.RdBu_r

        Notes
        -----
        Requires a FeatureDB to be attached.
        """
        if self.db is None:
            raise ValueError("FeatureDB required")
        db = gffutils.FeatureDB(self.db)

        def scored_feature_generator(d):
            for i in range(len(d)):
                try:
                    feature = db[d.ix[i]]
                except gffutils.FeatureNotFoundError:
                    raise gffutils.FeatureNotFoundError(d.ix[i])
                score = -10 * np.log10(d.padj[i])
                lfc = d.log2FoldChange[i]
                if np.isnan(lfc):
                    score = 0
                if lfc < 0:
                    score *= -1
                feature.score = str(score)
                feature = extend_fields(
                    gff2bed(gffutils.helpers.asinterval(feature)), 9)
                fields = feature.fields[:]
                fields[6] = fields[1]
                fields[7] = fields[2]
                fields.append(str(d.padj[i]))
                fields.append(str(d.pval[i]))
                fields.append('%.3f' % d.log2FoldChange[i])
                fields.append('%.3f' % d.baseMeanB[i])
                fields.append('%.3f' % d.baseMeanB[i])
                yield pybedtools.create_interval_from_list(fields)

        x = pybedtools.BedTool(scored_feature_generator(self)).saveas()
        norm = x.colormap_normalize()
        if cmap is None:
            cmap = cm.RdBu_r
        cmap = colormap_adjust.cmap_center_point_adjust(
            cmap, [norm.vmin, norm.vmax], 0)

        def score_zeroer(f):
            f.score = '0'
            return f
        return x.each(add_color, cmap=cmap, norm=norm)\
                .sort()\
                .each(score_zeroer)\
                .truncate_to_chrom(genome)\
                .saveas()

    def autosql_file(self):
        """
        Generate the autosql for DESeq results (to create bigBed)

        Returns a temp filename containing the autosql defining the extra
        fields.

        This for creating bigBed files from BED files created by
        colormapped_bed.  When a user clicks on a feature, the DESeq results
        will be reported.
        """
        fn = pybedtools.BedTool._tmp()

        AUTOSQL = dedent(
            """
            table example
            "output from DESeq"
            (
            string  chrom;  "chromosome"
            uint chromStart; "start coord"
            uint chromEnd; "stop coord"
            string name; "name of feature"
            uint score; "always zero"
            char[1] strand; "+ or - for strand"
            uint    thickStart; "Coding region start"
            uint    thickEnd;  "Coding region end"
            uint reserved; "color according to score"
            string padj; "DESeq adjusted p value"
            string pval; "DESeq raw p value"
            string logfoldchange; "DESeq log2 fold change"
            string basemeana; "DESeq baseMeanA"
            string basemeanb; "DESeq baseMeanB"
        )
        """)

        fout = open(fn, 'w')
        fout.write(AUTOSQL)
        fout.close()
        return fn


class DESeq2Results(DESeqResults):
    __doc__ = _base_doc % dedent(
        """
        Class for working with results from DESeq2.

        Just like a DifferentialExpressionResults object, but sets the
        pval_column, lfc_column, and mean_column to the names used in edgeR's
        output.
        """)
    pval_column = 'padj'
    lfc_column = 'log2FoldChange'
    mean_column = 'baseMean'


class LazyDict(object):
    def __init__(self, fn_dict, dbfn, index_file, extra=None, cls=DESeqResults,
                 modifier=None):
        """
        Dictionary-like object that lazily-loads ResultsTable objects.

        Parameters
        ----------
        fn_dict : dict
            Keys of `fn_dict` will be the keys of this LazyDict object.  Values
            should be filenames which will be loaded into ResultsTable object
            upon access for the first time.

        index_file : str
            Path to a file that contains one ID per line.  This file is used to
            ensure all ResultsTable objects are aligned to the same index.

        dbfn : str
            Filename to a gffutils database.  This enables gene info to be
            attached to the dataframe.

        extra : pandas.dataframe
            This dataframe hat will be merged into the data in each file.  This
            is useful for attaching things like gene lengths, alt names, etc.
            In order for it to work, this dataframe must be indexed the same
            way the ResultsTable files are indexed.

        cls : ResultsTable class or subclass
            Each filename in `fn_dict` will be converted using this class.

        modifier : callable
            Upon first access, each newly-constructed ResultsTable will first
            have the `extra` data attached, and then will be provided as this
            callable's only argument.  The callable can make any modifications
            to the ResultsTable, and return a new version that will be used in
            the future when the same key is accessed.  For example, exonic bp
            data can be provided as part of the `extra` dataframe, and then the
            `modifier` can be a function that adds an RPKM column.

        Notes
        -----
        When a key is provided for the first time, the workflow is
        ResultsTable(fn, **kwargs) -> attach `extra` -> send to `modifier` ->
        return extended and modified ResultsTable.  Subsequent access of the
        same key will immediately return the extended-and-modified
        ResultsTable.

        """
        self.fn_dict = fn_dict
        self._dict = {}
        self.dbfn = dbfn
        self.index = [i.strip() for i in open(index_file)]
        if extra is not None:
            self.extra = extra.ix[self.index]
        else:
            self.extra = extra
        self.modifier = modifier
        self._cls = cls

    def __getitem__(self, key):
        if key not in self._dict:
            fn = self.fn_dict[key]
            obj = self._cls(fn, db=self.dbfn)
            obj.data = obj.data.ix[self.index]
            if self.extra is not None:
                obj.data = pandas.merge(obj.data, self.extra, left_index=True,
                                        right_index=True)
            if self.modifier:
                obj = self.modifier(obj)

            self._dict[key] = obj

        return self._dict[key]

    def __repr__(self):
        s = "<%s> with possible keys\n:%s\n" \
            % (self.__class__.__name__, self.fn_dict.keys())
        s += "and existing keys:\n"
        s += repr(self._dict)
        return s

    def keys(self):
        return self.fn_dict.keys()

    def values(self):
        return [self._dict[key] for key in self.keys()]

    def items(self):
        return list((key, self._dict[key]) for key in self.keys())
