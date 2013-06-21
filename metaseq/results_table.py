import pybedtools
from textwrap import dedent
import numpy as np
import pandas
import gffutils
from gffutils.helpers import asinterval
from matplotlib import pyplot as plt


class ResultsTable(object):
    def __init__(self, data, db=None, import_kwargs=None):
        """
        Wrapper around a pandas.DataFrame that adds additional functionality.

        The underlying pandas.DataFrame is always available with the `data`
        attribute.

        Any attributes not explicitly in ResultsTable will be looked for in the
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
        return self.__class__(self.data.__getitem__(attr), **self._kwargs)

    def update(self, dataframe):
        """
        Updates the current data with a new dataframe.

        This extra step is required to get around the fancy pandas.DataFrame
        indexing (like .ix, .iloc, etc).
        """
        return self.__class__(dataframe, **self._kwargs)

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
                yield asinterval(self.db[i])
            except gffutils.FeatureNotFoundError:
                if ignore_unknown:
                    continue
                else:
                    raise gffutils.FeatureNotFoundError('%s not found' % i)

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

    def scatter(self, x, y, xfunc=None, yfunc=None, xscale=None, yscale=None,
                xlab=None, ylab=None, genes_to_highlight=None,
                label_genes=False,
                general_kwargs=dict(color="k", alpha=0.2, linewidths=0),
                marginal=False, offset_kwargs={}, label_kwargs=None, ax=None,
                one_to_one=None, callback=None, xlab_prefix=None,
                ylab_prefix=None, sizefunc=None):
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

        marginal : bool
            Toggles display of non-finite cases where `x` or `y` is
            nonzero but the other one is zero (and `xfunc` or `yfunc` are log)

        callback : callable
            Function to call upon clicking a point. Default is to print the
            gene name, but an example of another useful callback would be
            a mini-browser connected to a genomic_signal object from which the
            expression data were calculated.

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
        """
        _x = self.data[x]
        _y = self.data[y]

        # Construct defaults---------------------------------------------------
        def identity(x):
            return x

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
                                bbox=dict(facecolor='w', edgecolor='None',
                                          alpha=0.5))

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

        # Create an index of the 'non-interesting' features -- at least the
        # ones that are not highlighted in genes_to_highlight.
        allind = np.zeros_like(xi) == 0
        for i in genes_to_highlight:
            this_ind = i[0]
            if this_ind.dtype == 'bool':
                this_ind = np.where(this_ind)[0]
            allind[this_ind] = False

        # Plot everybody
        coll = ax.scatter(xi[allind], yi[allind], picker=5, **general_kwargs)
        coll.subdata = self
        coll.subind = allind

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
            coll = ax.scatter(xi[ind], yi[ind], picker=5, **updated_kwargs)
            coll.subdata = self
            coll.subind = ind

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
        if callback is not None:
            ax.figure.canvas.mpl_connect('pick_event', callback)

        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)

        return ax

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

    def genes_with_peak(self, peaks, transform_func=None,
                        intersect_kwargs=None, *args, **kwargs):
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
            a useful `transform_func`, and supply additional kwargs of
            `upstream=1000` and `downstream=0`.

        intersect_kwargs : dict
            kwargs passed to pybedtools.BedTool.intersect.
        """
        intersect_kwargs = intersect_kwargs or {}
        if not self._cached_features:
            self._cached_features = pybedtools\
                .BedTool(self.features())\
                .saveas()

        if transform_func:
            features = self._cached_features\
                .each(transform_func, *args, **kwargs)

        else:
            features = self._cached_features

        hits = [i.name for i in features.intersect(peaks, **intersect_kwargs)]
        return self.data.index.isin(hits)


class DESeqResults(ResultsTable):
    def __init__(self, data, db=None, header_check=True,
                 remove_deseq_extra=True, **kwargs):

        to_remove = [
            'not_aligned',
            'no_feature',
            'ambiguous',
            'alignment_not_unique',
            'too_low_aQual']

        import_kwargs = kwargs.pop('import_kwargs', {})
        if header_check and isinstance(data, basestring):
            comment_char = import_kwargs.get('comment', '#')
            for i, line in enumerate(open(data)):
                if line[0] != comment_char:
                    break
            import_kwargs['skiprows'] = i
        import_kwargs['na_values'] = ['nan']

        import_kwargs['index_col'] = import_kwargs.pop('index_col', 0)
        super(DESeqResults, self).__init__(
            data=data, db=db, import_kwargs=import_kwargs, **kwargs)
        self.data.log2FoldChange = self.data.log2FoldChange.astype(float)
        self.data.foldChange = self.data.foldChange.astype(float)
        if remove_deseq_extra:
            self.data = self.data[~self.data.index.isin(to_remove)]

    def changed(self, thresh=0.05, idx=True, col='padj'):
        """
        Changed features.

        {threshdoc}
        """
        ind = self.data[col] <= thresh
        if idx:
            return ind
        return self[ind]

    def unchanged(self, thresh=0.05, idx=True, col='padj'):
        """
        Changed features.

        {threshdoc}
        """
        ind = (self.data[col] > thresh) | np.isnan(self.data[col])
        if idx:
            return ind
        return self[ind]

    def enriched(self, thresh=0.05, idx=True, col='padj'):
        """
        Enriched features.

        {threshdoc}
        """
        return self.upregulated(thresh=thresh, idx=idx)

    def upregulated(self, thresh=0.05, idx=True, col='padj'):
        """
        Upregulated features.

        {threshdoc}
        """
        ind = (self.data[col] <= thresh) & (self.data['log2FoldChange'] > 0)
        if idx:
            return ind
        return self[ind]

    def downregulated(self, thresh=0.05, idx=True, col='padj'):
        """
        Downregulated features.

        {threshdoc}
        """
        ind = (self.data[col] <= thresh) & (self.data['log2FoldChange'] < 0)
        if idx:
            return ind
        return self[ind]

    def disenriched(self, thresh=0.05, idx=True, col='padj'):
        """
        Disenriched features.

        {threshdoc}
        """
        return self.downregulated(thresh=thresh, idx=idx)

    def _default_callback(self, event):
        for i in event.ind:
            print
            print event.artist.subdata.ix[event.artist.subind].ix[i]

    def ma_plot(self, thresh, col='padj', up_kwargs=None, dn_kwargs=None,
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

        col : str
            Column name that `thresh` will be applied to

        up_kwargs, dn_kwargs : None or dict
            Kwargs passed to matplotlib's scatter(), used for styling up/down
            regulated features (defined by `thresh` and `col`)

        zero_line : None or dict
            Kwargs passed to matplotlib.axhline(0).

        """
        genes_to_highlight = kwargs.pop('genes_to_highlight', [])
        genes_to_highlight.append(
            (self.upregulated(thresh, col=col),
             up_kwargs or dict(color='r')))
        genes_to_highlight.append(
            (self.downregulated(thresh, col=col),
             dn_kwargs or dict(color='b')))
        if zero_line is None:
            zero_line = {}
        x = self.baseMean
        y = self.log2FoldChange

        if 'xfunc' not in kwargs:
            xfunc = np.log1p
        x = xfunc(x)
        ax = self.scatter(
            x=x,
            y=y,
            genes_to_highlight=genes_to_highlight,
            **kwargs)
        if zero_line:
            ax.axhline(0, **zero_line)
        return ax

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
                feature = extend_fields(gff2bed(asinterval(feature)), 9)
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

    threshdoc = """
    Parameters
    ----------
    thresh : float
        Only features with <= `thresh` will be returned

    idx : bool
        If True, a boolean index will be returned.  If False, a new object will
        be returned that has been subsetted.

    col : str
        Name of the column to apply `thresh` to
    """
    enriched.__doc__ = enriched.__doc__.format(threshdoc=threshdoc)
    disenriched.__doc__ = disenriched.__doc__.format(threshdoc=threshdoc)
    upregulated.__doc__ = upregulated.__doc__.format(threshdoc=threshdoc)
    downregulated.__doc__ = downregulated.__doc__.format(threshdoc=threshdoc)


if __name__ == "__main__":
    import metaseq
    from matplotlib import pyplot as plt

    db = metaseq.example_filename('dmel-all-r5.33-cleaned.gff.db')
    import_kwargs = dict(comment='#')
    d = DESeqResults(
        metaseq.example_filename('rrp6-s2-polyA.final.summary'),
        db=db,
        import_kwargs=import_kwargs,
    )

    e = DESeqResults(
        metaseq.example_filename('rrp40-s2-polyA.final.summary'),
        db=db,
        import_kwargs=import_kwargs,
    )

    d = d.align_with(e)

    en = e.enriched()

    def sizefunc(x):
        return np.abs(x.log2FoldChange * 2)
    general_kwargs = dict(color='k', s=e.log2FoldChange ** 2)
    e.scatter(
        e.baseMeanA,
        e.baseMeanB,
        genes_to_highlight=[
            (en, dict(color='r', alpha=1.0)),
        ],
        general_kwargs=general_kwargs,
        xfunc=np.log1p,
        yfunc=np.log1p,
        one_to_one=dict(color='b', linestyle=':'),
        callback=e._default_callback,
    )

    e.ma_plot(
        0.05,
        callback=e._default_callback,
        zero_line=dict(color='y', linestyle=':'),
        general_kwargs=general_kwargs)
    plt.show()
