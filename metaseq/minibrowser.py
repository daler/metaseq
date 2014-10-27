"""
Module for spawning mini genome browsers using a plugin structure, making it
possible to build rather complex mini-browsers.  The goal is to point the
mini-browser to some data, and call its plot() method with a feature.  This
will spawn a new figure showing the data for that interval.

MiniBrowser classes are just a general way of mapping data-manipulation or
data-visualization methods to an Axes on which the data should be displayed.

To make a new subclass:

1. Create one or more methods that accept an Axes object and a pybedtools
   Interval object and return a feature.  The simplest do-nothing method would
   be::

        def my_panel(self, ax, feature)
            return feature

   A more useful method might be one that plots genomic signal over the
   region::

        def my_panel(self, ax, feature):
            # for simplicity, assume just use the first genomic_signal
            gs = self.genomic_signal_objs[0]
            x, y = gs.local_coverage(feature, bins=100)
            ax.plot(x, y, **kwargs)
            ax.axis('tight')
            return feature

2. Then, override the `panels()` method.  This method:

    * Creates Axes as needed; assumes that self.make_fig() has already been
      called so that self.fig is available.
    * Returns a list of (ax, method) tuples.  This list maps created Axes to
      methods that should operate on them (like `my_panel` method above).

    For example::

        def panels(self):
            ax = self.fig.add_subplot(111)
            return [(ax, self.my_panel)]


A figure is spawned by calling the `plot` method on a pybedtools genomic
interval, e.g.::

    s = SignalMiniBrowser(ip, control])
    s.plot(feature)

"""
from copy import copy
from matplotlib import pyplot as plt
from pybedtools.contrib.plotting import Track
import pybedtools
import gffutils
from gffutils.helpers import asinterval
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FormatStrFormatter
import _genomic_signal


class BaseMiniBrowser(object):
    """
    Base class for plotting a genomic region.

    This class is designed to be sub-classed, so it really just a shell with
    some example methods filled in.
    """
    _all_figures = []

    def __init__(self, genomic_signal_objs):
        self.genomic_signal_objs = genomic_signal_objs

    def close_all(self):
        """
        Close all figures spawned by this class.
        """
        for f in self._all_figures:
            plt.close(f)

    def plot(self, feature):
        """
        Spawns a new figure showing data for `feature`.

        :param feature: A `pybedtools.Interval` object

        Using the pybedtools.Interval `feature`, creates figure specified in
        :meth:`BaseMiniBrowser.make_fig` and plots data on panels according to
        `self.panels()`.
        """
        if isinstance(feature, gffutils.Feature):
            feature = asinterval(feature)
        self.make_fig()
        axes = []
        for ax, method in self.panels():
            feature = method(ax, feature)
            axes.append(ax)
        return axes

    def make_fig(self):
        """
        Figure constructor, called before `self.plot()`
        """
        self.fig = plt.figure(figsize=(8, 4))
        self._all_figures.append(self.fig)

    def panels(self):
        """
        Method to be overriden by subclasses.

        This method should create Axes objects on self.fig and register which
        methods operate on which axes.

        Each method's signature should be func(self, ax, feature), and it will
        take a pybedtools.Interval `feature` and plot it or something derived
        from it onto `ax`. It should return a pybedtools.Interval object,
        usually just the same one that was passed in.
        """
        ax = self.fig.add_subplot(111)
        return [(ax, self.example_panel)]

    def example_panel(self, ax, feature):
        """
        A example panel that just prints the text of the feature.
        """
        txt = '%s:%s-%s' % (feature.chrom, feature.start, feature.stop)
        ax.text(0.5, 0.5, txt, transform=ax.transAxes)
        return feature


class ChIPSeqMiniBrowser(BaseMiniBrowser):
    def __init__(self, ip, control, db=None, peaks=None,
                 local_coverage_kwargs=dict(stranded=False), ip_style=None,
                 control_style=None, peaks_style=None, max_bins=500):
        """
        Mini-browser to show IP and control signal, along with optional peaks
        and gene models.  See module-level documentation for details.

        Parameters
        ----------
        ip, control : genomic_signal object
            These can be any kind of genomic signal object.  If BamSignal, then
            the signal will be scaled by million mapped reads.

        db : filename or gffutils.FeatureDB
            Optional database of annotations.  If provided, gene models will be
            plotted on an additional axes.

        peaks : filename or pybedtools.BedTool object
            Optional file of called peaks.  If provided, will be plotted on an
            additional axes.

        local_coverage_kwargs : dict
            Kwargs to pass on to the local_coverage method of `ip` and
            `control`.  It is recommended that at least `stranded=False` is
            included, since the mini-browser always displays the plus strand.

        ip_style, control_style : dict
            Kwargs to pass to Axes.fill_between (e.g. colors, alpha, line
            styles)

        peaks_style : dict
            Kwargs to pass to PolyCollection (e.g., face_color, edgewidth)

        max_bins : int
            Maximum number of bins to use when getting genomic_signal.  This
            can be overridden by providing the `bins` kwarg in
            `local_coverage_kwargs`.

        Notes
        -----
        After creation, settings can be changed in between calls to the
        `plot()` method.

        For example, the self.settings dictionary, which
        contains arguments passed to the gffutils.contrib.plotting.Gene class,
        can be modified to specify which kinds of features should be plotted.
        The `gs` attribute contains the matplotlib.gridspec.GridSpec object,
        which can be used to modify the relative sizes of the subplots.


        """
        super(ChIPSeqMiniBrowser, self).__init__([ip, control])
        self.local_coverage_kwargs = local_coverage_kwargs or {}
        self.ip = ip
        self.control = control
        if isinstance(db, basestring):
            db = gffutils.FeatureDB(db)
        self.db = db
        self.ip_style = ip_style or {}
        self.peaks_style = peaks_style or {}
        self.control_style = control_style or {}

        if peaks is not None:
            peaks = pybedtools.BedTool(peaks)
        self.peaks = peaks
        self.max_bins = max_bins

        self.settings = {
            'transcripts': None,
            'cds': ['CDS'],
            'utrs': ['exon'],
            'color': '0.5',
            'featuretype': 'gene',
        }

    def panels(self):
        self.fig.set_facecolor('w')
        method_dispatch = {
            'ip': self.ip_panel,
            'control': self.control_panel,
            'peaks': self.peaks_panel,
            'genes': self.genes_panel}
        if self.db and self.peaks:
            gs = gridspec.GridSpec(4, 1, height_ratios=[1, 1, .15, .5])
            ip_ax = plt.subplot(gs[0])
            control_ax = plt.subplot(gs[1], sharex=ip_ax, sharey=ip_ax)
            peaks_ax = plt.subplot(gs[2], sharex=ip_ax)
            gene_ax = plt.subplot(gs[3], sharex=ip_ax)
            axes = {'ip': ip_ax,
                    'control': control_ax,
                    'peaks': peaks_ax,
                    'genes': gene_ax,
                    }

        elif self.db and self.peaks is None:
            gs = gridspec.GridSpec(3, 1, height_ratios=[1, 1, .5])
            ip_ax = plt.subplot(gs[0])
            control_ax = plt.subplot(gs[1], sharex=ip_ax, sharey=ip_ax)
            gene_ax = plt.subplot(gs[2], sharex=ip_ax)
            axes = {'ip': ip_ax,
                    'control': control_ax,
                    'genes': gene_ax,
                    }

        elif self.db is None and self.peaks:
            gs = gridspec.GridSpec(3, 1, height_ratios=[1, 1, .3])
            ip_ax = plt.subplot(gs[0])
            control_ax = plt.subplot(gs[1], sharex=ip_ax, sharey=ip_ax)
            peaks_ax = plt.subplot(gs[2], sharex=ip_ax)
            axes = {'ip': ip_ax,
                    'control': control_ax,
                    'peaks': peaks_ax,
                    }

        elif self.db is None and self.peaks is None:
            gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1])
            ip_ax = plt.subplot(gs[0])
            control_ax = plt.subplot(gs[1], sharex=ip_ax, sharey=ip_ax)
            axes = {'ip': ip_ax,
                    'control': control_ax,
                    }

        self.axes = axes
        self.gridspec = gs

        # Order shouldn't matter, but just in case...
        mapping = []
        for kind in ['ip', 'control', 'peaks', 'genes']:
            if kind in axes:
                mapping.append((axes[kind], method_dispatch[kind]))

        self._first_ax = mapping[0][0]
        self._last_ax = mapping[-1][0]

        return mapping

    def _bins(self, feature):
        return min(len(feature), self.max_bins)

    def _zoomed_feature(self, feature):
        extra = int(self.settings['zoom'] * len(feature) / 2)
        feature = copy(feature)
        feature.start -= extra
        feature.stop += extra
        return feature

    def ip_panel(self, ax, feature):
        bins = self._bins(feature)
        x, y = self.ip.local_coverage(
            feature, bins=bins, **self.local_coverage_kwargs)
        if isinstance(self.ip, _genomic_signal.BamSignal):
            y /= (self.ip.mapped_read_count() / 1e6)
        ax.fill_between(x, y, y2=0, **self.ip_style)
        ax.axis('tight')
        ax.spines['right'].set_color('None')
        ax.spines['top'].set_color('None')
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        return feature

    def control_panel(self, ax, feature):
        bins = self._bins(feature)
        x, y = self.control.local_coverage(
            feature, bins=bins, **self.local_coverage_kwargs)
        if isinstance(self.control, _genomic_signal.BamSignal):
            y /= (self.control.mapped_read_count() / 1e6)
        ax.fill_between(x, y, y2=0, **self.control_style)
        ax.axis('tight')
        ax.spines['right'].set_color('None')
        ax.spines['top'].set_color('None')
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        return feature

    def peaks_panel(self, ax, feature):
        hits = self.peaks.intersect([feature], u=True)
        track = Track(hits, **self.peaks_style)
        ax.add_collection(track)
        for x in (
            ax.get_yticklines() + ax.get_xticklines() + ax.get_yticklabels()
            + ax.get_xticklabels()
        ):
            x.set_visible(False)
        ax.set_ylabel('Peaks', rotation=0, horizontalalignment='right',
                      verticalalignment='center')
        ax.yaxis.set_ticks_position('none')
        ax.xaxis.set_ticks_position('none')
        ax.set_frame_on(False)

        return feature

    def plot(self, feature):
        self._current_feature = feature
        axes = super(ChIPSeqMiniBrowser, self).plot(feature)
        for ax in self.fig.axes:
            if ax is not self._last_ax:
                for txt in ax.get_xticklabels():
                    txt.set_visible(False)
        self.fig.subplots_adjust(wspace=0.05)
        self._first_ax.axis(xmin=feature.start, xmax=feature.stop)
        self._last_ax.set_xlabel(feature.chrom)
        self._last_ax.xaxis.set_major_formatter(FormatStrFormatter("%d"))
        return axes

    def coords(self):
        chrom = self._current_feature.chrom
        start, stop, _, _ = self._first_ax.axis()
        return '%s:%s-%s' % (chrom, int(start), int(stop))

    def genes_panel(self, ax, feature):
        """
        Plots gene models on an Axes.

        Queries the database

        :param ax: matplotlib.Axes object
        :param feature: pybedtools.Interval

        """
        from gffutils.contrib.plotting import Gene
        extent = [feature.start, feature.stop]
        nearby_genes = self.db.region(
            (feature.chrom, feature.start, feature.stop), featuretype='gene')
        ybase = 0
        ngenes = 0
        for nearby_gene in nearby_genes:
            # TODO: there should be a better way of identifying which gene is
            # the same as the feature requested.  Might need to expose an "ID"
            # kwarg.
            try:
                if nearby_gene['ID'][0] == feature['ID']:
                    color = '0.2'
                else:
                    color = '0.5'
            except KeyError:
                color = '0.5'
            ngenes += 1
            extent.extend([nearby_gene.start, nearby_gene.stop])
            gene_collection = Gene(
                self.db,
                nearby_gene,
                transcripts=None,
                cds=['CDS'],
                utrs=['exon'],
                ybase=ybase,
                color=color)
            gene_collection.name = nearby_gene.id
            gene_collection.add_to_ax(ax)
            ybase += gene_collection.max_y

        xmin = min(extent)
        xmax = max(extent)
        ymax = ngenes

        # 1% padding seems to work well
        padding = (xmax - xmin) * 0.01
        ax.axis('tight')

        # add lines indicating extent of current feature
        # vline_kwargs = dict(color='k', linestyle='--')
        # ax.axvline(feature.start, **vline_kwargs)
        # ax.axvline(feature.stop, **vline_kwargs)

        # Make a new feature to represent the region plus surrounding genes
        interval = pybedtools.create_interval_from_list(feature.fields)
        interval.strand = '.'
        for txt in ax.get_yticklabels():
            txt.set_visible(False)
        for tick in ax.get_yticklines():
            tick.set_visible(False)
        ax.set_ylabel('Genes')
        ax.spines['right'].set_color('None')
        ax.spines['left'].set_color('None')
        ax.spines['top'].set_color('None')
        ax.yaxis.set_ticks_position('none')
        ax.xaxis.set_ticks_position('bottom')
        ax.set_ylabel('Genes', rotation=0, horizontalalignment='right',
                      verticalalignment='center')

        return interval


class SignalMiniBrowser(BaseMiniBrowser):
    def __init__(self, genomic_signal_objs, local_coverage_kwargs=None,
                 plotting_kwargs=None):
        """
        Base class for plotting genomic signal.

        Plots genomic signal over a particular area of the genome.  Designed to
        be extended by subclasses, but can stand alone if all you want is
        a simple one-panel browser.

        :param genomic_signal_objs: list of genomic signal objects (e.g.,
            :class:`metaseq.genomic_signal.BamSignal` instances).
        :param local_coverage_kwargs: a dictionary of kwargs to send to each
            genomic signals' `local_coverage()` method, e.g.::

                local_coverage_kwargs = dict(fragment_size=200).
        :param plotting_kwargs:  a list of dictionaries, one for each genomic
            signals object, e.g,::

                plotting_kwargs = [dict(color='r', label='IP', dict(color='k',
                    label='input')]

        """
        super(SignalMiniBrowser, self).__init__(genomic_signal_objs)
        self.plotting_kwargs = \
            plotting_kwargs or [{} for i in genomic_signal_objs]
        self.local_coverage_kwargs = local_coverage_kwargs or {}

    def panels(self):
        """
        Add a single panel to the figure
        """
        ax = self.fig.add_subplot(111)
        return [(ax, self.signal_panel)]

    def signal_panel(self, ax, feature):
        """
        Plots each genomic signal as a line using the corresponding
        plotting_kwargs
        """
        for gs, kwargs in zip(self.genomic_signal_objs, self.plotting_kwargs):
            x, y = gs.local_coverage(feature, **self.local_coverage_kwargs)
            ax.plot(x, y, **kwargs)
        ax.axis('tight')
        return feature


class GeneModelMiniBrowser(SignalMiniBrowser):
    def __init__(self, genomic_signal_objs, db, **kwargs):
        """
        Mini-browser to show a signal panel on top and gene models on the
        bottom.

        :param genomic_signal_objs: a list of genomic_signals objects
        :param db: a `gffutils.FeatureDB`
        """
        super(GeneModelMiniBrowser, self).__init__(
            genomic_signal_objs, **kwargs)
        self.db = db

    def panels(self):
        """
        Add 2 panels to the figure, top for signal and bottom for gene models
        """
        ax1 = self.fig.add_subplot(211)
        ax2 = self.fig.add_subplot(212, sharex=ax1)
        return (ax2, self.gene_panel), (ax1, self.signal_panel)

    def gene_panel(self, ax, feature):
        """
        Plots gene models on an Axes.

        Queries the database

        :param ax: matplotlib.Axes object
        :param feature: pybedtools.Interval

        """
        from gffutils.contrib.plotting import Gene
        extent = [feature.start, feature.stop]
        nearby_genes = self.db.region(
            (feature.chrom, feature.start, feature.stop), featuretype='gene')
        ybase = 0
        ngenes = 0
        for nearby_gene in nearby_genes:
            ngenes += 1
            extent.extend([nearby_gene.start, nearby_gene.stop])
            gene_collection = Gene(
                self.db,
                nearby_gene,
                transcripts=['mRNA'],
                cds=['CDS'],
                utrs=['exon'],
                ybase=ybase,
                color="0.5", picker=5)
            gene_collection.name = nearby_gene.id
            gene_collection.add_to_ax(ax)
            ybase += gene_collection.max_y

        xmin = min(extent)
        xmax = max(extent)
        ymax = ngenes

        # 1% padding seems to work well
        padding = (xmax - xmin) * 0.01
        ax.axis('tight')

        # add lines indicating extent of current feature
        vline_kwargs = dict(color='k', linestyle='--')
        ax.axvline(feature.start, **vline_kwargs)
        ax.axvline(feature.stop, **vline_kwargs)

        # Make a new feature to represent the region plus surrounding genes
        interval = pybedtools.create_interval_from_list(feature.fields)
        interval.start = xmin - padding
        interval.stop = xmax + padding
        interval.strand = '.'
        return interval


class PeakMiniBrowser(SignalMiniBrowser):
    def __init__(self, genomic_signal_objs, bed, **kwargs):
        """
        Signal on the top panel, peaks in the bottom panel. `bed` is a filename
        or BedTool object (or list of these things) that will be used to make
        a pybedtools.contrib.plotting.Track.

        If genomic_signal_objs is None, then only show the peaks axes
        """
        super(PeakMiniBrowser, self).__init__(genomic_signal_objs, **kwargs)
        self.bed = bed

    def panels(self):
        ax1 = self.fig.add_subplot(211)
        ax2 = self.fig.add_subplot(212, sharex=ax1)
        return (ax1, self.signal_panel), (ax2, self.peak_panel)

    def peak_panel(self, ax, feature):
        bedtool = pybedtools.BedTool(self.bed)
        features = bedtool.intersect([feature], u=True)
        track = Track(features)
        ax.add_collection(track)
        # ax.axis('tight')
        return feature


if __name__ == "__main__":
    import metaseq
    import gffutils
    import pybedtools

    G = gffutils.FeatureDB(
        metaseq.example_filename('Homo_sapiens.GRCh37.66.cleaned.gtf.db'))

    ip = metaseq.genomic_signal(
        metaseq.example_filename('wgEncodeUwTfbsK562CtcfStdAlnRep1.bam'),
        'bam')
    inp = metaseq.genomic_signal(
        metaseq.example_filename('wgEncodeUwTfbsK562InputStdAlnRep1.bam'),
        'bam')
    peaks = pybedtools.BedTool(metaseq.example_filename(
        'wgEncodeUwTfbsK562CtcfStdPkRep1.narrowPeak.gz'))

    plotting_kwargs = [
        dict(color='r', label='IP'),
        dict(color='k', linestyle=':', label='input')]

    local_coverage_kwargs = dict(fragment_size=200)

    b = SignalMiniBrowser(
        [ip, inp],
        plotting_kwargs=plotting_kwargs,
        local_coverage_kwargs=local_coverage_kwargs)

    g = GeneModelMiniBrowser(
        [ip, inp],
        G,
        plotting_kwargs=plotting_kwargs,
        local_coverage_kwargs=local_coverage_kwargs)

    p = PeakMiniBrowser(
        [ip, inp],
        peaks,
        plotting_kwargs=plotting_kwargs,
        local_coverage_kwargs=local_coverage_kwargs)

    feature = peaks[3]

    b.plot(feature)
    g.plot(feature)
    p.plot(feature)
    plt.show()
