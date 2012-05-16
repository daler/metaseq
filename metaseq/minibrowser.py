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

from matplotlib import pyplot as plt
from pybedtools.contrib.plotting import Track
import pybedtools


class BaseMiniBrowser(object):
    """
    Base class for plotting genomic signal.

    This class is designed to be sub-classed, so it really just a shell with
    some example methods filled in.
    """
    def __init__(self, genomic_signal_objs):
        self.genomic_signal_objs = genomic_signal_objs

    def plot(self, feature):
        """
        Spawns a new figure showing data for `feature`.

        :param feature: A `pybedtools.Interval` object

        Using the pybedtools.Interval `feature`, creates figure specified in
        :meth:`BaseMiniBrowser.make_fig` and plots data on panels according to
        `self.panels()`.
        """
        self.make_fig()
        for ax, method in self.panels():
            feature = method(ax, feature)

    def make_fig(self):
        """
        Figure constructor, called before `self.plot()`
        """
        self.fig = plt.figure(figsize=(8, 4))

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


class SignalMiniBrowser(BaseMiniBrowser):
    def __init__(self, genomic_signal_objs, local_coverage_kwargs=None,
            plotting_kwargs=None):
        """
        Base class for plotting genomic signal over a particular area of the
        genome, designed to be extended by subclasses, but can stand alone if
        all you want is a simple one-panel browser.

        `genomic_signal_objs` is a list of genomic signal objects (e.g.,
        BamSignal or BigWigSignal).

        `local_coverage_kwargs` is a dictionary of kwargs to send to each
        genomic signals' local_coverage method, e.g.::

            local_coverage_kwargs = dict(fragment_size=200).

        `plotting_kwargs` is a list of dictionaries, one for each genomic
        signals object, e.g,::

            plotting_kwargs = [dict(color='r', label='IP', dict(color='k',
                label='input')]

        """
        super(SignalMiniBrowser, self).__init__(genomic_signal_objs)
        self.plotting_kwargs = plotting_kwargs or [{} for i in genomic_signals]
        self.local_coverage_kwargs = local_coverage_kwargs or {}

    def panels(self):
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
        Signal panel on top, gene models on the bottom. `db` is
        a gffutils.FeatureDB.
        """
        super(GeneModelMiniBrowser, self).__init__(
                genomic_signal_objs, **kwargs)
        self.db = db

    def panels(self):
        ax1 = self.fig.add_subplot(211)
        ax2 = self.fig.add_subplot(212, sharex=ax1)
        return (ax2, self.gene_panel), (ax1, self.signal_panel)

    def gene_panel(self, ax, feature):
        from gffutils.contrib.plotting import Gene
        extent = [feature.start, feature.stop]
        nearby_genes = self.db.overlapping_features(
                feature.chrom, feature.start, feature.stop, featuretype='gene')
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
        or BedTool object that will be used to make
        a pybedtools.contrib.plotting.Track.
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
        ax.axis('tight')
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

    b = SignalMiniBrowser([ip, inp],
            plotting_kwargs=plotting_kwargs,
            local_coverage_kwargs=local_coverage_kwargs)

    g = GeneModelMiniBrowser([ip, inp], G,
            plotting_kwargs=plotting_kwargs,
            local_coverage_kwargs=local_coverage_kwargs)

    p = PeakMiniBrowser([ip, inp], peaks,
            plotting_kwargs=plotting_kwargs,
            local_coverage_kwargs=local_coverage_kwargs)

    feature = peaks[3]

    b.plot(feature)
    g.plot(feature)
    p.plot(feature)
    plt.show()
