"""
This module integrates parts of metaseq that are useful for ChIP-seq analysis.
"""
import gffutils
from gffutils.helpers import asinterval
import metaseq
from metaseq.minibrowser import SignalMiniBrowser, GeneModelMiniBrowser
import numpy as np
from matplotlib import pyplot as plt


class Chipseq(object):
    """
    Class for visualizing and interactively exploring ChIP-seq data.

    Needs two BAM files (for IP and control) and a gffutils database
    filename in order to display gene models.

    Typical usage is to create a normalized array of signal over each
    feature with the `diff_array` method, and then plot with the `plot`
    method.

    The resulting figure has the matrix as a heatmap, the average signal
    over features, and a panel with points that can be zoomed and clicked,
    spawning a mini-browser window for the corresponding feature.

    Configuration can be done by adjusting the following attributes after
    creating a Chipseq instance:

        * `_strip_kwargs` (the style for the dots in the left panel)
        * `browser_plotting_kwargs` (style of signal lines in the mini-browser)

    Example usage:

        >>> dbfn=metaseq.example_filename(
        ...        'Homo_sapiens.GRCh37.66.cleaned.gtf.db')
        >>> C = Chipseq(
        ...        ip_bam=metaseq.example_filename(
        ...            'wgEncodeUwTfbsK562CtcfStdAlnRep1.bam'),
        ...        control_bam=metaseq.example_filename(
        ...            'wgEncodeUwTfbsK562InputStdAlnRep1.bam'),
        ...        dbfn=dbfn)
        >>> local_coverage_kwargs = dict(
        ...         fragment_size=200,
        ...         bins=100, chunksize=50, processes=6)

        >>> # make some features to use
        >>> G = gffutils.FeatureDB(dbfn)
        >>> genes = G.features_of_type('gene')
        >>> features = []
        >>> for i in range(100):
        ...    features.append(asinterval(genes.next()))

        >>> # x-axis for plots
        >>> x = np.arange(100)

        >>> # Create the array
        >>> C.diff_array(features=features, array_kwargs=local_coverage_kwargs)

        >>> # sort genes by
        >>> row_order = np.argsort(
        ...     plotutils.tip_zscores(C.diffed_array))[::-1]
        >>> C.plot(x=x, row_order=row_order)
        >>> plt.show()

    """
    def __init__(self, ip_bam, control_bam, dbfn):
        """
        Set up a :class:`Chipseq` object.


        :param ip_bam: filename of BAM file for ChIP data
        :param control_bam: filename of BAM file for control data
        :param dbfn: filename of gffutils database
        """
        self.ip = metaseq.genomic_signal(ip_bam, kind='bam')
        self.control = metaseq.genomic_signal(control_bam, kind='bam')
        self.dbfn = dbfn
        self.db = gffutils.FeatureDB(dbfn)
        self.ip_array = None
        self.control_array = None

        self._strip_kwargs = dict(color='.5', markeredgewidth=0, marker='o',
                linestyle='None', picker=5)
        self.browser_plotting_kwargs = [
                dict(color='r', label='IP'),
                dict(color='k', linestyle=':', label='control')]

    def diff_array(self, features, force=True, func=None,
            array_kwargs=dict()):
        """
        Scales the control and IP data to million mapped reads, then subtracts
        scaled control from scaled IP, applies `func(diffed)` to the diffed
        array, and finally sets `self.diffed_array` to be the result.

        Arrays `self.ip` and `self.control` are set as well, and if
        `force=False`, then previously-created arrays will be used instead of
        re-calculating new ones.  This is useful if you want to easily try
        multiple `func` functions without having to re-calculate the data.

        Another side-effect is that `self.features` is set so that it can be
        accesed by other methods.

        :param features: a list of pybedtools.Interval objects
        :param array_kwargs: extra keyword args passed to genomic_signal.array;
            typically this will include `bins`, `processes`, and `chunksize`
            arguments.
        :param func: a function to apply to the diffed arrays. By default
            this is :func:`metaseq.plotutils.nice_log`; another option might be
            `lambda x: x`, or `lambda x: 1e6*x`
        :param force: Force a re-calculation of the arrays; otherwise uses
            cached values
        """
        self.features = list(features)
        self.browser_local_coverage_kwargs = array_kwargs.copy()
        self.browser_local_coverage_kwargs.pop('processes')
        self.browser_local_coverage_kwargs.pop('chunksize')

        if (self.ip_array is None) or force:
            self.ip_array = self.ip.array(features, **array_kwargs)
            self.ip_array /= self.ip.million_mapped_reads()

        if (self.control_array is None) or force:
            self.control_array = self.control.array(features, **array_kwargs)
            self.control_array /= self.control.million_mapped_reads()

        if func is None:
            func = metaseq.plotutils.nice_log
        self.diffed_array = func(self.ip_array - self.control_array)

    def plot(self, x, row_order=None, imshow_kwargs=None):
        """
        Plot the scaled ChIP-seq data.

        :param x: X-axis to use (e.g, for TSS +/- 1kb with 100 bins, this would
            be `np.linspace(-1000, 1000, 100)`)
        :param row_order: Array-like object containing row order -- typically
            the result of an `np.argsort` call.
        """
        nrows = self.diffed_array.shape[0]
        if row_order is None:
            row_order = np.arange(nrows)
        extent = (min(x), max(x), 0, nrows)
        axes_info = metaseq.plotutils.matrix_and_line_shell(strip=True)
        fig, matrix_ax, line_ax, strip_ax, cbar_ax = axes_info
        _imshow_kwargs = dict(
                aspect='auto', extent=extent, interpolation='nearest')
        if imshow_kwargs:
            _imshow_kwargs.update(imshow_kwargs)

        mappable = matrix_ax.imshow(
                self.diffed_array[row_order],
                **_imshow_kwargs)
        plt.colorbar(mappable, cbar_ax)
        line_ax.plot(x, self.diffed_array.mean(axis=0))
        line, = strip_ax.plot(np.zeros((nrows,)), np.arange(nrows) + 0.5,
                **self._strip_kwargs)
        line.features = self.features
        line.ind = row_order

        matrix_ax.axis('tight')
        strip_ax.xaxis.set_visible(False)
        matrix_ax.yaxis.set_visible(False)
        matrix_ax.xaxis.set_visible(False)

        self.minibrowser = GeneModelMiniBrowser(
                [self.ip, self.control],
                db=self.db,
                plotting_kwargs=self.browser_plotting_kwargs,
                local_coverage_kwargs=self.browser_local_coverage_kwargs)

        fig.canvas.mpl_connect('pick_event', self.callback)

        self.axes = {
                'matrix_ax': matrix_ax,
                 'strip_ax': strip_ax,
                  'line_ax': line_ax,
                  'cbar_ax': cbar_ax
                }

    def callback(self, event):
        """
        Callback function to spawn a mini-browser when a feature is clicked.
        """
        artist = event.artist
        ind = artist.ind
        limit = 5
        browser = True
        if len(event.ind) > limit:
            print "more than %s genes selected; not spawning browsers" % limit
            browser = False
        for i in event.ind:
            feature = artist.features[ind[::-1][i]]
            print feature,
            if browser:
                self.minibrowser.plot(feature)


if __name__ == "__main__":

    dbfn=metaseq.example_filename('Homo_sapiens.GRCh37.66.cleaned.gtf.db')
    C = Chipseq(
            ip_bam=metaseq.example_filename(
                'wgEncodeUwTfbsK562CtcfStdAlnRep1.bam'),
            control_bam=metaseq.example_filename(
                'wgEncodeUwTfbsK562InputStdAlnRep1.bam'),
            dbfn=dbfn)

    local_coverage_kwargs = dict(fragment_size=200, bins=100, chunksize=50, processes=6)

    # make some features to use
    G = gffutils.FeatureDB(dbfn)
    genes = G.features_of_type('gene')
    features = []
    for i in range(1000):
        features.append(asinterval(genes.next()))

    # x-axis for plots
    x = np.arange(100)

    # Create the array
    C.diff_array(features=features, array_kwargs=local_coverage_kwargs)

    # sort genes by
    row_order = np.argsort(metaseq.plotutils.tip_zscores(C.diffed_array))[::-1]
    C.plot(x=x, row_order=row_order)
    plt.show()
