
.. _example_session:

Example session
===============
This example session will demonstrate reading in BAM files from a ChIP-seq
experiment and plotting the normalized ChIP signal for TSSs of selected genes
as a heatmap.

.. note::

    **Using this documentation**: To follow along:

    * first install `IPython <http://ipython.org/ipython-doc/dev/index.html>`_
    * Start IPython in pylab mode: ``ipython --pylab``
    * For each block of code below, use the `%cpaste` magic function to paste
      in the code, and end the pasted entry with a ``--`` on a line by itself.
      For example::

        In [1]: %cpaste
        Pasting code; enter '--' alone on the line to stop or use Ctrl-D.
        :
        :In [1]: import metaseq
        :
        :# ATF3 ChIP in K562 cells from ENCODE
        : In [2]: ip_bam = metaseq.example_filename(
        :   ...:     'wgEncodeHaibTfbsK562Atf3V0416101AlnRep1.bam')
        :
        :--

    Do not paste in lines containing output (e.g., ``Out [1]: 5345``) as these
    will result in a syntax error.

Example data
------------
Example data can be downloaded from ENCODE and GEO using the installed script,
:file:`download_metaseq_example_data.py`.  This will download several GB of
publicly available data and store it in the :mod:`metaseq` installation folder.

After downloading the data, it can be easily accessed with the
:func:`metaseq.example_filename` function:

.. ipython::

    In [1]: import metaseq

    # ATF3 ChIP in K562 cells from ENCODE
    In [2]: ip_bam = metaseq.example_filename(
       ...:     'wgEncodeHaibTfbsK562Atf3V0416101AlnRep1.bam')

    # Reverse-crosslinked chromatin as control
    In [3]: control_bam = metaseq.example_filename(
       ...:     'wgEncodeHaibTfbsK562RxlchV0416101AlnRep1.bam')


    # A gffutils database that was created by the data download script from
    # the Ensembl GTF annotations for GRCh37.66.
    In [4]: dbfn = metaseq.example_filename(
       ...:     'Homo_sapiens.GRCh37.66.cleaned.gtf.db')


Create a :class:`Chipseq` object
--------------------------------
The :class:`metaseq.integration.chipseq.Chipseq` class creates
:class:`metaseq.genomic_signal.BamSignal` objects and has some useful methods
to manipulate them that we'll take advantage of later (specificialy,
:meth:`metaseq.integration.chipseq.Chipseq.plot` and
:meth:`metaseq.integration.chipseq.Chipseq.diff_array`). So let's create one
using the example filenames:

.. ipython::

    In [1]: from metaseq.integration import chipseq

    In [1]: chip = chipseq.Chipseq(ip_bam=ip_bam, control_bam=control_bam,
       ...:    dbfn=dbfn)


Select features of interest
---------------------------
Let's create some features to look at. A common task for transcription factors
is to focus on transcription start sites (TSSs) of annotated genes.

We can use the prepared database to get individual isoforms for
all genes. Recall that GTF files only contain exon/CDS/start codon/stop codon features,
not whole genes or transcripts.  One step of the database creation with
:mod:`gffutils` is to collect all exons for a particular transcript and infer
the start/stop coords of the full transcript -- then do the same for all
transcripts of a gene to infer the full gene coords.  In the end, it allows us
to access the transcript coords easily in the :func:`feature_generator`
function below.

.. note::

    This documentation is run every time it is generated to ensure correctness.
    As a result, many of the example analyses are truncated (for example, here
    restricting to only chromsome 19) so that the doctests run quickly.  In
    actual analyses, you wouldn't need such truncation.

Here, we write a generator function that only returns transcripts on chr19, and
then use the :mod:`pybedtools.featurefuncs.TSS` function on each one to get
just the TSS +/- 1kb for each transcript.


.. ipython::

    In [1]: import gffutils

    # asinterval will convert a gffutils.Feature into a pybedtools.Interval
    In [1]: from gffutils.helpers import asinterval

    # Connect to the database
    In [1]: db = gffutils.FeatureDB(dbfn)

    In [1]: def feature_generator():
       ...:     for gene in db.features_of_type('gene', chrom='chr19'):
       ...:         for transcript in db.children(gene):
       ...:             yield asinterval(transcript)

    In [1]: import pybedtools

    In [1]: from pybedtools.featurefuncs import TSS

    # Saves a temp file of transcripts
    In [1]: transcripts = pybedtools.BedTool(feature_generator()).saveas()

    # Converts transcripts to TSSs +/- 1kb and saveas a temp file
    In [1]: tss_features = transcripts\
       ...:     .each(TSS, upstream=1000, downstream=1000)\
       ...:     .saveas()

    # How many TSSs are there on chr19?
    In [1]: len(tss_features)


Calculate normalized enrichment values
--------------------------------------
We now have an interesting set of features and an object (`chip`) that
encapulates data we want to look at for these features.

In the end we would like to have some normalized value that we can think of as
"enrichment".  In order to do this, we need to first correct for differences in
library sizes between the control and IP, and then we need to use the control
sample to correct the IP -- for example, open chromatin at promoters can lead
to strong signal in the control that is not specific to the IP.

The :meth:`metaseq.integration.chipseq.Chipseq.diff_array` method does all of
this.  Specifically, it:
 * takes an iterable of features and other configuration info
 * computes and bins the coverage across each feature
 * scales the coverage to reads per million mapped reads (RPMMR)
 * does this for IP and control BAM files
 * subtracts control from IP, resulting in a matrix of enrichment across each
   feature
 * log2-transforms this difference matrix, dealing with negative numbers
   appropriately
 * sets the :attr:`Chipseq.diffed_array` attribute for access to the newly
   created array.

Most of this behavior is configurable if the defaults aren't suitable.

`array_kwargs` are used to configure parallel processing. In this case,
8 processes will be used and each process will get 50 features to work on at
a time.  Each read will be extended in the 3' direction to a total of 300 bp,
and coverage will be binned into 100 bins.

.. ipython::

    In [1]: chip.diff_array(
       ...:     features=tss_features, array_kwargs=dict(processes=8, chunksize=50,
       ...:     bins=100, fragment_size=300))

    # The created array has a row for each feature and a column for each bin
    In [1]: chip.diffed_array.shape


Plotting
--------
Determine sort order
~~~~~~~~~~~~~~~~~~~~

The sort order of the rows can be important for interpretation. One method is
to sort by the TIP score (see (see Cheng et al. 2001, Bioinformatics
27(23):3221-3227)):

.. ipython::

    In [1]: import numpy as np

    In [1]: tip_order = np.argsort(
       ...:     metaseq.plotutils.tip_zscores(chip.diffed_array))


Another, albeit contrived, method is to sort by the position of the highest
value in the row:

.. ipython::

    In [1]: other_order = np.argsort(chip.diffed_array.argmax(axis=1))

Colormap
~~~~~~~~

Colormap choice is also important for interpretation.  Often, values vary so
widely that using a default linear colormap results in a washed-out heatmap.
One option is to log-transform the data.  Another is to use
a `matplotlib.colors.Normalize` instance, setting vmax to the 99th percentile.
Or use a `matplotlib.colors.LogNorm` instance instead.

:func:`metaseq.colormap_adjust.smart_colormap` function centers the colormap on
zero and shows positive and negative values with different hues but equivalent
saturation and value (see http://en.wikipedia.org/wiki/HSL_and_HSV).

But first, let's put the array on a log scale, but in such a way that the
negative numbers (that is, disenriched regions) stay negative.  The
`metaseq.plotutils.nice_log` is useful for this:

.. ipython::

    In [1]: from metaseq.colormap_adjust import smart_colormap

    In [1]: from metaseq.plotutils import nice_log

    In [1]: # make a copy of the diffed array

    In [1]: backup = chip.diffed_array.copy()

    In [1]: chip.diffed_array = nice_log(chip.diffed_array)

    In [1]: cmap = smart_colormap(chip.diffed_array.min(), chip.diffed_array.max())

X-axis
~~~~~~

Next, we need to construct a nice x-axis that makes sense for this plot, with
upstream coords as negative and downstream as positive. With that, we can use
the :meth:`metaseq.integration.chipseq.Chipseq.plot` method that shows a nice
multi-panel figure:

.. ipython::


    In [1]: x = np.linspace(-1000, 1000, 100)

    @savefig first_array.png width=4in
    In [1]: fig = chip.plot(x, row_order=tip_order,
       ...:     imshow_kwargs=dict(cmap=cmap))

    In [1]: import matplotlib.pyplot as plt

    In [1]: plt.show()

Axes labels, etc
~~~~~~~~~~~~~~~~
The :meth:`Chipseq.plot` method sets the :attr:`axes` attribute on the :class:`Chipseq`
object so we can access the :class:`Axes` objects for further tweaking.  Lets
add axes labels, a title, and a dashed line indicating the position of the TSS:

.. ipython::

    In [1]: chip.axes['line_ax'].set_xlabel('Distance from TSS (bp)');

    In [1]: chip.axes['line_ax'].set_ylabel('Average enrichment (RPMMR)');

    In [1]: chip.axes['line_ax'].axvline(0, color='k', linestyle='--');

    In [1]: chip.axes['strip_ax'].set_ylabel('TSSs');

    In [1]: chip.axes['cbar_ax'].set_ylabel('Enrichment (RPMMR)');

    In [1]: chip.axes['matrix_ax'].axvline(0, color='k', linestyle='--');

    In [1]: chip.axes['matrix_ax'].set_title('ATF3 binding profile over TSSs on chr19');

    @savefig second_array.png width=4in
    In [1]: plt.draw()


.. Or, sort by the position:
.. 
.. .. ipython::

..    @savefig first_array.png width=4in
..    In [1]: fig = chip.plot(x, row_order=other_order,
..       ...:     imshow_kwargs=dict(cmap=cmap))

Interactive exploration
-----------------------

The left-hand axes contains a point for each TSS.  Clicking on a dot will open
a mini-browser window.  Interactively, you can use the :mod:`matplotlib` zoom
tools to zoom in to the, say, top 10 genes:

.. ipython::

    # Can do this interactively, or set axes limits manually:
    In [1]: chip.axes['matrix_ax'].axis(ymax=10)

    @savefig third_array.png width=4in
    In [1]: plt.draw()

Clicking on a point -- say, the one for the 9th TSS -- spawns a new browser
window.  Alternatively, the same window can be spawned using the
:meth:`Chipseq.minibrowser` method using the feature:

.. ipython::

    In [1]: feature = tss_features[tip_order[::-1][8]]

    @savefig minibrowser.png width=5in
    In [1]: chip.minibrowser.plot(feature)

Currently, the mini-browser shows the extent of the actual feature as vertical
dashed lines.  For context, it also shows the full extent of any genes that
happen to overlap the feature.  This behavior is extremely customizable by
creating subclasses -- see the docs for the :mod:`metaseq.minibrowser` module
for more info.

Improving the analysis
----------------------

OK, so it looks like there's a peak on centered on the TSS in the matrix plots.
We don't yet have a negative control or anything else we can compare this peak
to.  How about transcription termination site, or TTS, as a comparison?  Let's
calculate another diffed array using TTS features instead of TSS (this will
also demonstrate creating a custom feature manipulator). First, the new set of
features:

.. ipython::

    In [1]: def TTS(feature, upstream, downstream):
       ...:     if feature.strand == '-':
       ...:         tts = feature.start
       ...:         start = tts - downstream
       ...:         stop = tts + upstream
       ...:     else:
       ...:         tts = feature.stop
       ...:         start = tts - upstream
       ...:         stop = tts + downstream
       ...:     start = max(start, 0)
       ...:     feature.start = start
       ...:     feature.stop = stop
       ...:     return feature

    In [1]: tts_features = transcripts\
       ...:     .each(TTS, upstream=1000, downstream=1000)\
       ...:     .saveas()


Now we just need to pass the new features to :class:`chip` to re-create a new
array.  Let's save the old array first though so we don't lose it:

.. ipython::

    In [1]: tss_array = chip.diffed_array.copy()

    In [1]: chip.diff_array(
       ...:     features=tts_features, array_kwargs=dict(processes=8, chunksize=50,
       ...:     bins=100, fragment_size=300))

    In [1]: tts_array = chip.diffed_array.copy()

And then we can plot the average of both matrices on the same plot:

.. ipython::

    In [1]: fig = plt.figure()

    In [1]: ax = fig.add_subplot(111)

    In [1]: ax.plot(x, tss_array.mean(axis=0), color='r', label='TSS');

    In [1]: ax.plot(x, tts_array.mean(axis=0), color='b', label='TTS');

    In [1]: ax.set_xlabel('Distance from TTS or TSS (bp)');

    In [1]: ax.set_title('Average ATF3 signal for transcripts on chr19');

    In [1]: ax.legend(loc='best');

    @savefig comparison_plot.png width=4in
    In [1]: ax.set_ylabel('RPMMR');

Another potential improvement is being more sophisticated about choosing TSSs
and TTSs.  Specifically, as things stand now, if all isoforms of a gene have
the same TSS then all of those TSSs will be included in the plot leading to
duplication and possible over-estimating of the enrichment.  We might want to
require that all TSSs be unique (some "genome algebra" with :mod:`pybedtools`
would be useful here).  There are many other improvements that can be done, as
informed by the biology of the system being studied.  :mod:`metaseq` attempts
to make these manipulations easy (or at least straightforward) to do.
