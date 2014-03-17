import numpy as np
import os
import metaseq

ip_filename = metaseq.helpers.example_filename(
    'wgEncodeHaibTfbsK562Atf3V0416101AlnRep1_chr17.bam')
input_filename = metaseq.helpers.example_filename(
    'wgEncodeHaibTfbsK562RxlchV0416101AlnRep1_chr17.bam')

ip_signal = metaseq.genomic_signal(ip_filename, 'bam')
input_signal = metaseq.genomic_signal(input_filename, 'bam')

# If you already have TSSs, skip this part.
import gffutils
db = gffutils.FeatureDB(
    metaseq.example_filename('Homo_sapiens.GRCh37.66_chr17.gtf.db'))

import pybedtools
from pybedtools.featurefuncs import TSS
from gffutils.helpers import asinterval


def tss_generator():
    for transcript in db.features_of_type('transcript'):
        yield TSS(asinterval(transcript), upstream=1000, downstream=1000)

if not os.path.exists('tsses.gtf'):
    tsses = pybedtools.BedTool(tss_generator()).saveas('tsses.gtf')
tsses = pybedtools.BedTool('tsses.gtf')

from metaseq import persistence
if not os.path.exists('example.npz'):
    ip_array = ip_signal.array(tsses, bins=100, processes=8)
    input_array = input_signal.array(tsses, bins=100, processes=8)
    ip_array /= ip_signal.mapped_read_count() / 1e6
    input_array /= input_signal.mapped_read_count() / 1e6
    persistence.save_features_and_arrays(
        features=tsses,
        arrays={'ip': ip_array, 'input': input_array},
        prefix='example',
        link_features=True,
        overwrite=True)

features, arrays = persistence.load_features_and_arrays(prefix='example')
normalized = arrays['ip'] - arrays['input']
ind = metaseq.plotutils.tip_zscores(normalized)
fig = metaseq.plotutils.imshow(
    normalized,
    vmin=5,
    vmax=99.,
    percentile=True,
    sort_by=ind,
    imshow_kwargs=dict(interpolation='bilinear'),
    line_kwargs=dict(color='k'),
    fill_kwargs=dict(color='k', alpha=0.4),
    x=np.linspace(-1000, 1000, 100),
    height_ratios=(2, 1, 1)
)

fig.array_axes.xaxis.set_visible(False)
fig.array_axes.set_ylabel('Transcripts on chr17')
fig.array_axes.axvline(0, color='k', linestyle='--')

fig.line_axes.set_xlabel('Distance from TSS')
fig.line_axes.axvline(0, color='k', linestyle='--')

from matplotlib import pyplot as plt
import matplotlib

d = metaseq.results_table.ResultsTable(
    metaseq.example_filename('GSM847566_SL2592.table'),
    import_kwargs=dict(index_col=0))

d = d.reindex_to(features, attribute='transcript_id')
import pandas
labels = pandas.qcut(d.fpkm, 4).labels
ulabels = sorted(list(set(labels)))
colors = matplotlib.cm.YlOrBr((np.array(ulabels) + 2) / 5.)
bottom_axes = plt.subplot(fig.gs[2, 0])
for q, color in zip(ulabels, colors):
    ind = labels == q
    print q, color
    metaseq.plotutils.ci_plot(
        np.linspace(-1000, 1000, 100),
        normalized[ind, :],
        ax=bottom_axes,
        line_kwargs=dict(color=color, label=q),
        fill_kwargs=dict(color=color, alpha=0.5),
    )

fig.line_axes.xaxis.set_visible(False)
bottom_axes.set_xlabel('Distance from TSS')
bottom_axes.legend(loc='best', fontsize=10)
fig.array_axes.set_ylabel('Transcripts')
fig.cax.set_ylabel('Enrichment')
fig.subplots_adjust(left=0.2)
bottom_axes.set_ylabel('Enrichment')
fig.line_axes.set_ylabel('Enrichment')

plt.show()
