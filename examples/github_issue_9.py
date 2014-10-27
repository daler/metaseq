import multiprocessing
from matplotlib import pyplot as plt
import matplotlib
import numpy as np
import metaseq
import pybedtools

# Use example data and generate some random features
gs = metaseq.genomic_signal(metaseq.example_filename('x.bam'), 'bam')

features = pybedtools.BedTool()\
    .window_maker(
        b=pybedtools.BedTool('chr2L 0 500000',
                             from_string=True).fn,
        w=1000)\
    .shuffle(seed=1, genome={'chr2L': (0, 5000000)})
genes = []
for i, f in enumerate(features):
    genes.append('gene_%s' % i)
genes = np.array(genes)
arr = gs.array(features, processes=multiprocessing.cpu_count(), bins=100)

# At this point, each item in `genes` corresponds to the same row in `arr`

ind, breaks = metaseq.plotutils.clustered_sortind(arr, k=5)

# Boundaries of clusters are provided in `breaks`.
# So the first cluster's original indices into `arr` are:
cluster_1_inds = ind[0:breaks[0]]

# Which means the genes in the first cluster are:
cluster_1_genes = genes[cluster_1_inds]

# More generally:
gene_clusters = []
start = 0
for b in breaks:
    gene_clusters.append(genes[ind[start:b]])
    start = b

# Plot everything
fig = plt.figure()
ax1 = fig.add_subplot(121)
ax1.matshow(arr, cmap=matplotlib.cm.hot)
ax1.set_title('unclustered')
ax1.axis('tight')

ax2 = fig.add_subplot(122)
ax2.matshow(arr[ind], cmap=matplotlib.cm.hot)
ax2.set_title('clustered, k=5')
for b in breaks:
    ax2.axhline(b, color='b')
ax2.axis('tight')


plt.show()
