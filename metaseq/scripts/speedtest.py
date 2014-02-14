#!/usr/bin/env python
usage="""
Benchmarking script to decide how many processors are optimal to use.

This will probably be a function of CPU speed and disk I/O throughput, so it's
best to run on each new hardware configuration.

Note that the actual time per 1k features on your own data will also be
a function of data density per feature (for BAM/BED/bigBed).

The result is a plot showing the actual time in seconds per 1000 features, over
a range of processors.  A dotted line for the optimal 1/x speedup.  The plot is
also saved as a PDF.

The default arguments work for data shipped with metaseq.
"""
import sys
import numpy as np
import time
from matplotlib import pyplot as plt
import pybedtools
import metaseq
import multiprocessing
import argparse
import os

ap = argparse.ArgumentParser(usage=usage)
ap.add_argument(
    '-n', '--nfeatures', type=int, default=1000,
    help='Number of random windows to make')
ap.add_argument(
    '-c', '--chunksize', type=int, default=100,
    help='Number of windows to send to each process')
ap.add_argument(
    '--start', type=int, default=10000,
    help='start coord of possible windows')
ap.add_argument(
    '--stop', type=int, default=5000000,
    help='stop coord of possible windows')
ap.add_argument('--chrom', default='chr2L', help='chromsome to make windows on')
ap.add_argument('--type', default='all', 
                help='Only use the specified file types.  Either "all" (default) or a comma-separated list of [bam, bigwig, bed, bigbed]')
ap.add_argument('--prefix', default=os.path.join(metaseq.data_dir(), 'x'), help='Prefix of filenames to use.  Expects files '
                'with this prefix, and the following suffixes: .bam, .bigwig, '
                '.bed.gz (should already be tabixed), .bigbed.  Default: %(default)s')
ap.add_argument('--plot-prefix', default='speedtest', help='Filename used to save the resulting plot. Default is %(default)s')
args = ap.parse_args()

requested = args.type.split(',')
allowed = ['bam', 'bed', 'bigwig', 'bigbed', 'all']
for req in requested:
    if req not in allowed:
        raise ValueError("%s not in %s" % (req, allowed))

plt.rcParams['font.size'] = 10

intervals = pybedtools.BedTool().window_maker(
    genome={args.chrom: (args.start, args.stop)}, n=args.nfeatures)\
    .shuffle(genome={args.chrom: (args.start, args.stop)}, seed=1)

size = (args.stop - args.start) / args.nfeatures
print size
print args.type

if args.type == 'all':
    requested = ['bigwig', 'bam', 'bigbed', 'bed']
signals = []
for req in requested:
    if req != 'bed':
        signals.append(metaseq.genomic_signal(args.prefix + '.' + req, req))
    else:
        signals.append(metaseq.genomic_signal(args.prefix + '.bed.gz', 'bed'))

d = {}
results = {}
processes = range(1, multiprocessing.cpu_count() + 1)
for x in signals:
    print x.kind
    sys.stdout.flush()
    times = []
    for p in processes:
        print p, 'processes'
        if p == 1:
            p = None
        t = []
        t0 = time.time()
        a = x.array(intervals, bins=100, processes=p,
                    chunksize=args.chunksize)
        elapsed = time.time() - t0
        times.append(elapsed)
    results[x.kind] = a.sum(axis=0)
    d[x.kind] = times

fig = plt.figure(figsize=(12, 6))
i = 1
for k, t in d.items():
    ax = fig.add_subplot(2, len(signals), i)
    t = np.array(t) / (len(intervals) / 1000.)
    ax.plot(processes, t, 'k.-', label=k)
    ax.plot(processes, 1 / np.array(processes, dtype=float) * t[0], 'k:')
    ax.set_title(k)
    plt.xlabel('# processes')
    plt.ylabel('time per 1k windows (s)')

    ax = fig.add_subplot(2, len(signals), i + len(signals))
    ax.plot(results[k], 'k')
    ax.set_ylabel('reads\nper bin')
    ax.set_xlabel('bins')
    i += 1
fig.subplots_adjust(left=0.1, right=0.95, wspace=0.5, hspace=0.5)
fig.suptitle('%s features, %sbp each, chunksize=%s' % (args.nfeatures, size, args.chunksize))
fig.savefig(
    'speedtest-%s_features-%s_bp_chunksize=%s.pdf' % (args.nfeatures, size, args.chunksize))
plt.show()
