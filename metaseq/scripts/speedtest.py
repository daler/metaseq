#!/usr/bin/python
"""
Benchmarking script to decide how many processors are optimal to use.

This will probably be a function of CPU speed and disk I/O throughput, so it's
best to run on each new hardware configuration.

Note that the actual time per 1k features on your own data will also be
a function of data density per feature (for BAM/BED/bigBed).
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

ap = argparse.ArgumentParser()
ap.add_argument('-n', type=int, help='Number of windows to make')
ap.add_argument('-c', type=int, help='Number of windows to send to each process')
ap.add_argument('--start', type=int, help='start coord')
ap.add_argument('--stop', type=int, help='stop coord')
ap.add_argument('--chrom', help='chromsome')
ap.add_argument('--prefix', help='Prefix of filenames to use.  Expects files '
                'with this prefix, and the following suffixes: .bam, .bigwig, '
                '.bed.gz (should already be tabixed), .bigbed')

args = ap.parse_args()

plt.rcParams['font.size'] = 10

intervals = pybedtools.BedTool().window_maker(
    genome={args.chrom: (args.start, args.stop)}, n=args.n)\
        .shuffle(genome={args.chrom: (args.start, args.stop)}, seed=1)

size = (args.stop - args.start) / args.n
print size

signals = [
    metaseq.genomic_signal(args.prefix + '.' + ext, ext) for ext in ['bam', 'bigwig', 'bigbed']
]

signals.append(metaseq.genomic_signal(args.prefix + '.bed.gz', 'bed'))

d = {}
results = {}
processes = range(1, multiprocessing.cpu_count() + 1)
for x in signals:
    print x.kind
    sys.stdout.flush()
    times = []
    for p in processes:
        t = []
        t0 = time.time()
        a = x.array(intervals, bins=100, processes=p,
                    chunksize=args.c)
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
    ax.plot(processes, 1/np.array(processes, dtype=float) * t[0], 'k:')
    ax.set_title(k)
    plt.xlabel('# processes')
    plt.ylabel('time per 1k windows (s)')

    ax = fig.add_subplot(2, len(signals), i + len(signals))
    ax.plot(results[k], 'k')
    ax.set_ylabel('reads\nper bin')
    ax.set_xlabel('bins')
    i += 1
fig.subplots_adjust(left=0.1, right=0.95, wspace=0.5, hspace=0.5)
fig.suptitle('%s features, %sbp each, chunksize=%s' % (args.n, size, args.c))
fig.savefig('speedtest-%s_features-%s_bp_chunksize=%s.pdf' % (args.n, size, args.c))
plt.show()
