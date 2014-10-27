#!/usr/bin/env python
usage="""
This is a benchmarking script for the metaseq Python package
(https://github.com/daler/metaseq).  It can help decide how many processors are
optimal to use for your particular hardware configuration.

The default arguments work for data shipped along with metaseq.  The default
arguments generate 1000 random features of about 5kb in size and then extract
signal data for BED, BAM, bigWig, and bigBed formats. This is done using from
1 to however many CPUs are availble.

The result is a plot showing the actual time in seconds per 1000 features, over
a range of processors.  A dotted line for the optimal 1/x speedup.  The plot is
also saved as a PDF.

Note that results are a function of CPU speed, disk I/O throughput, and nature
of the data (for BAM, BED, and bigBed formats), so your mileage may vary.  For
example, multiple CPUs will not help if there's only a small number of
features.
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
    help='Number of random windows to make.  The size of each feature is '
    '(start - stop) / nfeatures, so changing this parameter will affect '
    'the size of each feature.')
ap.add_argument(
    '-c', '--chunksize', type=int, default=100,
    help='Number of windows to send to each process')
ap.add_argument(
    '--start', type=int, default=10000,
    help='start coord of possible windows. See --nfeatures help regarding size.')
ap.add_argument(
    '--stop', type=int, default=5000000,
    help='stop coord of possible windows. See --nfeatures help regarding size.')
ap.add_argument(
    '--chrom', default='chr2L',
    help='chromsome to make windows on')
ap.add_argument(
    '--type', default='all',
    help='Only use the specified file types.  Either "all" '
    '(default) or a comma-separated list of [bam, bigwig, bed, '
    'bigbed].')
ap.add_argument(
    '--prefix', default=os.path.join(metaseq.data_dir(), 'x'),
    help='Prefix of filenames to use.  Expects files '
    'with this prefix, and the following suffixes: .bam, .bigwig, '
    '.bed.gz (should already be tabixed), .bigbed.  Default: %(default)s')
ap.add_argument(
    '--plot-prefix', default='./speedtest',
    help='Filename used to save the resulting plot. Default is %(default)s')
ap.add_argument(
    '--bins', default=100,
    help='Number of bins for each feature')
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
if args.type == 'all':
    requested = ['bigwig', 'bam', 'bigbed', 'bed']
signals = []
for req in requested:
    if req != 'bed':
        signals.append(metaseq.genomic_signal(args.prefix + '.' + req, req))
    else:
        signals.append(metaseq.genomic_signal(args.prefix + '.bed.gz', 'bed'))

files = '* ' + '\n  * '.join([i.fn for i in signals])
plot_filename = (args.plot_prefix + '-%s_features-%s_bp_chunksize=%s.pdf'
                 % (args.nfeatures, size, args.chunksize))
processes = range(1, multiprocessing.cpu_count() + 1)
max_proc = processes[-1]
print """
{usage}
Parameters
----------
This script will generate {args.nfeatures} random features, each about {size}
bp, from genomic coordinates {args.chrom}:{args.start}-{args.stop}.

Genomic signal will be extracted from these files:
  {files}

For each file, the signal for the {args.nfeatures} features will be extracted
in parallel and binned into {args.bins} bins, using from 1 to {max_proc} CPUs.
Each CPU will get as many as {args.chunksize} features at a time.

Plot will be saved as {plot_filename}.
""".format(**locals())



d = {}
results = {}
for x in signals:
    print '\n' + x.kind, '[%s]' % x.fn
    sys.stdout.flush()
    times = []
    print 'CPUs:',
    for p in processes:
        print p,
        sys.stdout.flush()
        if p == 1:
            p = None
        t = []
        t0 = time.time()
        a = x.array(intervals, bins=args.bins, processes=p,
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
fig.savefig(plot_filename)
plt.show()
