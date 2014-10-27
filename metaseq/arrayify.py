"""
Helper functions for converting genome-wide data into large NumPy arrays
"""
import os
import pandas
import pybedtools
import numpy as np
import _genomic_signal


class Binner(object):
    def __init__(self, genome, windowsize, chrom=None, window_cache_dir=".",
                 npz_dir='.', metric='mean0'):
        """
        Class to handle converting bigWig files into NumPy arrays.

        bigWigAverageOverBed needs to be available on the path.

        The arrays are saved to disk as .npz files, which can then be
        memory-mapped for fast, lightweight re-use.  Each .npz file contains
        the coordinates of the bin midpoints (x) and values in each bin (y).

        The class is designed to be set up once, then used many times on
        different bigWig files.

        General usage:

            >>> b = Binner('mm9', 1000, chrom='chr19')
            >>> b.to_npz('PolII.bigwig')
            >>> pol = np.load('PolII.bigwig.chr19.npz')

            Assuming matplotlib is installed,
            >>> from matplotlib import pyplot as plt
            >>> plt.plot(pol['x'], pol['y'])
            >>> plt.show()

        Parameters
        ----------
        genome : str
            Assembly name to use (e.g., hg19, dm3).  This is used for creating
            the right number of windows.

        windowsize : int
            Bp to use in each window

        chrom : None or str
            If None, all chromosomes will be used; otherwise you can specify
            a single chromosome here.

        window_cache_dir : str
            Path where BED files containing windowed chromsome coordinates will
            be stored.  These files are cached to avoid creating them every
            time, and have the filename pattern
            {window_cache_dir}/{chrom}.{windowsize}bp_windows.bed

        """
        self.chromsizes = pybedtools.chromsizes(genome)
        if chrom is None:
            self.chroms = sorted(self.chromsizes.keys())
        else:
            self.chroms = [chrom]
        self.windowsize = windowsize
        self.window_cache_dir = window_cache_dir

    def make_windows(self, chrom, force=False):
        chromsize = self.chromsizes[chrom][1]
        bed = pybedtools.BedTool(
            '{0} 0 {1}'.format(chrom, chromsize),
            from_string=True
        )
        output = os.path.join(
            self.window_cache_dir,
            '%s.%sbp_windows.bed' % (chrom, self.windowsize))
        if not os.path.exists(output) and not force:
            windows = pybedtools.BedTool()\
                .window_maker(
                    b=bed,
                    i='winnum',
                    w=self.windowsize,
                    output=output,
                )
        return output

    def to_npz(self, bigwig, metric='mean0', outdir=None):
        """
        Bin data for bigwig and save to disk.

        The .npz file will have the pattern
        {outdir}/{bigwig}.{chrom}.{windowsize}.{metric}.npz and will have two
        arrays, x (genomic coordinates of midpoints of each window) and
        y (metric for each window).  It can be loaded like this::

            d = np.load(filename, mmap_mode='r')

        bigwig : str or BigWigSignal object
            BigWig data that will be used to create the array

        metric : 'covered', 'sum', 'mean0', 'mean'
            Metric to store in array, as reported by bigWigAverageOverBed:
                * "covered": the number of bases covered by the bigWig.
                * "sum": sum of values over all bases covered
                * "mean0": average over bases with non-covered bases counted as
                  zeros
                * mean: average over just the covered bases

        outdir : str or None
            Where to store output filenames.  If None, store the file in the
            same directory as the bigwig file.

        """
        if isinstance(bigwig, _genomic_signal.BigWigSignal):
            bigwig = bigwig.fn

        if outdir is None:
            outdir = os.path.dirname(bigwig)

        basename = os.path.basename(bigwig)
        windowsize = self.windowsize

        outfiles = []
        for chrom in self.chroms:
            tmp_output = pybedtools.BedTool._tmp()
            windows = self.make_windows(chrom)

            outfile = os.path.join(
                outdir,
                '{basename}.{chrom}.{windowsize}.{metric}'.format(**locals())
                + '.npz')
            cmds = [
                'bigWigAverageOverBed',
                bigwig,
                windows,
                tmp_output]
            os.system(' '.join(cmds))
            names = ['name', 'size', 'covered', 'sum', 'mean0', 'mean']
            df = pandas.read_table(tmp_output, names=names)
            x = df.size.cumsum() - df.size / 2
            y = df[metric]
            np.savez(outfile, x=x, y=y)
            outfiles.append(outfile)
            del x, y, df
        return outfiles
