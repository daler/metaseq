"""
Various helper functions for the atf3_peaks example script -- mostly for
generating and manipulating features
"""
import pybedtools
from gffutils.helpers import asinterval
import ctcf_peaks_settings as settings
chromsizes = pybedtools.chromsizes('hg19')


def chromfilter(feature):
    """
    Only pass selected chromosomes
    """
    return feature.chrom in settings.CHROMS


def TSS(feature, upstream=1000, downstream=1000):
    """
    Transforms a pybedtools.Interval, `feature`, into a TSS extended by
    upstream/downstream, paying attention to strand and proximity to chromosome
    limits.

    Also edits the feature type to be "TSS"
    """
    chrom_size = chromsizes[feature.chrom][1]
    if feature.strand == '-':
        start = max(0, feature.stop - downstream)
        stop = min(feature.stop + upstream, chrom_size)
    else:
        start = max(0, feature.start - upstream)
        stop = min(feature.start + downstream, chrom_size)

    # Modify featuretype
    feature[2] = 'TSS'
    feature.start = start
    feature.stop = stop
    return asinterval(feature)


def gene_generator():
    """
    The database has inferred full gene models from the GTF, so we can simply
    iterate over them here.

    More complex generators can be created as well -- for example, one that
    only returns unique TSS sites from all isoforms of all genes.
    """
    for g in settings.G.features_of_type('gene'):
        if g.chrom not in settings.CHROMS:
            continue
        yield asinterval(g)



def intron_generator():
    """
    Construct intron features by subtracting all exons from all genes.
    """
    genes = pybedtools.BedTool(
            asinterval(g) for g in settings.G.features_of_type('gene')\
                    if g.chrom in settings.CHROMS)
    exons = pybedtools.BedTool(
            asinterval(e) for e in settings.G.features_of_type('exon')\
                    if e.chrom in settings.CHROMS)
    for feature in genes.subtract(exons).saveas():
        yield feature


def peak_filter(x):
    """
    ENCODE data's 8th field have -log10 pvals; use a reasonably stringent
    filter here
    """
    if float(x[7]) > 5:
        return True


def peak_extender(x):
    """
    ENCODE data peaks are quite narrow; here, extend each by some distance up
    and downstream
    """
    chrom_size = chromsizes[x.chrom][1]
    midpoint = x.start + (x.stop - x.start) / 2
    start = max(0, midpoint - settings.UPSTREAM)
    stop = min(midpoint + settings.DOWNSTREAM, chrom_size)
    x.start = start
    x.stop = stop
    return x

