"""
Settings for the ctcf_peaks example script
"""
import gffutils
import metaseq

UPSTREAM = 1000
DOWNSTREAM = 1000
BINS = 100
FRAGMENT_SIZE = 200
GENOME = 'hg19'
CHROMS = ['chr1', 'chr2']

gtfdb = metaseq.example_filename('Homo_sapiens.GRCh37.66.cleaned.gtf.db')
G = gffutils.FeatureDB(gtfdb)
