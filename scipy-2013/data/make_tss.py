import gffutils
import pybedtools
from pybedtools.featurefuncs import TSS

db = gffutils.FeatureDB('data/Mus_musculus.GRCm38.71.gtf.chr.db')

def tssgen():
    for i in db.features_of_type('gene'):
        yield TSS(gffutils.helpers.asinterval(i), upstream=1, downstream=1)

tsses = pybedtools.BedTool(tssgen())\
    .saveas()\
    .remove_invalid()\
    .sort().saveas('TSS.gff')
