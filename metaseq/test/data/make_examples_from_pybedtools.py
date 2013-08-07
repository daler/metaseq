"""
take advantage of the example data shipped with pybedtools, but strip out read
names to save (slightly) on file size
"""
import os
import pysam
import pybedtools
import metaseq

# new bam file without the read names
x = pysam.Samfile(pybedtools.example_filename('x.bam'))
s = pysam.Samfile('x.bam', template=x, mode='wb')
for i, r in enumerate(x):
    r.qname = str(i)
    r.cigarstring = '36M'
    s.write(r)
s.close()

os.system('samtools index x.bam')
bam = pybedtools.BedTool('x.bam')

pybedtools.contrib.bigwig.bam_to_bigwig(
    bam.fn, 'dm3', 'x.bigwig', scale=False)

bed = bam.bam_to_bed(output='x.bed', split=True)
pybedtools.contrib.bigbed.bigbed(bed.sort(), 'dm3', 'x.bigbed')

bed.tabix(in_place=True)
os.unlink(bed.fn)
