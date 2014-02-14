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
s = pysam.Samfile('x.bam.tmp', template=x, mode='wb')
for i, r in enumerate(x):
    r.qname = str(i)
    r.cigarstring = '36M'
    s.write(r)
s.close()

bam = pybedtools.BedTool('x.bam.tmp')
pybedtools.contrib.bigwig.bam_to_bigwig(
    bam.fn, 'dm3', 'x.bigwig', scale=False)

bed = bam.bam_to_bed(output='x.bed', split=True)
pybedtools.contrib.bigbed.bigbed(bed.sort(), 'dm3', 'x.bigbed')
bed.tabix(in_place=True)
os.unlink(bed.fn)

os.system(' '.join([
    'samtools view -h x.bam | awk -F "\t"',
    '\'{OFS="\t"; if (substr($1, 1, 1) == "@") print $0; else print $1,$2,$3,$4,$5,$6,$7,$8,$9,"*", "*",$12}\'',
    '| samtools view -Sb - > x2.bam'])
)

os.system('mv x2.bam x.bam')
os.system('samtools index x.bam')
os.system('rm x.bam.tmp')
