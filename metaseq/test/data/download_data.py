# Download large data for running metaseq tests.
#
# > 1.6 GB of downloads

import os
import hashlib

# md5 hex digests for example files
MD5 = """
47b33c89d7dcc58067550c7c28796794  http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwTfbs/wgEncodeUwTfbsK562CtcfStdAlnRep1.bam
fa991dcf9da229136ed9ecdb710017d4  http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwTfbs/wgEncodeUwTfbsK562CtcfStdAlnRep1.bam.bai
0f3c15c9399e988b5b4e44a1c6d84cae  http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwTfbs/wgEncodeUwTfbsK562CtcfStdHotspotsRep1.broadPeak.gz
e8332ac71efe7497fe88b557ce628c9c  http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwTfbs/wgEncodeUwTfbsK562CtcfStdPkRep1.narrowPeak.gz
0290f47c68f51e4cfb256fb752cf3859  http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwTfbs/wgEncodeUwTfbsK562CtcfStdRawRep1.bigWig
87bf94aee05213d3abd24e4256b14d29  http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwTfbs/wgEncodeUwTfbsK562InputStdAlnRep1.bam
abd48cc689cd6efc49d2e42249540391  http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwTfbs/wgEncodeUwTfbsK562InputStdAlnRep1.bam.bai
171d60dcfd75ec4e646532c6688a1ef3  http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwTfbs/wgEncodeUwTfbsK562InputStdRawRep1.bigWig
25e76f628088daabd296447d06abe16b  ftp://ftp.ensembl.org/pub/release-66/gtf/homo_sapiens/Homo_sapiens.GRCh37.66.gtf.gz
"""
MD5 = dict([(i.split()[1], i.split()[0]) for i in MD5.splitlines() if len(i) > 0])

header = "[metachip download]"

for full_path, md5 in MD5.items():
    fn = os.path.basename(full_path)
    if os.path.exists(fn):
        if hashlib.md5(open(fn).read()).hexdigest() == md5:
            print header, fn, "up to date"
            continue
        else:
            print header, fn, "md5 hash does not match; downloading..."
    else:
        print header, fn, "does not exist, downloading..."
    cmds = [
            'wget',
            full_path]
    os.system(' '.join(cmds))


