# Download large data for running metaseq tests.
#
# > 1.6 GB of downloads

import os
import hashlib
import gffutils

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

chroms_to_ignore = [
'MT', 'Y', 'GL000191.1', 'GL000192.1', 'GL000193.1', 'GL000194.1',
'GL000195.1', 'GL000199.1', 'GL000201.1', 'GL000204.1', 'GL000205.1',
'GL000209.1', 'GL000211.1', 'GL000212.1', 'GL000213.1', 'GL000216.1',
'GL000218.1', 'GL000219.1', 'GL000220.1', 'GL000222.1', 'GL000223.1',
'GL000224.1', 'GL000225.1', 'GL000228.1', 'GL000229.1', 'GL000230.1',
'GL000233.1', 'GL000236.1', 'GL000240.1', 'GL000241.1', 'GL000242.1',
'GL000243.1', 'GL000247.1', 'HG1000_2_PATCH', 'HG1032_PATCH',
'HG104_HG975_PATCH', 'HG115_PATCH', 'HG14_PATCH', 'HG183_PATCH', 'HG185_PATCH',
'HG186_PATCH', 'HG19_PATCH', 'HG243_PATCH', 'HG281_PATCH', 'HG480_HG481_PATCH',
'HG506_HG1000_1_PATCH', 'HG531_PATCH', 'HG536_PATCH', 'HG544_PATCH',
'HG686_PATCH', 'HG706_PATCH', 'HG730_PATCH', 'HG736_PATCH', 'HG745_PATCH',
'HG75_PATCH', 'HG79_PATCH', 'HG7_PATCH', 'HG858_PATCH', 'HG905_PATCH',
'HG946_PATCH', 'HG987_PATCH', 'HG989_PATCH', 'HG990_PATCH', 'HG991_PATCH',
'HG996_PATCH', 'HG998_1_PATCH', 'HG998_2_PATCH', 'HG999_1_PATCH',
'HG999_2_PATCH', 'HSCHR17_1', 'HSCHR4_1', 'HSCHR6_MHC_APD', 'HSCHR6_MHC_COX',
'HSCHR6_MHC_DBB', 'HSCHR6_MHC_MANN', 'HSCHR6_MHC_MCF', 'HSCHR6_MHC_QBL',
'HSCHR6_MHC_SSTO',]


gtf_fn, gtf_md5 = ('Homo_sapiens.GRCh37.66.gtf.gz', '25e76f628088daabd296447d06abe16b')
cleaned_fn, cleaned_md5 = ('Homo_sapiens.GRCh37.66.cleaned.gtf', '877208f1d322d0751600ccb1904b2d2d')
db_fn, db_md5 = ('Homo_sapiens.GRCh37.66.cleaned.gtf.db', '65d4014fee14198062e94bc3a9f18d18')

if os.path.exists(cleaned_fn) and hashlib.md5(open(cleaned_fn).read()).hexdigest() == cleaned_md5:
    print header, cleaned_fn, "up to date"
else:
    print header, "cleaning GTF..."
    gffutils.clean_gff(fn=gtf_fn, newfn=cleaned_fn, addchr=True, sanity_check=True, chroms_to_ignore=chroms_to_ignore)

if os.path.exists(db_fn) and hashlib.md5(open(db_fn).read()).hexdigest() == db_md5:
    print header, db_fn, "up to date"
else:
    gffutils.create_db(cleaned_fn, db_fn, verbose=True, force=True)
