#! /usr/bin/python

# Download large data for running metaseq tests.
#
# > 1.6 GB of downloads

import sys
import os
import hashlib
import gffutils
import pybedtools
import metaseq
import logging
logging.basicConfig(level=logging.DEBUG, format='[%(name)s] [%(asctime)s]: %(message)s')
logger = logging.getLogger('metaseq data download')

DATA_DIR = metaseq.data_dir()

# md5 hex digests for example files ===========================================

# ENCODE data
DATA = """
ff6979ace9befe82e71b6a05609d36e1  http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibTfbs/wgEncodeHaibTfbsK562RxlchV0416101AlnRep1.bam
ab2f3d2efd5a0281092e7ad542dfad36  http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibTfbs/wgEncodeHaibTfbsK562RxlchV0416101AlnRep1.bam.bai
fa20b05ea082dcb063463b73b6a5af2f  http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibTfbs/wgEncodeHaibTfbsK562RxlchV0416101RawRep1.bigWig
b0716bd81170efe1fd0a8e411fb669d8  http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibTfbs/wgEncodeHaibTfbsK562Atf3V0416101AlnRep1.bam
6e8f85d3ab428ef95e3382237b1b2419  http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibTfbs/wgEncodeHaibTfbsK562Atf3V0416101PkRep1.broadPeak.gz
cf869424dc915e59d9f1b3f73d720883  http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibTfbs/wgEncodeHaibTfbsK562Atf3V0416101AlnRep1.bam.bai
fb3b9dc8e85636a3a1226d22f8c1dbec  http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibTfbs/wgEncodeHaibTfbsK562Atf3V0416101RawRep1.bigWig
"""

# Ensembl annotations
DATA += """
25e76f628088daabd296447d06abe16b  ftp://ftp.ensembl.org/pub/release-66/gtf/homo_sapiens/Homo_sapiens.GRCh37.66.gtf.gz
"""

# Cufflinks results files from GSE33816, ATF3 samples
#
# GSM847565_SL2585 = uninduced rep 1
# GSM847566_SL2592 = induced rep 1
# GSM847567_SL4337 = uninduced rep 2
# GSM847568_SL4326 = induced rep 2
DATA += """
4d02dcbd813a538bbbcf27b3732ef7aa  ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/samples/GSM847nnn/GSM847568/GSM847568_SL4326.gtf.gz
a419d585a4a214885d6f249b5fc9a3a4  ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/samples/GSM847nnn/GSM847567/GSM847567_SL4337.gtf.gz
ec23fb05d3b46cae6a3fcd38fd64a8c5  ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/samples/GSM847nnn/GSM847566/GSM847566_SL2592.gtf.gz
45cffa476d6e74b144744ef515bb433e  ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/samples/GSM847nnn/GSM847565/GSM847565_SL2585.gtf.gz
"""

# bigWig files from GSE33816, ATF3 samples
DATA += """
6449cb8a49a78a45de654aa8af54c732  ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/samples/GSM847nnn/GSM847568/GSM847568_SL4326.bw
9e02208fec0f23b2859f44df4ad6d7af  ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/samples/GSM847nnn/GSM847567/GSM847567_SL4337.bw
82919ea67c564f6786e29a9150255d2d  ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/samples/GSM847nnn/GSM847566/GSM847566_SL2592.bw
b2e33ceb52bbd35629c0c666ad820ac7  ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/samples/GSM847nnn/GSM847565/GSM847565_SL2585.bw
"""

items = []
for i in DATA.splitlines():
    if (len(i) > 0) and not (i.startswith('#')):
        items.append(i.strip().split())

header = "[metaseq download]"

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

gtf_fn, gtf_md5 = (os.path.join(DATA_DIR, 'Homo_sapiens.GRCh37.66.gtf.gz'), '25e76f628088daabd296447d06abe16b')
cleaned_fn, cleaned_md5 = (os.path.join(DATA_DIR, 'Homo_sapiens.GRCh37.66.cleaned.gtf'), '6964313797754c68ea0e892abbfdc9d4')
db_fn, db_md5 = (os.path.join(DATA_DIR, 'Homo_sapiens.GRCh37.66.cleaned.gtf.db'), 'bf0a69d0787d01d0e3241ee23b1c66e3')

def _up_to_date(md5, fn):
    if os.path.exists(fn):
        logger.info('calculating md5 for %s...' % os.path.basename(fn))
        if hashlib.md5(open(fn).read()).hexdigest() == md5:
            logger.info('up to date')
            return True
        else:
            logger.info('md5sum does not match.')
            os.unlink(fn)
            return False


def _just_download():
    for md5, full_path in items:
        fn = os.path.join(DATA_DIR, os.path.basename(full_path))
        if not _up_to_date(md5, fn):
            logger.info('%s, downloading...' % fn)
            cmds = [
                    'wget',
                    full_path,
                    ]
            os.system(' '.join(cmds))


def _gffutils_prep():
    if not _up_to_date(cleaned_md5, cleaned_fn):
        logger.info("cleaning GTF...")
        gffutils.clean_gff(fn=gtf_fn, newfn=cleaned_fn, addchr=True, sanity_check=True, chroms_to_ignore=chroms_to_ignore)

    if not _up_to_date(db_md5, db_fn):
        os.unlink(db_fn)
        gffutils.create_db(cleaned_fn, db_fn, verbose=True, force=True)


def _cufflinks_conversion():
    # convert Cufflinks output GTF files into tables
    fns = [
        ('5c601c78e89e8d76fa998e2462ea718f', 'GSM847568_SL4326.gtf.gz'),
        ('9616609cd5c2ef012341471c07988e69', 'GSM847567_SL4337.gtf.gz'),
        ('b46eb07e01abc6c76c096896582f4a2d', 'GSM847566_SL2592.gtf.gz'),
        ('845cf7e32c61703781cf3316b9452029', 'GSM847565_SL2585.gtf.gz'),
        ]
    for md5, fn in fns:
        fn = os.path.join(DATA_DIR, fn)
        table = fn.replace('.gtf.gz', '.table')
        if not _up_to_date(md5, table):
            logger.info("Converting Cufflinks GTF %s to table" % fn)
            fout = open(table, 'w')
            fout.write('id\tscore\tfpkm\n')
            x = pybedtools.BedTool(fn)
            seen = set()
            for i in x:
                accession = i['transcript_id'].split('.')[0]
                if accession not in seen:
                    seen.update([accession])
                    fout.write('\t'.join([accession, i.score, i['FPKM']]) + '\n')
            fout.close()

def download_and_prep_data():
    _just_download()
    _gffutils_prep()
    _cufflinks_conversion()

if __name__ == "__main__":
    download_and_prep_data()
