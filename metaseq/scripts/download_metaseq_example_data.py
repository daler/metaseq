#!/usr/bin/env python
import os
import subprocess
import logging
import hashlib
import pybedtools
import gffutils
import metaseq

logging.basicConfig(
    level=logging.DEBUG, format='[%(name)s] [%(asctime)s]: %(message)s')
logger = logging.getLogger('metaseq data download')

hg19 = pybedtools.chromsizes('hg19')
genomes_file = pybedtools.chromsizes_to_file('hg19', 'hg19')

usage = """
Downloads data from UCSC, GEO, and Ensembl.
"""

import argparse
ap = argparse.ArgumentParser(usage=usage)
ap.add_argument(
    '--data-dir',
    default=metaseq.data_dir(),
    help='Location to store downloaded and prepped data.  '
    'Default is %(default)s')
args = ap.parse_args()

CHROM = 'chr17'
COORD = "%s:%s-%s" % (CHROM, 0, hg19[CHROM][-1])


def requirements_check():
    """
    Ensure we have programs needed to download/manipulate the data
    """
    required_programs = [
        ('samtools',
         'http://samtools.sourceforge.net/'),
        ('bedtools',
         'http://bedtools.readthedocs.org/en/latest/'),
        ('bigWigToBedGraph',
         'http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/'),
        ('bedGraphToBigWig',
         'http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/'),
    ]
    for req, url in required_programs:
        try:
            p = subprocess.Popen(
                [req], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            stdout, stderr = p.communicate()
        except OSError:
            raise ValueError("Please install %s (%s)" % (req, url))


def _up_to_date(md5, fn):
    """
    Make sure md5sum(fn) == md5, and if not, delete `fn`.
    """
    if os.path.exists(fn):
        if hashlib.md5(open(fn).read()).hexdigest() == md5:
            logger.info('md5sum match for %s' % fn)
            return True
        else:
            logger.info('wrong md5sum for %s' % fn)
            os.unlink(fn)

# Each list contains 3-tuples of (human-readable size, md5sum, filename).
# `filename` is a URL in some cases.
#
# ChIP-seq BAM files from ENCODE
bams = [
    (
        # Input chromatin BAM, K562 cells
        '13M',
        '0270b88c30339659b38f11cd32b7fb11',
        ('http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHai'
         'bTfbs/wgEncodeHaibTfbsK562RxlchV0416101AlnRep1.bam')
    ),


    (
        # ATF3 ChIP-seq BAM, K562 cells
        '37M',
        'fd5fee8c3db467238e1d11d15c36f3c0',
        ('http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHai'
         'bTfbs/wgEncodeHaibTfbsK562Atf3V0416101AlnRep1.bam')
    ),
]

# Corresponding .bai files for the above BAMs.  These will be created.
bais = [
    (
        '173K',
        '873b0d236c30e089fe346f3901502329',
        os.path.join(
            args.data_dir,
            'wgEncodeHaibTfbsK562Atf3V0416101AlnRep1_chr17.bam.bai')),

    (
        '161K',
        'f40de02458ce71db708480953a307507',
        os.path.join(
            args.data_dir,
            'wgEncodeHaibTfbsK562RxlchV0416101AlnRep1_chr17.bam.bai')),
]

# bigWig files from ENCODE.
bigwigs = [
    (
        # Input chromatin bigWig
        '97K',
        '9408d7d066c2f2c9b9f8f61a3985e289',
        ('http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHai'
         'bTfbs/wgEncodeHaibTfbsK562RxlchV0416101RawRep1.bigWig')
    ),
    (
        # ATF3 ChIP-seq bigWig
        '831K',
        '33d8a7b43fbecc2fb54944711cfe179a',
        ('http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHai'
         'bTfbs/wgEncodeHaibTfbsK562Atf3V0416101RawRep1.bigWig')
    ),
]

# hg19 GTF annotations from Ensembl
GTF = (
    '19M',
    '0612ecf8bc707f6d8145d85337b95cba',
    ('ftp://ftp.ensembl.org/pub/release-66/gtf/homo_sapiens/Homo_sapiens.GRCh3'
     '7.66.gtf.gz')
)

# gffutils database to be created
DB = (
    '55M',
    'e9e91b7231fc0c9305da1cfbc58e6bbf',
    os.path.join(
        args.data_dir,
        'Homo_sapiens.GRCh37.66_%s.gtf.db' % CHROM))

# Tables of cufflinks GTF files, from GEO entry GSM847566.  These are results
# from ATF3 knockdown experiments in K562 cells.  Only one replicate is
# included here.
cufflinks = [
    (
        # ATF3 knockdow,, replicate 1
        '6.9M',
        'ec23fb05d3b46cae6a3fcd38fd64a8c5',
        ('ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/samples/GSM847'
         'nnn/GSM847566/GSM847566_SL2592.gtf.gz')
    ),
    (
        # ATF3 no knockdown, replicate 1
        '6.9M',
        '45cffa476d6e74b144744ef515bb433e',
        ('ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/samples/GSM847'
         'nnn/GSM847565/GSM847565_SL2585.gtf.gz')
    ),
]

# The cufflinks GTF files will be converted into tab-delimited files that can
# be later loaded.
cufflinks_tables = [
    (
        '6.9M',
        'b46eb07e01abc6c76c096896582f4a2d',
        'GSM847566_SL2592.gtf.gz'),
    (
        '6.9M',
        '845cf7e32c61703781cf3316b9452029',
        'GSM847565_SL2585.gtf.gz'),
]


def logged_command(cmds):
    "helper function to log a command and then run it"
    logger.info(' '.join(cmds))
    os.system(' '.join(cmds))


def get_cufflinks():
    "Download cufflinks GTF files"
    for size, md5, url in cufflinks:
        cuff_gtf = os.path.join(args.data_dir, os.path.basename(url))
        if not _up_to_date(md5, cuff_gtf):
            cmds = [
                'wget', url, '-O', cuff_gtf]
            logged_command(cmds)


def get_bams():
    """
    Download BAM files if needed, extract only chr17 reads, and regenerate .bai
    """
    for size, md5, url in bams:
        bam = os.path.join(
            args.data_dir,
            os.path.basename(url).replace('.bam', '_%s.bam' % CHROM))
        if not _up_to_date(md5, bam):
            logger.info(
                'Downloading reads on chromosome %s from %s to %s'
                % (CHROM, url, bam))
            cmds = ['samtools', 'view', '-b', url, COORD, '>', bam]
            logged_command(cmds)
        bai = bam + '.bai'
        if not os.path.exists(bai):
            logger.info('indexing %s' % bam)
            logger.info(' '.join(cmds))
            cmds = [
                'samtools',
                'index',
                bam]
            logged_command(cmds)
        if os.path.exists(os.path.basename(url) + '.bai'):
            os.unlink(os.path.basename(url) + '.bai')

    for size, md5, fn in bais:
        if not _up_to_date(md5, fn):
            cmds = [
                'samtools', 'index', bai.replace('.bai', '')]
            logged_command(cmds)


def get_bigwigs():
    for size, md5, url in bigwigs:
        bigwig = os.path.join(
            args.data_dir,
            os.path.basename(url).replace('.bigWig', '_%s.bigWig' % CHROM))
        bedgraph = os.path.join(
            args.data_dir,
            os.path.basename(url).replace('.bigWig', '.bedGraph'))
        subset_bedgraph = bedgraph.replace('.bedGraph', '_%s.bedGraph' % CHROM)
        if not _up_to_date(md5, bigwig):
            logger.info('Downloading bigWig %s to %s' % (url, bedgraph))
            cmds = [
                'bigWigToBedGraph', url, bedgraph]
            logged_command(cmds)

            cmds = [
                'awk -F "\\t" \'{if ($1 == "%s") print $0}\'' % CHROM,
                bedgraph, '>', subset_bedgraph]
            logged_command(cmds)

            cmds = [
                'bedGraphToBigWig',
                subset_bedgraph,
                genomes_file,
                bigwig]
            logged_command(cmds)

        if os.path.exists(subset_bedgraph):
            os.unlink(subset_bedgraph)
        if os.path.exists(bedgraph):
            os.unlink(bedgraph)


def get_gtf():
    """
    Download GTF file from Ensembl, only keeping the chr17 entries.
    """
    size, md5, url = GTF
    full_gtf = os.path.join(args.data_dir, os.path.basename(url))
    subset_gtf = os.path.join(
        args.data_dir,
        os.path.basename(url).replace('.gtf.gz', '_%s.gtf' % CHROM))

    if not _up_to_date(md5, subset_gtf):
        cmds = [
            'wget', url, '-O', full_gtf]
        logged_command(cmds)

        cmds = [
            'zcat',
            full_gtf,
            '|', 'awk -F "\\t" \'{if ($1 == "%s") print $0}\''
            % CHROM.replace('chr', ''),
            '|', 'awk \'{print "chr"$0}\'', '>', subset_gtf]

        logged_command(cmds)


def make_db():
    """
    Create gffutils database
    """
    size, md5, fn = DB
    if not _up_to_date(md5, fn):
        gffutils.create_db(fn.replace('.db', ''), fn, verbose=True, force=True)


def cufflinks_conversion():
    """
    convert Cufflinks output GTF files into tables of score and FPKM.
    """
    for size, md5, fn in cufflinks_tables:
        fn = os.path.join(args.data_dir, fn)
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
                    fout.write(
                        '\t'.join([accession, i.score, i['FPKM']]) + '\n')
            fout.close()


if __name__ == "__main__":
    get_bams()
    get_bigwigs()
    get_gtf()
    make_db()
    get_cufflinks()
    cufflinks_conversion()
