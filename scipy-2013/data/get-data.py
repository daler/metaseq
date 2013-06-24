#!/usr/bin/python
"""
This script downloads a small-ish subset of mouse cell line ChIP-seq data from
ENCODE, genome annotations from Ensembl, and RNA-seq results from GEO.  The
data sets are then filtered to only keep data from chr7.

It will take 10 mins or so to complete, resulting in ~1GB or so of data.

External non-Python dependencies:

    wget
    samtools
    bigbedToBigWig
    R, with DESeq package installed from BioConductor

Python dependencies:

    pybedtools
    gffutils

"""
import os
import time
import pybedtools
from pybedtools.contrib import bigwig
from pybedtools.contrib import bigbed
import logging

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(message)s')

HISTONE_URL = ("http://hgdownload.cse.ucsc.edu/goldenPath/mm9/encodeDCC/"
               "wgEncodePsuHistone")
TF_URL = ("http://hgdownload.cse.ucsc.edu/goldenPath/mm9/encodeDCC/"
          "wgEncodePsuTfbs")
CHROM = "chr7"


def do_bam():
    """
    - download data only for CHROM (using samtools view)
    - create a new index for this subsetted data
    - create a bigWig file out of the subsetted data
    """

    bam_data = {
        # TF-like
        'G1E-ER_Gata1': (TF_URL, 'wgEncodePsuTfbsG1eer4Gata1aME0S129InputAlnRep1'),
        'G1E-ER_Input': (TF_URL, 'wgEncodePsuTfbsG1eer4InputME0S129InputAlnRep1'),
        #'G1E-ER_Tal1': (TF_URL, 'wgEncodePsuTfbsG1eer4e2Tal1ME0S129InputAlnRep1'),
        #'G1E-ER_Pol2': (TF_URL, 'wgEncodePsuTfbsG1eer4e2Pol24h8ME0S129InputAlnRep1'),
        #'G1E-ER_Ctcf': (TF_URL, 'wgEncodePsuTfbsG1eer4e2CtcfME0S129InputAlnRep1'),

        # Histone mods
        #'G1E-ER_H3K4me3': (HISTONE_URL, 'wgEncodePsuHistoneG1eer4e2H3k04me3ME0S129InputAlnRep1'),
        #'G1E-ER_H3K4me1': (HISTONE_URL, 'wgEncodePsuHistoneG1eer4e2H3k04me1ME0S129AlnRep1'),
        #'G1E-ER_H3K9me3': (HISTONE_URL, 'wgEncodePsuHistoneG1eer4e2H3k09me3ME0S129AlnRep1'),
        #'G1E-ER_H3K27me3': (HISTONE_URL, 'wgEncodePsuHistoneG1eer4e2H3k27me3ME0S129InputAlnRep1'),
        #'G1E-ER_H3K36me3': (HISTONE_URL, 'wgEncodePsuHistoneG1eer4e2H3k36me3BE0S129InputAlnRep1'),
        #'G1E-ER_histone_input': (HISTONE_URL, 'wgEncodePsuHistoneG1eer4e2InputME0S129InputAlnRep1'),

        #'G1E_Gata1': (TF_URL, 'wgEncodePsuTfbsG1eGata1aME0S129InputAlnRep1'),
        #'G1E_Input': (TF_URL, 'wgEncodePsuTfbsG1eInputME0S129InputAlnRep1'),
        #'G1E_H3K4me3': (HISTONE_URL, 'wgEncodePsuHistoneG1eH3k04me3ME0S129InputAlnRep1'),
        #'G1E_H3K4me1': (HISTONE_URL, 'wgEncodePsuHistoneG1eH3k04me1ME0S129InputAlnRep1'),
        #'G1E_histone_input': (HISTONE_URL, 'wgEncodePsuHistoneG1eer4e2InputME0S129InputAlnRep1'),
    }
    for label, bam in bam_data.items():
        _URL, bam = bam
        url = os.path.join(_URL, bam + '.bam')
        logging.info("""

    ===Downloading and indexing BAMs from ENCODE===
    URL  : %s
    label: %s
    chrom: %s
    """ % (url, label, CHROM))

        dest = label + '_%s.bam' % CHROM

        if os.path.exists(dest):
            logging.info("%s exists; skipping download" % dest)
        else:
            cmds = [
                'samtools',
                'view',
                '-b',
                '-F 0x04',
                url,
                '"%s"' % CHROM,
                '>', dest
            ]
            os.system(' '.join(cmds))

        # create a new index on the subsetted bam
        index = dest + '.bai'
        if os.path.exists(index):
            logging.info('%s exists; skipping indexing' % index)
        else:
            logging.info('indexing %s' % dest)
            cmds = [
                'samtools',
                'index',
                dest]
            os.system(' '.join(cmds))

        bw = label + '.bigwig'
        if os.path.exists(bw):
            logging.info('%s exists; skipping bigwig creation' % bw)
        else:
            logging.info("Making bigwig")
            bigwig.bam_to_bigwig(dest, genome='mm9', output=label + '.bigwig')


def do_gtf():
    """
    GTF file preparation
     - download full file from Ensembl
     - remove features from other chromosomes, and add "chr" to each chromosome
       to conform to UCSC names
     - create a gffutils db out of the cleaned file
    """
    # reads were aligned to mm9, so make sure we get the mm9 annotations
    gtf_url = ('ftp://ftp.ensembl.org/pub/release-67/gtf/mus_musculus/'
               'Mus_musculus.NCBIM37.67.gtf.gz')
    GTF = 'Mus_musculus.NCBIM37.67.gtf.gz'
    logging.info("""

    ===Downloading M. musculus GTF annotations from Ensembl===
    URL  : %s
    dest: %s
    """ % (gtf_url, GTF))
    cmds = [
        'wget',
        gtf_url,
        '-O', GTF
    ]
    if os.path.exists(GTF):
        logging.info("%s exists; skipping GTF download" % GTF)
    else:
        os.system(' '.join(cmds))

    logging.info("Filtering and trimming GTF")
    cleaned = 'Mus_musculus.NCBIM37.67_%s_cleaned.gtf' % CHROM
    if os.path.exists(cleaned):
        logging.info('%s exists; skipping GTF clean' % cleaned)
    else:
        def addchr(f):
            f[0] = 'chr' + f[0]
            return f

        g = pybedtools.BedTool(GTF)\
            .filter(lambda x: x[0] == CHROM.replace('chr', ''))\
            .each(addchr)\
            .saveas(cleaned)

    logging.info('Creating gffutils db')
    db = cleaned + '.db'
    if os.path.exists(db):
        logging.info("%s exists; skipping gffutils db creation" % db)
    else:
        os.system('gffutils-cli create %s' % cleaned)


peaks_data = {
    'G1E-ER_Gata1': 'wgEncodePsuTfbsG1eer4e2Gata1aME0S129InputPk.broadPeak.gz',
    'G1E-ER_Tal1': 'wgEncodePsuTfbsG1eer4e2Tal1ME0S129InputPk.broadPeak.gz',
}
def do_peaks():

    for label, bedgz in peaks_data.items():
        url = os.path.join(TF_URL, bedgz)
        logging.info("""

    ===Downloading and filtering peaks from ENCODE===
    URL  : %s
    label: %s
    chrom: %s
    """ % (url, label, CHROM))
        dest = label + '_%s.broadPeak' % CHROM
        if os.path.exists(dest):
            logging.info("%s exists; skipping peaks download and filtering" % dest)
        else:
            cmds = [
                'wget',
                url,
                '-O',
                dest + '.tmp']
            os.system(' '.join(cmds))

            logging.info("Keeping only features on %s" % CHROM)
            x = pybedtools.BedTool(dest + '.tmp')\
                .filter(lambda x: x.chrom == 'chr7')\
                .saveas(dest)
            logging.info("Removing tempfile")
            os.unlink(dest + '.tmp')


def do_bigbed():
    def convert(f):
        return pybedtools.create_interval_from_list([
            f.chrom,
            str(f.start),
            str(f.stop)])
    for label in peaks_data.keys():
        x = pybedtools.BedTool(label + '_' + CHROM + '.broadPeak')
        outfn = x.fn.replace('.broadPeak', '.bigbed')
        x = x.each(convert).sort()
        bigbed.bigbed(x.fn, genome='mm9', output=outfn)
        logging.info('wrote %s' % outfn)


def do_counts():
    """
    Run DESeq on data from PMID:23390196
    GEO accession: GSE43041

    WT and Ldb1 -/- KO mouse ES cells RNA-seq (HiSeq 2000, Bowtie 1.0 to mm9.
    """
    logging.info("""
     ===Downloading counts data from GSE43041 and running DESeq===
     """)
    final_result = 'deseq_results.txt'
    if os.path.exists(final_result):
        logging.info("%s exists; skipping" % final_result)
        return

    logging.info("Downloading count data")
    cmds = [
        'wget',
        'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE43nnn/GSE43041/suppl/GSE43041_Flk1%2B_WT_VS_KO_countable.txt.gz',
        '-O',
        'counts.tmp.txt.gz',
    ]
    os.system(' '.join(cmds))
    os.system('gunzip counts.tmp.txt.gz')

    logging.info("Trimming and filtering")
    os.system('cut -f 1,3,4,5,6,7,8,9,10 counts.tmp.txt | tail -n +2 > counts.txt')

    with open('DESeq.R', 'w') as fout:
        fout.write("""
        library(DESeq)
        counts <- read.table(
            'counts.txt',
            sep='\t',
            row.names=1,
            col.names=c('ID', 'WT', 'WT', 'WT', 'WT', 'KO', 'KO', 'KO', 'KO')
        )

        cds <- newCountDataSet(
            counts,
            conditions=c('WT', 'WT', 'WT', 'WT', 'KO', 'KO', 'KO', 'KO')
            )

        # Optional independent filtering
        #rs <- rowSums(counts(cds))
        #theta <- 0.4
        #use <- (rs > quantile(rs, probs=theta))
        #cds <- cds[use, ]

        cds <- estimateSizeFactors(cds)
        cds <- estimateDispersions(cds)
        res <- nbinomTest(cds, "WT", "KO")
        write.table(res, file='deseq_results.txt', row.names=FALSE, sep='\\t')
        """)


    logging.info("Running R script (DESeq.R) on counts; see DESeq.Rout for log")
    os.system('R CMD BATCH DESeq.R')


def do_tss():
    """
    after the gffutils db is created, make some TSS
    """
    logging.info("creating file of TSSs")
    import gffutils
    from pybedtools.featurefuncs import TSS
    if os.path.exists('TSS.gff'):
        logging.info("TSS.gff exists; skipping")
        return

    db = gffutils.FeatureDB('Mus_musculus.NCBIM37.67_chr7_cleaned.gtf.db')

    def tssgen():
        for i in db.features_of_type('gene'):
            yield TSS(gffutils.helpers.asinterval(i), upstream=1, downstream=1)

    tsses = pybedtools.BedTool(tssgen())\
        .saveas()\
        .remove_invalid()\
        .sort().saveas('TSS.gff')
    logging.info('created TSS.gff')

do_bam()
do_gtf()
do_peaks()
do_counts()
do_bigbed()
do_tss()
logging.info('Done.')
