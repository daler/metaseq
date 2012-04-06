# Download ENCODE data from UCSC used for testing metaseq.
#
# These data are from CTCF ChIP-seq in K562 cells.  For input and IP there are
# .bam, .bam.bai, and .bigWig files.  There are also called peaks and hotspots
# files.
#
# Total data size is ~1.6 GB

CMD="wget"
HOST="http://hgdownload.cse.ucsc.edu"

# 3.2 MB broadPeak file for CTCF hotspots
${CMD} ${HOST}/goldenPath/hg19/encodeDCC/wgEncodeUwTfbs/wgEncodeUwTfbsK562CtcfStdHotspotsRep1.broadPeak.gz

# 842 KB narrowPeak file for CTCF peaks
${CMD} ${HOST}/goldenPath/hg19/encodeDCC/wgEncodeUwTfbs/wgEncodeUwTfbsK562CtcfStdPkRep1.narrowPeak.gz

# 5.8 MB BAM index for input
${CMD} ${HOST}/goldenPath/hg19/encodeDCC/wgEncodeUwTfbs/wgEncodeUwTfbsK562InputStdAlnRep1.bam.bai

# 5.8 MB BAM index for CTCF ChIP
${CMD} ${HOST}/goldenPath/hg19/encodeDCC/wgEncodeUwTfbs/wgEncodeUwTfbsK562CtcfStdAlnRep1.bam.bai

# 128 MB bigWig for CTCF ChIP
${CMD} ${HOST}/goldenPath/hg19/encodeDCC/wgEncodeUwTfbs/wgEncodeUwTfbsK562CtcfStdRawRep1.bigWig

# 153 MB bigWig for input
${CMD} ${HOST}/goldenPath/hg19/encodeDCC/wgEncodeUwTfbs/wgEncodeUwTfbsK562InputStdRawRep1.bigWig

# 675 MB BAM for CTCF ChIP
${CMD} ${HOST}/goldenPath/hg19/encodeDCC/wgEncodeUwTfbs/wgEncodeUwTfbsK562CtcfStdAlnRep1.bam

# 656 MB BAM for input
${CMD} ${HOST}/goldenPath/hg19/encodeDCC/wgEncodeUwTfbs/wgEncodeUwTfbsK562InputStdAlnRep1.bam
