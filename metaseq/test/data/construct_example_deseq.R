set.seed(0)
library( DESeq )
fbgns = read.table('fbgns.txt')
fbgns = fbgns$V1
exampleFile = system.file( "extra/TagSeqExample.tab", package="DESeq" )
countsTable <- read.delim( exampleFile, header=TRUE, stringsAsFactors=TRUE )
rownames( countsTable ) <- countsTable$gene
countsTable <- countsTable[ , -1 ]
conds <- c( "T", "T", "T", "T2", "N", "N" )
cds <- newCountDataSet( countsTable, conds )
cds <- cds[ ,-1 ]
cds <- estimateSizeFactors( cds )
cds <- estimateVarianceFunctions( cds )
res <- nbinomTest( cds, "N", "T" )
res.5000 <- res[sample(nrow(res), 5000),]
res.5000$id <- sample(fbgns, 5000)
write.table(res.5000, file="ex.deseq", sep="\t", row.names=F)
