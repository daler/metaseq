
        library(DESeq)
        counts <- read.table(
            'counts.txt',
            sep='	',
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
        write.table(res, file='deseq_results.txt', row.names=FALSE, sep='\t')
        