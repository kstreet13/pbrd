sce <- readRDS('data/sce.rds')

# Pbrd is the experimental group and prd is the reference control group.

# We want to know what genes are altered in their expression in Pbrd in
# comparison to Prd, especially in the cells expressing EGFP and tdTomato.

# EGFP is listed as "EGFP_bGHpA" and TdTomato as "ThreePrimeTdTomato".


# cluster 1+13 and 6+8 are the two sets of clusters expressing EGFP and TdTomato

# do DE analysis (Pbrd vs. prd) within those pairs of clusters

sce$Sample <- factor(sce$Sample, levels = c('prd_1','pbrd_1'))


require(Seurat)
seu <- as.Seurat(sce)
Idents(seu) <- 'clus'
seu <- NormalizeData(seu)

de1 <- FindMarkers(seu, ident.1 = 'pbrd_1', group.by = 'Sample', subset.ident = '1')
de13 <- FindMarkers(seu, ident.1 = 'pbrd_1', group.by = 'Sample', subset.ident = '13')
de6 <- FindMarkers(seu, ident.1 = 'pbrd_1', group.by = 'Sample', subset.ident = '6')
de8 <- FindMarkers(seu, ident.1 = 'pbrd_1', group.by = 'Sample', subset.ident = '8')



# volcano plots
volcano_plot <- function(fc, pval, max = 20, ...){
    pval <- -log10(pval)
    pval[pval > max] <- max
    plot(fc, pval, ...)
}

png(filename = '~/Desktop/volcano_1.png', width = 800, height = 800)
volcano_plot(de1$avg_log2FC, de1$p_val_adj, 
             main = "Cluster 1", xlab = "log FC (Pbrd - Prd)", 
             ylab = "FDR-adjusted p-value", pch=16)
abline(h = -log10(.05), col=2, lty=2)
dev.off()

png(filename = '~/Desktop/volcano_13.png', width = 800, height = 800)
volcano_plot(de13$avg_log2FC, de13$p_val_adj, 
             main = "Cluster 13", xlab = "log FC (Pbrd - Prd)", 
             ylab = "FDR-adjusted p-value", pch=16)
abline(h = -log10(.05), col=2, lty=2)
dev.off()

png(filename = '~/Desktop/volcano_6.png', width = 800, height = 800)
volcano_plot(de6$avg_log2FC, de6$p_val_adj, 
             main = "Cluster 6", xlab = "log FC (Pbrd - Prd)", 
             ylab = "FDR-adjusted p-value", pch=16)
abline(h = -log10(.05), col=2, lty=2)
dev.off()




############

# other methods, too sensitive



idx1 <- which(sce$clus %in% c(1,13))

#glm(assay(sce,'counts')[1,idx1] ~ sce$Sample[idx1], family = poisson, offset = sce$nUMI[idx1])

pvals <- sapply(1:nrow(sce), function(i){
    fit <- glm(assay(sce,'counts')[i,idx1] ~ as.numeric(sce$Sample[idx1]), family = "poisson", offset = sce$nUMI[idx1])
    return(summary(fit)$coefficient[2,4])
})



require(DESeq2)
ds <- DESeqDataSetFromMatrix(countData = assay(sce,'counts')[,idx1],
                              colData = colData(sce)[idx1,],
                              design= ~ Sample)
res <- DESeq(ds)
tab <- results(res, name=resultsNames(res)[2])
rownames(tab)[which(tab$padj < .01)]




