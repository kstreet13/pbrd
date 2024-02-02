sce <- readRDS('data/sce.rds')

# want full dataset, in case marker genes not included in top 5000
require(DropletUtils)
sceFull <- read10xCounts(c('data/PRD_wk1_filtered_feature_bc_matrix','data/PBRD_wk1_filtered_feature_bc_matrix'),
                     sample.names = c('prd_1','pbrd_1'))
genes <- read.table('data/m33_gencode_Ens110_annotation.csv', sep = '\t', header = TRUE)
require(dplyr)
rowData(sceFull) <- left_join(as.data.frame(rowData(sceFull)), genes, join_by(ID == EnsemblID))
rm(genes)
sceFull$nUMI <- colSums(assay(sceFull,'counts'))
sceFull$nGene <- colSums(assay(sceFull,'counts') > 0)
norm <- log1p(1e4*t(t(assay(sceFull,'counts')) / sceFull$nUMI))
assay(sceFull,'logcounts') <- norm
rm(norm)

# From Brandon:
# Look at endothelial cell markers and HSC markers
# General vascular endothelial: CD34, PECAM (CD31*), CDH5 (VE-Cad*)
# Arterial blood endothelial: GJA4, HEY1, MECOM, SOX17, JAG1
# Venous blood endothelial: EFNB2, ACKR1, ICAM1, ERG2*, LRG1, SELP
# Lymphatic endothelial: PROX1, LYVE1, PDPN, NR2F2, VEGFR3*
# Hematopoeitic stem cell: SCA1, KIT (C-KIT*), SLAMF1, FLK2*
# * = not found
# +++ tracing markers: EGFP_bGHpA, ThreePrimeTdTomato

upgenes <- sapply(unique(sce$clus), function(clID){
    idx <- which(sce$clus == clID)
    m1 <- rowMeans(assay(sce,'logcounts')[,idx])
    m2 <- rowMeans(assay(sce,'logcounts')[,-idx])
    return(which.max(m1-m2))
})
upgenes <- unique(rownames(sce)[upgenes])

markers <- c("CD34", "PECAM1", "CDH5",
             "GJA4", "HEY1", "MECOM", "SOX17", "JAG1",
             "EFNB2", "ACKR1", "ICAM1", "LRG1", "SELP",
             "PROX1", "LYVE1", "PDPN", "NR2F2",
             "SCA1", "KIT", "SLAMF1",
             "EGFP_bGHpA", "ThreePrimeTdTomato")
markers <- sapply(markers, function(m){
    grep(paste0('^', m, '$'), rowData(sceFull)$Symbol, ignore.case = TRUE)
})
markers <- do.call(c, markers)
markers <- unique(rownames(sceFull)[markers])
markers <- markers[! markers %in% upgenes] # give priority to the unsupervised method
#

genes <- match(c(upgenes, markers), rownames(sceFull))
genecol <- c(rep(4,length(upgenes)), rep(6,length(markers)-2), 3, 2)

#

means <- sapply(unique(sce$clus), function(clID){
    sapply(genes, function(gID){
        mean(assay(sceFull,'logcounts')[gID, which(sce$clus == clID)])
    })
})
rownames(means) <- rowData(sceFull)$Symbol[genes]
pcts <- sapply(unique(sce$clus), function(clID){
    sapply(genes, function(gID){
        mean((assay(sceFull,'logcounts')[gID, which(sce$clus == clID)] > 0))
    })
})
rownames(pcts) <- rowData(sceFull)$Symbol[genes]

hc.genes <- hclust(dist(means))
hc.clus <- hclust(dist(t(means)))

means <- means[hc.genes$order, hc.clus$order]
pcts <- pcts[hc.genes$order, hc.clus$order]
genecol <- genecol[hc.genes$order]

png(filename = '~/Desktop/dots.png', width = 800, height = 1500, res=130)
plot(c(1,ncol(means)), c(1,nrow(means)), col='white', asp=1, axes=FALSE, xlab='Cluster', ylab='', main='Potential Marker Genes')
#abline(v = 1:ncol(means), col = 'lightgrey', lty = 2)
#abline(h = 1:nrow(means), col = 'lightgrey', lty = 2)
points(rep(1:ncol(means), each = nrow(means)),
       rep(1:nrow(means), times = ncol(means)),
       col = alpha(rep(genecol, times = ncol(means)), alpha=as.numeric(pcts)),
       cex = 1.8*sqrt(as.numeric(means)),
       pch = 16)
axis(1, at=1:ncol(means), labels = unique(sce$clus)[hc.clus$order], cex.axis=.7)
axis(2, at=1:nrow(means), labels = rownames(means), las=1, cex.axis=.6)
dev.off()


