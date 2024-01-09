sce <- readRDS('data/sce.rds')

upgenes <- sapply(unique(sce$clus), function(clID){
    idx <- which(sce$clus == clID)
    m1 <- rowMeans(assay(sce,'logcounts')[,idx])
    m2 <- rowMeans(assay(sce,'logcounts')[,-idx])
    return(which.max(m1-m2))
})
genes <- c(upgenes, which(rownames(sce) == "EGFP_bGHpA"), 
           which(rownames(sce) == "ThreePrimeTdTomato"))
genecol <- c(rep(4,length(upgenes)), 3, 2)

means <- sapply(unique(sce$clus), function(clID){
    sapply(genes, function(gID){
        mean(assay(sce,'logcounts')[gID, which(sce$clus == clID)])
    })
})
rownames(means) <- rowData(sce)$Symbol[genes]
pcts <- sapply(unique(sce$clus), function(clID){
    sapply(genes, function(gID){
        mean((assay(sce,'logcounts')[gID, which(sce$clus == clID)] > 0))
    })
})
rownames(pcts) <- rowData(sce)$Symbol[genes]

hc.genes <- hclust(dist(means))
hc.clus <- hclust(dist(t(means)))

means <- means[hc.genes$order, hc.clus$order]
pcts <- pcts[hc.genes$order, hc.clus$order]
genecol <- genecol[hc.genes$order]

#png(filename = '~/Desktop/dots.png', width = 1000, height = 1000, res=130)
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
#dev.off()


