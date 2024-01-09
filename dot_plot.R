
upgenes <- sapply(unique(sce$clus), function(clID){
    idx <- which(sce$clus == clID)
    m1 <- rowMeans(assay(sce,'logcounts')[,idx])
    m2 <- rowMeans(assay(sce,'logcounts')[,-idx])
    return(which.max(m1-m2))
})

means <- sapply(unique(sce$clus), function(clID){
    sapply(upgenes, function(gID){
        mean(assay(sce,'logcounts')[gID, which(sce$clus == clID)])
    })
})
rownames(means) <- rowData(sce)$Symbol[upgenes]
pcts <- sapply(unique(sce$clus), function(clID){
    sapply(upgenes, function(gID){
        mean((assay(sce,'logcounts')[gID, which(sce$clus == clID)] > 0))
    })
})
rownames(pcts) <- rowData(sce)$Symbol[upgenes]

hc.genes <- hclust(dist(means))
hc.clus <- hclust(dist(t(means)))

means <- means[hc.genes$order, hc.clus$order]
pcts <- pcts[hc.genes$order, hc.clus$order]

png(filename = '~/Desktop/dots.png', width = 1000, height = 1000, res=130)
plot(c(1,ncol(means)), c(1,nrow(means)), col='white', asp=1, axes=FALSE, xlab='Cluster', ylab='', main='Potential Marker Genes')
points(rep(1:nrow(means), each = ncol(means)),
       rep(1:ncol(means), times = nrow(means)),
       col = alpha(4, alpha=as.numeric(pcts)),
       cex = 1.8*sqrt(as.numeric(means)),
       pch = 16)
axis(1, at=1:ncol(means), labels = unique(sce$clus)[hc.clus$order], cex.axis=.75)
axis(2, at=1:nrow(means), labels = rownames(means), las=1, cex.axis=.75)
dev.off()


