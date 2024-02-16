source('setup_2.R')
set.seed(1)

seu <- as.Seurat(sce)
seu <- FindNeighbors(seu, reduction = "mnn", dims = 1:25)

params <- data.frame(
    alg = rep(2:4, each=6),
    res = rep(c(.5,.6,.7,.8,.9,1), times=3)
)

clusMat <- sapply(1:nrow(params), function(i){
    seu <- FindClusters(seu, algorithm = params$alg[i], resolution = params$res[i])
    return(as.numeric(seu$seurat_clusters))
})
for(k in 4:10){
    km <- kmeans(reducedDim(sce,'mnn')[,1:25], centers = k)
    clusMat <- cbind(clusMat, km$cluster)
}
rm(km)
hc <- hclust(dist(reducedDim(sce,'mnn')[,1:25]))
for(k in 4:10){
    clusMat <- cbind(clusMat, cutree(hc, k = k))
}
rm(hc)




big <- apply(clusMat,1,function(x){
    apply(clusMat,1,function(y){
        mean(x==y)
    })
})
h <- hclust(as.dist(1-big))



# which clustering is most consistent with overall agreement?
diag(big) <- NA
within <- apply(clusMat,2,function(cl){
    means <- sapply(unique(cl), function(clID){
        idx <- which(cl == clID)
        return(c(mean(big[idx,idx], na.rm=TRUE), length(idx)))
    })
    return(sum(means[1,]*means[2,]) / length(cl))
})
between <- apply(clusMat,2,function(cl){
    means <- sapply(unique(cl), function(clID){
        idx <- which(cl == clID)
        return(c(mean(big[idx,-idx], na.rm=TRUE), length(idx)))
    })
    return(sum(means[1,]*means[2,]) / length(cl))
})
diag(big) <- 1

require(cluster)
d.all <- dist(reducedDim(sce,'mnn')[,1:25])
sil <- apply(clusMat,2,function(cl){
    s <- silhouette(as.integer(cl), d.all)
    return(mean(s[,'sil_width']))
})
rm(d.all)

require(dbscan)
nn <- kNN(reducedDim(sce,'mnn')[,1:25], k = 20)$id
locAg <- apply(clusMat,2,function(cl){
    nnClus <- matrix(cl[nn], ncol=20)
    return(mean(cl == nnClus))
})
rm(nn)


# plot(c(1,ncol(clusMat)), 0:1, col='white')
# points(within, col=4)
# points(between, col=2)
# abline(v = which.max(within - between))
# points(sil, col=3)
# points(locAg, col = 5)


# choose the winner
clus <- factor(clusMat[,which.max(within - between + locAg + sil)])




table(clus, sce$Sample)


# umap plot colored by cluster


sce$clus <- clus


# barplots showing breakdown of each cluster by sample
# table(clus, sce$Sample)


# save it, so I don't have to keep re-running all this
# saveRDS(sce, file='data/sce.rds')
# don't do this, just load the version from OneDrive

