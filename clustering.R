# source('setup_2.R')

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
png(filename = '~/Desktop/coclus.png', width = 1000, height = 1000)
image(big[h$order,h$order])
dev.off()



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



plot(c(1,ncol(clusMat)), 0:1, col='white')
points(within, col=4)
points(between, col=2)
abline(v = which.max(within - between))
points(sil, col=3)
points(locAg, col = 5)



clus <- factor(clusMat[,which.max(within - between + locAg + sil)])


table(clus, sce$Sample)

png(filename = '~/Desktop/umap.png', width = 1000, height = 1000, res=100)
plot(reducedDim(sce,'umap'),asp=1, col=colorby(clus), pch=16, main='UMAP - pbrd_1 and prd_1')

pal <- colorby(factor(1:lenu(clus)))
centers <- t(sapply(levels(clus), function(clID){
    colMeans(reducedDim(sce,'umap')[which(clus==clID),])
}))
legend('topright', legend=levels(clus), pch=16, col=pal, bty='n')
points(centers,pch=1,cex=2.5)
points(centers,pch=16,cex=2.5, col=1)
text(centers, labels = levels(clus), col = pal, font=2)
dev.off()


sce$clus <- clus



table(clus, sce$Sample)
png(filename = '~/Desktop/cluster_by_samp.png', width = 1000, height = 1000)
layout(matrix(1:24, ncol=4))
par(mar=c(3,3,3,1))
for(i in 1:22){
    barplot(table(clus, sce$Sample)[i,], col=pal[i], 
            main=paste('Cluster',i),
            ylim = c(0,max(table(clus,sce$Sample))))
}
dev.off()
layout(1)
par(mar=c(5,4,4,2)+.1)




