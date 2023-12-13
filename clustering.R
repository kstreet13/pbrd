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
for(k in 3:8){
    km <- kmeans(reducedDim(sce,'mnn'), centers = k)
    clusMat <- cbind(clusMat, km$cluster)
}



big <- apply(clusMat,1,function(x){
    apply(clusMat,1,function(y){
        mean(x==y)
    })
})
h <- hclust(as.dist(1-big))
image(big[h$order,h$order])

clus <- factor(cutree(h, k = 17))

table(clus, sce$Sample)
plot(reducedDim(sce,'umap'),asp=1, col=colorby(clus))

pal <- colorby(factor(1:lenu(clus)))
centers <- t(sapply(levels(clus), function(clID){
    colMeans(reducedDim(sce,'umap')[which(clus==clID),])
}))
legend('right', legend=levels(clus), pch=16, col=pal, bty='n')
points(centers,pch=1,cex=2.5)
points(centers,pch=16,cex=2.5, col='grey80')
text(centers, labels = levels(clus), col = pal, font=2)

