# some of these won't work just by loading data/sce.rds, will need artifacts from clustering.R
# source('clustering.R')

sce <- readRDS('data/sce.rds')

# FROM SETUP

# umap plot, colored by sample
# png(filename = '~/Desktop/umap_samp.png', width = 1000, height = 1000, res=100)
ind <- sample(ncol(sce))
plot(reducedDim(sce,'umap')[ind,],asp=1, col=colorby(sce$Sample[ind], alpha=.5))
layout(matrix(1:2,nrow=1))
plot(reducedDim(sce,'umap'),asp=1, col='grey90')
points(reducedDim(sce,'umap')[sce$Sample=='pbrd_1', ], col=alpha(brewer.pal(9,'Set1')[1], alpha=.5))
plot(reducedDim(sce,'umap'),asp=1, col='grey90')
points(reducedDim(sce,'umap')[sce$Sample=='prd_1', ], col=alpha(brewer.pal(9,'Set1')[2], alpha=.5))
layout(1)
# dev.off()


# feature plot
colormap <- function(counts, max = 3871, relative = FALSE, hue=.56, sat=.9){
    stopifnot(all(counts%%1==0))
    cc <- rep('grey80', length(counts))
    #pal <- colorRampPalette(colors = c('grey80', brewer.pal(11,'Spectral')[6:11]))(100)[51:100]
    if(relative){
        max <- max(c(counts,2))
    }
    b <- seq(log1p(1), log1p(max), length.out = 31)
    nz <- which(counts > 0)
    pal <- hsv(hue, sat, seq(1,0, length.out=30))
    cc[nz] <- pal[cut(log1p(counts[nz]), breaks = b, include.lowest = TRUE)]
    return(cc)
}
feature_plot <- function(sce, gene, reduction = 'umap', hue = .56, sat = .9, ...){
    g.idx <- ifelse(is.numeric(gene), gene, which(rownames(sce)==gene))
    counts <- assay(sce,'counts')[g.idx,]
    
    nz <- which(counts > 0)
    gn <- rownames(sce)[g.idx]
    
    plot(range(reducedDim(sce,reduction)[,1]), range(reducedDim(sce,reduction)[,2]), col = 'white',
         xlab='',ylab='', main=gn, asp = 1, ...)
    points(reducedDim(sce,reduction)[-nz,], col = alpha('grey80',alpha=.5), cex = .5)
    points(reducedDim(sce,reduction)[nz,], col = colormap(counts[nz], relative = FALSE, hue = hue, sat = sat), cex = .5)
}

feature_plot(sce,'ThreePrimeTdTomato', hue = .99)
feature_plot(sce,'EGFP_bGHpA', hue = .35)


# FROM CLUSTERING

# co-clustering matrix
#png(filename = '~/Desktop/coclus.png', width = 1000, height = 1000)
image(big[h$order,h$order])
#dev.off()


# evaluate clusterings
plot(c(1,ncol(clusMat)), 0:1, col='white')
points(within, col=4)
points(between, col=2)
abline(v = which.max(within - between + locAg + sil))
points(sil, col=3)
points(locAg, col = 5)


# umap plot, colored by cluster
# png(filename = '~/Desktop/umap.png', width = 1000, height = 1000, res=100)
plot(reducedDim(sce,'umap'),asp=1, col=colorby(clus), pch=16, main='UMAP - pbrd_1 and prd_1')
# with cluster labels on top
pal <- colorby(factor(1:lenu(clus)))
centers <- t(sapply(levels(clus), function(clID){
    colMeans(reducedDim(sce,'umap')[which(clus==clID),])
}))
legend('topright', legend=levels(clus), pch=16, col=pal, bty='n')
points(centers,pch=1,cex=2.5)
points(centers,pch=16,cex=2.5, col=1)
text(centers, labels = levels(clus), col = pal, font=2)
# dev.off()


# barplots showing breakdown of each cluster by sample
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
