# attempt to annotated PBRD data from public data (Jacob, et al. Developmental Cell, 2023)
require(zellkonverter)
ref <- readH5AD('~/Projects/pbrd/data/Jacob2023/Adata-object_All.h5')
# assays: 
# zscore is based on X
# X is cut off at 10 (max)
# normalized is log1p(1e6 * CPM)

# possible cluster labels:
# HC_named, subclustering, subclustering_grouped

sce <- readRDS('data/sce.rds')
clus <- sce$clus
rm(sce)
require(DropletUtils)
sceFull <- read10xCounts(c('data/PRD_wk1_filtered_feature_bc_matrix','data/PBRD_wk1_filtered_feature_bc_matrix'),
                         sample.names = c('prd_1','pbrd_1'))

# check overlap
mean(rownames(ref) %in% rowData(sceFull)$Symbol)
mean(rowData(sceFull)$Symbol %in% rownames(ref))

require(SingleR)
assay(sceFull,'normalized') <- log1p(t(t(assay(sceFull,'counts'))/colSums(assay(sceFull,'counts'))) * 10000)
rownames(sceFull) <- rowData(sceFull)$Symbol




# by clusters
pred.clus <- SingleR(test = sceFull, ref = ref, 
                    assay.type.test = 'normalized', assay.type.ref = 'normalized',
                    labels = ref$subclustering_grouped,
                    clusters = clus)

# by cells
pred.cell <- SingleR(test = sceFull, ref = ref, 
                    assay.type.test = 'normalized', assay.type.ref = 'normalized',
                    labels = ref$subclustering_grouped)
pred.cell <- readRDS('data/pred_labs_cell.rds')


clus.pred.lab <- pred.clus$labels[clus]
cell.pred.lab <- pred.cell$labels
mean(clus.pred.lab == cell.pred.lab)


plot(reducedDim(sce,'umap'),asp=1, col=colorby(sce$clus), pch=16)

plot(reducedDim(sce,'umap'),asp=1, col=c(1,3)[1+(clus.pred.lab == cell.pred.lab)], pch=16)

plot(reducedDim(sce,'umap'),asp=1, col=colorby(pred.clus$delta.next[clus]), pch=16)
plot(reducedDim(sce,'umap'),asp=1, col=colorby(pred.cell$delta.next), pch=16)


# cluster level confidence estimates
entropy <- function(x, base = exp(1)){
    p <- table(x)
    p <- p / length(x)
    p <- p[p > 0]
    -sum(p * log(p, base = base))
}
df <- data.frame(
    clus.delta = pred.clus$delta.next,
    lab.agreement = sapply(1:20, function(cl){ mean((clus.pred.lab == cell.pred.lab)[clus==cl]) }),
    cell.delta = sapply(1:20, function(cl){ mean(pred.cell$delta.next[clus==cl]) }),
    cell.entropy = sapply(1:20, function(cl){ entropy(cell.pred.lab[clus==cl]) })
)

most.likely <- sapply(1:20, function(cl){ names(which.max(table(cell.pred.lab[clus==cl]))) })
most.likely.pct <- sapply(1:20, function(cl){ max(table(cell.pred.lab[clus==cl])) / sum(clus==cl) })


png(filename = '~/Desktop/umap.png', width = 2000, height = 1000, res=100)
layout(matrix(1:2, nrow=1))
plot(reducedDim(sce,'umap'),asp=1, col=colorby(sce$clus), pch=16, main='UMAP - pbrd_1 and prd_1')
# with cluster labels on top
pal <- colorby(factor(1:lenu(sce$clus)))
centers <- t(sapply(levels(sce$clus), function(clID){
    colMeans(reducedDim(sce,'umap')[which(sce$clus==clID),])
}))
#legend('topright', legend=levels(sce$clus), pch=16, col=pal, bty='n')
points(centers,pch=1,cex=2.5)
points(centers,pch=16,cex=2.5, col=1)
text(centers, labels = levels(sce$clus), col = pal, font=2)

plot.new()
legend('left', col = pal, pch=16, bty='n', cex=1.25,
       legend = paste0(1:20, ': ',most.likely, ' (', format(100*most.likely.pct, digits=3, trim=TRUE), '%)'))
dev.off()




