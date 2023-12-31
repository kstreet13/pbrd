
# setup
require(DropletUtils)
sce <- read10xCounts(c('data/PRD_wk1_filtered_feature_bc_matrix','data/PBRD_wk1_filtered_feature_bc_matrix'),
                     sample.names = c('prd_1','pbrd_1'))
genes <- read.table('data/m33_gencode_Ens110_annotation.csv', sep = '\t', header = TRUE)
require(dplyr)
rowData(sce) <- left_join(as.data.frame(rowData(sce)), genes, join_by(ID == EnsemblID))
rm(genes)
sce$nUMI <- colSums(assay(sce,'counts'))
sce$nGene <- colSums(assay(sce,'counts') > 0)

# gene filtering
require(scry)
sce <- devianceFeatureSelection(sce, fam = 'binomial', batch = factor(sce$Sample))
#plot(sort(rowData(sce)$binomial_deviance, decreasing = TRUE))
#plot(log1p(sort(rowData(sce)$binomial_deviance, decreasing = TRUE)))
sce <- devianceFeatureSelection(sce, fam = 'binomial', nkeep = 5000, batch = factor(sce$Sample))

# initial dimred (~PCA)
#sce <- nullResiduals(sce, fam = 'binomial', type = 'deviance', batch = factor(sce$Sample))
res <- matrix(0.0, nrow = nrow(sce), ncol = ncol(sce))
for(b in unique(sce$Sample)){
    idx <- which(sce$Sample == b)
    res[,idx] <- as.matrix(nullResiduals(assay(sce,'counts')[,idx], fam = 'binomial', type = 'deviance'))
}
assay(sce,'binomial_deviance_residuals', withDimnames=FALSE) <- res
rm(res)


pca <- BiocSingular::runPCA(t(assay(sce,'binomial_deviance_residuals')), rank = 50)
reducedDim(sce, 'pca') <- pca$x
rm(pca)
# plot(pca$sdev^2) # 20 looks like enough, maybe 25

# UMAP
require(Seurat)
norm <- log1p(1e4*t(t(assay(sce,'counts')) / colSums(assay(sce,'counts'))))
assay(sce,'logcounts') <- norm
rm(norm)
seu <- as.Seurat(sce)
seu <- RunUMAP(seu, reduction = 'pca', dims = 1:25)
reducedDim(sce,'umap') <- seu@reductions$umap@cell.embeddings
rm(seu)

ind <- sample(ncol(sce))
plot(reducedDim(sce,'umap')[ind,],asp=1, col=colorby(sce$Sample[ind], alpha=.5))
layout(matrix(1:2,nrow=1))
plot(reducedDim(sce,'umap'),asp=1, col='grey90')
points(reducedDim(sce,'umap')[sce$Sample=='pbrd_1', ], col=alpha(brewer.pal(9,'Set1')[1], alpha=.5))
plot(reducedDim(sce,'umap'),asp=1, col='grey90')
points(reducedDim(sce,'umap')[sce$Sample=='prd_1', ], col=alpha(brewer.pal(9,'Set1')[2], alpha=.5))
layout(1)


# some batch effects present in UMAP, will try again with MNN

