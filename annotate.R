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



