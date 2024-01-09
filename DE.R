sce <- readRDS('data/sce.rds')

# Pbrd is the experimental group and prd is the reference control group.

# We want to know what genes are altered in their expression in Pbrd in
# comparison to Prd, especially in the cells expressing EGFP and tdTomato.

# EGFP is listed as "EGFP_bGHpA" and TdTomato as "ThreePrimeTdTomato".


# cluster 1+13 and 6+8 are the two sets of clusters expressing EGFP and TdTomato

# do DE analysis (Pbrd vs. prd) within those pairs of clusters

sce$Sample <- factor(sce$Sample, levels = c('prd_1','pbrd_1'))


idx1 <- which(sce$clus %in% c(1,13))

#glm(assay(sce,'counts')[1,idx1] ~ sce$Sample[idx1], family = poisson, offset = sce$nUMI[idx1])

pvals <- sapply(1:nrow(sce), function(i){
    fit <- glm(assay(sce,'counts')[i,idx1] ~ sce$Sample[idx1], family = poisson)
    return(summary(fit)$coefficient[2,4])
})





