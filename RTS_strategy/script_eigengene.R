library(WGCNA)
library(parallel)
library(doParallel)
ncores = parallel::detectCores()
cl = makePSOCKcluster(ncores)
registerDoParallel(cl)

ecostab = readRDS("data/metabolites_ecostab.rds")

rts = readRDS("data/rts.rds")
selected_features = readRDS("data/selected_features.rds")

################### Permutation
i <<- 1
j <<- 1
perm = lapply(1:6, function(x) {
  
  dummy <<- rep(1, length(selected_features[[i]]))
  i <<- i + 1
  n <<- 1
  perm = lapply(1:1000, function(x) {
    
    print(j)
    set.seed(n)
    random_features = as.data.frame(rts[, 
                                        sample(1:1680, length(dummy))])
    names(dummy) = names(random_features)
    eg = moduleEigengenes(random_features, dummy)
    n <<- n + 1
    cor(as.numeric(ecostab[j, ]), eg$eigengenes[, 1], method = "spearman")
    
  })
  j <<- j + 1
  perm
})
saveRDS(perm, "permuted_correlations_RTS.rds")

# Compute the correlations and their pvalue by manually changing the metabolite (met = 1 for Alcohol
# met = 2 for Aldehydes etc...)
met = 1

gene_set = as.data.frame(rts[, selected_features[[met]]])
dummy = rep(1, length(selected_features[[met]]))
names(dummy) = colnames(gene_set)
eg = moduleEigengenes(gene_set, dummy)$eigengenes
true_cor = cor(as.numeric(ecostab[met, ]), 
                   eg[, 1], method = "spearman")
true_cor
random_cor = unlist(perm[[met]])
range(random_cor)
hist(random_cor)
length(random_cor[random_cor >= as.numeric(true_cor)])/1000 # The true correlations can be positive or negative, >= and <= are set accoringly


