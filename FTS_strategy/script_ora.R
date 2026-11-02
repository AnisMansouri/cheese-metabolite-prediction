library(clusterProfiler)

selected_features = unlist(readRDS("data/selected_features.rds"))
selected_features = unique(selected_features)

# Set universe
fts = readRDS("data/ecostab_norm.rds")
universe = colnames(fts)

ora = enrichKEGG(gene = selected_features, organism = "ko",
                 universe = universe)
write.csv(ora@result, "ora_FTS.csv")




