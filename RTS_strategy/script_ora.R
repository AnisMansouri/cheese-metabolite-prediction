library(clusterProfiler)

selected_features = unlist(readRDS("data/selected_features.rds"))
selected_features = unique(selected_features)

universe = colnames(readRDS("data/rts.rds"))

ora = enrichKEGG(gene = selected_features, organism = "ko",
                 universe = universe)
write.csv(ora@result, "ora_RTS.csv")

