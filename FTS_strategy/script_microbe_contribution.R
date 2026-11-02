library(stringr)


transcriptomics_ecostab = readRDS("data/transcriptomics_ecostab.rds")
fts = readRDS("data/ecostab_norm.rds")


selected_features = readRDS("data/selected_features.rds")
selected_features = unique(unlist(selected_features))

# % of contribution of each microbe
selection_inter_fts = lapply(transcriptomics_ecostab, function(x) {
  step_1 = intersect(rownames(x), colnames(fts))
  length(intersect(step_1, selected_features))/length(step_1)
  })

#### 
prop_fts = unlist(selection_inter_fts)/max(unlist(selection_inter_fts))
colors = c("orange", "orange", "orange", "orange",
           "orange", "orange", "pink", "pink", "pink") 

barplot(unlist(selection_inter_fts)*100, 
        names.arg = gsub("_normalized", "", str_to_title(names(selection_inter_fts))),
        xlab = "Microbes", ylab = "Contribution %", 
        col = colors)
