library(stringr)


transcriptomics_ecostab = readRDS("data/transcriptomics_ecostab.rds")
rts = readRDS("data/rts.rds")


selected_features = readRDS("data/selected_features.rds")
selected_features = unique(unlist(selected_features))

# % of contribution of each microbe

selection_inter_rts = lapply(transcriptomics_ecostab, function(x) {
  step_1 = intersect(rownames(x), colnames(rts))
  length(intersect(step_1, selected_features))/length(step_1)
  })



# Plot
colors = c("orange", "orange", "orange", "orange",
           "orange", "orange", "pink", "pink", "pink") 

barplot(unlist(selection_inter_rts)*100, 
        names.arg = gsub("_normalized", "", 
                         str_to_title(names(selection_inter_rts))),
        xlab = "Microbes", ylab = "Contribution %", 
        col = colors)
