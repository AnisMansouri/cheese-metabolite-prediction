
metabolites_ecostab = readRDS("data/metabolites_ecostab.rds")
metabolites_execo = readRDS("data/metabolites_execo.rds")

# ============================================================
# 1. Run k-means on TRAINING data only
# ============================================================
set.seed(1)
clustering_ecostab = apply(
  metabolites_ecostab,
  1,
  function(x) {
    kmeans(
      x,
      centers = 2,
      nstart = 25
    )
  }
)

# ============================================================
# 2. Extract Low/High centres and thresholds
# ============================================================

cluster_info = t(
  sapply(
    clustering_ecostab,
    function(km) {
      
      centers = sort(as.numeric(km$centers))
      
      low_center  = centers[1]
      high_center = centers[2]
      
      threshold = mean(c(
        low_center,
        high_center
      ))
      
      c(
        Low_center = low_center,
        High_center = high_center,
        Threshold = threshold
      )
    }
  )
)

cluster_info = as.data.frame(cluster_info)

cluster_info$Metabolite = rownames(metabolites_ecostab)

cluster_info = cluster_info[
  ,
  c(
    "Metabolite",
    "Low_center",
    "High_center",
    "Threshold"
  )
]

# ============================================================
# 3. Classify TRAINING samples using training thresholds
# ============================================================

ecostab_classes = lapply(
  seq_len(nrow(metabolites_ecostab)),
  function(i) {
    
    threshold = cluster_info$Threshold[i]
    
    ifelse(
      metabolites_ecostab[i, ] < threshold,
      "Low",
      "High"
    )
  }
)

names(ecostab_classes) = rownames(metabolites_ecostab)

# ============================================================
# 4. Classify TEST samples using the SAME training thresholds
# ============================================================

execo_classes = lapply(
  seq_len(nrow(metabolites_execo)),
  function(i) {
    
    threshold = cluster_info$Threshold[i]
    
    ifelse(
      metabolites_execo[i, ] < threshold,
      "Low",
      "High"
    )
  }
)

names(execo_classes) <- rownames(metabolites_execo)

# ============================================================
# 5. Inspect results
# ============================================================

print(cluster_info)

cat("\nTraining class distribution:\n")
print(lapply(ecostab_classes, table))

cat("\nTest class distribution:\n")
print(lapply(execo_classes, table))

# ============================================================
# 6. Save outputs
# ============================================================

saveRDS(
  ecostab_classes,
  "data/ecostab_classes.rds"
)

saveRDS(
  execo_classes,
  "data/execo_classes.rds"
)

write.csv(
  round(cluster_info[, -1]),
  "metabolite_kmeans_thresholds.csv",
  quote = FALSE
)

  
