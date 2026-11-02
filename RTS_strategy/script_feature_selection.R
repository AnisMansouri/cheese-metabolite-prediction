library(caret)
library(glmnet)

# Training set
normalized_counts = readRDS("data/ecostab_norm.rds")

# Test set
#testset = readRDS("data/execo_norm.rds")
#testset = as.data.frame(testset)
#testset = replace(testset, is.na(testset), 0)
#saveRDS(testset, "data/testset.rds")

# Building of RTS
#common_genes = intersect(colnames(normalized_counts), colnames(testset))
#rts = normalized_counts[, common_genes]
#saveRDS(rts, "data/rts.rds")

rts = readRDS("data/rts.rds")

# Loading of metabolites' classes
ecostab_classes = readRDS("data/ecostab_classes.rds")

# Train the EN models and select the best parameters alpha and lambda
elnet = lapply(ecostab_classes, function(x) {
    train(rts, x, method = "glmnet", metric = "Kappa",
          trControl = trainControl("LOOCV", classProbs = TRUE, savePredictions = "final"),
          tuneLength = 10)})


# Feature selection
met = 1  
selected_features = lapply(elnet, function(x) {
    
    en_exact = glmnet(
      x = rts,
      y = ecostab_classes[[met]],
      family = "binomial",
      alpha = x$bestTune$alpha,
      lambda = x$bestTune$lambda
    )
    
    # Features' coefficients
    coefs = as.matrix(coef(en_exact))
    
    # Remove the intercept
    coefs = as.data.frame(coefs[-1, ])
    features = rownames(coefs)[coefs[, 1] != 0]
    met <<- met + 1
    features
    
  })

saveRDS(selected_features, "data/selected_features.rds")  

