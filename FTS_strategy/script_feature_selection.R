library(caret)
library(glmnet)

# Training set
fts = readRDS("data/ecostab_norm.rds")

# Loading of metabolites' classes
ecostab_classes = readRDS("data/ecostab_classes.rds")

# Train the EN models and select the best parameters alpha and lambda
elnet = lapply(ecostab_classes, function(x) {
    train(fts, x, method = "glmnet", metric = "Kappa",
          trControl = trainControl("LOOCV", classProbs = TRUE, savePredictions = "final"),
          tuneLength = 10)})


# Feature selection
met = 1  
selected_features = lapply(elnet, function(x) {
    
    en_exact = glmnet(
      x = fts,
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

