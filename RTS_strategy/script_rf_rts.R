library(caret)
library(randomForest)
library(ROCR)
library(parallel)
library(doParallel)
library(foreach)
library(dplyr)


# ============================================================
# 1. Load data
# ============================================================

# RTS training set
rts = readRDS("data/rts.rds")

# Independent test set
testset = readRDS("data/testset.rds")

# Metabolite classes
ecostab_classes = readRDS("data/ecostab_classes.rds")
execo_classes   = readRDS("data/execo_classes.rds")

# EN-selected features
selected_features = readRDS("data/selected_features.rds")


# ============================================================
# 2. Metabolite names
# ============================================================

metabolite_names = c(
  "Alcohols",
  "Aldehydes",
  "Alkanes",
  "Ketones",
  "Esters",
  "Sulphur compounds"
)


# ============================================================
# 3. mtry grids
# ============================================================

mtry_list = list(
  
  # Alcohols
  c(1, 2, 3, 5, 7, 10, 14, 20, 30, 50),
  
  # Aldehydes
  c(1, 2, 3, 5, 7, 10, 14, 20, 30, 50),
  
  # Alkanes
  c(1, 2, 3, 5, 7, 10, 15, 20, 30, 50),
  
  # Ketones
  1:10,
  
  # Esters
  c(1, 2, 3, 5, 7, 10, 13, 20, 30, 50),
  
  # Sulphur compounds
  c(1, 2, 3, 5, 7, 10, 15, 20, 30)
)

names(mtry_list) = metabolite_names


# Remove mtry values larger than the number of selected features
mtry_list = lapply(
  seq_along(mtry_list),
  function(i) {
    
    p = length(selected_features[[i]])
    
    sort(
      unique(
        mtry_list[[i]][mtry_list[[i]] <= p]
      )
    )
  }
)

names(mtry_list) = metabolite_names


# ============================================================
# 4. Function to generate reproducible caret seeds
# ============================================================

make_caret_seeds = function(
    n_resamples,
    n_models,
    base_seed
) {
  
  set.seed(base_seed)
  
  seeds <- vector(
    mode = "list",
    length = n_resamples + 1
  )
  
  # Seeds for every mtry model within every LOOCV resample
  for (i in seq_len(n_resamples)) {
    
    seeds[[i]] <- sample.int(
      1000000,
      size = n_models,
      replace = FALSE
    )
  }
  
  # Seed for final RF model
  seeds[[n_resamples + 1]] <- sample.int(
    1000000,
    size = 1
  )
  
  seeds
}


# ============================================================
# 5. Basic consistency checks
# ============================================================

stopifnot(
  length(selected_features) == 6,
  length(ecostab_classes) == 6,
  length(execo_classes) == 6,
  length(mtry_list) == 6
)

stopifnot(
  all(
    sapply(
      ecostab_classes,
      length
    ) == nrow(rts)
  )
)

stopifnot(
  all(
    sapply(
      execo_classes,
      length
    ) == nrow(testset)
  )
)


# Check features are present in RTS and test set
for (mt in seq_along(selected_features)) {
  
  missing_train = setdiff(
    selected_features[[mt]],
    colnames(rts)
  )
  
  missing_test <- setdiff(
    selected_features[[mt]],
    colnames(testset)
  )
  
  if (length(missing_train) > 0) {
    
    stop(
      paste(
        "Training predictors missing for",
        metabolite_names[mt],
        ":",
        paste(missing_train, collapse = ", ")
      )
    )
  }
  
  if (length(missing_test) > 0) {
    
    stop(
      paste(
        "Test predictors missing for",
        metabolite_names[mt],
        ":",
        paste(missing_test, collapse = ", ")
      )
    )
  }
}


# ============================================================
# 6. Parallel backend
# ============================================================

ncores = suppressWarnings(
  as.integer(
    Sys.getenv("SLURM_CPUS_PER_TASK")
  )
)

if (is.na(ncores) || ncores < 1) {
  
  ncores = suppressWarnings(
    as.integer(
      Sys.getenv("SLURM_CPUS_ON_NODE")
    )
  )
}

# Fallback outside SLURM
if (is.na(ncores) || ncores < 1) {
  
  ncores = parallel::detectCores()
}


cat(
  "Using",
  ncores,
  "CPUs\n"
)


# Prevent nested multithreading
Sys.setenv(
  OMP_NUM_THREADS = 1,
  OPENBLAS_NUM_THREADS = 1,
  MKL_NUM_THREADS = 1
)


cl = makePSOCKcluster(ncores)

registerDoParallel(cl)


cat(
  "Registered parallel workers:",
  foreach::getDoParWorkers(),
  "\n"
)


# ============================================================
# 7. Number of RF repetitions
# ============================================================

n_iterations = 1000


# ============================================================
# 8. Run RTS Random Forest analysis
# ============================================================

results_rts = tryCatch({
  
  lapply(
    seq_len(n_iterations),
    function(iter) {
      
      cat(
        "\n==============================\n",
        "Iteration:", iter, "\n",
        "==============================\n"
      )
      
      
      metabolite_results = lapply(
        seq_along(selected_features),
        function(mt) {
          
          cat(
            metabolite_names[mt],
            "\n"
          )
          
          
          # ==================================================
          # EN-selected predictors
          # ==================================================
          
          features = selected_features[[mt]]
          
          
          x_train = rts[
            ,
            features,
            drop = FALSE
          ]
          
          
          x_test = testset[
            ,
            features,
            drop = FALSE
          ]
          
          
          # ==================================================
          # Outcomes
          # ==================================================
          
          y_train = factor(
            ecostab_classes[[mt]],
            levels = c("Low", "High")
          )
          
          
          y_test = factor(
            execo_classes[[mt]],
            levels = c("Low", "High")
          )
          
          
          if (any(is.na(y_train))) {
            
            stop(
              paste(
                "Unexpected training class label for",
                metabolite_names[mt]
              )
            )
          }
          
          
          if (any(is.na(y_test))) {
            
            stop(
              paste(
                "Unexpected test class label for",
                metabolite_names[mt]
              )
            )
          }
          
          
          # Training data must contain both classes
          if (length(unique(y_train)) < 2) {
            
            stop(
              paste(
                "Training data contain only one class for",
                metabolite_names[mt]
              )
            )
          }
          
          
          # ==================================================
          # mtry grid
          # ==================================================
          
          tune_grid = expand.grid(
            mtry = mtry_list[[mt]]
          )
          
          
          # ==================================================
          # Reproducible seeds
          # ==================================================
          
          base_seed = iter * 100 + mt
          
          
          caret_seeds = make_caret_seeds(
            n_resamples = nrow(rts),
            n_models = nrow(tune_grid),
            base_seed = base_seed
          )
          
          
          set.seed(base_seed)
          
          
          # ==================================================
          # Train Random Forest
          #
          # mtry is tuned using LOOCV on RTS.
          # The final model is then fitted using all RTS
          # training observations.
          # ==================================================
          
          rf_model = train(
            
            x = x_train,
            y = y_train,
            
            method = "rf",
            
            metric = "Kappa",
            
            importance = TRUE,
            
            tuneGrid = tune_grid,
            
            ntree = 1000,
            
            trControl = trainControl(
              
              method = "LOOCV",
              
              classProbs = TRUE,
              
              savePredictions = "final",
              
              allowParallel = TRUE,
              
              seeds = caret_seeds
            )
          )
          
          
          # ==================================================
          # Best mtry
          # ==================================================
          
          best_mtry = rf_model$bestTune$mtry
          
          
          # ==================================================
          # Independent TEST-set probabilities
          # ==================================================
          
          prob = predict(
            rf_model,
            newdata = x_test,
            type = "prob"
          )
          
          
          # ==================================================
          # Independent TEST-set class predictions
          # ==================================================
          
          pred_class = predict(
            rf_model,
            newdata = x_test,
            type = "raw"
          )
          
          
          pred_class = factor(
            pred_class,
            levels = c("Low", "High")
          )
          
          
          # ==================================================
          # AUC
          #
          # AUC requires both classes in the independent
          # test set.
          # ==================================================
          
          if (length(unique(y_test)) == 2) {
            
            roc_obj = ROCR::prediction(
              predictions = prob[, "Low"],
              labels = as.integer(
                y_test == "Low"
              )
            )
            
            
            AUC = ROCR::performance(
              roc_obj,
              measure = "auc"
            )@y.values[[1]]
            
          } else {
            
            AUC = NA_real_
          }
          
          
          # ==================================================
          # Confusion matrix
          #
          # Low = positive class
          # ==================================================
          
          conf = confusionMatrix(
            data = pred_class,
            reference = y_test,
            positive = "Low",
            mode = "everything"
          )
          
          
          # ==================================================
          # Extract performance metrics
          #
          # Some metrics can be NA when the independent
          # test set contains only one class.
          # ==================================================
          
          F1 = unname(
            conf$byClass["F1"]
          )
          
          Sensitivity = unname(
            conf$byClass["Sensitivity"]
          )
          
          Specificity = unname(
            conf$byClass["Specificity"]
          )
          
          Balanced_Accuracy = unname(
            conf$byClass["Balanced Accuracy"]
          )
          
          
          # ==================================================
          # Metrics for this RF run
          # ==================================================
          
          metrics = data.frame(
            
            iteration = iter,
            
            metabolite = metabolite_names[mt],
            
            seed = base_seed,
            
            n_features = length(features),
            
            best_mtry = best_mtry,
            
            AUC = AUC,
            
            F1 = F1,
            
            Sensitivity = Sensitivity,
            
            Specificity = Specificity,
            
            Balanced_Accuracy = Balanced_Accuracy
          )
          
          
          # ==================================================
          # Store independent test predictions
          # ==================================================
          
          test_predictions = data.frame(
            
            iteration = iter,
            
            metabolite = metabolite_names[mt],
            
            best_mtry = best_mtry,
            
            rowIndex = seq_len(nrow(testset)),
            
            observed = y_test,
            
            predicted = pred_class,
            
            prob_Low = prob[, "Low"],
            
            prob_High = prob[, "High"]
          )
          
          
          # ==================================================
          # Return compact results
          #
          # Do not store 6000 complete RF models because this
          # would consume a large amount of RAM.
          # ==================================================
          
          list(
            
            metrics = metrics,
            
            predictions = test_predictions
          )
        }
      )
      
      
      names(
        metabolite_results
      ) = metabolite_names
      
      
      metabolite_results
    }
  )
  
}, finally = {
  
  stopCluster(cl)
  
  registerDoSEQ()
  
  cat(
    "\nParallel cluster stopped.\n"
  )
})


# ============================================================
# 9. Extract all performance metrics
# ============================================================

metrics_rts = do.call(
  rbind,
  lapply(
    results_rts,
    function(iter_result) {
      
      do.call(
        rbind,
        lapply(
          iter_result,
          function(x) {
            x$metrics
          }
        )
      )
    }
  )
)


rownames(metrics_rts) = NULL


# Inspect
head(metrics_rts)


# ============================================================
# 10. Extract all independent test-set predictions
# ============================================================

predictions_rts = do.call(
  rbind,
  lapply(
    results_rts,
    function(iter_result) {
      
      do.call(
        rbind,
        lapply(
          iter_result,
          function(x) {
            x$predictions
          }
        )
      )
    }
  )
)


rownames(predictions_rts) = NULL


# Inspect
head(predictions_rts)


# ============================================================
# 11. Save raw results
# ============================================================

write.csv(
  metrics_rts,
  "RF_RTS_1000_metrics.csv",
  row.names = FALSE
)


saveRDS(
  metrics_rts,
  "RF_RTS_1000_metrics.rds"
)


write.csv(
  predictions_rts,
  "RF_RTS_1000_test_predictions.csv",
  row.names = FALSE
)


saveRDS(
  results_rts,
  "RF_RTS_1000_results.rds"
)


# ============================================================
# 12. Summary over the 1000 RF repetitions
# ============================================================

summary_rts = metrics_rts %>%
  
  group_by(metabolite) %>%
  
  summarise(
    
    # AUC
    AUC_mean = mean(
      AUC,
      na.rm = TRUE
    ),
    
    AUC_sd = sd(
      AUC,
      na.rm = TRUE
    ),
    
    # F1
    F1_mean = mean(
      F1,
      na.rm = TRUE
    ),
    
    F1_sd = sd(
      F1,
      na.rm = TRUE
    ),
    
    # Sensitivity
    Sensitivity_mean = mean(
      Sensitivity,
      na.rm = TRUE
    ),
    
    Sensitivity_sd = sd(
      Sensitivity,
      na.rm = TRUE
    ),
    
    # Specificity
    Specificity_mean = mean(
      Specificity,
      na.rm = TRUE
    ),
    
    Specificity_sd = sd(
      Specificity,
      na.rm = TRUE
    ),
    
    # Balanced Accuracy
    Balanced_Accuracy_mean = mean(
      Balanced_Accuracy,
      na.rm = TRUE
    ),
    
    Balanced_Accuracy_sd = sd(
      Balanced_Accuracy,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )


print(summary_rts)


write.csv(
  summary_rts,
  "RF_RTS_1000_summary.csv",
  row.names = FALSE
)


saveRDS(
  summary_rts,
  "RF_RTS_1000_summary.rds"
)