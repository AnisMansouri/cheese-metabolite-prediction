library(caret)
library(randomForest)
library(ROCR)
library(parallel)
library(doParallel)
library(foreach)


fts = readRDS("data/ecostab_norm.rds")
ecostab_classes = readRDS("data/ecostab_classes.rds")
selected_features = readRDS("data/selected_features.rds")


# ============================================================
# 1. Metabolite names
# ============================================================

metabolite_names <- c(
  "Alcohols",
  "Aldehydes",
  "Alkanes",
  "Ketones",
  "Esters",
  "Sulphur compounds"
)


# ============================================================
# 2. Build mtry grids automatically from the number
#    of EN-selected predictors
# ============================================================

make_mtry_grid <- function(p) {
  
  if (p <= 10) {
    return(1:p)
  }
  
  candidates <- c(
    1, 2, 3, 5, 7, 10,
    round(sqrt(p)),
    15, 20, 30, 50, 75, 100
  )
  
  candidates <- candidates[candidates <= p]
  
  sort(unique(candidates))
}


mtry_list <- lapply(
  selected_features,
  function(features) {
    make_mtry_grid(length(features))
  }
)

names(mtry_list) <- metabolite_names


# Inspect grids
mtry_list


# ============================================================
# 3. Function for caret seeds
#
# For LOOCV with 34 observations:
#
# seeds[[1]]  ... seeds[[34]]
#       -> seeds for tuning models
#
# seeds[[35]]
#       -> seed for final RF model
# ============================================================

make_caret_seeds <- function(n_resamples,
                             n_models,
                             base_seed) {
  
  set.seed(base_seed)
  
  seeds <- vector(
    mode = "list",
    length = n_resamples + 1
  )
  
  for (i in seq_len(n_resamples)) {
    
    seeds[[i]] <- sample.int(
      1000000,
      size = n_models,
      replace = FALSE
    )
  }
  
  # final model
  seeds[[n_resamples + 1]] <- sample.int(
    1000000,
    size = 1
  )
  
  seeds
}


# ============================================================
# 4. Basic checks
# ============================================================

stopifnot(
  length(selected_features) == 6,
  length(ecostab_classes) == 6
)

stopifnot(
  all(
    sapply(
      ecostab_classes,
      length
    ) == nrow(fts)
  )
)


# Check features are present in FTS
for (mt in seq_along(selected_features)) {
  
  missing_features <- setdiff(
    selected_features[[mt]],
    colnames(fts)
  )
  
  if (length(missing_features) > 0) {
    
    stop(
      paste(
        "Missing predictors for",
        metabolite_names[mt],
        ":",
        paste(missing_features, collapse = ", ")
      )
    )
  }
}


# ============================================================
# 5. SLURM parallel backend
# ============================================================

ncores <- suppressWarnings(
  as.integer(
    Sys.getenv("SLURM_CPUS_PER_TASK")
  )
)

if (is.na(ncores) || ncores < 1) {
  
  ncores <- suppressWarnings(
    as.integer(
      Sys.getenv("SLURM_CPUS_ON_NODE")
    )
  )
}

# fallback outside SLURM
if (is.na(ncores) || ncores < 1) {
  
  ncores <- parallel::detectCores()
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


cl <- makePSOCKcluster(ncores)

registerDoParallel(cl)


cat(
  "Registered parallel workers:",
  foreach::getDoParWorkers(),
  "\n"
)


# ============================================================
# 6. Number of RF repetitions
# ============================================================

n_iterations <- 1000


# ============================================================
# 7. Run RF analysis
# ============================================================

results_fts <- tryCatch({
  
  lapply(
    seq_len(n_iterations),
    function(iter) {
      
      cat(
        "\n==============================\n",
        "Iteration:", iter, "\n",
        "==============================\n"
      )
      
      
      metabolite_results <- lapply(
        seq_along(selected_features),
        function(mt) {
          
          cat(
            metabolite_names[mt],
            "\n"
          )
          
          
          # ==================================================
          # Predictors selected by Elastic Net
          # ==================================================
          
          features <- selected_features[[mt]]
          
          x_train <- fts[
            ,
            features,
            drop = FALSE
          ]
          
          
          # ==================================================
          # Outcome
          #
          # ecostab_classes contains character vectors,
          # therefore convert explicitly to factor.
          # ==================================================
          
          y_train <- factor(
            ecostab_classes[[mt]],
            levels = c("Low", "High")
          )
          
          
          if (any(is.na(y_train))) {
            
            stop(
              paste(
                "Unexpected class label for",
                metabolite_names[mt]
              )
            )
          }
          
          
          # ==================================================
          # mtry grid
          # ==================================================
          
          tune_grid <- expand.grid(
            mtry = mtry_list[[mt]]
          )
          
          
          # ==================================================
          # Reproducible seeds
          #
          # Unique seed for every:
          #
          # iteration x metabolite
          # ==================================================
          
          base_seed <- iter * 100 + mt
          
          
          caret_seeds <- make_caret_seeds(
            n_resamples = nrow(fts),
            n_models = nrow(tune_grid),
            base_seed = base_seed
          )
          
          
          set.seed(base_seed)
          
          
          # ==================================================
          # Train Random Forest
          #
          # mtry is selected according to Kappa
          # using LOOCV.
          #
          # savePredictions = "final" means:
          #
          # rf_model$pred contains the LOO predictions
          # corresponding to the BEST mtry only.
          # ==================================================
          
          rf_model <- train(
            
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
          
          best_mtry <- rf_model$bestTune$mtry
          
          
          # ==================================================
          # LOO-CV predictions from BEST mtry
          # ==================================================
          
          cv_pred <- rf_model$pred
          
          
          # Confirm that only the best mtry is present
          if ("mtry" %in% colnames(cv_pred)) {
            
            if (
              !all(
                cv_pred$mtry == best_mtry
              )
            ) {
              
              stop(
                "rf_model$pred contains predictions from more than the best mtry"
              )
            }
          }
          
          
          # ==================================================
          # Ensure identical factor levels
          # ==================================================
          
          cv_pred$obs <- factor(
            cv_pred$obs,
            levels = c("Low", "High")
          )
          
          cv_pred$pred <- factor(
            cv_pred$pred,
            levels = c("Low", "High")
          )
          
          
          # ==================================================
          # AUC
          #
          # Low = positive class
          #
          # Each probability is produced for the observation
          # left out during LOO-CV.
          # ==================================================
          
          roc_obj <- ROCR::prediction(
            predictions = cv_pred$Low,
            labels = as.integer(
              cv_pred$obs == "Low"
            )
          )
          
          
          AUC <- ROCR::performance(
            roc_obj,
            measure = "auc"
          )@y.values[[1]]
          
          
          # ==================================================
          # Confusion matrix from pooled LOO predictions
          # ==================================================
          
          conf <- confusionMatrix(
            
            data = cv_pred$pred,
            
            reference = cv_pred$obs,
            
            positive = "Low",
            
            mode = "everything"
          )
          
          
          # ==================================================
          # Extract performance metrics
          # ==================================================
          
          F1 <- unname(
            conf$byClass["F1"]
          )
          
          Sensitivity <- unname(
            conf$byClass["Sensitivity"]
          )
          
          Specificity <- unname(
            conf$byClass["Specificity"]
          )
          
          Balanced_Accuracy <- unname(
            conf$byClass["Balanced Accuracy"]
          )
          
          
          # ==================================================
          # Metrics for this RF replicate
          # ==================================================
          
          metrics <- data.frame(
            
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
          # Store LOO predictions too
          # ==================================================
          
          cv_predictions <- data.frame(
            
            iteration = iter,
            
            metabolite = metabolite_names[mt],
            
            best_mtry = best_mtry,
            
            rowIndex = cv_pred$rowIndex,
            
            observed = cv_pred$obs,
            
            predicted = cv_pred$pred,
            
            prob_Low = cv_pred$Low,
            
            prob_High = cv_pred$High
          )
          
          
          # ==================================================
          # Return
          #
          # Deliberately do NOT store rf_model:
          # 6000 complete 1000-tree models would consume
          # a very large amount of RAM.
          # ==================================================
          
          list(
            
            metrics = metrics,
            
            predictions = cv_predictions
          )
        }
      )
      
      
      names(
        metabolite_results
      ) <- metabolite_names
      
      
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
# 8. Extract all performance metrics
# ============================================================

metrics_fts <- do.call(
  rbind,
  lapply(
    results_fts,
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


rownames(metrics_fts) <- NULL


# Inspect
head(metrics_fts)


# ============================================================
# 9. Extract all LOO-CV predictions
# ============================================================

predictions_fts <- do.call(
  rbind,
  lapply(
    results_fts,
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


rownames(predictions_fts) <- NULL


# Inspect
head(predictions_fts)


# ============================================================
# 10. Save raw results
# ============================================================

write.csv(
  metrics_fts,
  "RF_FTS_1000_metrics.csv",
  row.names = FALSE
)

saveRDS(
  metrics_fts,
  "RF_FTS_1000_metrics.rds"
)

write.csv(
  predictions_fts,
  "RF_FTS_1000_LOOCV_predictions.csv",
  row.names = FALSE
)


saveRDS(
  results_fts,
  "RF_FTS_1000_results.rds"
)


# ============================================================
# 11. Summary over the 1000 RF repetitions
# ============================================================

library(dplyr)


summary_fts <- metrics_fts %>%
  
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


print(summary_fts)


write.csv(
  summary_fts,
  "RF_FTS_1000_summary.csv",
  row.names = FALSE
)

saveRDS(
  summary_fts,
  "RF_FTS_1000_summary.rds"
)
