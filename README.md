# Predicting cheese metabolite levels from microbial gene expression

This repository contains the R scripts and supporting data used to investigate associations between microbial gene expression and the levels of six classes of volatile metabolites during cheese ripening. The main analysis combines **Elastic Net (EN) feature selection** with **Random Forest (RF) classification**, followed by biological interpretation of the selected KEGG Orthology (KO) features.

The repository accompanies the manuscript *Estimation of metabolite levels in cheese from microbial gene expression* (publication details and DOI to be added).

## Study overview

The analysis considers six metabolite classes:

- Alcohols
- Aldehydes
- Alkanes
- Ketones
- Esters
- Sulphur compounds

Each class is treated as a separate **Low/High binary classification task**. The training experiment (**EcoStab**) contains 34 samples. A second experiment (**ExEco**, 24 samples) is used for external evaluation in the reduced-feature strategy.

Two complementary analysis strategies are provided:

| Strategy | Predictor space | Model evaluation | Main purpose |
| --- | --- | --- | --- |
| **FTS — Full Training Set** | All available training-set KO predictors (4,245 in the supplied matrix) | Leave-one-out cross-validation (LOOCV) on the 34 training samples | Identify a broad set of potentially informative genes for biological interpretation |
| **RTS — Reduced Training Set** | KOs shared between the EcoStab and ExEco datasets (1,680 in the supplied reduced matrix) | Evaluation on 24 independent ExEco samples | Assess generalization to a separate experiment using a compatible feature space |

## Repository structure

```text
.
├── processing/
│   ├── training_set/
│   │   ├── script_processing.R
│   │   ├── annotation/
│   │   └── [EcoStab count, annotation, and sample metadata files]
│   ├── test_set/
│   │   ├── script_processing.R
│   │   └── [ExEco count, annotation, and sample metadata files]
│   └── binarization/
│       ├── script_binarization.R
│       └── data/
├── FTS_strategy/
│   ├── data/
│   ├── script_feature_selection.R
│   ├── script_rf_fts.R
│   ├── script_eigengene.R
│   ├── script_ora.R
│   ├── script_microbe_contribution.R
│   ├── script_growth_aroma.R
│   └── readme.txt
└── RTS_strategy/
    ├── data/
    ├── script_feature_selection.R
    ├── script_rf_rts.R
    ├── script_eigengene.R
    ├── script_ora.R
    ├── script_microbe_contribution.R
    └── readme.txt
```

The `data/` directories under the two strategy folders contain precomputed RDS objects, including normalized expression matrices, metabolite class labels, and selected feature sets. These allow users to start from intermediate data without rerunning every upstream processing step. Each strategy also has its own short `readme.txt` listing its analysis scripts.

## Methods

### 1. Expression preprocessing and binary labels

The scripts in `processing/training_set/` and `processing/test_set/` process sequencing count matrices and KO annotations to construct expression matrices. The preprocessing includes normalization and aggregation of KO-annotated features.

`processing/binarization/script_binarization.R` constructs Low/High metabolite labels using two-cluster k-means on the **EcoStab training measurements**; thresholds derived from the training experiment are then applied to the ExEco data. The six metabolite classes are modeled separately. In the supplied analysis, alcohols and aldehydes have identical binary labels and consequently represent the same classification target, despite differences in their continuous measurements.

### 2. Elastic Net feature selection

For each metabolite class, EN hyperparameters (`alpha` and `lambda`) are tuned with `caret::train(method = "glmnet")`, using LOOCV on the 34 training samples and maximizing **Cohen's Kappa**. The selected EN model is then refitted on all 34 training samples; genes with nonzero coefficients define the corresponding fixed feature signature.

The FTS and RTS strategies perform this procedure independently in their respective predictor spaces:

- `FTS_strategy/script_feature_selection.R`
- `RTS_strategy/script_feature_selection.R`

The feature lists are saved as `data/selected_features.rds` within each strategy folder.

### 3. Random Forest classification

Each RF contains **1,000 decision trees**. For each metabolite class, the `mtry` parameter (number of candidate predictors considered at each split) is optimized using LOOCV on the 34 training samples, with Cohen's Kappa as the selection criterion. This involves training on 33 samples and predicting the remaining sample in each fold.

- **FTS:** the fixed EN-selected signature is used for RF LOOCV; pooled held-out predictions are used to calculate classification metrics.
- **RTS:** RF hyperparameters are tuned on the EcoStab training samples, then the fitted model is evaluated on the independent ExEco dataset using the common KO features.

RF training, tuning and evaluation are repeated **1,000 times** with different random seeds to characterize variability associated with the stochastic learning algorithm. Predictions and performance metrics are recorded for each repetition and summarized across repetitions. These repetitions do **not** constitute 1,000 independent biological validation datasets.

The RF scripts export per-repetition metrics, prediction tables, and summary files (for example, `RF_FTS_1000_metrics.csv` and `RF_RTS_1000_metrics.csv`) to their respective working directories.

### 4. Biological interpretation and additional analyses

| Script | Purpose |
| --- | --- |
| `script_eigengene.R` | Compute eigengenes for EN-selected genes and compare their correlations with continuous metabolite measurements against randomly sampled gene sets. |
| `script_ora.R` | Perform KO-based over-representation analysis using `clusterProfiler`. |
| `script_microbe_contribution.R` | Estimate microbial species' contributions to selected KO features. |
| `FTS_strategy/script_growth_aroma.R` | Examine associations between microbial growth/abundance and metabolite classes and produce a correlation plot. |

## Requirements

The analysis is written in **R**. The scripts use the following packages:

```r
install.packages(c(
  "caret", "glmnet", "randomForest", "ROCR",
  "doParallel", "foreach", "dplyr", "stringr", "corrplot",
  "xlsx", "WGCNA"
))

# Bioconductor packages:
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("edgeR", "clusterProfiler"))
```

`parallel` is included with base R. The `xlsx` package may require a working Java installation. Package versions used to generate the published results have **not yet been pinned**; adding an `renv.lock` file or an archived `sessionInfo()` is recommended for full reproducibility.

## Running the analysis

**Working-directory convention:** The scripts use relative paths such as `data/ecostab_norm.rds`. Run each script from its own folder, rather than from the repository root.

### Reproduce the analyses from the included intermediate data

From the repository root:

```bash
cd FTS_strategy
Rscript script_feature_selection.R
Rscript script_rf_fts.R
# Additional analyses as needed:
Rscript script_ora.R
Rscript script_microbe_contribution.R

cd ../RTS_strategy
Rscript script_feature_selection.R
Rscript script_rf_rts.R
# Additional analyses as needed:
Rscript script_ora.R
Rscript script_microbe_contribution.R
```

The RF scripts perform 1,000 repetitions and can require substantial compute time and memory. They use parallel processing; inspect their cluster and core settings before running on shared computing resources.

**Manual steps:** The two `script_eigengene.R` files currently require manual selection of the metabolite index (`met`) to inspect results for all classes. They do not yet produce one consolidated, fully automated results table. `FTS_strategy/script_growth_aroma.R` can be run separately from within `FTS_strategy/` to generate its correlation plot.

## Citation

If you use this code, please cite the accompanying manuscript once its bibliographic details are available:

> *Estimation of metabolite levels in cheese from microbial gene expression*. [Authors, journal, year, DOI — to be completed.]

## Data provenance and reuse

The repository includes count matrices, annotations, metabolite data and/or intermediate processed matrices supplied as part of the research workflow. **Data redistribution permissions, third-party attribution, and an explicit software license should be confirmed before making the repository public.** No reuse license is asserted by this README.
