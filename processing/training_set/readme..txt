script_processing.R performs:

1. Filtering and normalization of the raw read counts counts_matrix_2020-10-01.csv.
2. Merging of orthologous genes and creation of the data set ecostab_norm.rds to be used for downstream analysis.

20230419_EcoStab_Aroma_class_table.csv is the latest and updated volatile compounds data sent by our colleagues from INRAE. It is converted to .rds file called metabolites_ecostab.rds used for the downstream analysis.

Aroma_mapping.csv contains the classifications of the detected compounds into 6 classes.

annotation folder contains genes mapping to GO and KEGG IDs, as shared by our colleagues from INRAE.
