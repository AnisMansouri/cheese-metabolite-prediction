script_processing.R performs:

1. Normalization of the raw read counts raw_counts.csv.
2. Merging of orthologous genes and creation of the data set execo_norm.rds to be used for downstream analysis.

Execo_Aroma_class_table.csv is the latest and updated volatile compounds data shared by our colleagues from INRAE. It is converted to .rds file called metabolites_execo.rds used for the downstream analysis.

Aroma_mapping.csv contains the classifications of the detected compounds into 6 classes.

Execo_Geneid_annotations.csv contains genes mapping to GO and KEGG IDs, as shared by our colleagues from INRAE.
