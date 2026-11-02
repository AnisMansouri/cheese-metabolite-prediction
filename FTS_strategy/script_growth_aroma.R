library(corrplot)
library(xlsx)

growth = read.csv("data/growth.csv", row.names = 1, check.names = FALSE)
aroma = read.csv("data/EcoStab_Aroma_class_table.csv", sep = ";", row.names = 1,
                  check.names = FALSE)


# Experimental unit: the same replicate in the same group across four days.
# E.g., CT7_a, CT14_a, CT24_a, CT31_a are unit CT_a.
ids = colnames(growth)
unit = sub("^(CT|DT|GT)(7|14|24|31)_([abc])$", "\\1_\\3", ids)


# Rank-transform x and y and remove the experimental-unit means.
# Approximate p-value uses n - number of units - 1 degrees of freedom.
# This estimates association *within units*.
rank_rmcor = function(x, y, unit) {
  x = rank(as.numeric(x), ties.method = "average")
  y = rank(as.numeric(y), ties.method = "average")
  x = x - ave(x, unit, FUN = mean)
  y = y - ave(y, unit, FUN = mean)
  if (sum(x^2) == 0 || sum(y^2) == 0) return(c(r = NA, p = NA))
  r = sum(x * y) / sqrt(sum(x^2) * sum(y^2))
  df = length(x) - length(unique(unit)) - 1
  t_stat = r * sqrt(df / max(1e-14, 1 - r^2))
  c(r = r, p = 2 * pt(-abs(t_stat), df = df))
}

r_mat = matrix(NA_real_, nrow(growth), nrow(aroma),
                dimnames = list(rownames(growth), rownames(aroma)))
p_mat = r_mat

for (i in seq_len(nrow(growth))) {
  for (j in seq_len(nrow(aroma))) {
    s = rank_rmcor(as.numeric(growth[i, ]),
                    as.numeric(aroma[j, ]), unit)
    r_mat[i, j] = unname(s["r"])
    p_mat[i, j] = unname(s["p"])
  }
}

# BH correction is global: nine microorganisms x six classes = 54 tests.
p_fdr = matrix(p.adjust(as.vector(p_mat), method = "BH"),
                nrow = nrow(p_mat), ncol = ncol(p_mat),
                dimnames = dimnames(p_mat))

# Ga is the new name of the species 'Aa'
rownames(r_mat)[5] <- "Ga"
rownames(p_fdr)[5] <- "Ga"
row_order = c(5, 7, 9, 6, 4, 8, 2, 3, 1)
r_mat = r_mat[row_order, , drop = FALSE]
p_fdr = p_fdr[row_order, , drop = FALSE]
colnames(r_mat)[6] = "Sulphur compounds"
colnames(p_fdr)[6] = "Sulphur compounds"


corrplot(corr = r_mat, method = "circle", order = "original",
         col = colorRampPalette(c("blue", "white", "red"))(200),
         p.mat = p_fdr, sig.level = 0.05,
         insig = "pch", pch = 4,
         tl.col = "black", cl.cex = 0.7, tl.cex = 0.5)


#write.xlsx(r_mat, "growth_metabolites_correlations.xlsx", sheetName = "Correlation", append = FALSE)
#write.xlsx(p_fdr, "growth_metabolites_correlations.xlsx", sheetName = "FDR", append = TRUE)
