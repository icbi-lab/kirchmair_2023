
# R script to plot heatmaps of the NGSCheckMate results

library(datamisc)

cor <- read.delim("data/rnaseq/MEM/samplecheck/output_corr_matrix.txt", row.names = 1)
dimnames(cor) <- dimnames(cor) %L>% gsub(x = ., ".markdup.sorted", "") %L>% make.names()
cor <- as.matrix(cor)
diag(cor) <- NA
ids <- c("TN1","TCM1","TEM1","TSCM1","TN2","TCM2","TEM2","TSCM2","TN3","TCM3","TEM3","TSCM3")
cxheatmap(cor[ids,ids], cluster_rows = FALSE, cluster_cols = FALSE, na_col = "grey35",
          heatpal = circlize::colorRamp2(c(0,0.5,0.8), c("#ffffff", "#edebf5", "#3700ff")) ) %>%
  saveplot("data/rnaseq/MEM/samplecheck/samplecheck_heatmap_mem.png")


cor <- read.delim("data/rnaseq/EXH/samplecheck/output_corr_matrix.txt", row.names = 1)
dimnames(cor) <- dimnames(cor) %L>% gsub(x = ., ".markdup.sorted", "") %L>% make.names()
cor <- as.matrix(cor)
diag(cor) <- NA
ids <- c("TEFF1a","TEFF1b","TEX1a","TEX1b","TEFF2","TEX2","TEFF3","TEX3")
cxheatmap(cor[ids,ids], cluster_rows = FALSE, cluster_cols = FALSE, na_col = "grey35",
          heatpal = circlize::colorRamp2(c(0,0.5,0.8), c("#ffffff", "#edebf5", "#3700ff")) ) %>%
  saveplot("data/rnaseq/EXH/samplecheck/samplecheck_heatmap_exh.png")
