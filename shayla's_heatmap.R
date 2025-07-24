library(GEOquery)
library(limma)

gse <- getGEO("GSE93272", GSEMatrix = TRUE, AnnotGPL = TRUE)
gse_data <- gse[[1]]
exprs_data <- exprs(gse_data)
#View(exprs_data)
pheno_data <- pData(gse_data)
feature_data <- fData(gse_data)

View(pheno_data)

write.csv(exprs_data, "~/Downloads/Desktop/omics_project/week1/GSE272/datasets/GSE93272_expression_matrix.csv")
write.csv(pheno_data, "~/Downloads/Desktop/omics_project/week1/GSE272/datasets/GSE93272_sample_metadata.csv")
write.csv(feature_data, "~/Downloads/Desktop/omics_project/week1/GSE272/datasets/GSE93272_probe_annotation.csv")


# Remove probes with too many NAs or zero variance
exprs_data <- exprs_data[complete.cases(exprs_data), ]
exprs_data <- exprs_data[apply(exprs_data, 1, var) > 0.01, ]  # change threshold as needed
# Filter out low expression probes (e.g., below 25th percentile)

exprs_data <- log2(exprs_data + 1)

threshold <- quantile(rowMeans(exprs_data), 0.25)
exprs_data <- exprs_data[rowMeans(exprs_data) > threshold, ]



# Ensure annotation includes gene symbols
feature_data <- fData(gse_data)
head(feature_data$`Gene symbol`)  # Confirm it has gene symbols

# Combine expression data and annotation
exprs_annotated <- data.frame(probe_id = rownames(exprs_data),
                              exprs_data,
                              stringsAsFactors = FALSE)

# Add gene symbols
exprs_annotated$gene_symbol <- feature_data[rownames(exprs_data), "Gene symbol"]
exprs_annotated$gene_symbol <- sapply(strsplit(exprs_annotated$gene_symbol, "///", fixed = TRUE), `[`, 1)

write.csv(exprs_annotated, "~/Downloads/Desktop/omics_project/week1/GSE272/datasets/GSE272_Exprs_annotated.csv")


female_samples <- pheno_data$`gender:ch1` == "F"
exprs_female <- exprs_annotated[, female_samples]
pheno_female <- pheno_data[female_samples, ]


exprs_female_matrix <- exprs_female[, !(colnames(exprs_female) %in% c("probe_id", "gene_symbol"))]



group_female <- factor(pheno_female$`disease state:ch1`)
design_female <- model.matrix(~ 0 + group_female)
colnames(design_female) <- make.names(colnames(design_female))
contrast.matrix <- makeContrasts(RA_vs_Control = group_femaleRA - group_femalehealthy.control, levels = design_female)




fit <- lmFit(exprs_female_matrix, design_female)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

results <- topTable(fit2, adjust = "fdr", number = Inf)

results$probe_id <- rownames(results)
results_annotated <- merge(results, exprs_annotated[, c("probe_id", "gene_symbol")], by = "probe_id", all.x = TRUE)


#significant_DEGs <- results_annotated[results_annotated$adj.P.Val < 0.05 & abs(results_annotated$logFC) > 0.1, ]
write.csv(results_annotated, "~/Downloads/Desktop/omics_project/week1/GSE272/results_as_csv_272/OVERALL_FEMALE_272.csv", row.names = FALSE)




# Define thresholds
logFC_threshold <- 0.2
adjP_threshold <- 0.05

results_annotated$Regulation <- ifelse(
  results_annotated$adj.P.Val < adjP_threshold & results_annotated$logFC > logFC_threshold, "Upregulated",
  ifelse(results_annotated$adj.P.Val < adjP_threshold & results_annotated$logFC < -logFC_threshold, "Downregulated", "Not Significant")
)

table(results_annotated$Regulation)

write.csv(results_annotated, "~/Downloads/Desktop/omics_project/week1/GSE272/results_as_csv_272/OVERALL_FEMALE_272_with_regulation.csv", row.names = FALSE)
significant_genes <- results_annotated[results_annotated$Regulation != "Not Significant", ]


###### heatmap



library(pheatmap)

top_DEGs <- significant_genes[order(significant_genes$adj.P.Val), ]
top_DEGs <- top_DEGs[!duplicated(top_DEGs$gene_symbol), ]
top_DEGs <- top_DEGs[!is.na(top_DEGs$gene_symbol) & top_DEGs$gene_symbol != "", ]
top50_genes <- head(top_DEGs, 50)
heatmap_data <- exprs_female_matrix[top50_genes$probe_id, ]
rownames(heatmap_data) <- top50_genes$gene_symbol
heatmap_scaled <- t(scale(t(heatmap_data)))
annotation_col <- data.frame(Disease = group_female)
rownames(annotation_col) <- colnames(heatmap_scaled)

#pheatmap
pheatmap(
  heatmap_scaled,
  annotation_col = annotation_col,
  fontsize_row = 7,
  fontsize_col = 7,
  main = "Top 50 DEGs: RA vs Control (Females)",
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  border_color = NA
)






