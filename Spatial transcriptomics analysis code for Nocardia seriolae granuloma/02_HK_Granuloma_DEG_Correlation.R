# ==============================================================================
# Head Kidney Granuloma - 2. Differential Gene Expression & Correlation
# ==============================================================================

# --- 1. Find all marker genes for each cluster ---
# Identify markers for every cluster compared to all remaining cells
TS.markers <- FindAllMarkers(
  TS_fliter, 
  only.pos = TRUE,       # Only return positive markers (upregulated)
  min.pct = 0.1,         # Genes must be detected in at least 10% of cells in either of the two populations
  test.use = "wilcox",   # Statistical test to use
  logfc.threshold = 0.1  # Filter out genes with log2FC < 0.1
)

# Filter for significant markers based on adjusted p-value
significant.markers <- TS.markers[TS.markers$p_val_adj < 0.1, ]

# Extract the top 20 marker genes per cluster based on log2FC
top20_markers <- TS.markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)

# --- 2. Find markers between specific clusters ---
# Example: Inner E-Mac vs Outer E-Mac
OE_Mac_vs_IE_Mac.markers <- FindMarkers(
  TS_fliter, 
  ident.1 = "Inner E-Mac", 
  ident.2 = "Outer E-Mac", 
  min.pct = 0.1
)

# --- 3. Calculate Average Gene Expression per Cluster ---
# Average expression using normalized data
cluster_avg_exp <- AverageExpression(TS_fliter, assays = "SCT", layer = "data")

# Average expression using raw counts
cluster_avg_exp_counts <- AverageExpression(TS_fliter, assays = "SCT", layer = "counts")

# Log1p transformation of the average counts for better distribution/comparison
log_avg_counts <- log1p(cluster_avg_exp_counts$SCT)

# --- 4. Subset Average Expression Matrix for Significant Markers ---
# Ensure gene names match between the significant markers and the expression matrix
common_genes <- intersect(rownames(significant.markers), rownames(cluster_avg_exp$SCT))

# Retain only the matching genes
if (length(common_genes) > 0) {
  diff_gene_exp <- cluster_avg_exp$SCT[common_genes, ]
} else {
  print("No matching genes found, please check the data.")
}

# --- 5. Cluster Correlation Heatmap ---
# Calculate the Pearson correlation matrix between clusters
cor_matrix <- cor(as.matrix(cluster_avg_exp$SCT))

library(pheatmap)

# Define color gradient: lightblue -> white -> red (100 color scale)
my_colors <- colorRampPalette(c("lightblue", "white", "red"))(100)

# Define color breaks to map white to the middle of your expected correlation range
# Assuming correlations range roughly from 0.5 to 1.0 based on your original breaks
my_breaks <- seq(0.5, 1, length.out = 101)

# Generate the heatmap
p_cor_heatmap <- pheatmap(
  cor_matrix,
  color = my_colors,           
  breaks = my_breaks,          
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "Cluster Correlation Heatmap (SCT Data)"
)

# Print the heatmap to the plots pane
print(p_cor_heatmap)