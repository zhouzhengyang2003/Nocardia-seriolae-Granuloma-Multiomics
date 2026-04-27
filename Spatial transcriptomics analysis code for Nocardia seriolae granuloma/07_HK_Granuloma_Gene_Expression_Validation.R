# ==============================================================================
# Head Kidney Granuloma - 7. Gene Expression Validation & Statistical Plotting
# ==============================================================================
suppressMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# Set base saving directory
save_dir <- "./results/Gene_Expression_Plots/"
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Part 1: Boxplot & Stacked Violin for Profibrotic Signals (All Clusters)
# ==============================================================================
cat("Generating plots for all selected clusters...\n")

target_genes <- c("LOC119888537", "bmp1a", "p4ha1b", "loxa", "LOC119897617", "LOC119886711")
target_clusters <- c("Eryt", "Eryt-Gran", "Gran", "Non-G Mac", "B", "T-B-APCs", "Outer E-Mac", "Inner E-Mac", "Fibro-Endo")

DefaultAssay(TS_fliter) <- "SCT"
valid_genes <- intersect(target_genes, rownames(TS_fliter))
cells <- Cells(TS_fliter)

# Extract Expression Matrix: SCT@data (log1p corrected counts)
expr_mat <- LayerData(TS_fliter, assay="SCT", layer="data")[valid_genes, cells, drop=FALSE]

# Melt into a long dataframe
df_all <- as.data.frame(t(expr_mat)) |>
  mutate(ident = Idents(TS_fliter)[cells]) |>
  filter(ident %in% target_clusters) |>
  pivot_longer(all_of(valid_genes), names_to = "gene", values_to = "expr") |>
  mutate(
    ident = factor(ident, levels = target_clusters),
    gene  = factor(gene, levels = valid_genes)
  )

# Failsafe: Ensure we are using log-normalized data, not raw counts
stopifnot(max(df_all$expr, na.rm = TRUE) < 20)

# --- 1A. Grouped Boxplot ---
dodge <- position_dodge(width = 0.78)

p_box_all <- ggplot(df_all, aes(x = ident, y = expr, fill = gene)) +
  geom_boxplot(position = dodge, width = 0.64, outlier.shape = NA, color = "grey25", linewidth = 0.5) +
  geom_point(aes(color = gene), position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.78),
             alpha = 0.40, size = 1.3, shape = 15, stroke = 0) +
  scale_fill_brewer(palette = "Set3") +
  scale_color_brewer(palette = "Set3", guide = "none") +
  labs(x = NULL, y = "SCT data (log1p corrected counts)", fill = "Gene") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x      = element_text(angle = 20, hjust = 1, face = "bold"),
    axis.text.y      = element_text(face = "bold"),
    legend.position  = "top",
    legend.title     = element_text(face = "bold"),
    panel.border     = element_rect(color = "grey85", fill = NA),
    axis.line        = element_line(color = "black")
  ) +
  coord_cartesian(ylim = c(0, max(df_all$expr, na.rm = TRUE) * 1.05))

ggsave(file.path(save_dir, "1_Genes_byGroup_Boxplot.pdf"), p_box_all, width = 8.5, height = 5)

# --- 1B. Stacked Facet Violin Plot ---
p_vln_stacked <- ggplot(df_all, aes(x = ident, y = expr, fill = ident)) +
  geom_violin(scale = "width", adjust = 1.5, trim = TRUE, color = "black", linewidth = 0.3) +
  stat_summary(fun = mean, geom = "point", size = 1.2, color = "white") +
  stat_summary(fun = mean, geom = "point", size = 0.6, color = "black") +
  facet_grid(gene ~ ., scales = "free_y") +
  scale_fill_brewer(palette = "Set3") + 
  scale_y_continuous(n.breaks = 3) + 
  labs(x = NULL, y = "Expression Level (SCT)") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold", color = "black"),
    axis.text.y      = element_text(size = 9, color = "black"), 
    axis.ticks.y     = element_line(color = "black"),
    strip.text.y     = element_text(angle = 0, face = "bold.italic", size = 11), 
    strip.background = element_blank(), 
    legend.position  = "none", 
    panel.spacing    = unit(0.1, "lines"), 
    axis.line        = element_line(color = "black")
  )

ggsave(file.path(save_dir, "2_Genes_byGroup_StackedVlnPlot.pdf"), p_vln_stacked, width = 8.5, height = 5)


# ==============================================================================
# Part 2: FindMarkers Validation for Core Genes
# ==============================================================================
cat("\nExecuting FindMarkers for structural/paracrine genes...\n")

# Test 1: Outer E-Mac vs Non-G Mac (Validating epithelioid transition)
res_vs_Mac <- FindMarkers(TS_fliter, ident.1 = "Outer E-Mac", ident.2 = "Non-G Mac", 
                          features = valid_genes, logfc.threshold = 0, min.pct = 0)
cat("=== Outer E-Mac vs Non-G Mac Results ===\n")
print(res_vs_Mac)

# Test 2: Outer E-Mac vs Fibro-Endo (Validating macrophage origin of enzymes)
res_vs_Fibro <- FindMarkers(TS_fliter, ident.1 = "Outer E-Mac", ident.2 = "Fibro-Endo", 
                            features = valid_genes, logfc.threshold = 0, min.pct = 0)
cat("=== Outer E-Mac vs Fibro-Endo Results ===\n")
print(res_vs_Fibro)


# ==============================================================================
# Part 3: Macrophage Subset Plotting with Custom FDR Stars
# ==============================================================================
cat("\nGenerating significant annotated plots for Macrophage subsets...\n")

sub_genes <- c("spaca4l", "LOC119906807", "wt1b", "gata6", "LOC119894306", "LOC119914475", "aldh1a2", "gprc5c")
sub_clusters <- c("Non-G Mac", "Outer E-Mac", "Inner E-Mac")
baseline <- "Non-G Mac"

DefaultAssay(object_sub) <- "SCT"
valid_sub_genes <- intersect(sub_genes, rownames(object_sub))
cells_sub <- Cells(object_sub)

expr_mat_sub <- LayerData(object_sub, assay="SCT", layer="data")[valid_sub_genes, cells_sub, drop=FALSE]

df_sub <- as.data.frame(t(expr_mat_sub)) |>
  mutate(ident = Idents(object_sub)[cells_sub]) |>
  filter(ident %in% sub_clusters) |>
  pivot_longer(all_of(valid_sub_genes), names_to = "gene", values_to = "expr") |>
  mutate(
    ident = factor(ident, levels = sub_clusters),
    gene  = factor(gene, levels = valid_sub_genes)
  )

# Failsafe
stopifnot(max(df_sub$expr, na.rm = TRUE) < 20)

# --- Statistical Calculation & Annotation ---
pair_tests <- df_sub |>
  group_by(gene) |>
  do({
    dd <- .
    base <- dd$expr[dd$ident == baseline]
    out  <- dd$expr[dd$ident == "Outer E-Mac"]
    inn  <- dd$expr[dd$ident == "Inner E-Mac"]
    tibble(
      gene = unique(dd$gene),
      ident = c("Outer E-Mac", "Inner E-Mac"), 
      contrast = c("Outer vs Non-G", "Inner vs Non-G"),
      p = c(wilcox.test(out, base, exact = FALSE)$p.value,
            wilcox.test(inn, base, exact = FALSE)$p.value)
    )
  }) |>
  ungroup() |>
  mutate(p.adj = p.adjust(p, method = "BH"))

p_to_star <- function(p) {
  ifelse(p < 1e-4, "****",
         ifelse(p < 1e-3, "***",
                ifelse(p < 0.01, "**",
                       ifelse(p < 0.05, "*", "ns"))))
}

pair_tests <- pair_tests |> mutate(label = p_to_star(p.adj))

# Calculate Y-position for stars
y_span <- diff(range(df_sub$expr, na.rm = TRUE))
y_bump <- ifelse(is.finite(y_span) && y_span > 0, 0.06 * y_span, 0.3)

ypos <- df_sub |>
  group_by(gene, ident) |>
  summarise(y = max(expr, na.rm = TRUE) + y_bump, .groups = "drop")

annot_df <- pair_tests |> left_join(ypos, by = c("gene", "ident"))

# --- 3A. Annotated Boxplot ---
p_box_sub <- ggplot(df_sub, aes(x = ident, y = expr, fill = gene)) +
  geom_boxplot(position = dodge, width = 0.64, outlier.shape = NA, color = "grey25", linewidth = 0.5) +
  geom_point(aes(color = gene), position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.78),
             alpha = 0.40, size = 1.3, shape = 15, stroke = 0) +
  scale_fill_brewer(palette = "Paired") +
  scale_color_brewer(palette = "Paired", guide = "none") +
  labs(x = NULL, y = "Log-Normalized Expression", fill = "Gene") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x      = element_text(angle = 20, hjust = 1, face = "bold"),
    axis.text.y      = element_text(face = "bold"),
    legend.position  = "top",
    legend.title     = element_text(face = "bold"),
    panel.border     = element_rect(color = "grey85", fill = NA),
    axis.line        = element_line(color = "black")
  ) +
  coord_cartesian(ylim = c(0, max(df_sub$expr, na.rm = TRUE) * 1.15)) +
  geom_text(data = annot_df, aes(x = ident, y = y, label = label, color = gene),
            position = dodge, vjust = 0, size = 4.25, show.legend = FALSE)

ggsave(file.path(save_dir, "3_Sub_Annotated_Boxplot.pdf"), p_box_sub, width = 6, height = 4)


# --- 3B. Annotated Violin Plot ---
p_violin_sub <- ggplot(df_sub, aes(x = ident, y = expr, fill = gene)) +
  geom_violin(position = dodge, width = 0.78, scale = "width", trim = TRUE, color = "grey25", linewidth = 0.5) +
  geom_point(aes(color = gene), position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.78),
             alpha = 0.40, size = 1.3, shape = 15, stroke = 0) +
  scale_fill_brewer(palette = "Set3") +
  scale_color_brewer(palette = "Set3", guide = "none") +
  labs(x = NULL, y = "Log-Normalized Expression", fill = "Gene") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x      = element_text(angle = 20, hjust = 1, face = "bold"),
    axis.text.y      = element_text(face = "bold"),
    legend.position  = "top",
    legend.title     = element_text(face = "bold"),
    panel.border     = element_rect(color = "grey85", fill = NA),
    axis.line        = element_line(color = "black")
  ) +
  coord_cartesian(ylim = c(0, max(df_sub$expr, na.rm = TRUE) * 1.15)) +
  geom_text(data = annot_df, aes(x = ident, y = y, label = label, color = gene),
            position = dodge, vjust = 0, size = 4.25, show.legend = FALSE)

ggsave(file.path(save_dir, "4_Sub_Annotated_Violin.pdf"), p_violin_sub, width = 6, height = 4)

cat("✅ All statistical plots generated and saved successfully!\n")