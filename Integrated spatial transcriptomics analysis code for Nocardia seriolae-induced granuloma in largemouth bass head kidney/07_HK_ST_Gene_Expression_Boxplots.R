# ==============================================================================
# Head Kidney Granuloma - 7. Gene Expression Boxplots and Violin Plots
# ==============================================================================
#
# Purpose:
#   Extract SCT-normalized expression values for selected representative genes
#   and visualize their expression across annotated spatial domains.
#
# Required input objects:
#   TS_fliter:
#     Infected-only head kidney spatial transcriptomics object with annotated
#     domains, including Eryt, Gran, Non-G Mac, Outer E-Mac, Inner E-Mac and
#     Fibro-Endo.
#
#   object_sub:
#     Macrophage/E-Mac subset object containing:
#       Mac (Homeostatic), Mac (Infection-associated),
#       Outer E-Mac and Inner E-Mac.
#
# Main outputs generated in memory:
#   expression_plot_results
#   marker_test_results
#
# Notes:
#   This script does not save figures to disk. Users can export plots as needed.
# ==============================================================================


# ==============================================================================
# --- 0. Load R packages ---
# ==============================================================================

suppressMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(RColorBrewer)
})

options(future.globals.maxSize = 9e10)
options(Seurat.object.assay.version = "v5")


# ==============================================================================
# --- 1. Example input object loading ---
# ==============================================================================

project_dir <- "path/to/project"

infected_only_object_path <- file.path(
  project_dir,
  "processed_objects",
  "head_kidney_infected_only_ST_annotated.rds"
)

macrophage_subset_object_path <- file.path(
  project_dir,
  "processed_objects",
  "head_kidney_macrophage_EMac_subset.rds"
)

if (!exists("TS_fliter")) {
  TS_fliter <- readRDS(infected_only_object_path)
}

if (!exists("object_sub")) {
  object_sub <- readRDS(macrophage_subset_object_path)
}

DefaultAssay(TS_fliter) <- "SCT"
DefaultAssay(object_sub) <- "SCT"

message("TS_fliter identities:")
print(table(Idents(TS_fliter)))

message("object_sub identities:")
print(table(Idents(object_sub)))


# ==============================================================================
# --- 2. Helper functions ---
# ==============================================================================

extract_sct_expression_long <- function(
    object,
    genes,
    groups_keep,
    assay = "SCT",
    layer = "data"
) {
  
  DefaultAssay(object) <- assay
  
  genes_present <- intersect(genes, rownames(object))
  genes_missing <- setdiff(genes, genes_present)
  
  if (length(genes_missing) > 0) {
    warning(
      "The following genes were not found and will be skipped: ",
      paste(genes_missing, collapse = ", ")
    )
  }
  
  if (length(genes_present) == 0) {
    stop("None of the selected genes were found in the object.")
  }
  
  cells <- Cells(object)
  
  expr_mat <- LayerData(
    object = object,
    assay = assay,
    layer = layer
  )[genes_present, cells, drop = FALSE]
  
  ident_vec <- as.character(Idents(object)[cells])
  names(ident_vec) <- cells
  
  df <- as.data.frame(as.matrix(t(expr_mat))) %>%
    rownames_to_column("cell_id") %>%
    mutate(
      ident = ident_vec[cell_id]
    ) %>%
    filter(ident %in% groups_keep) %>%
    pivot_longer(
      cols = all_of(genes_present),
      names_to = "gene",
      values_to = "expr"
    ) %>%
    mutate(
      ident = factor(ident, levels = groups_keep),
      gene = factor(gene, levels = genes_present)
    )
  
  if (max(df$expr, na.rm = TRUE) >= 20) {
    warning(
      "Maximum expression value is >= 20. Please confirm that SCT data layer was used."
    )
  }
  
  return(df)
}


run_gene_pairwise_tests <- function(
    df,
    baseline,
    comparison_groups,
    p_adjust_method = "BH"
) {
  
  safe_wilcox <- function(x, y) {
    x <- x[!is.na(x)]
    y <- y[!is.na(y)]
    
    if (length(x) == 0 || length(y) == 0) {
      return(NA_real_)
    }
    
    wilcox.test(x, y, exact = FALSE)$p.value
  }
  
  p_to_star <- function(p) {
    ifelse(
      is.na(p), NA_character_,
      ifelse(
        p < 1e-4, "****",
        ifelse(
          p < 1e-3, "***",
          ifelse(
            p < 0.01, "**",
            ifelse(p < 0.05, "*", "ns")
          )
        )
      )
    )
  }
  
  pair_tests <- expand_grid(
    gene = unique(as.character(df$gene)),
    comparison_group = comparison_groups
  ) %>%
    rowwise() %>%
    mutate(
      baseline_group = baseline,
      contrast = paste0(comparison_group, " vs ", baseline_group),
      p = safe_wilcox(
        df$expr[
          as.character(df$gene) == gene &
            as.character(df$ident) == comparison_group
        ],
        df$expr[
          as.character(df$gene) == gene &
            as.character(df$ident) == baseline_group
        ]
      )
    ) %>%
    ungroup() %>%
    mutate(
      p.adj = p.adjust(p, method = p_adjust_method),
      p.adjust.method = p_adjust_method,
      significance_label = p_to_star(p.adj)
    )
  
  return(pair_tests)
}


make_gene_boxplot <- function(
    df,
    pair_tests = NULL,
    y_label = "SCT data (log1p corrected counts)",
    fill_palette = "Set3",
    add_points = TRUE
) {
  
  dodge <- position_dodge(width = 0.78)
  
  p <- ggplot(
    df,
    aes(x = ident, y = expr, fill = gene)
  ) +
    geom_boxplot(
      position = dodge,
      width = 0.64,
      outlier.shape = NA,
      color = "grey25",
      linewidth = 0.5
    ) +
    scale_fill_brewer(palette = fill_palette) +
    labs(
      x = NULL,
      y = y_label,
      fill = "Gene"
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(
        angle = 20,
        hjust = 1,
        face = "bold"
      ),
      axis.text.y = element_text(face = "bold"),
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      panel.border = element_rect(
        color = "grey85",
        fill = NA
      ),
      axis.line = element_line(color = "black")
    ) +
    coord_cartesian(
      ylim = c(0, max(df$expr, na.rm = TRUE) * 1.15)
    )
  
  if (add_points) {
    p <- p +
      geom_point(
        aes(color = gene),
        position = position_jitterdodge(
          jitter.width = 0.12,
          dodge.width = 0.78
        ),
        alpha = 0.40,
        size = 1.3,
        shape = 15,
        stroke = 0
      ) +
      scale_color_brewer(
        palette = fill_palette,
        guide = "none"
      )
  }
  
  if (!is.null(pair_tests)) {
    
    y_span <- diff(range(df$expr, na.rm = TRUE))
    y_bump <- ifelse(
      is.finite(y_span) && y_span > 0,
      0.06 * y_span,
      0.3
    )
    
    ypos <- df %>%
      group_by(gene, ident) %>%
      summarise(
        y = max(expr, na.rm = TRUE) + y_bump,
        .groups = "drop"
      )
    
    annot_df <- pair_tests %>%
      rename(ident = comparison_group) %>%
      left_join(
        ypos,
        by = c("gene", "ident")
      ) %>%
      mutate(
        gene = factor(gene, levels = levels(df$gene)),
        ident = factor(ident, levels = levels(df$ident))
      )
    
    p <- p +
      geom_text(
        data = annot_df,
        aes(
          x = ident,
          y = y,
          label = significance_label,
          color = gene
        ),
        position = dodge,
        vjust = 0,
        size = 4.25,
        show.legend = FALSE
      )
  }
  
  return(p)
}


make_gene_violinplot <- function(
    df,
    pair_tests = NULL,
    y_label = "Log-Normalized Expression",
    fill_palette = "Set3"
) {
  
  dodge <- position_dodge(width = 0.78)
  
  p <- ggplot(
    df,
    aes(x = ident, y = expr, fill = gene)
  ) +
    geom_violin(
      position = dodge,
      width = 0.78,
      scale = "width",
      trim = TRUE,
      color = "grey25",
      linewidth = 0.5
    ) +
    geom_point(
      aes(color = gene),
      position = position_jitterdodge(
        jitter.width = 0.12,
        dodge.width = 0.78
      ),
      alpha = 0.40,
      size = 1.3,
      shape = 15,
      stroke = 0
    ) +
    scale_fill_brewer(palette = fill_palette) +
    scale_color_brewer(
      palette = fill_palette,
      guide = "none"
    ) +
    labs(
      x = NULL,
      y = y_label,
      fill = "Gene"
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(
        angle = 20,
        hjust = 1,
        face = "bold"
      ),
      axis.text.y = element_text(face = "bold"),
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      panel.border = element_rect(
        color = "grey85",
        fill = NA
      ),
      axis.line = element_line(color = "black")
    ) +
    coord_cartesian(
      ylim = c(0, max(df$expr, na.rm = TRUE) * 1.15)
    )
  
  if (!is.null(pair_tests)) {
    
    y_span <- diff(range(df$expr, na.rm = TRUE))
    y_bump <- ifelse(
      is.finite(y_span) && y_span > 0,
      0.06 * y_span,
      0.3
    )
    
    ypos <- df %>%
      group_by(gene, ident) %>%
      summarise(
        y = max(expr, na.rm = TRUE) + y_bump,
        .groups = "drop"
      )
    
    annot_df <- pair_tests %>%
      rename(ident = comparison_group) %>%
      left_join(
        ypos,
        by = c("gene", "ident")
      ) %>%
      mutate(
        gene = factor(gene, levels = levels(df$gene)),
        ident = factor(ident, levels = levels(df$ident))
      )
    
    p <- p +
      geom_text(
        data = annot_df,
        aes(
          x = ident,
          y = y,
          label = significance_label,
          color = gene
        ),
        position = dodge,
        vjust = 0,
        size = 4.25,
        show.legend = FALSE
      )
  }
  
  return(p)
}


make_stacked_violinplot <- function(
    df,
    y_label = "Expression Level (SCT)",
    fill_palette = "Set3"
) {
  
  p <- ggplot(
    df,
    aes(x = ident, y = expr, fill = ident)
  ) +
    geom_violin(
      scale = "width",
      adjust = 1.5,
      trim = TRUE,
      color = "black",
      linewidth = 0.3
    ) +
    stat_summary(
      fun = mean,
      geom = "point",
      size = 1.2,
      color = "white"
    ) +
    stat_summary(
      fun = mean,
      geom = "point",
      size = 0.6,
      color = "black"
    ) +
    facet_grid(
      gene ~ .,
      scales = "free_y"
    ) +
    scale_fill_brewer(palette = fill_palette) +
    scale_y_continuous(n.breaks = 3) +
    labs(
      x = NULL,
      y = y_label
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        face = "bold",
        color = "black"
      ),
      axis.text.y = element_text(
        size = 9,
        color = "black"
      ),
      axis.ticks.y = element_line(color = "black"),
      strip.text.y = element_text(
        angle = 0,
        face = "bold.italic",
        size = 11
      ),
      strip.background = element_blank(),
      legend.position = "none",
      panel.spacing = unit(0.1, "lines"),
      axis.line = element_line(color = "black")
    )
  
  return(p)
}


# ==============================================================================
# --- 3. Gene sets for expression visualization ---
# ==============================================================================

gene_sets <- list(
  
  "Profibrotic_EMac_factors" = c(
    "LOC119888537", "bmp1a", "p4ha1b",
    "loxa", "LOC119897617", "LOC119886711"
  ),
  
  "Inflammatory_markers" = c(
    "LOC119885585", "LOC119888605", "LOC119912424",
    "LOC119908110", "LOC119890600", "ctsz",
    "tnfa", "m17"
  ),
  
  "Immune_homeostasis_markers" = c(
    "LOC119910049", "LOC119884197", "si-dkey-243k1.3",
    "ywhag1", "LOC119893531", "hic1", "junbb", "LOC119918304"
  ),
  
  "Transporter_channel_markers" = c(
    "LOC119901983", "rhbg", "LOC119904410",
    "slc6a14", "LOC119885976", "LOC119905556",
    "slc5a8l", "slc28a1"
  ),
  
  "Metabolic_adaptation_markers" = c(
    "LOC119910845", "aass", "fabp4a", "LOC119910978",
    "cdab", "nudt4b", "LOC119912177", "LOC119900507"
  )
)


# ==============================================================================
# --- 4. Expression visualization across infected-only spatial domains ---
# ==============================================================================

groups_infected_only <- c(
  "Eryt",
  "Eryt-Gran",
  "Gran",
  "Non-G Mac",
  "B",
  "T-B-APCs",
  "Outer E-Mac",
  "Inner E-Mac",
  "Fibro-Endo"
)

expression_plot_results <- list()

df_profibrotic_all_domains <- extract_sct_expression_long(
  object = TS_fliter,
  genes = gene_sets$Profibrotic_EMac_factors,
  groups_keep = groups_infected_only,
  assay = "SCT",
  layer = "data"
)

expression_plot_results$profibrotic_boxplot_all_domains <- make_gene_boxplot(
  df = df_profibrotic_all_domains,
  pair_tests = NULL,
  y_label = "SCT data (log1p corrected counts)",
  fill_palette = "Set3"
)

expression_plot_results$profibrotic_stacked_violin_all_domains <- make_stacked_violinplot(
  df = df_profibrotic_all_domains,
  y_label = "Expression Level (SCT)",
  fill_palette = "Set3"
)

print(expression_plot_results$profibrotic_boxplot_all_domains)
print(expression_plot_results$profibrotic_stacked_violin_all_domains)


# ==============================================================================
# --- 5. Pairwise Seurat marker tests for profibrotic E-Mac factors ---
# ==============================================================================

target_genes <- gene_sets$Profibrotic_EMac_factors

marker_test_results <- list()

marker_test_results$Outer_EMac_vs_NonG_Mac <- FindMarkers(
  object = TS_fliter,
  ident.1 = "Outer E-Mac",
  ident.2 = "Non-G Mac",
  features = target_genes,
  logfc.threshold = 0,
  min.pct = 0
)

marker_test_results$Outer_EMac_vs_Fibro_Endo <- FindMarkers(
  object = TS_fliter,
  ident.1 = "Outer E-Mac",
  ident.2 = "Fibro-Endo",
  features = target_genes,
  logfc.threshold = 0,
  min.pct = 0
)

print(marker_test_results$Outer_EMac_vs_NonG_Mac)
print(marker_test_results$Outer_EMac_vs_Fibro_Endo)


# ==============================================================================
# --- 6. Expression visualization within macrophage and E-Mac domains ---
# ==============================================================================

groups_macrophage_emac <- c(
  "Mac (Homeostatic)",
  "Mac (Infection-associated)",
  "Outer E-Mac",
  "Inner E-Mac"
)

baseline_group <- "Mac (Infection-associated)"

comparison_groups <- c(
  "Outer E-Mac",
  "Inner E-Mac"
)

selected_gene_set <- "Metabolic_adaptation_markers"

df_macrophage_expression <- extract_sct_expression_long(
  object = object_sub,
  genes = gene_sets[[selected_gene_set]],
  groups_keep = groups_macrophage_emac,
  assay = "SCT",
  layer = "data"
)

pair_tests_macrophage <- run_gene_pairwise_tests(
  df = df_macrophage_expression,
  baseline = baseline_group,
  comparison_groups = comparison_groups,
  p_adjust_method = "BH"
)

expression_plot_results$macrophage_boxplot <- make_gene_boxplot(
  df = df_macrophage_expression,
  pair_tests = pair_tests_macrophage,
  y_label = "Log-Normalized Expression",
  fill_palette = "Paired"
)

expression_plot_results$macrophage_violinplot <- make_gene_violinplot(
  df = df_macrophage_expression,
  pair_tests = pair_tests_macrophage,
  y_label = "Log-Normalized Expression",
  fill_palette = "Set3"
)

print(pair_tests_macrophage)
print(expression_plot_results$macrophage_boxplot)
print(expression_plot_results$macrophage_violinplot)


# ==============================================================================
# --- 7. Loop over all macrophage/E-Mac marker gene sets ---
# ==============================================================================

macrophage_gene_sets <- c(
  "Inflammatory_markers",
  "Immune_homeostasis_markers",
  "Transporter_channel_markers",
  "Metabolic_adaptation_markers"
)

macrophage_expression_results <- list()

for (gene_set_name in macrophage_gene_sets) {
  
  message("Processing macrophage/E-Mac gene set: ", gene_set_name)
  
  df_tmp <- extract_sct_expression_long(
    object = object_sub,
    genes = gene_sets[[gene_set_name]],
    groups_keep = groups_macrophage_emac,
    assay = "SCT",
    layer = "data"
  )
  
  pair_tests_tmp <- run_gene_pairwise_tests(
    df = df_tmp,
    baseline = baseline_group,
    comparison_groups = comparison_groups,
    p_adjust_method = "BH"
  )
  
  macrophage_expression_results[[gene_set_name]] <- list(
    expression_data = df_tmp,
    pair_tests = pair_tests_tmp,
    boxplot = make_gene_boxplot(
      df = df_tmp,
      pair_tests = pair_tests_tmp,
      y_label = "Log-Normalized Expression",
      fill_palette = "Paired"
    ),
    violinplot = make_gene_violinplot(
      df = df_tmp,
      pair_tests = pair_tests_tmp,
      y_label = "Log-Normalized Expression",
      fill_palette = "Set3"
    )
  )
}


# Example:
print(macrophage_expression_results$Metabolic_adaptation_markers$boxplot)


# ==============================================================================
# --- 8. Session information ---
# ==============================================================================

sessionInfo()