# ==============================================================================
# Public Human TB Granuloma - Spatial Transcriptomics Analysis
# ==============================================================================
#
# Purpose:
#   Reanalyze public human lung Mycobacterium tuberculosis granuloma
#   spatial transcriptomics data and project teleost E-Mac-associated
#   gene signatures onto human TB granuloma sections.
#
# Main steps:
#   1. Load human TB granuloma Visium spatial transcriptomics data.
#   2. Perform basic QC visualization.
#   3. Filter low-quality spatial spots.
#   4. Normalize data using SCTransform.
#   5. Perform PCA, UMAP and clustering.
#   6. Visualize clusters in UMAP and spatial coordinates.
#   7. Identify cluster marker genes.
#   8. Project representative teleost E-Mac orthologous genes onto human TB data.
#
# Required input:
#   A 10x Genomics Visium-format spatial transcriptomics dataset.
#   The input directory should contain the feature-barcode matrix and spatial folder.
#
# Notes:
#   This script does not save figures or tables to disk.
#   Users can export plots and results as needed.
# ==============================================================================


# ==============================================================================
# --- 0. Load R packages ---
# ==============================================================================

suppressMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
  library(scales)
  library(viridis)
  library(ggrepel)
  library(patchwork)
  library(cowplot)
  library(scico)
})

options(future.globals.maxSize = 9e10)
options(Seurat.object.assay.version = "v5")


# ==============================================================================
# --- 1. Example input path and data loading ---
# ==============================================================================

# Example project directory.
# Users should modify this path according to their local file organization.
project_dir <- "path/to/project"

# Example input path for the public human TB granuloma Visium dataset.
# This directory should contain the 10x spatial transcriptomics output,
# including the matrix files and the spatial folder.
human_tb_data_dir <- file.path(
  project_dir,
  "public_data",
  "human_TB_granuloma",
  "Visium"
)

# Load human TB granuloma spatial data.
human_tb <- Load10X_Spatial(
  data.dir = human_tb_data_dir
)

human_tb$orig.ident <- "Human_TB_Granuloma"

message("Input object:")
print(human_tb)


# ==============================================================================
# --- 2. Basic QC before filtering ---
# ==============================================================================

# Human mitochondrial genes are usually prefixed with "MT-".
human_tb[["percent.mt"]] <- PercentageFeatureSet(
  object = human_tb,
  pattern = "^MT-"
)

qc_features <- c(
  "nFeature_Spatial",
  "nCount_Spatial",
  "percent.mt"
)

p_qc_vln_unfiltered <- VlnPlot(
  object = human_tb,
  features = qc_features,
  ncol = 3,
  pt.size = 0.1
)

print(p_qc_vln_unfiltered)


qc_spatial_plots_unfiltered <- list()

for (feature in qc_features) {
  
  p <- SpatialFeaturePlot(
    object = human_tb,
    features = feature,
    image.alpha = 0.1,
    pt.size.factor = 2.0,
    crop = FALSE
  ) +
    theme(legend.position = "right") +
    scale_fill_viridis_c(
      option = "H",
      name = feature
    ) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 16
      )
    ) +
    ggtitle(paste("Human TB:", feature))
  
  qc_spatial_plots_unfiltered[[feature]] <- p
}

print(qc_spatial_plots_unfiltered$nFeature_Spatial)
print(qc_spatial_plots_unfiltered$nCount_Spatial)
print(qc_spatial_plots_unfiltered$percent.mt)


# ==============================================================================
# --- 3. Spot filtering ---
# ==============================================================================

# Filtering criteria used in the exploratory analysis:
#   nFeature_Spatial > 300
#   nFeature_Spatial < 9000
#   percent.mt < 20
human_tb <- subset(
  x = human_tb,
  subset = nFeature_Spatial > 300 &
    nFeature_Spatial < 9000 &
    percent.mt < 20
)

message("After filtering:")
print(human_tb)


# ==============================================================================
# --- 4. Basic QC after filtering ---
# ==============================================================================

p_qc_vln_filtered <- VlnPlot(
  object = human_tb,
  features = qc_features,
  ncol = 3,
  pt.size = 0.1
)

print(p_qc_vln_filtered)


qc_spatial_plots_filtered <- list()

for (feature in qc_features) {
  
  p <- SpatialFeaturePlot(
    object = human_tb,
    features = feature,
    image.alpha = 0.1,
    pt.size.factor = 2.0,
    crop = FALSE
  ) +
    theme(legend.position = "right") +
    scale_fill_viridis_c(
      option = "H",
      name = feature
    ) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 16
      )
    ) +
    ggtitle(paste("Human TB filtered:", feature))
  
  qc_spatial_plots_filtered[[feature]] <- p
}

print(qc_spatial_plots_filtered$nFeature_Spatial)
print(qc_spatial_plots_filtered$nCount_Spatial)
print(qc_spatial_plots_filtered$percent.mt)


# ==============================================================================
# --- 5. SCTransform normalization, dimensional reduction and clustering ---
# ==============================================================================

human_tb <- SCTransform(
  object = human_tb,
  assay = "Spatial",
  verbose = FALSE
)

human_tb <- RunPCA(
  object = human_tb,
  assay = "SCT",
  verbose = FALSE
)

human_tb <- FindNeighbors(
  object = human_tb,
  reduction = "pca",
  dims = 1:30
)

human_tb <- FindClusters(
  object = human_tb,
  verbose = FALSE,
  resolution = 0.6
)

human_tb <- RunUMAP(
  object = human_tb,
  reduction = "pca",
  dims = 1:30
)


# ==============================================================================
# --- 6. Color scheme setup ---
# ==============================================================================

base_colors_paired <- brewer.pal(
  min(12, brewer.pal.info["Paired", "maxcolors"]),
  "Paired"
)

base_colors_dark2 <- brewer.pal(
  min(8, brewer.pal.info["Dark2", "maxcolors"]),
  "Dark2"
)

base_colors_set3 <- brewer.pal(
  min(12, brewer.pal.info["Set3", "maxcolors"]),
  "Set3"
)

base_colors_set1 <- brewer.pal(
  min(9, brewer.pal.info["Set1", "maxcolors"]),
  "Set1"
)

all_available_colors <- c(
  base_colors_paired,
  base_colors_dark2,
  base_colors_set3,
  base_colors_set1,
  colors()[grep("blue|green|orange|purple|yellow|pink|brown", colors())]
)


# ==============================================================================
# --- 7. UMAP visualization ---
# ==============================================================================

umap_df <- human_tb@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  cbind(
    seurat_clusters = human_tb@active.ident
  )

seurat_clusterspos <- umap_df %>%
  group_by(seurat_clusters) %>%
  summarise(
    umap_1 = median(umap_1),
    umap_2 = median(umap_2),
    .groups = "drop"
  )

p_umap_base <- DimPlot(
  object = human_tb,
  reduction = "umap",
  label = FALSE,
  pt.size = 0.8
) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(
      face = "bold",
      hjust = 0.03
    )
  ) +
  NoLegend()

p_umap_clusters <- p_umap_base +
  geom_label_repel(
    aes(
      x = umap_1,
      y = umap_2,
      label = seurat_clusters
    ),
    data = seurat_clusterspos,
    fontface = "bold",
    box.padding = 0.5
  ) +
  scale_color_manual(values = all_available_colors) +
  ggtitle("Human TB: UMAP clusters")

print(p_umap_clusters)


# ==============================================================================
# --- 8. Spatial cluster visualization ---
# ==============================================================================

p_spatial_clusters <- SpatialDimPlot(
  object = human_tb,
  label = FALSE,
  image.alpha = 0.8,
  pt.size.factor = 2.0,
  crop = FALSE
) +
  scale_fill_manual(values = all_available_colors) +
  guides(
    fill = guide_legend(
      override.aes = list(size = 4)
    )
  ) +
  theme(legend.title = element_blank()) +
  ggtitle("Human TB: spatial clusters")

print(p_spatial_clusters)


# ==============================================================================
# --- 9. Cluster marker gene identification ---
# ==============================================================================

human_tb_markers <- FindAllMarkers(
  object = human_tb,
  only.pos = TRUE,
  min.pct = 0.1,
  test.use = "wilcox",
  logfc.threshold = 0.1
)

significant_markers <- human_tb_markers %>%
  filter(p_val_adj < 0.1)

top20_markers <- significant_markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 20,
    with_ties = FALSE
  ) %>%
  ungroup()

message("Number of marker genes identified:")
print(nrow(human_tb_markers))

message("Number of significant marker genes with adjusted P < 0.1:")
print(nrow(significant_markers))

message("Top 20 marker genes per cluster:")
print(top20_markers)


# ==============================================================================
# --- 10. Teleost E-Mac signature projection onto human TB granuloma ---
# ==============================================================================

# Representative human orthologous genes corresponding to selected teleost
# E-Mac-associated features.
# Users may replace this list with the final ortholog table used in the manuscript.
ortholog_genes <- c(
  "WT1",
  "GJA1",
  "CDH1",
  "KRT8",
  "MMP9",
  "SPP1",
  "ALDH1A2"
)

valid_ortholog_genes <- ortholog_genes[
  ortholog_genes %in% rownames(human_tb)
]

missing_ortholog_genes <- setdiff(
  ortholog_genes,
  valid_ortholog_genes
)

if (length(missing_ortholog_genes) > 0) {
  warning(
    "The following orthologous genes were not found and will be skipped: ",
    paste(missing_ortholog_genes, collapse = ", ")
  )
}

if (length(valid_ortholog_genes) == 0) {
  warning("No valid orthologous genes were found for module score projection.")
} else {
  
  emac_signature_list <- list(
    Teleost_EMac_Signature = valid_ortholog_genes
  )
  
  human_tb <- AddModuleScore(
    object = human_tb,
    features = emac_signature_list,
    name = "Teleost_EMac_Score",
    assay = "SCT",
    seed = 1
  )
  
  p_teleost_emac_projection <- SpatialFeaturePlot(
    object = human_tb,
    features = "Teleost_EMac_Score1",
    pt.size.factor = 1.5,
    alpha = c(0.1, 1),
    crop = FALSE
  ) +
    scale_fill_viridis_c(option = "plasma") +
    ggtitle("Teleost E-Mac signature projected onto human TB granuloma")
  
  print(p_teleost_emac_projection)
}


# ==============================================================================
# --- 11. Session information ---
# ==============================================================================

sessionInfo()