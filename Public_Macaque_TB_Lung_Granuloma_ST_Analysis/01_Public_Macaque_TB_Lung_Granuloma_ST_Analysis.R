# ==============================================================================
# Public Macaque TB Granuloma - Spatial Transcriptomics Analysis
# ==============================================================================
#
# Purpose:
#   Reanalyze public spatial transcriptomics data from macaque lung
#   Mycobacterium tuberculosis granulomas.
#
# Main steps:
#   1. Load public macaque TB granuloma Visium data.
#   2. Perform basic QC visualization.
#   3. Filter low-quality spatial spots.
#   4. Normalize data using SCTransform.
#   5. Perform PCA, UMAP and clustering.
#   6. Visualize clusters in UMAP and spatial coordinates.
#   7. Identify cluster marker genes.
#
# Required input:
#   sample_data_trans:
#     A Seurat spatial object containing the public macaque TB granuloma
#     Visium dataset.
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
# --- 1. Example input object loading ---
# ==============================================================================

# Example project directory.
# Users should modify this path according to their local file organization.
project_dir <- "path/to/project"

# Example public-data object path.
# The object should be prepared from the public macaque TB granuloma Visium data.
public_macaque_tb_object_path <- file.path(
  project_dir,
  "public_data",
  "macaque_TB_granuloma",
  "TB_Granuloma_Visium.RData"
)

# Load object if it is not already present in the R session.
# The loaded object is expected to be named sample_data_trans.
if (!exists("sample_data_trans")) {
  load(public_macaque_tb_object_path)
}

if (!exists("sample_data_trans")) {
  stop(
    "Object 'sample_data_trans' was not found. ",
    "Please load the public macaque TB granuloma spatial Seurat object first."
  )
}

message("Input object:")
print(sample_data_trans)


# ==============================================================================
# --- 2. Basic QC visualization before filtering ---
# ==============================================================================

qc_features <- c(
  "nFeature_Spatial",
  "nCount_Spatial"
)

p_qc_vln_unfiltered <- VlnPlot(
  object = sample_data_trans,
  features = qc_features,
  ncol = 2,
  pt.size = 0.1
)

print(p_qc_vln_unfiltered)


qc_spatial_plots_unfiltered <- list()

for (feature in qc_features) {
  
  p <- SpatialFeaturePlot(
    object = sample_data_trans,
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
    ggtitle(paste("Macaque TB:", feature))
  
  qc_spatial_plots_unfiltered[[feature]] <- p
}

print(qc_spatial_plots_unfiltered$nFeature_Spatial)
print(qc_spatial_plots_unfiltered$nCount_Spatial)


# ==============================================================================
# --- 3. Spot filtering ---
# ==============================================================================

# Original analysis used nFeature_Spatial > 300.
sample_data_trans <- subset(
  x = sample_data_trans,
  subset = nFeature_Spatial > 300
)

message("After filtering:")
print(sample_data_trans)


# ==============================================================================
# --- 4. Basic QC visualization after filtering ---
# ==============================================================================

p_qc_vln_filtered <- VlnPlot(
  object = sample_data_trans,
  features = qc_features,
  ncol = 2,
  pt.size = 0.1
)

print(p_qc_vln_filtered)


qc_spatial_plots_filtered <- list()

for (feature in qc_features) {
  
  p <- SpatialFeaturePlot(
    object = sample_data_trans,
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
    ggtitle(paste("Macaque TB filtered:", feature))
  
  qc_spatial_plots_filtered[[feature]] <- p
}

print(qc_spatial_plots_filtered$nFeature_Spatial)
print(qc_spatial_plots_filtered$nCount_Spatial)


# ==============================================================================
# --- 5. SCTransform normalization, dimensional reduction and clustering ---
# ==============================================================================

sample_data_trans <- SCTransform(
  object = sample_data_trans,
  assay = "Spatial",
  verbose = FALSE
)

sample_data_trans <- RunPCA(
  object = sample_data_trans,
  assay = "SCT",
  verbose = FALSE
)

sample_data_trans <- FindNeighbors(
  object = sample_data_trans,
  reduction = "pca",
  dims = 1:20
)

sample_data_trans <- FindClusters(
  object = sample_data_trans,
  verbose = FALSE,
  resolution = 0.6
)

sample_data_trans <- RunUMAP(
  object = sample_data_trans,
  reduction = "pca",
  dims = 1:20
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

umap_df <- sample_data_trans@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  cbind(
    seurat_clusters = sample_data_trans@active.ident
  )

seurat_clusterspos <- umap_df %>%
  group_by(seurat_clusters) %>%
  summarise(
    umap_1 = median(umap_1),
    umap_2 = median(umap_2),
    .groups = "drop"
  )

p_umap_base <- DimPlot(
  object = sample_data_trans,
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
  scale_color_manual(values = all_available_colors)

print(p_umap_clusters)


# ==============================================================================
# --- 8. Spatial cluster visualization ---
# ==============================================================================

p_spatial_clusters <- SpatialDimPlot(
  object = sample_data_trans,
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
  ggtitle("Macaque TB: spatial clusters")

print(p_spatial_clusters)


# ==============================================================================
# --- 9. Cluster marker gene identification ---
# ==============================================================================

sample_data_trans_markers <- FindAllMarkers(
  object = sample_data_trans,
  only.pos = TRUE,
  min.pct = 0.1,
  test.use = "wilcox",
  logfc.threshold = 0.1
)

significant_markers <- sample_data_trans_markers %>%
  filter(p_val_adj < 0.1)

message("Number of marker genes identified:")
print(nrow(sample_data_trans_markers))

message("Number of significant marker genes with adjusted P < 0.1:")
print(nrow(significant_markers))


# ==============================================================================
# --- 10. Session information ---
# ==============================================================================

sessionInfo()