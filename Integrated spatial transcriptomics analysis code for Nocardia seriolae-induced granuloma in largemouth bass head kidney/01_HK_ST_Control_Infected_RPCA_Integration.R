# ==============================================================================
# Head Kidney Granuloma - Control and Infected Spatial Transcriptomics Integration
# RPCA-based mild integration of BMKMANU S3000 spatial transcriptomics data
#
# Purpose:
#   Integrate control and Nocardia seriolae-infected largemouth bass head kidney
#   spatial transcriptomics data while preserving infection-associated spatial
#   domains.
#
# Platform:
#   BMKMANU S3000 spatial transcriptomics
# ==============================================================================


# ==============================================================================
# --- 0. Load R packages ---
# ==============================================================================

suppressMessages({
  library(Seurat)
  library(scCustomize)
  library(ggplot2)
  library(tidydr)
  library(dplyr)
  library(patchwork)
  library(ggrepel)
  library(ggsci)
  library(viridis)
  library(ComplexHeatmap)
  library(PRECAST)
  library(scMayoMap)
  library(harmony)
  library(RColorBrewer)
  library(tidyr)
  library(scales)
  library(gghalves)
  library(ggpubr)
  library(scico)
  library(Matrix)
})

options(future.globals.maxSize = 9e10)
options(Seurat.object.assay.version = "v5")


# ==============================================================================
# --- 1. Set example input paths and parameters ---
# ==============================================================================

# Example project directory.
# Users should modify this path according to their local file organization.
project_dir <- "path/to/project"

# Custom function for importing BMKMANU S3000 spatial transcriptomics data.
# The script CreateBmkObject_v2.R should be provided with the code repository.
source(file.path(project_dir, "scripts", "CreateBmkObject_v2.R"))

# Example input paths for BMKMANU S3000 L4-level data.
control_matrix_path <- file.path(
  project_dir,
  "data",
  "spatial",
  "head_kidney_control_TSDZ",
  "L4_heAuto"
)

control_image_path <- file.path(
  control_matrix_path,
  "he_roi_small.png"
)

infected_matrix_path <- file.path(
  project_dir,
  "data",
  "spatial",
  "head_kidney_infected_TS",
  "L4_heAuto"
)

infected_image_path <- file.path(
  infected_matrix_path,
  "he_roi_small.png"
)

# Analysis parameters
min_features_import <- 100
min_features_qc     <- 250
min_cells_gene      <- 3
nfeatures_integrate <- 3000
npcs_use            <- 30
dims_use            <- 1:20
cluster_resolution  <- 0.60


# ==============================================================================
# --- 2. Gene lists for filtering ---
# ==============================================================================

mito_genes <- c(
  "ND1", "ND2", "COX1", "COX2", "ATP8", "ATP6",
  "COX3", "ND3", "ND4L", "ND4", "ND5", "ND6", "CYTB"
)

# Ribosomal genes can be provided as a separate text file in the repository.
# One gene symbol per line is recommended.
ribo_gene_file <- file.path(
  project_dir,
  "gene_lists",
  "largemouth_bass_ribosomal_genes.txt"
)

if (file.exists(ribo_gene_file)) {
  ribo_genes <- readLines(ribo_gene_file)
} else {
  ribo_genes <- character(0)
  warning("Ribosomal gene file not found. Proceeding without ribosomal gene removal.")
}


# ==============================================================================
# --- 3. Import BMKMANU S3000 spatial transcriptomics data ---
# ==============================================================================

TSDZ_raw <- CreateBmkObject(
  matrix_path  = control_matrix_path,
  png_path     = control_image_path,
  min.cells    = 0,
  min.features = min_features_import,
  type         = "S3000"
)

TSDZ_raw$condition  <- "TSDZ"
TSDZ_raw$orig.ident <- "Control"


TS_raw <- CreateBmkObject(
  matrix_path  = infected_matrix_path,
  png_path     = infected_image_path,
  min.cells    = 0,
  min.features = min_features_import,
  type         = "S3000"
)

TS_raw$condition  <- "TS"
TS_raw$orig.ident <- "Treat"


# ==============================================================================
# --- 4. Merge samples before gene filtering ---
# ==============================================================================

# Merge immediately after import with min.cells = 0.
# This helps preserve sample-specific genes, including infection-associated genes.
seu_int_rpca <- merge(TSDZ_raw, y = TS_raw)

rm(TSDZ_raw, TS_raw)
gc()

message(
  "Merged object: ",
  nrow(seu_int_rpca), " genes; ",
  ncol(seu_int_rpca), " spatial units."
)


# ==============================================================================
# --- 5. Unified QC and gene filtering ---
# ==============================================================================

seu_int_rpca[["percent.mt"]] <- PercentageFeatureSet(
  object = seu_int_rpca,
  features = mito_genes
)

# Filter low-quality spatial units
seu_int_rpca <- subset(
  seu_int_rpca,
  subset = nFeature_Spatial > min_features_qc
)

# Remove ribosomal genes
if (length(ribo_genes) > 0) {
  genes_no_ribo <- setdiff(rownames(seu_int_rpca), ribo_genes)
  seu_int_rpca <- seu_int_rpca[genes_no_ribo, ]
}

# Remove very lowly expressed genes across Seurat v5 count layers
counts_layers <- grep(
  pattern = "^counts",
  x = Layers(seu_int_rpca[["Spatial"]]),
  value = TRUE
)

total_spots_expressing <- Reduce(
  f = `+`,
  x = lapply(counts_layers, function(layer_i) {
    Matrix::rowSums(
      LayerData(seu_int_rpca, assay = "Spatial", layer = layer_i) > 0
    )
  })
)

genes_valid <- names(total_spots_expressing[total_spots_expressing >= min_cells_gene])
seu_int_rpca <- seu_int_rpca[genes_valid, ]

message(
  "After QC: ",
  nrow(seu_int_rpca), " genes; ",
  ncol(seu_int_rpca), " spatial units."
)


# ==============================================================================
# --- 6. Prepare objects for SCT-based RPCA integration ---
# ==============================================================================

obj.list <- SplitObject(seu_int_rpca, split.by = "condition")

rm(seu_int_rpca)
gc()

obj.list <- lapply(obj.list, function(x) {
  x <- SCTransform(
    object = x,
    assay = "Spatial",
    vst.flavor = "v2",
    method = "glmGamPoi",
    return.only.var.genes = FALSE,
    verbose = TRUE
  )
  return(x)
})

obj.list <- lapply(obj.list, function(x) {
  x <- RunPCA(
    object = x,
    assay = "SCT",
    npcs = npcs_use,
    verbose = FALSE
  )
  return(x)
})


# ==============================================================================
# --- 7. SCT-based RPCA integration ---
# ==============================================================================

features <- SelectIntegrationFeatures(
  object.list = obj.list,
  nfeatures = nfeatures_integrate
)

obj.list <- PrepSCTIntegration(
  object.list = obj.list,
  anchor.features = features
)

anchors_rpca <- FindIntegrationAnchors(
  object.list = obj.list,
  normalization.method = "SCT",
  anchor.features = features,
  reduction = "rpca",
  dims = 1:npcs_use,
  verbose = TRUE
)

seu_int_rpca <- IntegrateData(
  anchorset = anchors_rpca,
  normalization.method = "SCT",
  dims = 1:npcs_use,
  verbose = TRUE
)

rm(obj.list, anchors_rpca)
gc()


# ==============================================================================
# --- 8. Dimensional reduction and clustering ---
# ==============================================================================

DefaultAssay(seu_int_rpca) <- "integrated"

seu_int_rpca <- RunPCA(
  object = seu_int_rpca,
  verbose = FALSE
)

seu_int_rpca <- RunUMAP(
  object = seu_int_rpca,
  reduction = "pca",
  dims = dims_use
)

seu_int_rpca <- FindNeighbors(
  object = seu_int_rpca,
  reduction = "pca",
  dims = dims_use
)

seu_int_rpca <- FindClusters(
  object = seu_int_rpca,
  resolution = cluster_resolution
)


# ==============================================================================
# --- 9. Color scheme setup ---
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

color_prgn <- colorRampPalette(
  rev(brewer.pal(n = 11, name = "PRGn"))
)(100)


# ==============================================================================
# --- 10. UMAP and spatial visualization ---
# ==============================================================================

# UMAP split by sample
p_umap_split <- DimPlot(
  object = seu_int_rpca,
  reduction = "umap",
  label = TRUE,
  pt.size = 0.3,
  raster = FALSE,
  split.by = "orig.ident"
) &
  theme(legend.title = element_blank()) &
  scale_color_manual(values = all_available_colors)

print(p_umap_split)


# UMAP colored by integrated clusters
p_umap <- DimPlot(
  object = seu_int_rpca,
  reduction = "umap",
  label = TRUE,
  pt.size = 0.3,
  raster = FALSE
) &
  theme(legend.title = element_blank()) &
  scale_color_manual(values = all_available_colors)

print(p_umap)


# UMAP colored by condition
p_umap_condition <- DimPlot(
  object = seu_int_rpca,
  reduction = "umap",
  group.by = "condition",
  pt.size = 0.3,
  raster = TRUE
)

print(p_umap_condition)


# Spatial visualization of integrated clusters
p_spatial_clusters <- SpatialDimPlot(
  object = seu_int_rpca,
  label = FALSE,
  label.size = 6,
  image.alpha = 0,
  pt.size.factor = 4.0,
  crop = FALSE
) &
  guides(fill = guide_legend(override.aes = list(size = 4))) &
  theme(legend.title = element_blank()) &
  scale_fill_manual(values = all_available_colors)

print(p_spatial_clusters)


# ==============================================================================
# --- 11. Marker gene DotPlot for domain annotation ---
# ==============================================================================

DefaultAssay(seu_int_rpca) <- "SCT"

marker_list <- list(
  "Vascular niche" = c(
    "col1a1b", "col6a1", "sparc",
    "sox18", "LOC119912382", "dll4",
    "pecam1", "vwf", "plvapb",
    "acta2", "tagln"
  ),
  
  "B" = c(
    "LOC119890508", "LOC119890540", "LOC119890555", "LOC119890570",
    "ebf1a", "sox4b", "LOC119882634", "cd79a",
    "LOC119884751", "LOC119897215"
  ),
  
  "T" = c(
    "LOC119903450", "LOC119910797", "LOC119883393", "LOC119883434",
    "LOC119914246", "LOC119897409", "LOC119897410",
    "cd8a", "cd8b", "lck", "LOC119900445"
  ),
  
  "T-B interaction zone" = c(
    "LOC119891695", "cxcl19", "LOC119888586", "LOC119894313"
  ),
  
  "Antigen Processing & Presentation" = c(
    "LOC119888605", "LOC119917015", "LOC119912422",
    "LOC119884279", "LOC119918458",
    "ifi30", "cd83", "cd40"
  ),
  
  "Mac" = c(
    "csf1rb", "marco", "mrc1b", "LOC119891768",
    "si-ch211-212k18.7", "LOC119894959", "LOC119883615",
    "LOC119885585", "LOC119882213", "LOC119890208"
  ),
  
  "Eryt" = c(
    "zgc-163057", "LOC119891533",
    "slc4a1a", "cahz",
    "alas2", "fech", "slc25a37",
    "sptb", "klf1", "epor"
  ),
  
  "Gran" = c(
    "LOC119914447", "LOC119900902", "LOC119889900",
    "LOC119894083", "zgc-92027", "mmp9",
    "LOC119903646", "LOC119904387",
    "LOC119900132", "LOC119891034"
  ),
  
  "E-Mac" = c(
    "LOC119917376", "jupa", "LOC119889888", "dspa",
    "LOC119902058", "LOC119917243", "dsc2l",
    "LOC119900663", "cldn11a",
    "LOC119908668", "LOC119908847"
  )
)

p_dotplot <- DotPlot(
  object = seu_int_rpca,
  features = marker_list
) +
  scale_color_gradientn(colors = color_prgn) +
  labs(
    x = "Cell Clusters",
    y = "Marker Genes",
    color = "Avg Expression",
    size = "Percent Expressed"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 9,
      face = "bold"
    ),
    axis.text.y = element_text(
      size = 9,
      face = "italic"
    ),
    axis.title = element_text(
      size = 10,
      face = "bold"
    ),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    panel.grid = element_blank(),
    legend.position = "right"
  )

print(p_dotplot)


# ==============================================================================
# --- 12. Cell/domain composition calculation ---
# ==============================================================================

cell_table <- table(
  cluster = Idents(seu_int_rpca),
  sample = seu_int_rpca$orig.ident
)

cell_ratio <- prop.table(cell_table, margin = 2) %>%
  as.data.frame()

colnames(cell_ratio) <- c("Cluster", "Sample", "Ratio")

cell_counts <- as.data.frame(cell_table)
colnames(cell_counts) <- c("Cluster", "Sample", "Count")

cell_ratio <- left_join(
  cell_ratio,
  cell_counts,
  by = c("Cluster", "Sample")
)

cell_ratio$Cluster <- factor(cell_ratio$Cluster)

cell_ratio <- cell_ratio %>%
  arrange(Sample, Cluster)

p_cellratio <- ggplot(cell_ratio) +
  geom_bar(
    aes(x = Sample, y = Ratio, fill = Cluster),
    stat = "identity",
    width = 0.7,
    size = 0.5,
    colour = "#222222"
  ) +
  theme_classic() +
  labs(
    x = "Sample",
    y = "Ratio"
  ) +
  scale_fill_manual(values = all_available_colors) +
  theme(
    panel.border = element_rect(
      fill = NA,
      color = "black",
      size = 0.5,
      linetype = "solid"
    )
  )

print(p_cellratio)


# ==============================================================================
# --- 13. Session information ---
# ==============================================================================

sessionInfo()