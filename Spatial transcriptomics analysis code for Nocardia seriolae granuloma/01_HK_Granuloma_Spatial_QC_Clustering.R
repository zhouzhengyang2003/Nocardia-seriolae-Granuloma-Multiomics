# ==============================================================================
# Head Kidney Granuloma - 1. Data Import, QC, Normalization, and Clustering
# ==============================================================================

# --- 0. Load R packages ---
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
})

# Set global data path
data_path <- "E:/zhengyangzhou_workking/TS-s.BSTViewer_project/subdata/L4_heAuto"

# --- 1. Data Import ---
source("C:/Users/26579/Desktop/CreateBmkObject_v2.R")

TS <- CreateBmkObject(
  matrix_path = data_path,
  png_path = file.path(data_path, "he_roi_small.png"),
  min.cells = 3,
  min.features = 100,
  type = "S3000"
)
TS$condition <- "TS"
TS$orig.ident <- 'Treat'

# --- 2. Gene Annotation & Ribosomal Gene Removal ---
# Mitochondrial and Ribosomal gene lists
mito_genes <- c("ND1", "ND2", "COX1", "COX2", "ATP8", "ATP6", 
                "COX3", "ND3", "ND4L", "ND4", "ND5", "ND6", "CYTB")

ribo_genes <- c("mrpl39", "rpl31", "mrps9", "rps20", "rpl7", "rpl21",
                "mrpl15", "rpl22l1", "LOC119898710", "LOC119898703", "mrpl45", "mrpl4", "rpl23", "rpsa", "mrps34", "mrps7",
                "rps2", "LOC119888983", "mrpl12", "rps11", "rpl19", "rpl27", "rpl3", "mrps15", "rps19", "mrpl17", "LOC119897373", "mrpl51",
                "rps9", "rpl18", "dap3", "rps27.2", "LOC119897812", "LOC119897846", "mrpl32", "LOC119898192", "LOC119898191", "rps3a", "rps6", "rpl9",
                "rps15a", "mrps28", "LOC119900488", "mrpl47", "LOC119900749", "mrpl54", "LOC119901212", "rpl36", "LOC119901312", "LOC119901314",
                "rps8a", "mrps14", "LOC119902226", "rps15", "mrpl34", "rpl5a", "rpl37a", "mrpl37", "LOC119903249", "mrpl30", "rpl8", "rps21",
                "mrps25", "rpl32", "rpl10a", "rpl10", "mrpl49", "LOC119904990", "rps18", "rpl11", "rps5", "mrps21", "rps27.1", "LOC119907920",
                "mrpl24", "mrpl9", "rplp1", "rpl30", "mrpl53", "mrpl3", "rpl14", "mrpl36", "rpl15", "mrps18b", "rpl7l1", "faua", "rps29", "mrps18a",
                "rps12", "mrps10", "mrpl19", "rps7", "rpl13a", "mrpl35", "mrpl33", "rps27a", "rpl34", "mrpl52", "mrpl48", "rps4x", "LOC119914721",
                "LOC119915000", "mrps18c", "mrps24", "mrps2", "rpl7a", "rpl28", "rpl37", "LOC119914137", "LOC119915072", "mrpl27", "mrpl58",
                "mrpl38", "LOC119916482", "rpl38", "mrpl43", "mrpl10", "mrpl2", "mrpl14", "mrps6", "rps24", "mrpl57", "mrps5", "mrps33",
                "mrpl42", "rps16", "LOC119918593", "LOC119919179", "mrps23", "mrps22", "rps25", "mrpl44", "mrpl28", "rps3", "LOC119882934",
                "rpl24", "rpl23a", "LOC119884849", "rps23", "rpl22", "rps10", "rplp2", "rpl29", "mrps16", "mrpl20", "LOC119886813", "rplp0",
                "mrps26", "rpl26", "mrpl18", "mrpl11", "rpl39", "rpl36a", "rpl12", "LOC119890837", "mrpl16", "LOC119889915", "rpl13",
                "mrpl46", "mrps11", "rpl35", "mrpl41", "mrpl1", "mrps30", "rpl6", "mrpl40", "mrps36", "mrps27", "rps28", "LOC119894123",
                "rpl27a", "rps13", "mrps35", "rpl4", "rps27l", "mrpl23", "LOC119894724", "rpl18a", "rps14", "mrpl22", "mrps17", "mrps31")  

# Remove ribosomal genes and calculate mitochondrial percentage
TS <- subset(TS, features = setdiff(rownames(TS), ribo_genes))
TS[['percent.mt']] <- PercentageFeatureSet(TS, features = mito_genes)

# --- 3. Quality Control (QC) & Filtering ---
print("Starting QC for TS...")

# 3.1 Pre-filter QC
p_vln_unfilter <- VlnPlot_scCustom(TS, group.by = "orig.ident", num_columns = 3, pt.size = 0.05, 
                                   features = c("nCount_Spatial","nFeature_Spatial","percent.mt"), assay = "Spatial")
print(p_vln_unfilter)

# 3.2 Execute filtering (Keep spots with nFeature_Spatial > 250)
TS_fliter <- subset(TS, nFeature_Spatial > 250)

# 3.3 Post-filter QC (VlnPlot)
p_vln_filter <- VlnPlot_scCustom(TS_fliter, group.by = "orig.ident", num_columns = 3, pt.size = 0.05, 
                                 features = c("nCount_Spatial","nFeature_Spatial","percent.mt"), assay = "Spatial")
print(p_vln_filter)

# 3.4 Post-filter QC (SpatialPlot)
qc_features_to_plot <- c("nCount_Spatial", "nFeature_Spatial", "percent.mt")
for (feature in qc_features_to_plot) {
  p_spatial_qc <- SpatialFeaturePlot(TS_fliter, features = feature, image.alpha = 0.1, pt.size.factor = 4, crop = FALSE) + 
    theme(legend.position = "right") + 
    scale_fill_viridis(option = "H")
  print(p_spatial_qc)
}

# --- 4. Color Scheme Setup ---
base_colors_paired <- brewer.pal(min(12, brewer.pal.info["Paired", "maxcolors"]), "Paired")
base_colors_dark2 <- brewer.pal(min(8, brewer.pal.info["Dark2", "maxcolors"]), "Dark2")
base_colors_set3 <- brewer.pal(min(12, brewer.pal.info["Set3", "maxcolors"]), "Set3")
base_colors_set1 <- brewer.pal(min(9, brewer.pal.info["Set1", "maxcolors"]), "Set1")
all_available_colors <- c(base_colors_paired, base_colors_dark2, base_colors_set3, base_colors_set1,  
                          colors()[grep("blue|green|orange|purple|yellow|pink|brown", colors())])
color_prgn <- colorRampPalette(rev(brewer.pal(n = 11, name = "PRGn")))(100)

# --- 5. SCTransform Normalization & Clustering ---
options(future.globals.maxSize = 9e10)
options(Seurat.object.assay.version = "v5")

TS_fliter <- SCTransform(TS_fliter, assay = "Spatial", verbose = TRUE, method = "glmGamPoi")
print("SCTransform complete. Starting downstream analysis...")

TS_fliter <- RunPCA(TS_fliter, assay = "SCT", npcs = 30, verbose = TRUE)
TS_fliter <- RunUMAP(TS_fliter, reduction = "pca", dims = 1:20)
TS_fliter <- FindNeighbors(TS_fliter, reduction = "pca", dims = 1:20)
TS_fliter <- FindClusters(TS_fliter, verbose = TRUE, resolution = 0.95)

# 5.1 UMAP Visualization
umap_df <- TS_fliter@reductions$umap@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(seurat_clusters = TS_fliter@active.ident) 

seurat_clusterspos <- umap_df %>%
  group_by(seurat_clusters) %>%
  summarise(umap_1 = median(umap_1), umap_2 = median(umap_2))

p_umap <- DimPlot(TS_fliter, reduction = "umap", label = FALSE, pt.size = 0.8) + 
  theme_dr(xlength = 0.2, ylength = 0.2, arrow = arrow(length = unit(0.2, "inches"), type = "closed")) +
  theme(panel.grid = element_blank(), axis.title = element_text(face = 2, hjust = 0.03)) +
  NoLegend() +
  geom_label_repel(aes(x = umap_1, y = umap_2, label = seurat_clusters), 
                   data = seurat_clusterspos, fontface = "bold", box.padding = 0.5) +
  scale_color_manual(values = all_available_colors)

print(p_umap)

# 5.2 Spatial Cluster Visualization
p_spatial_clusters <- SpatialDimPlot(TS_fliter, label = FALSE, label.size = 3, pt.size.factor = 4.0, crop = FALSE) +
  scale_fill_manual(values = all_available_colors) +
  guides(fill = guide_legend(override.aes = list(size = 2))) +
  theme(legend.title = element_blank())

print(p_spatial_clusters)

# --- 6. Marker Gene DotPlot ---
marker_list <- list(
  "Eryt" = c("zgc-163057", "LOC119891494", "slc4a1a", "cahz", "LOC119912177", "alas2", "tfr1a", "si-ch1073-184j22.1"),
  "Gran" = c("LOC119914447", "LOC119900902", "LOC119889900", "mmp9", "LOC119896105", "mmp25b", "LOC119897780", "si-ch211-284o19.8", "LOC119894963", "LOC119903646", "LOC119904387"),
  "Mac" = c("csf1rb", "marco", "mrc1a", "LOC119891768", "grna", "si-ch211-212k18.7", "LOC119883615", "LOC119915885", "LOC119910309"),
  "B" = c("cd79a", "LOC119884751", "LOC119890508", "LOC119890540", "LOC119890555", "LOC119890570"),
  "T" = c("LOC119910797", "LOC119883393", "LOC119883434", "LOC119914246", "LOC119897409", "LOC119897410"), 
  "E-Mac" = c("LOC119917376", "jupa", "LOC119889888", "dspa", "LOC119902058", "LOC119917243", "dsc2l", "LOC119900663", "cldn11a", "LOC119908668", "LOC119908847"),
  "Fibro-Endo" = c("col1a1b", "col1a2", "col1a1a", "col6a1", "pcolcea", "sparc")
)

p_dotplot <- DotPlot(object = TS_fliter, features = marker_list) +
  scale_color_gradientn(colors = color_prgn) +  
  labs(x = "Cell Clusters", y = "Marker Genes", color = "Avg Expression", size = "Percent Expressed") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9, face = "bold"),
    axis.text.y = element_text(size = 9, face = "italic"),
    axis.title = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    panel.grid = element_blank(),
    legend.position = "right"
  )

print(p_dotplot)

# --- 7. Cell Proportion Calculation & Pie Chart ---
cell_table <- table(Idents(TS_fliter), TS_fliter$orig.ident)
Cellratio <- as.data.frame(prop.table(cell_table, margin = 2))
colnames(Cellratio) <- c("CellType", "Sample", "Freq")

cell_counts <- as.data.frame(cell_table)
colnames(cell_counts) <- c("CellType", "Sample", "Count")
Cellratio <- left_join(Cellratio, cell_counts, by = c("CellType", "Sample"))

Cellratio$CellType <- factor(Cellratio$CellType)
Cellratio <- Cellratio %>% arrange(Sample, CellType)

# Construct legend labels
Cellratio$CellTypeLabel <- paste0(
  Cellratio$CellType, " (", sprintf("%.1f%%", Cellratio$Freq * 100), ", ", Cellratio$Count, ")"
)
# Convert to factor to prevent reordering in plot
Cellratio$CellTypeLabel <- factor(Cellratio$CellTypeLabel, levels = unique(Cellratio$CellTypeLabel))

p_pie <- ggplot(Cellratio, aes(x = "", y = Freq, fill = CellTypeLabel)) +   
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  facet_wrap(~ Sample, nrow = 1) +
  scale_fill_manual(values = all_available_colors) +
  theme_void() +
  labs(fill = "Cell Type (Ratio, Count)", title = "Cell Composition per Sample") +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9),
    plot.title = element_text(hjust = 0.5)
  )

print(p_pie)

# --- 8. Save Data Objects ---
# Change paths accordingly if you still wish to save the final Seurat objects
# saveRDS(object = TS, file = "TS_nomt_dims20_250.rds", compress = FALSE) 
# saveRDS(object = TS_fliter, file = "TS_filter_nomt_dims20_250_0.95.rds", compress = FALSE)