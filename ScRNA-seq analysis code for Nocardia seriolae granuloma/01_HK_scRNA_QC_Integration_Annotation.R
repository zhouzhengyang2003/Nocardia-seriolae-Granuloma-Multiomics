# ==============================================================================
# Head Kidney Granuloma - 8. scRNA-seq Processing, Integration, and Annotation
# ==============================================================================

# --- 0. Load Required Packages ---
suppressMessages({
  library(Seurat)
  library(tidyverse)
  library(scCustomize)
  library(cowplot)
  library(DoubletFinder)
  library(harmony)
  library(R.utils)
  library(ggrepel)
  library(ggsci)
  library(viridis)
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(scales)
  library(gghalves)
  library(ggpubr)
})

# Define Base Paths (Simulated relative paths for portability)
input_dir <- "./data/SC/TS/"
out_dir <- "./results/scRNA_Analysis/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 1. Batch Load scRNA Data & Initialize Seurat Objects
# ==============================================================================
cat("--- 1. Loading scRNA Data ---\n")
dir_names <- list.files(input_dir)
scRNAlist <- list()

for(i in seq_along(dir_names)){
  counts <- Read10X(data.dir = file.path(input_dir, dir_names[i]))
  scRNAlist[[i]] <- CreateSeuratObject(counts, project = dir_names[i], min.cells = 3, min.features = 100)
}

# Print Basic Stats
for (i in seq_along(scRNAlist)) {
  counts_matrix <- GetAssayData(scRNAlist[[i]], assay = "RNA", slot = "counts")
  total_counts <- Matrix::colSums(counts_matrix)
  detected_genes <- Matrix::colSums(counts_matrix > 0)
  
  cat(sprintf("Sample %d (%s):\n", i, dir_names[i]))
  cat(sprintf("  Mean reads per cell: %d\n", round(mean(total_counts))))
  cat(sprintf("  Median genes per cell: %d\n\n", median(detected_genes)))
}

# ==============================================================================
# 2. Quality Control (QC): Ribosomal & Mitochondrial Gene Filtering
# ==============================================================================
cat("--- 2. Calculating QC Metrics & Filtering ---\n")
mito_genes <- c("gene-ND1", "gene-ND2", "gene-COX1", "gene-COX2", "gene-ATP8", "gene-ATP6", 
                "gene-COX3", "gene-ND3", "gene-ND4L", "gene-ND4", "gene-ND5", "gene-ND6", "gene-CYTB")

ribo_genes <- c("gene-mrpl39", "gene-rpl31", "gene-mrps9", "gene-rps20", "gene-rpl7", "gene-rpl21",
                "gene-mrpl15", "gene-rpl22l1", "gene-LOC119898710", "gene-LOC119898703", "gene-mrpl45",
                "gene-mrpl4", "gene-rpl23", "gene-rpsa", "gene-mrps34", "gene-mrps7", "gene-rps2",
                "gene-LOC119888983", "gene-mrpl12", "gene-rps11", "gene-rpl19", "gene-rpl27",
                "gene-rpl3", "gene-mrps15", "gene-rps19", "gene-mrpl17", "gene-LOC119897373",
                "gene-mrpl51", "gene-rps9", "gene-rpl18", "gene-dap3", "gene-rps27.2",
                "gene-LOC119897812", "gene-LOC119897846", "gene-mrpl32", "gene-LOC119898192",
                "gene-LOC119898191", "gene-rps3a", "gene-rps6", "gene-rpl9", "gene-rps15a",
                "gene-mrps28", "gene-LOC119900488", "gene-mrpl47", "gene-LOC119900749",
                "gene-mrpl54", "gene-LOC119901212", "gene-rpl36", "gene-LOC119901312",
                "gene-LOC119901314", "gene-rps8a", "gene-mrps14", "gene-LOC119902226",
                "gene-rps15", "gene-mrpl34", "gene-rpl5a", "gene-rpl37a", "gene-mrpl37",
                "gene-LOC119903249", "gene-mrpl30", "gene-rpl8", "gene-rps21", "gene-mrps25",
                "gene-rpl32", "gene-rpl10a", "gene-rpl10", "gene-mrpl49", "gene-LOC119904990",
                "gene-rps18", "gene-rpl11", "gene-rps5", "gene-mrps21", "gene-rps27.1",
                "gene-LOC119907920", "gene-mrpl24", "gene-mrpl9", "gene-rplp1", "gene-rpl30",
                "gene-mrpl53", "gene-mrpl3", "gene-rpl14", "gene-mrpl36", "gene-rpl15",
                "gene-mrps18b", "gene-rpl7l1", "gene-faua", "gene-rps29", "gene-mrps18a",
                "gene-rps12", "gene-mrps10", "gene-mrpl19", "gene-rps7", "gene-rpl13a",
                "gene-mrpl35", "gene-mrpl33", "gene-rps27a", "gene-rpl34", "gene-mrpl52",
                "gene-mrpl48", "gene-rps4x", "gene-LOC119914721", "gene-LOC119915000",
                "gene-mrps18c", "gene-mrps24", "gene-mrps2", "gene-rpl7a", "gene-rpl28",
                "gene-rpl37", "gene-LOC119914137", "gene-LOC119915072", "gene-mrpl27",
                "gene-mrpl58", "gene-mrpl38", "gene-LOC119916482", "gene-rpl38",
                "gene-mrpl43", "gene-mrpl10", "gene-mrpl2", "gene-mrpl14", "gene-mrps6",
                "gene-rps24", "gene-mrpl57", "gene-mrps5", "gene-mrps33", "gene-mrpl42",
                "gene-rps16", "gene-LOC119918593", "gene-LOC119919179", "gene-mrps23",
                "gene-mrps22", "gene-rps25", "gene-mrpl44", "gene-mrpl28", "gene-rps3",
                "gene-LOC119882934", "gene-rpl24", "gene-rpl23a", "gene-LOC119884849",
                "gene-rps23", "gene-rpl22", "gene-rps10", "gene-rplp2", "gene-rpl29",
                "gene-mrps16", "gene-mrpl20", "gene-LOC119886813", "gene-rplp0",
                "gene-mrps26", "gene-rpl26", "gene-mrpl18", "gene-mrpl11", "gene-rpl39",
                "gene-rpl36a", "gene-rpl12", "gene-LOC119890837", "gene-mrpl16",
                "gene-LOC119889915", "gene-rpl13", "gene-mrpl46", "gene-mrps11",
                "gene-rpl35", "gene-mrpl41", "gene-mrpl1", "gene-mrps30", "gene-rpl6",
                "gene-mrpl40", "gene-mrps36", "gene-mrps27", "gene-rps28", "gene-LOC119894123",
                "gene-rpl27a", "gene-rps13", "gene-mrps35", "gene-rpl4", "gene-rps27l",
                "gene-mrpl23", "gene-LOC119894724", "gene-rpl18a", "gene-rps14",
                "gene-mrpl22", "gene-mrps17", "gene-mrps31") 

# Remove Ribosomal Genes & Calculate Mito Ratio
for(i in seq_along(scRNAlist)){
  sc <- scRNAlist[[i]]
  sc <- subset(sc, features = setdiff(rownames(sc), ribo_genes))
  sc[['percent.mt']] <- PercentageFeatureSet(sc, features = mito_genes)
  scRNAlist[[i]] <- sc
}

# Pre-filter Violin Plots
violin_before <- lapply(scRNAlist, function(x) {
  VlnPlot(x, layer = 'counts', features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
})
p_vln_before <- wrap_plots(violin_before, ncol = 1)
ggsave(file.path(out_dir, "1_QC_BeforeFilter_Vln.pdf"), p_vln_before, width = 6, height = 12)

# Execute Filtering
scRNAlist_filtered <- lapply(scRNAlist, function(x) {
  subset(x, subset = nFeature_RNA > 300 & nFeature_RNA < 5500 & percent.mt < 30 & nCount_RNA > 100)
})

# Post-filter Violin Plots
violin_after <- lapply(scRNAlist_filtered, function(x) {
  VlnPlot(x, layer = 'counts', features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
})
p_vln_after <- wrap_plots(violin_after, ncol = 1)
ggsave(file.path(out_dir, "2_QC_AfterFilter_Vln.pdf"), p_vln_after, width = 6, height = 12)

# ==============================================================================
# 3. DoubletFinder Processing
# ==============================================================================
cat("--- 3. Running DoubletFinder ---\n")
scRNA_rm_doublet <- list()

for(i in seq_along(scRNAlist_filtered)){
  sc <- scRNAlist_filtered[[i]]
  
  # Standard Pre-processing
  sc <- NormalizeData(sc) %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData() %>% RunPCA(verbose = FALSE) %>% RunUMAP(dims = 1:20, verbose = FALSE)
  
  # DoubletFinder Parameter Sweep
  sweep.res.list <- paramSweep(sc, PCs = 1:20, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  best_pk <- as.numeric(as.character(bcmvn$pK[bcmvn$BCmetric == max(bcmvn$BCmetric)]))
  nExp_poi <- round(((ncol(sc) / 1000) * 0.008) * ncol(sc))
  
  # Run DoubletFinder
  sc <- doubletFinder(sc, PCs = 1:20, pN = 0.25, pK = best_pk, nExp = nExp_poi, reuse.pANN = NULL, sct = FALSE)
  
  # Visualization
  df_classification_col <- names(sc@meta.data)[grep("DF.classifications", names(sc@meta.data))]
  p_df <- DimPlot(sc, group.by = df_classification_col, raster=FALSE) + ggtitle(paste0("DoubletFinder - ", Project(sc)))
  ggsave(file.path(out_dir, paste0("3_DoubletFinder_Check_", Project(sc), ".pdf")), p_df, width = 8, height = 6)
  
  # Extract Singlets ONLY
  singlet_cells_ids <- rownames(sc@meta.data)[sc@meta.data[[df_classification_col]] == "Singlet"]
  sc.filtered <- sc[, singlet_cells_ids]
  
  scRNA_rm_doublet[[i]] <- sc.filtered
  cat(sprintf("Sample %s processed: Removed %d doublets.\n", Project(sc), ncol(sc) - ncol(sc.filtered)))
}

names(scRNA_rm_doublet) <- names(scRNAlist_filtered)
saveRDS(scRNA_rm_doublet, file = file.path(out_dir, "scRNA_list_filtered_rmDoublets.rds"))

# ==============================================================================
# 4. Merging & Harmony Integration (Seurat V5)
# ==============================================================================
cat("--- 4. Merging & Harmony Integration ---\n")
scRNA.merged <- merge(scRNA_rm_doublet[[1]], y = scRNA_rm_doublet[2:length(scRNA_rm_doublet)], add.cell.ids = names(scRNA_rm_doublet))

scRNA.merged <- NormalizeData(scRNA.merged, normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000, verbose = FALSE) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE)

# Harmony Integration
scRNA.integrated <- IntegrateLayers(
  object = scRNA.merged,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony",
  group.by = "orig.ident",
  verbose = TRUE
)

scRNA.integrated[["RNA"]] <- JoinLayers(scRNA.integrated[["RNA"]])

# Clustering & UMAP
scRNA_harmony <- scRNA.integrated %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20, verbose = TRUE) %>%
  FindClusters(resolution = 0.1, verbose = TRUE) %>%
  RunUMAP(dims = 1:20, reduction = "harmony", verbose = TRUE)

saveRDS(scRNA_harmony, file = file.path(out_dir, "scRNA_harmony_Unannotated.rds"))

# ==============================================================================
# 5. Cell Type Annotation
# ==============================================================================
cat("--- 5. Annotating Cell Types ---\n")

# Identify Markers
scRNA_harmony.markers <- FindAllMarkers(scRNA_harmony, only.pos = TRUE, min.pct = 0.1, test.use = "wilcox", logfc.threshold = 0.2)
write.csv(scRNA_harmony.markers[scRNA_harmony.markers$p_val_adj < 0.1, ], file = file.path(out_dir, "significant_markers_Unannotated.csv"))

# Rename Clusters
new.cluster.ids <- c("B", "Gran", "Mac", "Eryt", "Eryt", "B", "HSPCs", "T", "B", "DC", "Endo")
names(new.cluster.ids) <- levels(scRNA_harmony)
scRNA_harmony <- RenameIdents(scRNA_harmony, new.cluster.ids)

# Reorder Idents
desired_order <- c("HSPCs", "Eryt", "Gran", "Mac", "DC", "B", "T", "Endo")
scRNA_harmony <- SetIdent(scRNA_harmony, value = factor(Idents(scRNA_harmony), levels = desired_order))
scRNA_harmony$annotated_cluster <- Idents(scRNA_harmony)

# ==============================================================================
# 6. Visualization & Statistics
# ==============================================================================
cat("--- 6. Generating Final Visualizations ---\n")

celltype_colors <- c(
  "HSPCs" = "#A0522D", "Eryt" = "#FFB6C1", "Gran"  = "#FF8C00", "Mac"   = "#228B22",
  "DC"    = "#98FB98", "B"    = "#FFFFCC", "T"     = "#A6CEE3", "Endo"  = "cyan"
)

# --- UMAP Visualizations ---
umap_df <- scRNA_harmony@reductions$umap@cell.embeddings %>%
  as.data.frame() %>% cbind(seurat_clusters = scRNA_harmony@active.ident)

seurat_clusterspos <- umap_df %>% group_by(seurat_clusters) %>%
  summarise(umap_1 = median(umap_1), umap_2 = median(umap_2))

p_umap_base <- DimPlot(scRNA_harmony, reduction = "umap", label = FALSE, pt.size = 0.9) +
  theme_dr(xlength = 0.2, ylength = 0.2, arrow = arrow(length = unit(0.2, "inches"), type = "closed")) +
  theme(panel.grid = element_blank(), axis.title = element_text(face = 2, hjust = 0.03)) +
  NoLegend() + labs(title = NULL) +
  geom_label_repel(aes(x = umap_1, y = umap_2, label = seurat_clusters), data = seurat_clusterspos, fontface = "bold", box.padding = 0.5)

p_umap_annotated <- p_umap_base + scale_color_manual(values = celltype_colors)
ggsave(file.path(out_dir, "4_Annotated_UMAP.pdf"), p_umap_annotated, width = 8, height = 8)

p_umap_split <- p_umap_base + facet_wrap(~ scRNA_harmony$orig.ident) + scale_color_manual(values = celltype_colors)
ggsave(file.path(out_dir, "5_Annotated_UMAP_Split.pdf"), p_umap_split, width = 16, height = 9)

# --- DotPlot Visualization ---
marker_list <- list(
  "HSPCs" = c("gene-tal1","gene-gata2b","gene-meis1b","gene-gfi1b","gene-zeb2a","gene-fli1","gene-lmo2"),
  "Eryt"  = c("gene-zgc-163057","gene-LOC119891494","gene-slc4a1a","gene-cahz","gene-LOC119912177","gene-alas2","gene-tfr1a","gene-si-ch1073-184j22.1"),
  "Gran"  = c("gene-LOC119914447", "gene-LOC119900902", "gene-LOC119889900", "gene-mmp9", "gene-LOC119896105","gene-mmp25b", "gene-LOC119897780", "gene-si-ch211-284o19.8", "gene-LOC119894963","gene-LOC119903646","gene-LOC119904387"),
  "Mac"   = c("gene-csf1rb", "gene-marco" ,"gene-LOC119901899", "gene-si-ch211-212k18.7", "gene-LOC119883615", "gene-LOC119910309","gene-lgals2a"),
  "DC"    = c("gene-spi1a",  "gene-LOC119897155", "gene-LOC119914170","gene-mfge8b", "gene-cd83", "gene-irf8", "gene-zbtb46", "gene-lgmn","gene-ifi30", "gene-LOC119888605", "gene-LOC119917015", "gene-LOC119912422","gene-LOC119912429", "gene-LOC119912010"),
  "B"     = c("gene-cd79a","gene-LOC119884751","gene-LOC119895379", "gene-pax5","gene-sox4a",  "gene-rag1", "gene-rag2","gene-LOC119887406", "gene-LOC119890508", "gene-LOC119890540", "gene-LOC119890555", "gene-LOC119890570"),
  "T"     = c("gene-LOC119883393","gene-LOC119883434","gene-LOC119914246","gene-LOC119897409","gene-LOC119897410","gene-cd8a","gene-cd8b"),
  "Endo"  = c("gene-eng","gene-cdh5","gene-kdrl","gene-gpr182","gene-sox18","gene-id1")
)

p_dotplot <- DotPlot(scRNA_harmony, features = marker_list) + 
  scale_color_viridis(option = "viridis", direction = -1) +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9, face = "bold"),
        axis.text.y = element_text(size = 9, face = "italic"), axis.title = element_text(size = 10, face = "bold"),
        legend.title = element_text(size = 10), legend.text = element_text(size = 8),
        panel.grid = element_blank(), legend.position = "right")

ggsave(file.path(out_dir, "6_Annotated_DotPlot.pdf"), p_dotplot, width = 16, height = 4)

# --- Cell Proportions Barplot ---
cell_table <- table(Idents(scRNA_harmony), scRNA_harmony$orig.ident)
Cellratio <- as.data.frame(prop.table(cell_table, margin = 2))
colnames(Cellratio) <- c("CellType", "Sample", "Ratio")

p_bar <- ggplot(Cellratio) + 
  geom_bar(aes(x = Sample, y = Ratio, fill = CellType), stat = "identity", width = 0.7, linewidth = 0.5, colour = '#222222') + 
  theme_classic() +
  labs(x = 'Sample', y = 'Ratio') +
  scale_fill_manual(values = celltype_colors) +
  theme(panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5, linetype = "solid"))

ggsave(file.path(out_dir, "7_Cell_Proportions.pdf"), p_bar, width = 4, height = 6)

# ==============================================================================
# 7. Save Final Annotated Object
# ==============================================================================
saveRDS(scRNA_harmony, file = file.path(out_dir, "scRNA_harmony_Annotated_Final.rds"))
cat("✅ Pipeline Completed Successfully!\n")