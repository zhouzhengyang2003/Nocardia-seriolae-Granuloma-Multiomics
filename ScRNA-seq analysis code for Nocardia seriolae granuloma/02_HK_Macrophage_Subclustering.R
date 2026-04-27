# ==============================================================================
# Head Kidney Granuloma - 02. Macrophage Sub-clustering & Annotation
# ==============================================================================

suppressMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(RColorBrewer)
})

# Define base output directory
out_dir <- "./results/Macrophage_Subclustering/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 1. Subset Macrophages & Re-cluster
# ==============================================================================
cat("--- 1. Subsetting Macrophages & Re-clustering ---\n")

# Assuming 'scRNA_harmony' is loaded from the previous pipeline
Mac_scRNA_harmony <- subset(scRNA_harmony, idents = "Mac")

# Re-run the standard pipeline on the subset
Mac_scRNA_harmony <- NormalizeData(Mac_scRNA_harmony, normalization.method = "LogNormalize", scale.factor = 1e4) 
Mac_scRNA_harmony <- FindVariableFeatures(Mac_scRNA_harmony, selection.method = 'vst', nfeatures = 1500)
Mac_scRNA_harmony <- ScaleData(Mac_scRNA_harmony, verbose = TRUE)
Mac_scRNA_harmony <- RunPCA(Mac_scRNA_harmony, dims = 1:15) 
Mac_scRNA_harmony <- FindNeighbors(Mac_scRNA_harmony, dims = 1:10)
Mac_scRNA_harmony <- FindClusters(Mac_scRNA_harmony, resolution = 0.3)
Mac_scRNA_harmony <- RunUMAP(Mac_scRNA_harmony, dims = 1:10)

# Define Color Palette
base_colors_paired <- brewer.pal(min(12, brewer.pal.info["Paired", "maxcolors"]), "Paired")
base_colors_dark2 <- brewer.pal(min(8, brewer.pal.info["Dark2", "maxcolors"]), "Dark2")
base_colors_set3 <- brewer.pal(min(12, brewer.pal.info["Set3", "maxcolors"]), "Set3")
base_colors_set1 <- brewer.pal(min(9, brewer.pal.info["Set1", "maxcolors"]), "Set1")
cluster_colors <- c(base_colors_paired, base_colors_dark2, base_colors_set3, base_colors_set1,  
                    colors()[grep("blue|green|orange|purple|yellow|pink|brown", colors())])

# ==============================================================================
# 2. UMAP Visualization (Unannotated)
# ==============================================================================
cat("--- 2. Generating Initial UMAPs ---\n")

umap_df <- Mac_scRNA_harmony@reductions$umap@cell.embeddings %>%  
  as.data.frame() %>% cbind(seurat_clusters = Mac_scRNA_harmony@active.ident)

seurat_clusterspos <- umap_df %>% group_by(seurat_clusters) %>%
  summarise(umap_1 = median(umap_1), umap_2 = median(umap_2))

p_umap_base <- DimPlot(Mac_scRNA_harmony, reduction = "umap", label = FALSE, pt.size = 1.5) +
  theme_dr(xlength = 0.2, ylength = 0.2, arrow = arrow(length = unit(0.2, "inches"), type = "closed")) +
  theme(panel.grid = element_blank(), axis.title = element_text(face = 2, hjust = 0.03)) +
  NoLegend() +
  geom_label_repel(aes(x = umap_1, y = umap_2, label = seurat_clusters), data = seurat_clusterspos, fontface = "bold", box.padding = 0.5, size = 5) +
  scale_color_manual(values = cluster_colors)

ggsave(file.path(out_dir, "1_UMAP_Unannotated.pdf"), p_umap_base, width = 9, height = 9)

p_umap_split <- p_umap_base + facet_wrap(~ Mac_scRNA_harmony$orig.ident)
ggsave(file.path(out_dir, "2_UMAP_Unannotated_Split.pdf"), p_umap_split, width = 16, height = 9)

# ==============================================================================
# 3. Find Markers & Calculate Average Expression
# ==============================================================================
cat("--- 3. Finding Markers ---\n")

Mac_markers <- FindAllMarkers(Mac_scRNA_harmony, only.pos = TRUE, min.pct = 0.1, test.use = "wilcox", logfc.threshold = 0.2)
significant_markers <- Mac_markers[Mac_markers$p_val_adj < 0.1,]
write.csv(significant_markers, file = file.path(out_dir, "Markers_Unannotated_res0.3.csv"))

top20_markers <- Mac_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.csv(top20_markers, file = file.path(out_dir, "Markers_Top20_Unannotated_res0.3.csv"), row.names = TRUE)

cluster_avg_exp <- AverageExpression(Mac_scRNA_harmony, slot = "data")
write.csv(cluster_avg_exp$RNA, file = file.path(out_dir, "Average_Expression_RNA.csv"))

# ==============================================================================
# 4. Cell Type Annotation
# ==============================================================================
cat("--- 4. Annotating Sub-clusters ---\n")

new.cluster.ids <- c(
  "Microenvironment-Regulating Mac", "GSTM3+ CTH+ Mac", "Homeostatic Mac",  
  "Proliferating Mac2", "Proliferating Mac1", "Mac/Gran-like",  
  "Mac/Eryt-like",  "MHCII+ Mac"
)

names(new.cluster.ids) <- levels(Mac_scRNA_harmony)
Mac_scRNA_harmony <- RenameIdents(Mac_scRNA_harmony, new.cluster.ids)

# Reorder Idents based on your preference
desired_order <- c("5", "8", "3", "7", "6", "1", "2", "0", "9", "4") # Update this array if you need the named idents instead
# Mac_scRNA_harmony <- SetIdent(Mac_scRNA_harmony, value = factor(Idents(Mac_scRNA_harmony), levels = desired_order))

# ==============================================================================
# 5. Cell Proportions (Bar & Pie Charts)
# ==============================================================================
cat("--- 5. Generating Proportion Plots ---\n")

cell_table <- table(Idents(Mac_scRNA_harmony), Mac_scRNA_harmony$orig.ident)
Cellratio <- as.data.frame(prop.table(cell_table, margin = 2))
colnames(Cellratio) <- c("CellType", "Sample", "Freq")

cell_counts <- as.data.frame(cell_table)
colnames(cell_counts) <- c("CellType", "Sample", "Count")
Cellratio <- left_join(Cellratio, cell_counts, by = c("CellType", "Sample"))

# Stacked Bar Plot
p_bar <- ggplot(Cellratio) + 
  geom_bar(aes(x = Sample, y = Freq, fill = CellType), stat = "identity", width = 0.7, size = 0.5, colour = '#222222') + 
  theme_classic() + labs(x = 'Sample', y = 'Ratio') + scale_fill_manual(values = cluster_colors) +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"))

ggsave(file.path(out_dir, "3_Proportions_Barplot.pdf"), p_bar, width = 4, height = 6)

# Pie Chart with Text Labels
Cellratio$TextLabel <- paste0(sprintf("%.1f%%", Cellratio$Freq * 100), "\n(", Cellratio$Count, ")")
Cellratio$TextLabel <- ifelse(Cellratio$Freq > 0.02, Cellratio$TextLabel, "")

p_pie <- ggplot(Cellratio, aes(x = "", y = Freq, fill = CellType)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") + facet_wrap(~ Sample, nrow = 1) +
  geom_text(aes(label = TextLabel), position = position_stack(vjust = 0.5), size = 3, fontface = "bold") +
  scale_fill_manual(values = cluster_colors) + 
  theme_void() + labs(fill = "Cell Type", title = "Macrophage Composition") +
  theme(strip.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 12), legend.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold"))

ggsave(file.path(out_dir, "4_Proportions_Piechart.pdf"), p_pie, width = 9, height = 5)

# ==============================================================================
# 6. Functional DotPlots (Automated Helper Function)
# ==============================================================================
cat("--- 6. Generating Functional DotPlots ---\n")

plot_mac_dotplot <- function(features_list, filename, width=7.5, height=3.5) {
  p <- DotPlot(Mac_scRNA_harmony, features = features_list) + RotatedAxis() +
    scale_color_gradientn(colors = c("lightblue", "white", "red")) +
    theme_classic(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9, face = "bold"),
      axis.text.y = element_text(size = 9, face = "italic"),
      axis.title = element_text(size = 10, face = "bold"),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      panel.grid = element_blank(), 
      legend.position = "right"
    )
  ggsave(file.path(out_dir, filename), p, width = width, height = height)
}

# Define and Plot Module Lists
list_junctions <- list(`AJ` = c("gene-LOC119917376","gene-cdh11","gene-jupa", "gene-ctnna1", "gene-ctnnd1"), `DSM` = c("gene-dsc2l", "gene-LOC119902058", "gene-LOC119905441", "gene-LOC119889888", "gene-dspa", "gene-LOC119917243"), `TJ` = c("gene-cldn5a", "gene-LOC119900663", "gene-cldn11a", "gene-tjp1a", "gene-cgnl1"), `GJ` = c("gene-LOC119908847", "gene-LOC119908668"))
plot_mac_dotplot(list_junctions, "DotPlot_Cell_Junctions.pdf")

list_cyto <- list(`Intermediate Filament` = c("gene-LOC119905733", "gene-LOC119886877", "gene-krt15", "gene-LOC119886876", "gene-krt94", "gene-tgm1l1"), `Actin` = c("gene-LOC119885789", "gene-pof1b", "gene-crip1", "gene-LOC119901854"), `Membrane Coupling` = c("gene-ahnak", "gene-LOC119888099", "gene-anxa2a", "gene-LOC119890527", "gene-LOC119898827"), `Calcium` = c("gene-s100a10a", "gene-chp1", "gene-cast", "gene-capns1b", "gene-LOC119912040", "gene-LOC119915637"), `Rho GTPase` = c("gene-plekhg5b", "gene-rnd3a", "gene-cdc42ep3", "gene-ppp1r27b"))
plot_mac_dotplot(list_cyto, "DotPlot_Cytoskeleton.pdf")

list_ecm <- list(`Proteases` = c("gene-mmp2", "gene-mmp15b", "gene-mmp19"), `Inhibitors` = c("gene-serpine1", "gene-serpine2", "gene-LOC119908472", "gene-LOC119906979"), `Assembly` = c("gene-loxa", "gene-bmp1a", "gene-LOC119897617", "gene-sparc", "gene-tnfaip6"), `Adhesion` = c("gene-cd44b", "gene-fn1a", "gene-fn1b", "gene-LOC119897849", "gene-cavin1b", "gene-LOC119907497", "gene-LOC119918071", "gene-LOC119893200", "gene-LOC119895007", "gene-LOC119899685"), `Integrins` = c("gene-itgav", "gene-LOC119892793", "gene-LOC119915730"))
plot_mac_dotplot(list_ecm, "DotPlot_ECM.pdf")

list_upstream <- list(`Microdomain` = c("gene-spaca4l", "gene-LOC119906807"), `Drivers` = c("gene-wt1b", "gene-aldh1a2", "gene-gata6","gene-LOC119894306", "gene-LOC119914475"))
plot_mac_dotplot(list_upstream, "DotPlot_Upstream_Regulators.pdf", width=4)

list_imm <- list(`Complement` = c("gene-LOC119885585","gene-c1qb","gene-LOC119885788","gene-LOC119897430","gene-LOC119892056", "gene-LOC119915106","gene-LOC119896603"), `MHCII` = c("gene-LOC119888605","gene-LOC119917015","gene-LOC119912422","gene-LOC119912423", "gene-LOC119912429","gene-LOC119912010"), `Cathepsins` = c("gene-ctsba","gene-ctsc","gene-LOC119890600","gene-LOC119890271","gene-LOC119882211", "gene-ctsk", "gene-ctsl.1","gene-LOC119892265","gene-ctss1", "gene-ctss2.1","gene-LOC119882210","gene-zgc-110239", "gene-ctsz", "gene-napsa"), `Inflammation` = c("gene-LOC119897691","gene-tnfa", "gene-m17", "gene-LOC119888199","gene-LOC119888219", "gene-LOC119890263","gene-il1b"))
plot_mac_dotplot(list_imm, "DotPlot_Immune_Response.pdf", width=8)

list_reg <- list(`Inhibition` = c("gene-LOC119910060", "gene-LOC119910049","gene-cst3", "gene-LOC119882640", "gene-cast"), `Regulation` = c("gene-LOC119884194", "gene-nfkbiaa"), `Stress` = c("gene-LOC119910277", "gene-ywhag1","gene-LOC119891659", "gene-pkig", "gene-tsc22d1", "gene-LOC119893531", "gene-LOC119893593"), `Transcription` = c("gene-bhlhe40", "gene-hic1"))
plot_mac_dotplot(list_reg, "DotPlot_Immune_Regulation.pdf", width=8)

list_oxphos <- list("COX" = c("gene-COX1", "gene-COX2", "gene-COX3", "gene-LOC119905811", "gene-LOC119890004", "gene-LOC119917832", "gene-LOC119892218", "gene-LOC119893803", "gene-LOC119897473", "gene-LOC119911148", "gene-LOC119885605", "gene-LOC119916265", "gene-LOC119896019", "gene-LOC119892403", "gene-LOC119893283", "gene-LOC119910933"), "ATPase" = c("gene-ATP6", "gene-LOC119891187", "gene-atp5f1c", "gene-atp5f1d", "gene-atp5f1e", "gene-atp5po", "gene-LOC119890133", "gene-atp5mc1", "gene-atp5pd", "gene-LOC119899395", "gene-LOC119914361", "gene-atp5mf", "gene-LOC119901952", "gene-atp5pf", "gene-atp5mj", "gene-atp5md"), "Reductase" = c("gene-LOC119900817", "gene-LOC119889032", "gene-CYTB"))
plot_mac_dotplot(list_oxphos, "DotPlot_OxPhos.pdf", width=8)

list_metab <- list("Amino Acid" = c("gene-LOC119892499", "gene-aass", "gene-LOC119890054"), "Lipid" = c("gene-LOC119903459", "gene-LOC119898927", "gene-fabp4a", "gene-plin2"), "Nucleotide" = c("gene-cdab", "gene-nudt4b"), "Transport" = c("gene-LOC119893713", "gene-LOC119910978"), "Homeostasis" = c("gene-LOC119912177"), "Adaptation" = c("gene-ndrg1a", "gene-nupr1b", "gene-egln3", "gene-LOC119900507"))
plot_mac_dotplot(list_metab, "DotPlot_Metabolism.pdf", width=8)

list_trans <- list("Ions" = c("gene-LOC119901983", "gene-rhbg", "gene-LOC119904410", "gene-atp1b1a", "gene-LOC119912124", "gene-slc12a7b", "gene-LOC119908031", "gene-chp1", "gene-LOC119918691", "gene-slc39a8"), "Amino/Neuro" = c("gene-slc6a8", "gene-LOC119905556"), "Vitamins" = c("gene-slc23a2"), "Nucleosides" = c("gene-slc28a1"), "ATP" = c("gene-si-dkey-251i10.1"), "Unknown" = c("gene-slc10a3", "gene-slc25a55a"))
plot_mac_dotplot(list_trans, "DotPlot_Transporters.pdf", width=8)

list_fibro <- list(`Pro-fibrotic` = c("gene-LOC119888537", "gene-bmp1a","gene-p4ha1b","gene-loxa","gene-LOC119897617","gene-LOC119886711"))
plot_mac_dotplot(list_fibro, "DotPlot_ProFibrotic.pdf", width=8)

# ==============================================================================
# 7. DoHeatmap for Curated Markers
# ==============================================================================
cat("--- 7. Generating Heatmap ---\n")

heatmap_markers <- unique(c("gene-LOC119890038", "gene-pdia6", "gene-calr", "gene-hspa5","gene-si-dkeyp-118h9.7", "gene-LOC119891167","gene-LOC119882585","gene-LOC119891533","gene-LOC119882652", "gene-cdk1", "gene-mki67","gene-ube2c", "gene-LOC119900344", "gene-tubb2b","gene-cenpw", "gene-cdca2","gene-kif22","gene-top2a", "gene-LOC119913275", "gene-si-ch73-281n10.2", "gene-pclaf","gene-LOC119889212", "gene-LOC119899816","gene-dtymk", "gene-LOC119915246", "gene-LOC119900609", "gene-cct7", "gene-LOC119907150","gene-LOC119900894", "gene-LOC119895682","gene-hsp90ab1", "gene-slc4a1a", "gene-sptb","gene-cahz", "gene-zgc-163057", "gene-ank1a", "gene-LOC119897155", "gene-LOC119914170", "gene-LOC119912414", "gene-LOC119888605","gene-LOC119917015", "gene-LOC119912422", "gene-LOC119912423","gene-LOC119912424","gene-LOC119912429","gene-LOC119912432", "gene-ctsba", "gene-LOC119890271", "gene-LOC119885585", "gene-c1qb","gene-LOC119885788", "gene-bpifcl", "gene-LOC119915906", "gene-mcl1b", "gene-LOC119915880","gene-nfkbiaa", "gene-LOC119910049",  "gene-cst3", "gene-LOC119890527",  "gene-LOC119898827", "gene-LOC119917376", "gene-jupa", "gene-cldn11a","gene-LOC119900663", "gene-LOC119902058", "gene-dsc2l",  "gene-LOC119889888","gene-dspa", "gene-LOC119908668", "gene-LOC119908847", "gene-LOC119886876", "gene-LOC119905733","gene-krt15", "gene-krt94",  "gene-LOC119886877", "gene-LOC119914447", "gene-LOC119900902","gene-LOC119890363", "gene-LOC119894963", "gene-mmp9"))

Mac_scRNA_harmony <- ScaleData(Mac_scRNA_harmony, features = heatmap_markers, assay = "RNA")

p_heatmap <- DoHeatmap(Mac_scRNA_harmony, features = heatmap_markers, assay = 'RNA', group.colors = cluster_colors, label = TRUE) + 
  scale_fill_gradientn(colors = c("lightblue", "white",  "red"), name = "Scaled Expression") + 
  theme(axis.text.y = element_text(size = 10, face = "italic"), legend.title = element_text(size = 10), legend.text = element_text(size = 10))

ggsave(file.path(out_dir, "5_Heatmap_CuratedMarkers.pdf"), p_heatmap, width = 14, height = 16)

# ==============================================================================
# 8. Save Final Object
# ==============================================================================
saveRDS(Mac_scRNA_harmony, file = file.path(out_dir, "Mac_scRNA_harmony_Final.rds"))
cat("✅ Macrophage Sub-clustering Pipeline Complete!\n")