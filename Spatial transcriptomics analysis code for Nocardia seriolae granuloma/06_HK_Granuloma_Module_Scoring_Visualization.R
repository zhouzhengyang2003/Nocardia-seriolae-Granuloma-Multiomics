# ==============================================================================
# Head Kidney Granuloma - 6. Module Scoring & Comprehensive Visualization
# ==============================================================================
suppressMessages({
  library(Seurat)
  library(ggplot2)
  library(ggpubr)
  library(gghalves)
  library(viridis)
  library(RColorBrewer)
  library(dplyr)
})

DefaultAssay(TS_fliter) <- "SCT"

# Set base saving directory
save_base <- "./results/Module_Scores/E-Mac/"
dir.create(save_base, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 1. Consolidate All Signatures into ONE Master List
# ==============================================================================
master_signature_list <- list(
  # 1. Cell Junctions
  `Cell Junctions` = c("LOC119917376", "jupa", "ctnna1", "ctnnd1", "dsc2l", "LOC119902058", "LOC119905441", "LOC119889888", "dspa", "LOC119917243", "cldn5a", "LOC119900663", "cldn11a", "tjp1a", "cgnl1", "LOC119908847", "LOC119908668"),
  
  # 2. Cytoskeleton Regulation
  `Cytoskeleton Regulation` = c("LOC119905733", "LOC119886877", "krt15", "LOC119886876", "krt94", "krt99", "eppk1", "LOC119912524", "LOC119886962", "tgm1l1", "LOC119885789", "LOC119909726", "LOC119895369", "pof1b","LOC119901854", "crip1","crip2","nudcd2", "ahnak", "LOC119909740", "ezrb", "ehd2b", "stoml3b", "LOC119914141", "mpp6b", "anxa1a", "LOC119888099", "anxa2a", "LOC119890527", "LOC119898827", "anxa11a", "anxa13", "s100a10a", "si-ch211-105c13.3", "s100v2", "calm1b","capns1b", "LOC119912040", "LOC119915637", "LOC119886455", "efnb2a","plekhg5b", "limk2", "cdc42ep3","ppp2cab", "ppp1r27b", "rnd3a", "plekhn1","LOC119890724"),
  
  # 3. Extracellular Matrix (ECM)
  `Extracellular Matrix` = c("mmp2", "mmp15b", "mmp19", "adam8a", "adam11", "LOC119886402", "serpine1", "serpine2", "LOC119908472", "LOC119906979", "timp2a", "timp2b", "LOC119918277", "ero1a","p4hb","cd44b", "tnfaip6", "fn1a", "fn1b", "LOC119897849", "cav1", "cavin1b", "LOC119907497", "LOC119918071", "cd81a", "LOC119893200", "LOC119895007","LOC119899685", "vwa11", "anos1b", "aplp2", "prom2", "LOC119886359", "itga3b", "itga5", "LOC119902203", "itgav", "LOC119901714", "LOC119892793", "LOC119915730", "itgb6"),
  
  # 4. Potential Upstream Regulators
  `Upstream Regulators` = c("spaca4l", "LOC119906807", "zanl","wt1b", "gata6","LOC119894306", "LOC119914475","aldh1a2", "gprc5c", "LOC119884110","crabp2a"),
  
  # 5. Inflammation-Related
  `Inflammation` = c("LOC119885585","c1qb","LOC119885788","LOC119897430","LOC119892056", "LOC119915106","LOC119896603", "LOC119888605","LOC119917015","LOC119912422","LOC119912429","LOC119912010", "ctsba","ctsc","LOC119890600","LOC119890271","ctsh","ctsk","LOC119882211","ctsl.1","LOC119892265","ctss1", "ctss2.1","LOC119882210","ctsz", "napsa", "tnfa", "m17","LOC119894803","LOC119890263","il1b"),
  
  # 6. Immune Regulation & Homeostasis
  `Immune Regulation` = c("LOC119910060","LOC119910049","cst3", "LOC119882640","cast", "LOC119884194","LOC119884197","nfkbiaa","il10","tgfb1a", "mcl1b","bcl2l1","si-dkey-243k1.3","babam2","LOC119910277","ywhag1","LOC119891659","rcan1b", "pkig","tsc22d1","tsc22d2","LOC119915309","LOC119893531","LOC119893593", "bhlhe40","hic1","LOC119899642","junbb"),
  
  # 7. Oxidative Phosphorylation
  `Oxidative Phosphorylation` = c("COX1", "COX2", "COX3", "LOC119905811", "LOC119890004", "LOC119917832", "LOC119892218", "LOC119893803", "LOC119897473", "LOC119907810", "LOC119911148", "LOC119885605", "LOC119881863", "LOC119882732", "LOC119909213", "LOC119916265", "LOC119896019", "LOC119892403", "LOC119887585", "LOC119893283", "LOC119910933", "ATP6","ATP8", "LOC119891187", "atp5fa1", "atp5f1b", "atp5f1c", "atp5f1d", "atp5f1e", "atp5po", "atp5pb", "LOC119915245", "LOC119890133", "LOC119903517", "atp5mc3a", "atp5mc1", "atp5pd", "LOC119899395", "LOC119914361", "atp5mf", "LOC119881880", "atp5l", "LOC119901952", "atp5pf", "atp5mj", "atp5md", "ND1", "ND2", "ND3", "ND4", "ND5", "ND6", "LOC119900817", "LOC119889032", "CYTB"),
  
  # 8. Metabolism
  `Metabolism` = c("LOC119910845","LOC119892499", "aass", "LOC119918948", "LOC119903459", "LOC119898927", "fabp4a", "LOC119897230", "LOC119918650", "LOC119918651", "plin2", "cdab", "nudt4b", "LOC119893713","LOC119910978", "LOC119912177", "ndrg1a", "nupr1b", "egln3", "LOC119900507", "flt1","ace","si-ch211-202p1.5","vegfab"),
  
  # 9. Channel Proteins
  `Channel Proteins` = c("LOC119901983",  "rhbg", "LOC119904410", "atp1b1a",  "LOC119912124",  "kcnk5a", "slc12a7b", "LOC119908031", "slc9a3r1b", "chp1", "LOC119918691", "slc39a8", "slc6a14", "LOC119918306","slc1a5", "slc1a3a", "slc15a2", "LOC119885976", "slc6a8", "LOC119905556", "slc5a8l", "LOC119882223", "slc23a2", "slc28a1", "si-dkey-251i10.1", "slc10a3", "slc25a55a"),
  
  # 10. Fibroblasts & Macrophages
  `Fibroblasts_Macrophages` = c("col1a1a", "col1a1b", "col1a2", "LOC119904231", "col5a1","col5a2a", "col6a1","col6a2", "LOC119904002","LOC119903447","fn1b", "mmp2", "sparc", "LOC119882173", "dcn", "LOC119882174", "pcolcea", "ccn1", "LOC119909085", "fstl1b", "ntd5", "id3", "rbp4", "loxa", "bmp1a", "LOC119897617", "LOC119886711", "LOC119888537")
)

# ==============================================================================
# 2. Add Module Scores
# ==============================================================================
# This calculates scores for ALL signatures at once. 
# They will be named Score1, Score2, ..., Score10
TS_fliter <- AddModuleScore(
  object = TS_fliter,
  features = master_signature_list,
  name = "Score"  
)

# Extract generated column names to verify
score_cols <- colnames(TS_fliter@meta.data)[grepl("^Score\\d+$", colnames(TS_fliter@meta.data))]
print(score_cols)

# Extract Celltype Names
celltype_names <- names(master_signature_list)

# Define Custom Colors for Plots
color_wp <- colorRampPalette(c("white", "purple"))(100)

# Pre-subset object_sub for efficiency (used in Raincloud plot later)
cell_groups_to_compare <- c("Mac", "Outer E-Mac", "Inner E-Mac") # Note: Adjusted to 'Mac' based on previous script naming
cells_to_plot <- WhichCells(TS_fliter, idents = cell_groups_to_compare)
object_sub <- subset(TS_fliter, cells = cells_to_plot)

# ==============================================================================
# 3. Automation Loop for All Visualizations
# ==============================================================================
for (i in seq_along(celltype_names)) {
  celltype <- celltype_names[i]
  module_score_col <- paste0("Score", i)  
  
  save_dir <- file.path(save_base, celltype)
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  
  cat(sprintf("Processing Signature %d: %s\n", i, celltype))
  
  # ---------------------------------------------------------
  # Plot 1: Spatial Feature Plot (Transparent)
  # ---------------------------------------------------------
  p1 <- SpatialFeaturePlot(TS_fliter, features = module_score_col, image.alpha = 0.0, pt.size.factor = 4, ncol = 1) +
    ggtitle(paste0(celltype, " (Signature Score - Transparent)")) +
    scale_fill_gradientn(colors = color_wp, limits = c(0, NA), oob = scales::squish) 
  
  ggsave(filename = file.path(save_dir, paste0(celltype, "_Spatial_Transparent.pdf")), plot = p1, width = 12, height = 12)
  ggsave(filename = file.path(save_dir, paste0(celltype, "_Spatial_Transparent.jpg")), plot = p1, width = 12, height = 12)
  
  # ---------------------------------------------------------
  # Plot 2: Spatial Feature Plot (Default)
  # ---------------------------------------------------------
  p2 <- SpatialFeaturePlot(TS_fliter, features = module_score_col, pt.size.factor = 4, ncol = 1) +
    ggtitle(paste0(celltype, " (Signature Score - Default)")) +
    scale_fill_distiller(palette = "Blues", direction = 1, name = "Score") 
  
  ggsave(filename = file.path(save_dir, paste0(celltype, "_Spatial_Default.pdf")), plot = p2, width = 12, height = 12)
  ggsave(filename = file.path(save_dir, paste0(celltype, "_Spatial_Default.jpg")), plot = p2, width = 12, height = 12)
  
  # ---------------------------------------------------------
  # Plot 3: Standard VlnPlot
  # ---------------------------------------------------------
  p3 <- VlnPlot(TS_fliter, features = module_score_col, pt.size = 0.01, cols = cluster_colors) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  
  ggsave(filename = file.path(save_dir, paste0(celltype, "_VlnPlot_Standard.pdf")), plot = p3, width = 8, height = 8)
  
  # ---------------------------------------------------------
  # Data Prep for Advanced Violin/Raincloud Plots
  # ---------------------------------------------------------
  my_comparisons <- list(
    c("Outer E-Mac", "Mac"),
    c("Inner E-Mac", "Mac")
  )
  
  plot_data <- FetchData(TS_fliter, vars = c("ident", module_score_col))
  
  # Dynamic Y-axis calculation
  base_plot_for_range <- ggplot(plot_data, aes(x = ident, y = .data[[module_score_col]])) + geom_violin()
  y_range <- ggplot_build(base_plot_for_range)$layout$panel_params[[1]]$y.range
  y_max <- y_range[2]
  
  label_y1 <- y_max * 1.05
  label_y2 <- y_max * 1.15
  y_limit_total <- y_max * 1.25
  
  # ---------------------------------------------------------
  # Plot 4: Enhanced Violin Plot
  # ---------------------------------------------------------
  p_violin_enhanced <- ggplot(plot_data, aes(x = ident, y = .data[[module_score_col]], fill = ident)) +
    geom_violin(trim = FALSE, scale = "width", alpha = 0.8) +
    geom_jitter(shape = 16, position = position_jitter(0.15), size = 0.5, alpha = 0.4) +
    scale_fill_manual(values = cluster_colors) + 
    theme_classic(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
          axis.text.y = element_text(size = 12), axis.title = element_text(size = 14, face = "bold"),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5), legend.position = "none") +
    labs(title = paste(celltype, "Signature Score"), x = "Cell Identity", y = "Signature Score") +
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif", label.y = c(label_y1, label_y2), size = 5) +
    coord_cartesian(ylim = c(NA, y_limit_total))
  
  ggsave(filename = file.path(save_dir, paste0(celltype, "_VlnPlot_Enhanced.pdf")), plot = p_violin_enhanced, width = 5, height = 5)
  
  # ---------------------------------------------------------
  # Plot 5: General Raincloud Plot
  # ---------------------------------------------------------
  p_raincloud <- ggplot(plot_data, aes(x = ident, y = .data[[module_score_col]], fill = ident)) +
    geom_half_violin(side = "r", trim = FALSE, scale = "width", alpha = 0.7) +
    geom_half_boxplot(side = "r", width = 0.2, outlier.shape = NA) +
    geom_half_point(side = "l", width = 0.2, height = 0, size = 1.2, alpha = 0.3) +
    scale_fill_manual(values = cluster_colors) +
    theme_classic(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
          axis.text.y = element_text(size = 12), axis.title = element_text(size = 14, face = "bold"),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5), legend.position = "none") +
    labs(title = paste(celltype, "(Rain Cloud)"), x = "Cell Identity", y = "Signature Score") +
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif", label.y = c(label_y1, label_y2), size = 5) +
    coord_cartesian(ylim = c(NA, y_limit_total), clip = "off")
  
  ggsave(filename = file.path(save_dir, paste0(celltype, "_RainCloudPlot.pdf")), plot = p_raincloud, width = 5, height = 6)
  
  # ---------------------------------------------------------
  # Plot 6: Subset Raincloud Plot (Macrophages Only)
  # ---------------------------------------------------------
  plot_data_sub <- FetchData(object_sub, vars = c("ident", module_score_col))
  
  base_plot_sub <- ggplot(plot_data_sub, aes(x = ident, y = .data[[module_score_col]])) + geom_violin()
  y_range_sub <- ggplot_build(base_plot_sub)$layout$panel_params[[1]]$y.range
  y_max_sub <- y_range_sub[2]
  
  label_y1_sub <- y_max_sub * 1.00
  label_y2_sub <- y_max_sub * 1.12
  label_y3_sub <- y_max_sub * 1.3 
  y_limit_total_sub <- y_max_sub * 1.40 
  
  my_comparisons_sub <- list(
    c("Mac", "Outer E-Mac"),
    c("Outer E-Mac", "Inner E-Mac"),
    c("Mac", "Inner E-Mac")
  )
  
  p_raincloud_sub <- ggplot(plot_data_sub, aes(x = ident, y = .data[[module_score_col]], fill = ident)) +
    geom_half_violin(side = "r", trim = FALSE, scale = "width", alpha = 0.7) +
    geom_half_boxplot(side = "r", width = 0.2, outlier.shape = 16, outlier.size = 0, alpha = 0.3) +
    geom_half_point(side = "l", width = 0.2,  size = 0.6, alpha = 0.3) +
    scale_fill_manual(values = cluster_colors) +
    theme_classic(base_size = 14) +
    theme(axis.text.x = element_text(size = 12, face = "bold"),
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 14, face = "bold"),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          legend.position = "none") +
    labs(x = "Cell Identity", y = "Score", title = paste(celltype, "Subset")) +
    stat_compare_means(comparisons = my_comparisons_sub, method = "wilcox.test", label = "p.signif", 
                       label.y = c(label_y1_sub, label_y2_sub, label_y3_sub), size = 6) +
    coord_cartesian(ylim = c(NA, y_limit_total_sub), clip = "off")
  
  ggsave(filename = file.path(save_dir, paste0(celltype, "_RainCloudPlot_Subset.pdf")), plot = p_raincloud_sub, width = 5.25, height = 3.5)
  
  # ---------------------------------------------------------
  # Plot 7: UMAP Feature Plot
  # ---------------------------------------------------------
  p7 <- FeaturePlot(TS_fliter, features = module_score_col, reduction = "umap") +
    scale_color_distiller(palette = "Blues", direction = 1, limits = c(0, NA), oob = scales::squish) +
    theme_classic() +
    theme_minimal(base_size = 14) +
    theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
          legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, face = "bold"), panel.grid = element_blank()) 
  
  ggsave(filename = file.path(save_dir, paste0(celltype, "_FeaturePlot_UMAP.pdf")), plot = p7, width = 10, height = 10)
  ggsave(filename = file.path(save_dir, paste0(celltype, "_FeaturePlot_UMAP.jpg")), plot = p7, width = 10, height = 10)
}

cat("\n✅ All module scores calculated and plots generated successfully!\n")