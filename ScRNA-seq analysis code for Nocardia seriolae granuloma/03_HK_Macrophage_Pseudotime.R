# ==============================================================================
# Head Kidney Granuloma - 03. Macrophage Pseudotime Analysis (Monocle3)
# ==============================================================================

suppressMessages({
  library(Seurat)
  library(SeuratWrappers)
  library(monocle3)
  library(dplyr)
  library(ggplot2)
})

# Define base output directory
base_out_dir <- "./results/Pseudotime_Analysis/"
dir.create(base_out_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 1. Subset Data & Convert Seurat to Monocle3 CellDataSet (CDS)
# ==============================================================================
cat("--- 1. Converting Seurat to Monocle3 CDS ---\n")

# Assuming 'Mac_scRNA_harmony' is loaded and annotated from the previous script
# Extract specific clusters representing the developmental/infection trajectory
target_clusters <- c("0", "1", "2", "4", "6", "9") # Update these if you changed to named idents earlier
Mac_infected <- subset(Mac_scRNA_harmony, idents = target_clusters)

# Convert to Monocle3 CDS
cds <- as.cell_data_set(Mac_infected)

# Sync Seurat UMAP coordinates to Monocle3
cds@int_colData$reducedDims$UMAP <- Mac_infected@reductions$umap@cell.embeddings

# ==============================================================================
# 2. Trajectory Inference & Pseudotime Ordering
# ==============================================================================
cat("--- 2. Learning Trajectory Graph ---\n")

# Cluster cells and learn the principal graph
cds <- cluster_cells(cds)
cds <- learn_graph(cds, use_partition = FALSE, learn_graph_control = list(ncenter = 200, minimal_branch_len = 8))

# Order cells (Set the root/starting point of the trajectory)
# Here, we set cells from Cluster "6" (e.g., HSPCs or Proliferating Macs) as the root
root_candidates <- intersect(colnames(Mac_infected)[Idents(Mac_infected) == "6"], colnames(cds))
cds <- order_cells(cds, root_cells = root_candidates[6]) # Picks a representative cell as the origin

# Visualize and Save the Trajectory
p_traj <- plot_cells(cds,
                     color_cells_by = "pseudotime",
                     label_groups_by_cluster = TRUE,
                     label_leaves = TRUE,
                     label_branch_points = TRUE,
                     graph_label_size = 3) +
  scale_color_gradientn(colors = c("lightblue", "white", "salmon", "red")) +
  theme(panel.grid = element_blank(), axis.title = element_text()) +
  NoLegend()

ggsave(file.path(base_out_dir, "1_Pseudotime_Trajectory_Map.pdf"), p_traj, width = 4, height = 4)

# ==============================================================================
# 3. Identify Pseudotime-Dependent Genes (graph_test)
# ==============================================================================
cat("--- 3. Identifying Trajectory-Correlated Genes ---\n")

deg_pseudotime <- graph_test(cds, neighbor_graph = "principal_graph", cores = 8)

# Filter for highly significant genes changing along the path
top_deg <- deg_pseudotime %>%
  arrange(q_value) %>%
  filter(q_value < 1e-4)   

write.csv(top_deg, file = file.path(base_out_dir, "2_Significant_Pseudotime_Genes.csv"))

cat(sprintf("Found %d significant pseudotime-dependent genes.\n", nrow(top_deg)))

# ==============================================================================
# 4. Define Helper Function for Single Gene Trend Plotting
# ==============================================================================
# Refactored to return a ggplot object for easier saving
plot_gene_trend_raw <- function(gene, cds, log_y = FALSE) {
  
  expr <- as.numeric(counts(cds[gene, ]))
  pt <- pseudotime(cds)
  
  # Remove NAs (cells not assigned a pseudotime)
  keep <- !is.na(pt)
  expr <- expr[keep]
  pt <- pt[keep]
  
  df <- data.frame(expression = expr, pseudotime = pt)
  
  p <- ggplot(df, aes(x = pseudotime, y = expression, color = pseudotime)) +
    geom_point(size = 2, alpha = 0.8) +
    geom_smooth(method = "loess", color = "black", size = 1.2, se = FALSE) +
    scale_color_gradientn(colors = c("lightblue", "white", "salmon", "red")) +
    labs(title = sub("gene-", "", gene), x = "Pseudotime", 
         y = ifelse(log_y, "log10(Expression + 1)", "Raw Expression (counts)")) +
    theme_classic(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", size = 0.6),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
    )
  
  if (log_y) {
    p <- p + scale_y_continuous(trans = scales::log10_trans())
  }
  
  return(p)
}

# ==============================================================================
# 5. Batch Process Curated Functional Signatures
# ==============================================================================
cat("--- 4. Batch Plotting Gene Expression Trends ---\n")

master_marker_list <- list(
  `Protease Inhibition` = c("LOC119910060","LOC119910049","cst3","LOC119882640","cast"),
  `Immune Regulation` = c("LOC119884194","LOC119884197","nfkbiaa","il10","tgfb1a"),
  `Stress Adaptation` = c("bcl2l1","si-dkey-243k1.3","babam2","LOC119910277","ywhag1","LOC119891659","rcan1b", "pkig","tsc22d1","tsc22d2","LOC119915309","LOC119893531","LOC119893593"),
  `Transcriptional regulation` = c("bhlhe40","hic1","LOC119899642","junbb","LOC119918304"),
  
  `Cytochrome c oxidase` = c("COX1", "COX2", "COX3", "LOC119905811", "LOC119890004", "LOC119917832", "LOC119892218", "LOC119911148", "LOC119885605", "LOC119882732", "LOC119916265", "LOC119896019", "LOC119892403", "LOC119893283"), 
  `F-type Atpase` = c("ATP6", "LOC119891187", "atp5fa1", "atp5f1b", "atp5f1c", "atp5f1d", "atp5f1e", "atp5po", "atp5pb", "LOC119915245", "LOC119890133", "LOC119903517", "atp5mc3a", "atp5mc1", "atp5pd", "LOC119899395", "LOC119914361", "atp5mf", "LOC119881880", "atp5l", "LOC119901952", "atp5pf", "atp5mj", "atp5md"),
  `NADH dehydrogenase` = c("ND1", "ND2", "ND3", "ND4", "ND5", "ND6"),
  `Cytochrome c reductase` = c("LOC119900817", "LOC119889032", "CYTB"),
  
  `Glycogen metabolism` = c("LOC119910845"),
  `Amino acid metabolism` = c("LOC119892499", "aass", "LOC119918948", "LOC119890054"),
  `Lipid metabolism` = c("LOC119903459", "LOC119898927", "fabp4a", "LOC119897230", "LOC119918650", "LOC119918651", "plin2"),
  `Nucleotide metabolism` = c("cdab"),
  `Inositol pyrophosphate metabolism` = c("nudt4b"),
  `Phosphocreatine transport` = c("LOC119893713","LOC119910978"),
  `Acid-Base and CO2 Homeostasis` = c("LOC119912177"),
  `Microenvironment Adaptation` = c("ndrg1a", "nupr1b", "egln3", "LOC119900507", "flt1","ace","si-ch211-202p1.5","vegfab"),
  
  `Water and inorganic ions` = c("LOC119901983", "rhbg", "LOC119904410", "atp1b1a", "LOC119912124", "kcnk5a", "slc12a7b", "LOC119908031", "slc9a3r1b", "chp1", "LOC119918691", "slc39a8"),
  `Amino Acids Peptides and Neurotransmitters` = c("slc6a14", "LOC119918306","slc1a5", "slc1a3a", "slc15a2", "LOC119885976", "slc6a8", "LOC119905556" ), 
  `Organic Acids` = c("slc5a8l", "LOC119882223"),
  `Vitamins` = c("slc23a2"),
  `Nucleosides` = c("slc28a1"),
  `ADP ATP` = c("si-dkey-251i10.1"),
  `Unknown` = c("slc10a3", "slc25a55a")
)

# Loop through categories and genes
for (category in names(master_marker_list)) {
  
  cat_dir <- file.path(base_out_dir, "Gene_Trends", category)
  dir.create(cat_dir, recursive = TRUE, showWarnings = FALSE)
  
  genes_to_plot <- master_marker_list[[category]]
  
  for (g in genes_to_plot) {
    # Ensure correct gene naming format required by the Monocle object
    current_gene <- paste0("gene-", g)
    
    if (current_gene %in% rownames(cds)) {
      message("  -> Plotting ", current_gene, " for category [", category, "]")
      
      # Generate the plot object
      p_trend <- plot_gene_trend_raw(gene = current_gene, cds = cds, log_y = TRUE)
      
      # Save PDF and PNG
      ggsave(filename = file.path(cat_dir, paste0("Trend_", g, ".pdf")), plot = p_trend, width = 8, height = 4)
      ggsave(filename = file.path(cat_dir, paste0("Trend_", g, ".jpg")), plot = p_trend, width = 8, height = 4, dpi = 300)
      
    } else {
      message("  [Warning] Gene '", current_gene, "' not found in CDS. Skipping...")
    }
  }
}

cat("\n✅ Pseudotime Analysis and Visualizations Completed!\n")