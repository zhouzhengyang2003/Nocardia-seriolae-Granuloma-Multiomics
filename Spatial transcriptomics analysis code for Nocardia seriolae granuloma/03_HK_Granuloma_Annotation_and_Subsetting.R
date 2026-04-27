# ==============================================================================
# Head Kidney Granuloma - 3. Cluster Annotation & Visualization (Merged Clusters)
# ==============================================================================

# --- 1. Rename Clusters ---
new_cluster_names <- c(
  "0" = "Mac",
  "1" = "T-B-APCs",
  "2" = "Eryt",
  "3" = "Gran",
  "4" = "Mac",
  "5" = "Gran",
  "6" = "Outer E-Mac", # Corrected from "E-Mac-Outer" to match desired_order
  "7" = "Eryt",
  "8" = "Eryt-Gran",
  "9" = "Eryt",
  "10" = "Fibro-Endo",
  "11" = "B",
  "12" = "Outer E-Mac",
  "13" = "Outer E-Mac", 
  "14" = "Inner E-Mac", 
  "15" = "B"
)

TS_fliter <- RenameIdents(TS_fliter, new_cluster_names)

# --- 2. Set Desired Order for Clusters ---
# Note: Ensure these names exactly match the output of new_cluster_names
desired_order <- c(
  "Eryt",                   
  "Eryt-Gran", 
  "Gran",                   
  "Mac",       # Using "Mac" instead of "Non-G Mac" as defined in step 1, change if needed
  "B",   
  "T-B-APCs",  
  "Outer E-Mac", 
  "Inner E-Mac", 
  "Fibro-Endo"                
)

TS_fliter <- SetIdent(TS_fliter, value = factor(Idents(TS_fliter), levels = desired_order))

# Assign renamed identities to a new metadata column
TS_fliter$annotated_cluster <- Idents(TS_fliter)

# --- 3. Define Cell Type Colors ---
cluster_colors <- c(
  "Eryt" = "pink",            
  "Eryt-Gran" = "#FF9933",        
  "Gran" = "#FFFFCC",  
  "Mac" = "#1F78B4",       
  "B" = "red",      
  "T-B-APCs" = "#A6CEE3", 
  "Outer E-Mac" = "cyan", 
  "Inner E-Mac" = "#4169E1", 
  "Fibro-Endo" = "#000080" 
)

# Ensure cluster_colors are sorted according to desired_order
cluster_colors <- cluster_colors[desired_order]
print(levels(TS_fliter))         # Verify cell type groups
print(names(cluster_colors))     # Verify color names match exactly

# --- 4. UMAP Visualization ---
library(ggrepel)
library(dplyr)

# Extract UMAP coordinates and merge with cluster IDs
umap_df <- TS_fliter@reductions$umap@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(seurat_clusters = TS_fliter@active.ident) 

# Calculate the center position for each cluster label
seurat_clusterspos <- umap_df %>%
  group_by(seurat_clusters) %>%
  summarise(
    umap_1 = median(umap_1),
    umap_2 = median(umap_2)
  )

p_umap <- DimPlot(TS_fliter, reduction = "umap", label = FALSE, pt.size = 0.8) + 
  theme_dr(xlength = 0.2, ylength = 0.2, arrow = arrow(length = unit(0.2, "inches"), type = "closed")) +
  theme(panel.grid = element_blank(), axis.title = element_text(face = 2, hjust = 0.03)) +
  NoLegend() +
  geom_label_repel(
    aes(x = umap_1, y = umap_2, label = seurat_clusters), 
    data = seurat_clusterspos,
    fontface = "bold", box.padding = 0.5, size = 6
  ) +
  scale_color_manual(values = cluster_colors)

print(p_umap)

# --- 5. Spatial Cluster Visualization ---
p_spatial_clusters <- SpatialDimPlot(
  TS_fliter,
  label = FALSE,
  label.size = 6,
  pt.size.factor = 4.0,
  image.alpha = 0.0,  
  crop = FALSE
) +
  scale_fill_manual(values = cluster_colors) +
  guides(fill = guide_legend(override.aes = list(size = 8))) +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 12) 
  )

print(p_spatial_clusters)

# --- 6. Marker Gene DotPlot ---
library(RColorBrewer)

marker_list <- list(
  "Eryt" = c("zgc-163057", "LOC119891494", "slc4a1a", "cahz", "LOC119912177", "alas2", "tfr1a", "sptb", "klf1"),
  "Gran" = c("LOC119914447", "LOC119900902", "LOC119889900", "mmp9", "LOC119896105", "mmp25b", "LOC119897780", "si-ch211-284o19.8", "LOC119894963", "LOC119903646", "LOC119904387"),
  "Mac" = c("csf1rb", "marco", "mrc1a", "LOC119891768", "grna", "si-ch211-212k18.7", "LOC119883615", "LOC119915885", "LOC119910309"),
  "B" = c("cd79a", "LOC119884751", "LOC119890508", "LOC119890540", "LOC119890555", "LOC119890570"),
  "T" = c("LOC119910797", "LOC119883393", "LOC119883434", "LOC119914246", "LOC119897409", "LOC119897410"), 
  "E-Mac" = c("LOC119917376", "jupa", "LOC119889888", "dspa", "LOC119902058", "LOC119917243", "dsc2l", "LOC119900663", "cldn11a", "LOC119908668", "LOC119908847"),
  "Fibro-Endo" = c("col1a1b", "col1a2", "col1a1a", "col6a1", "pcolcea", "sparc")
)

# PRGn Palette (Colorblind friendly)
color_prgn <- colorRampPalette(rev(brewer.pal(n = 11, name = "PRGn")))(100)

p_dotplot <- DotPlot(object = TS_fliter, features = marker_list) +
  scale_color_gradientn(colors = color_prgn) +  
  labs(x = "Cell Clusters", y = "Marker Genes", color = "Avg Expression", size = "Percent Expressed") +
  theme_minimal(base_size = 10) +
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

# --- 7. Cell Proportions Pie Chart ---
# Calculate ratio
cell_table <- table(Idents(TS_fliter), TS_fliter$orig.ident)
Cellratio <- as.data.frame(prop.table(cell_table, margin = 2))
colnames(Cellratio) <- c("CellType", "Sample", "Freq")

# Calculate count
cell_counts <- as.data.frame(cell_table)
colnames(cell_counts) <- c("CellType", "Sample", "Count")

# Merge dataframes
Cellratio <- merge(Cellratio, cell_counts, by = c("CellType", "Sample"))

# Create legend labels
Cellratio$LegendLabel <- paste0(
  Cellratio$CellType, " (", sprintf("%.1f%%", Cellratio$Freq * 100), ", ", Cellratio$Count, ")"
)

# Ensure Factor Levels strictly follow desired_order
Cellratio$CellType <- factor(Cellratio$CellType, levels = desired_order)
Cellratio <- Cellratio[order(Cellratio$CellType), ]
Cellratio$LegendLabel <- factor(Cellratio$LegendLabel, levels = unique(Cellratio$LegendLabel))

p_pie <- ggplot(Cellratio, aes(x = "", y = Freq, fill = CellType)) + 
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  facet_wrap(~ Sample, nrow = 1) +
  scale_fill_manual(
    values = cluster_colors,          
    breaks = Cellratio$CellType,      
    labels = Cellratio$LegendLabel    
  ) +
  theme_void() +
  labs(fill = "Cell Type (Ratio, Count)", title = "Cell Composition per Sample") +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9),
    plot.title = element_text(hjust = 0.5)
  )

print(p_pie)

# --- 8. Subsetting Macrophages for Downstream Analysis ---
cell_groups_to_compare <- c(
  "Mac",               # Replaced "Non-G Mac" with "Mac" based on your idents above
  "Outer E-Mac",       # Outer Epithelioid Macrophages
  "Inner E-Mac"        # Inner Epithelioid Macrophages
)

# Extract subset containing only these specified macrophage populations
cells_to_plot <- WhichCells(TS_fliter, idents = cell_groups_to_compare)
object_sub <- subset(TS_fliter, cells = cells_to_plot)

# Note: You can save 'TS_fliter' and 'object_sub' using saveRDS() here if needed.