# ==============================================================================
# Head Kidney Granuloma - 4. Functional Enrichment Analysis (GO & KEGG)
# ==============================================================================

# ==============================================================================
# 0. Load Required Packages (Ensure 'bass_db' is loaded via AnnotationHub prior)
# ==============================================================================
suppressMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(clusterProfiler)
  library(GO.db)
})

# ==============================================================================
# 1. Core Parameters & Path Configuration
# ==============================================================================
# Path to your significant DEGs table
path_TS <- "./data/TS-MacvsE_Mac-significant.markers.csv"

# Base output directory
base_out_dir <- "./results/Enrichment_Analysis/"

# Read the master DEG table
df_TS <- read.csv(path_TS, row.names = 1) 

# Define parameters for the batch loop
metrics <- c("p_val_adj")       # Could also include "avg_log2FC" later
top_ns <- seq(1000, 1000, by = 100) # From 1000 to 1000, step 100

# ==============================================================================
# 2. Execute Automated Dual-loop Processing
# ==============================================================================

for (metric in metrics) {
  for (top_n in top_ns) {
    
    cat("\n=======================================================\n")
    cat(sprintf("🚀 Running batch task: Sorting Metric = %s | Top Genes = %d\n", metric, top_n))
    cat("=======================================================\n")
    
    # ---------------------------------------------------------
    # Step 1: Sort by metric and subset top N genes
    # ---------------------------------------------------------
    if (metric == "p_val_adj") {
      # For p-value, smaller is better (ascending)
      df_sorted <- df_TS[order(df_TS$p_val_adj, decreasing = FALSE), ]
    } else if (metric == "avg_log2FC") {
      # For log2FC, larger is better (descending)
      df_sorted <- df_TS[order(df_TS$avg_log2FC, decreasing = TRUE), ]
    }
    
    # Safely extract genes (prevents error if requested top_n > total rows)
    actual_top_n <- min(top_n, nrow(df_sorted))
    target_genes <- rownames(df_sorted)[1:actual_top_n]
    
    # ---------------------------------------------------------
    # Step 2: Create a dedicated directory for current parameters
    # ---------------------------------------------------------
    folder_name <- paste0("TS(", metric, ",", top_n, ")")
    current_out_dir <- file.path(base_out_dir, folder_name)
    
    if (!dir.exists(current_out_dir)) {
      dir.create(current_out_dir, recursive = TRUE)
    }
    
    # Save the current Gene List
    write.csv(data.frame(Gene_Name = target_genes), 
              file = file.path(current_out_dir, paste0("GeneList_", metric, "_", top_n, ".csv")), 
              row.names = FALSE)
    
    # ---------------------------------------------------------
    # Step 3: Gene ID Conversion (SYMBOL to ENTREZID)
    # ---------------------------------------------------------
    gene_convert <- bitr(target_genes, 
                         fromType = "SYMBOL", 
                         toType = "ENTREZID", 
                         OrgDb = bass_db)
    target_entrez_ids <- gene_convert$ENTREZID
    
    if (length(target_entrez_ids) == 0) {
      cat("⚠️ ID conversion failed. Skipping this iteration.\n")
      next
    }
    
    # ---------------------------------------------------------
    # Step 4: GO Analysis (BP, CC, MF separately + Network Plot)
    # ---------------------------------------------------------
    ontologies <- c("BP", "CC", "MF")
    ont_titles <- c("Biological Process", "Cellular Component", "Molecular Function")
    names(ont_titles) <- ontologies
    
    for (ont in ontologies) {
      # Use tryCatch to prevent loop termination if a specific ontology yields no results
      tryCatch({
        ego <- enrichGO(gene = target_entrez_ids, OrgDb = bass_db, ont = ont, 
                        pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
        
        if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
          # Save result table
          write.csv(as.data.frame(ego), file.path(current_out_dir, paste0("GO_Enrichment_UP_", ont, ".csv")))
          
          # Plot standard dotplot
          p_go <- dotplot(ego, showCategory = 15, title = paste0("GO Enrichment: ", ont_titles[ont]))
          print(p_go)
          
          # Plot network diagram
          p_net <- cnetplot(ego, showCategory = 10, circular = FALSE, colorEdge = TRUE) + 
            ggtitle(paste0("Gene-Pathway Network (UP-regulated ", ont, ")"))
          print(p_net)
        }
      }, error = function(e) { cat(sprintf("  - Skipped: %s category failed to plot or yielded no enrichment\n", ont)) })
    }
    
    # ---------------------------------------------------------
    # Step 5: Combined GO Analysis (ALL)
    # ---------------------------------------------------------
    tryCatch({
      ego_all <- enrichGO(gene = target_entrez_ids, OrgDb = bass_db, ont = "ALL", 
                          pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
      
      if (!is.null(ego_all) && nrow(as.data.frame(ego_all)) > 0) {
        p_all <- dotplot(ego_all, showCategory = 12, split = "ONTOLOGY", 
                         title = paste0("Conserved E-Mac Signature (", metric, " Top", top_n, ")")) + 
          facet_grid(ONTOLOGY ~ ., scale = "free", space = "free_y") +
          scale_fill_gradientn(colors = c("#4A1486", "#9E9AC8", "#67A9CF", "#016C59")) +
          theme_bw(base_size = 14) + 
          theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
            axis.text.y = element_text(size = 12, color = "black"),
            axis.text.x = element_text(size = 12, color = "black"),
            axis.title = element_text(face = "bold", size = 14),
            strip.text.y = element_text(face = "bold", size = 12, angle = 270),
            strip.background = element_rect(fill = "#e9ecef", color = "black", linewidth = 1),
            legend.title = element_text(face = "bold"),
            panel.grid.major.y = element_line(color = "grey80", linetype = "dashed"),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "black", linewidth = 1)
          )
        
        print(p_all)
      }
    }, error = function(e) { cat("  - Skipped: Combined GO plot failed\n") })
    
    # ---------------------------------------------------------
    # Step 6: KEGG Analysis
    # ---------------------------------------------------------
    tryCatch({
      cat("🔍 Running KEGG Enrichment (Organism: msam)...\n")
      
      # 1. Initial run without cutoff to verify server fetch
      ekegg_raw <- enrichKEGG(gene = target_entrez_ids, 
                              organism = 'msam', 
                              keyType = 'ncbi-geneid', 
                              pvalueCutoff = 1, 
                              qvalueCutoff = 1)
      
      if (is.null(ekegg_raw) || nrow(as.data.frame(ekegg_raw)) == 0) {
        cat("❌ KEGG returned empty results. Please verify: 1. Connection; 2. Organism code 'msam'; 3. Entrez IDs.\n")
      } else {
        # 2. Check for significant pathways
        df_kegg <- as.data.frame(ekegg_raw)
        sig_count <- sum(df_kegg$pvalue < 0.05)
        cat(sprintf("💡 Found %d raw pathways, %d of which have pvalue < 0.05\n", nrow(df_kegg), sig_count))
        
        if (sig_count > 0) {
          # Save raw KEGG table
          write.csv(df_kegg, file.path(current_out_dir, "KEGG_Enrichment_UP_ALL.csv"))
          
          # 3. Rerun with strict thresholds for plotting
          ekegg_sig <- enrichKEGG(gene = target_entrez_ids, 
                                  organism = 'msam', 
                                  keyType = 'ncbi-geneid', 
                                  pvalueCutoff = 0.05)
          
          # Aesthetic Dotplot
          p_kegg <- dotplot(ekegg_sig, showCategory = 20, 
                            title = paste0("KEGG Enrichment (", metric, " Top", top_n, ")")) +
            scale_fill_gradientn(colors = c("#D73027", "#FDAE61", "#87CEFA", "#313695" )) + 
            scale_y_discrete(labels = function(x) {
              if (requireNamespace("stringr", quietly = TRUE)) {
                stringr::str_wrap(x, width = 45)
              } else {
                sapply(x, function(s) paste(strwrap(s, width = 45), collapse = "\n"))
              }
            }) +
            theme_bw(base_size = 14) +
            theme(
              plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 15)),
              axis.text.y = element_text(size = 11, color = "black"),
              axis.text.x = element_text(size = 12, color = "black"),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_line(color = "grey90", linetype = "dashed"),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
            )
          
          print(p_kegg)
          cat("🎨 KEGG plot rendered successfully!\n")
          
        } else {
          cat("⚠️ KEGG pathways found, but none met the pvalue < 0.05 threshold. Skipping plot.\n")
        }
      }
    }, error = function(e) { 
      cat("❌ KEGG Error: ", e$message, "\n") 
    })
    
    cat("✅ Batch iteration complete! Data safely written to disk.\n")
  }
}

cat("\n🎉🎉🎉 All tasks completed successfully!\n")