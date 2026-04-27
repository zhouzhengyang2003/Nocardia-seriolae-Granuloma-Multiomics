# ==============================================================================
# Head Kidney Granuloma - 5. Marker Gene Dotplots
# ==============================================================================
library(Seurat)
library(ggplot2)
library(viridis)
library(RColorBrewer)

# Define Base Output Directory
dotplot_out_dir <- "./results/Dotplots/"
dir.create(dotplot_out_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Helper Function to Generate Dotplots
# This function eliminates repetitive code
# ==============================================================================
create_custom_dotplot <- function(seurat_obj, features, color_scale, limits = NULL) {
  
  # Basic DotPlot
  p <- DotPlot(object = seurat_obj, features = features) +
    scale_size_continuous(range = c(-1.1, 8), limits = c(0, 100), breaks = c(25, 50, 75)) +
    labs(x = NULL, y = NULL, color = NULL, size = NULL) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9, face = "bold"),
      axis.text.y = element_text(size = 9, face = "italic"),
      axis.title = element_text(size = 10, face = "bold"),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      panel.grid = element_blank(), # Remove all grid lines
      legend.position = "right"
    )
  
  # Apply color scale (viridis or custom gradient)
  if (is.character(color_scale) && length(color_scale) == 1) {
    # If it's a string, assume it's a viridis option name
    if (is.null(limits)) {
      p <- p + scale_color_viridis(option = color_scale)
    } else {
      p <- p + scale_color_viridis(option = color_scale, limits = limits, oob = scales::squish)
    }
  } else {
    # If it's an array of colors, use gradientn
    if (is.null(limits)) {
      p <- p + scale_color_gradientn(colors = color_scale)
    } else {
      p <- p + scale_color_gradientn(colors = color_scale, limits = limits, oob = scales::squish)
    }
  }
  
  return(p)
}


# ==============================================================================
# Plot 1: Cell Junctions
# ==============================================================================
cat("Generating Plot 1: Cell Junctions...\n")
marker_list_1 <- list(
  `AJ` = c("LOC119917376","cdh11", "jupa", "ctnna1", "ctnnd1"),
  `DSM` = c("dsc2l", "LOC119902058", "LOC119905441", "LOC119889888", "dspa", "LOC119917243"),
  `TJ` = c("cldn5a", "LOC119900663", "cldn11a", "tjp1a", "cgnl1"),
  `GJ` = c("LOC119908847", "LOC119908668")
)

# All Clusters
p1_all <- create_custom_dotplot(TS_fliter, marker_list_1, "magma")
ggsave(file.path(dotplot_out_dir, "1_Dotplot_All_CellJunctions.pdf"), p1_all, width=7.5, height=3.1)

# Macrophage Subsets Only
p1_sub <- create_custom_dotplot(object_sub, marker_list_1, "magma", limits=c(-1.2, 1.0))
ggsave(file.path(dotplot_out_dir, "1_Dotplot_Mac_CellJunctions.pdf"), p1_sub, width=7.5, height=2.2)


# ==============================================================================
# Plot 2: Cytoskeleton Regulation
# ==============================================================================
cat("Generating Plot 2: Cytoskeleton...\n")
marker_list_2 <- list(
  `Keratin Cytoskeleton` = c("LOC119905733", "LOC119886877", "krt15", "LOC119886876", "krt94", "krt99", "eppk1", "LOC119912524", "LOC119886962", "tgm1l1"),
  `Actin Cytoskeleton & Organization` = c("LOC119885789", "LOC119909726", "LOC119895369", "pof1b","LOC119901854", "crip1","crip2","nudcd2"),
  `Membrane-Cytoskeleton Coupling & Repair` = c("ahnak", "LOC119909740", "ezrb", "ehd2b", "stoml3b", "LOC119914141", "mpp6b", "anxa1a", "LOC119888099", "anxa2a", "LOC119890527", "LOC119898827", "anxa11a", "anxa13"),
  `Calcium Signaling` = c("s100a10a", "si-ch211-105c13.3", "s100v2", "calm1b","capns1b", "LOC119912040", "LOC119915637", "LOC119886455"),
  `Rho GTPase Signaling` = c("efnb2a","plekhg5b", "limk2", "cdc42ep3","ppp2cab", "ppp1r27b", "rnd3a", "plekhn1","LOC119890724")
)

p2_all <- create_custom_dotplot(TS_fliter, marker_list_2, "viridis")
ggsave(file.path(dotplot_out_dir, "2_Dotplot_All_Cytoskeleton.pdf"), p2_all, width=14, height=3.1)

p2_sub <- create_custom_dotplot(object_sub, marker_list_2, "viridis", limits=c(-1.2, 1.0))
ggsave(file.path(dotplot_out_dir, "2_Dotplot_Mac_Cytoskeleton.pdf"), p2_sub, width=14, height=2.2)


# ==============================================================================
# Plot 3: Extracellular Matrix (ECM)
# ==============================================================================
cat("Generating Plot 3: ECM...\n")
marker_list_3 <- list(
  `ECM Proteases` = c("mmp2", "mmp15b", "mmp19", "adam8a", "adam11", "LOC119886402"),
  `ECM Proteolysis Inhibitors` = c("serpine1", "serpine2", "LOC119908472", "LOC119906979", "timp2a", "timp2b", "LOC119918277"),
  `Cell–Matrix Adhesion & Transmembrane Anchorage` = c("ero1a","p4hb","cd44b", "tnfaip6", "fn1a", "fn1b", "LOC119897849", "cav1", "cavin1b", "LOC119907497", "LOC119918071", "cd81a", "LOC119893200", "LOC119895007", "LOC119899685", "vwa11", "anos1b", "aplp2", "prom2", "LOC119886359"),
  `Integrins` = c( "itga3b", "itga5", "LOC119902203", "itgav",   "LOC119901714", "LOC119892793",  "LOC119915730", "itgb6")
)

p3_all <- create_custom_dotplot(TS_fliter, marker_list_3, "mako")
ggsave(file.path(dotplot_out_dir, "3_Dotplot_All_ECM.pdf"), p3_all, width=12, height=3.1)

p3_sub <- create_custom_dotplot(object_sub, marker_list_3, "mako", limits=c(-1.2, 1.0))
ggsave(file.path(dotplot_out_dir, "3_Dotplot_Mac_ECM.pdf"), p3_sub, width=12, height=2.2)


# ==============================================================================
# Plot 4: Potential Upstream Regulators
# ==============================================================================
cat("Generating Plot 4: Upstream Regulators...\n")
marker_list_4 <- list(
  `Potential membrane microdomain organizers` = c( "spaca4l", "LOC119906807", "zanl"),
  `potential upstream regulators` = c( "wt1b", "gata6","LOC119894306", "LOC119914475"),
  `RA-associated signaling` = c( "aldh1a2", "gprc5c", "LOC119884110","crabp2a")
)

p4_all <- create_custom_dotplot(TS_fliter, marker_list_4, "rocket")
ggsave(file.path(dotplot_out_dir, "4_Dotplot_All_UpstreamRegulators.pdf"), p4_all, width=5.5, height=3.1)

p4_sub <- create_custom_dotplot(object_sub, marker_list_4, "rocket", limits=c(-1.2, 1.0))
ggsave(file.path(dotplot_out_dir, "4_Dotplot_Mac_UpstreamRegulators.pdf"), p4_sub, width=5.5, height=2.2)


# ==============================================================================
# Plot 5: Inflammation-Related
# ==============================================================================
cat("Generating Plot 5: Inflammation...\n")
marker_list_5 <- list(
  `Complements` = c("LOC119885585","c1qb","LOC119885788","LOC119897430","LOC119892056", "LOC119915106","LOC119896603"),
  `MHCII` = c("LOC119888605","LOC119917015","LOC119912422","LOC119912429","LOC119912010"),
  `Proteases` = c("ctsba","ctsc","LOC119890600","LOC119890271","ctsh","ctsl.1","LOC119892265","ctss1", "ctss2.1","LOC119882210","ctsz", "napsa"), 
  `Immune effectors` = c("tnfa", "m17","LOC119894803","LOC119890263","il1b")
)

p5_all <- create_custom_dotplot(TS_fliter, marker_list_5, "plasma")
ggsave(file.path(dotplot_out_dir, "5_Dotplot_All_Inflammation.pdf"), p5_all, width=12, height=3.1)

p5_sub <- create_custom_dotplot(object_sub, marker_list_5, "plasma", limits=c(-1.2, 1.0))
ggsave(file.path(dotplot_out_dir, "5_Dotplot_Mac_Inflammation.pdf"), p5_sub, width=12, height=2.2)


# ==============================================================================
# Plot 6: Immune Regulation & Homeostasis
# ==============================================================================
cat("Generating Plot 6: Immune Regulation & Homeostasis...\n")
marker_list_6 <- list(
  `Protease Inhibition` = c("LOC119910060","LOC119910049","cst3", "LOC119882640","cast"),
  `Immune Regulation` = c("LOC119884194","LOC119884197","nfkbiaa","il10","tgfb1a"),
  `Stress Adaptation, Survival & Genome Stability` = c("mcl1b","bcl2l1","si-dkey-243k1.3","babam2","LOC119910277","ywhag1", "LOC119891659","rcan1b", "pkig","tsc22d1","tsc22d2","LOC119915309","LOC119893531","LOC119893593"),
  `Transcriptional regulation` = c("bhlhe40","hic1","LOC119899642","junbb","LOC119918304")
)

p6_all <- create_custom_dotplot(TS_fliter, marker_list_6, "mako")
ggsave(file.path(dotplot_out_dir, "6_Dotplot_All_ImmuneRegulation.pdf"), p6_all, width=12, height=3.1)

p6_sub <- create_custom_dotplot(object_sub, marker_list_6, "mako", limits=c(-1.2, 1.0))
ggsave(file.path(dotplot_out_dir, "6_Dotplot_Mac_ImmuneRegulation.pdf"), p6_sub, width=12, height=2.2)


# ==============================================================================
# Plot 7: Oxidative Phosphorylation
# ==============================================================================
cat("Generating Plot 7: Oxidative Phosphorylation...\n")
marker_list_7 <- list(
  "Cytochrome c oxidase" = c("COX1", "COX2", "COX3",  "LOC119905811", "LOC119890004", "LOC119917832", "LOC119892218", "LOC119893803", "LOC119897473", "LOC119907810", "LOC119911148", "LOC119885605", "LOC119881863", "LOC119882732", "LOC119909213",  "LOC119916265", "LOC119896019", "LOC119892403", "LOC119887585", "LOC119893283", "LOC119910933"),
  "F-type Atpase" = c("ATP6","ATP8", "LOC119891187", "atp5fa1", "atp5f1b", "atp5f1c", "atp5f1d", "atp5f1e", "atp5po", "atp5pb", "LOC119915245", "LOC119890133", "LOC119903517", "atp5mc3a", "atp5mc1", "atp5pd", "LOC119899395", "LOC119914361", "atp5mf", "LOC119881880", "atp5l", "LOC119901952", "atp5pf",  "atp5mj", "atp5md"),
  "NADH dehydrogenase" = c("ND1", "ND2", "ND3", "ND4", "ND5", "ND6"),
  "Cytochrome c reductase" = c("LOC119900817", "LOC119889032", "CYTB")
)

color_YlGnBu <- colorRampPalette((brewer.pal(n = 9, name = "YlGnBu")))(100)

p7_all <- create_custom_dotplot(TS_fliter, marker_list_7, color_YlGnBu)
ggsave(file.path(dotplot_out_dir, "7_Dotplot_All_OxPhos.pdf"), p7_all, width=14, height=3.1)

p7_sub <- create_custom_dotplot(object_sub, marker_list_7, color_YlGnBu, limits=c(-1.2, 1.3))
ggsave(file.path(dotplot_out_dir, "7_Dotplot_Mac_OxPhos.pdf"), p7_sub, width=14, height=2.2)


# ==============================================================================
# Plot 8: Metabolism
# ==============================================================================
cat("Generating Plot 8: Metabolism...\n")
marker_list_8 <- list(
  "Glycogen metabolism" = c("LOC119910845"),
  "Amino acid metabolism" = c("LOC119892499", "aass", "LOC119918948", "LOC119890054"),
  "Lipid metabolism" = c("LOC119903459", "LOC119898927", "fabp4a", "LOC119897230", "LOC119918650", "LOC119918651", "plin2"),
  "Nucleotide metabolism" = c("cdab"),
  "Inositol pyrophosphate metabolism" = c("nudt4b"),
  "Phosphocreatine transport" = c("LOC119893713","LOC119910978"),
  "Acid–Base and CO₂ Homeostasis" = c("LOC119912177"),
  "Microenvironment Adaptation" = c("ndrg1a", "nupr1b", "egln3", "LOC119900507", "flt1","ace","si-ch211-202p1.5","vegfab")
)

color_PuOr <- colorRampPalette((brewer.pal(n = 9, name = "PuOr")))(100)

p8_all <- create_custom_dotplot(TS_fliter, marker_list_8, color_PuOr)
ggsave(file.path(dotplot_out_dir, "8_Dotplot_All_Metabolism.pdf"), p8_all, width=9, height=3.1)

p8_sub <- create_custom_dotplot(object_sub, marker_list_8, color_PuOr, limits=c(-1.2, 1.0))
ggsave(file.path(dotplot_out_dir, "8_Dotplot_Mac_Metabolism.pdf"), p8_sub, width=9, height=2.2)


# ==============================================================================
# Plot 9: Channel Proteins
# ==============================================================================
cat("Generating Plot 9: Channel Proteins...\n")
marker_list_9 <- list(
  "Water and inorganic ions" = c("LOC119901983",  "rhbg", "LOC119904410", "atp1b1a",  "LOC119912124",  "kcnk5a", "slc12a7b", "LOC119908031", "slc9a3r1b", "chp1", "LOC119918691", "slc39a8"),
  "Amino Acids, Peptides & Neurotransmitters" = c("slc6a14", "LOC119918306","slc1a5", "slc1a3a", "slc15a2", "LOC119885976", "slc6a8", "LOC119905556" ),
  "Organic Acids" = c("slc5a8l", "LOC119882223"),
  "Vitamins" = c("slc23a2"),
  "Nucleosides" = c("slc28a1"),
  "ADP/ATP" = c("si-dkey-251i10.1"),
  "Unknown" = c("slc10a3", "slc25a55a")
)

color_RdYlBu <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100)

p9_all <- create_custom_dotplot(TS_fliter, marker_list_9, color_RdYlBu)
ggsave(file.path(dotplot_out_dir, "9_Dotplot_All_ChannelProteins.pdf"), p9_all, width=10, height=3.1)

p9_sub <- create_custom_dotplot(object_sub, marker_list_9, color_RdYlBu, limits=c(-1.2, 1.0))
ggsave(file.path(dotplot_out_dir, "9_Dotplot_Mac_ChannelProteins.pdf"), p9_sub, width=10, height=2.2)


# ==============================================================================
# Plot 10: Fibroblasts & Macrophages
# ==============================================================================
cat("Generating Plot 10: Fibroblasts & Macrophages...\n")
marker_list_10 <- list(
  `Structural ECM Proteins` = c("col1a1a",  "col1a1b", "col1a2", "LOC119904231", "col5a1","col5a2a", "col6a1","col6a2", "LOC119904002","LOC119903447","fn1b"),
  `ECM Modulators & Remodeling Enzymes` = c("mmp2", "sparc", "LOC119882173", "dcn", "LOC119882174", "pcolcea"),
  `Cellular Signaling Factors` = c("ccn1", "LOC119909085", "fstl1b", "ntd5", "id3", "rbp4"),
  `E-mac-related Paracrine & Remodeling Factors` = c("LOC119888537","bmp1a","p4ha1b", "loxa",  "LOC119897617", "LOC119886711")
)

color_Blue <- colorRampPalette((brewer.pal(n = 11, name = "RdYlBu")))(100)

p10_all <- create_custom_dotplot(TS_fliter, marker_list_10, color_Blue, limits=c(-2.3, 2.3))
ggsave(file.path(dotplot_out_dir, "10_Dotplot_All_FibroMacs.pdf"), p10_all, width=10, height=3.1)

# Note: You didn't include the Macrophage subset code for Plot 10 in your original snippet,
# but if you need it, you can just call the function again:
# p10_sub <- create_custom_dotplot(object_sub, marker_list_10, color_Blue, limits=c(-2.3, 2.3))
# ggsave(file.path(dotplot_out_dir, "10_Dotplot_Mac_FibroMacs.pdf"), p10_sub, width=10, height=2.2)

cat("\n✅ All Dotplots generated successfully!\n")