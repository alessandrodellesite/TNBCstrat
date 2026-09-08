# RNAseq data loading and preprocessing

library(limma)
library(circlize)
library(ComplexHeatmap)
library(RColorBrewer)
library(readxl)

rna <- readRDS("/mnt/petasan_ccb/alessandro/SCANB/rna_logtransformed.rds")
gene_mads <- apply(rna, 1, mad)
ordered_mads <- order(gene_mads, decreasing = TRUE)
top_3000_indices <- ordered_mads[1:3000]
rna <- rna[top_3000_indices, ]
dt_matrix <- as.matrix(rna)

cluster_data <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/singleomic/nmf_rna/rna_nmf_clusters.csv")
cluster_data <- cluster_data[match(colnames(dt_matrix), cluster_data$SampleID), ]

groups <- factor(cluster_data$Cluster) 
design <- model.matrix(~0 + groups)
colnames(design) <- c("C1", "C2", "C3")
fit <- lmFit(dt_matrix, design) 
cont <- makeContrasts(
  C1_vs_others = C1 - (C2 + C3)/2,  
  C2_vs_others = C2 - (C1 + C3)/2,
  C3_vs_others = C3 - (C1 + C2)/2,
  levels = design
)
fit2 <- contrasts.fit(fit, cont)
fit2 <- eBayes(fit2, trend = TRUE) 

get_cluster_genes <- function(fit_obj, coef_name, n_genes = 500) {
  top_tab <- topTable(fit_obj, coef = coef_name, number = n_genes, p.value = 0.05, sort.by = "logFC")
  return(rownames(top_tab)) 
}
cluster1_genes <- get_cluster_genes(fit2, "C1_vs_others")
cluster2_genes <- get_cluster_genes(fit2, "C2_vs_others")
cluster3_genes <- get_cluster_genes(fit2, "C3_vs_others")

library(pheatmap)
c1_markers <- get_cluster_genes(fit2, "C1_vs_others", n_genes = 50)
c2_markers <- get_cluster_genes(fit2, "C2_vs_others", n_genes = 50)
c3_markers <- get_cluster_genes(fit2, "C3_vs_others", n_genes = 50)
ordered_genes <- c(c1_markers, c2_markers, c3_markers)

# --- NEW: relabel old cluster IDs so sorting gives visual order 3,1,2 ---
# old "3" -> new "1", old "1" -> new "2", old "2" -> new "3"
relabel_map <- c("3" = "1", "1" = "2", "2" = "3")
cluster_data$Cluster_relabelled <- factor(relabel_map[as.character(cluster_data$Cluster)],
                                           levels = c("1", "2", "3"))
sample_order <- order(cluster_data$Cluster_relabelled)

plot_matrix <- dt_matrix[ordered_genes, sample_order]
plot_matrix <- t(scale(t(plot_matrix)))

annotation_col <- data.frame(RNAseq = cluster_data$Cluster_relabelled[sample_order])
rownames(annotation_col) <- colnames(plot_matrix)

methyl_probes      <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/singleomic/nmf_meth/methyl_nmf_clusters.csv")
methyl_probes_4    <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/singleomic/nmf_meth/methyl_nmf_clusters_4.csv")
cnv_data_clusters  <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/singleomic/nmf_cnv/cnv_nmf_clusters.csv")

dt_metadata <- read_excel("/mnt/petasan_ccb/juanra/SCANB/RNAseq/metadata/ids_cohorts_match.xlsx", sheet = "1a SCAN-B discovery")

colonna_id_excel <- "PD_ID" 
metadata_matched <- dt_metadata[match(cluster_data$SampleID[sample_order], dt_metadata[[colonna_id_excel]]), ]

cols_to_keep <- c("TMB", "TILs", "PAM50_Basal_NCN", "PAM50_NCN", 
                   "TNBCtype4_n235_notPreCentered", "TNBCtype6_n235_notPreCentered", 
                   "CibersortX.Tcell", "CibersortX.endothelial", "CibersortX.Bcell", 
                   "CibersortX.stroma", "CibersortX.macrophage", "CibersortX.epithelial", 
                   "ASCAT_PLOIDY", "ASCAT_TUM_FRAC")
meta_sub <- metadata_matched[, cols_to_keep]

numeric_cols <- c("TMB", "TILs", "CibersortX.Tcell", "CibersortX.endothelial", 
                   "CibersortX.Bcell", "CibersortX.stroma", "CibersortX.macrophage", 
                   "CibersortX.epithelial", "ASCAT_PLOIDY", "ASCAT_TUM_FRAC")
meta_sub[numeric_cols] <- lapply(meta_sub[numeric_cols], function(x) as.numeric(as.character(x)))

cat_cols <- setdiff(cols_to_keep, numeric_cols)
meta_sub[cat_cols] <- lapply(meta_sub[cat_cols], factor)

names(meta_sub)[names(meta_sub) == "TNBCtype4_n235_notPreCentered"] <- "TNBCtype4"
names(meta_sub)[names(meta_sub) == "TNBCtype6_n235_notPreCentered"] <- "TNBCtype6"

annotation_col <- data.frame(
  RNAseq = cluster_data$Cluster_relabelled[sample_order],
  meta_sub
)
rownames(annotation_col) <- colnames(plot_matrix)

desired_order <- c("RNAseq", "TNBCtype4", "TNBCtype6", "PAM50_NCN", "PAM50_Basal_NCN",
                    "TMB", "TILs", 
                    "CibersortX.epithelial", "CibersortX.macrophage", "CibersortX.stroma",
                    "CibersortX.Bcell", "CibersortX.endothelial", "CibersortX.Tcell",
                    "ASCAT_PLOIDY", "ASCAT_TUM_FRAC")
annotation_col <- annotation_col[, desired_order]

# --- NEW: Cluster colors follow the ORIGINAL cluster identity, keyed by the NEW label ---
# old "1" (tomato3) is now new "2"; old "2" (#0984E3) is now new "3"; old "3" (#00B894) is now new "1"
ann_colors = list(
  RNAseq = c("1" = "#00B894", "2" = "tomato3", "3" = "#0984E3"),
  TNBCtype4 = c("BL1" = "#A29BFE", "BL2" = "#74B9FF", 
                "M" = "#55E6C1", "LAR" = "#FDCB6E", "NA" = "#B2BEC3"),
  TNBCtype6 = c("BL1" = "#A29BFE", "BL2" = "#74B9FF", "M" = "#55E6C1", 
                "LAR" = "#FDCB6E", "NA" = "#B2BEC3", "IM" = "#006266", 
                "MSL" = "#FFEAA7", "UNS" = "#FFADAD"),
  PAM50_NCN = c(Basal = "#D63031", Her2 = "#A29BFE", LumB = "#0984E3", 
                LumA = "#74B9FF", Normal = "#00B894", unclassified = "#B2BEC3"),
  PAM50_Basal_NCN = c(Basal = "#D63031", nonBasal = "#E0E0E0"),
  
  TMB   = colorRamp2(c(0, max(annotation_col$TMB, na.rm = TRUE)), c("#F5F6FA", "#079992")),
  TILs  = colorRamp2(c(0, max(annotation_col$TILs, na.rm = TRUE)), c("#F5F6FA", "#6C5CE7")),
  
  CibersortX.epithelial  = colorRamp2(c(0, max(annotation_col$CibersortX.epithelial, na.rm = TRUE)), c("#F5F6FA", "#0984E3")),
  CibersortX.macrophage  = colorRamp2(c(0, max(annotation_col$CibersortX.macrophage, na.rm = TRUE)), c("#F5F6FA", "#E17055")),
  CibersortX.stroma      = colorRamp2(c(0, max(annotation_col$CibersortX.stroma, na.rm = TRUE)), c("#F5F6FA", "#00B894")),
  CibersortX.Bcell       = colorRamp2(c(0, max(annotation_col$CibersortX.Bcell, na.rm = TRUE)), c("#F5F6FA", "#6C5CE7")),
  CibersortX.endothelial = colorRamp2(c(0, max(annotation_col$CibersortX.endothelial, na.rm = TRUE)), c("#F5F6FA", "#FDCB6E")),
  CibersortX.Tcell       = colorRamp2(c(0, max(annotation_col$CibersortX.Tcell, na.rm = TRUE)), c("#F5F6FA", "#D63031")),
  
  ASCAT_PLOIDY   = colorRamp2(c(0, max(annotation_col$ASCAT_PLOIDY, na.rm = TRUE)), c("#F5F6FA", "#273C75")),
  ASCAT_TUM_FRAC = colorRamp2(c(0, max(annotation_col$ASCAT_TUM_FRAC, na.rm = TRUE)), c("#F5F6FA", "#079992"))
)

col_ann <- HeatmapAnnotation(
  df = annotation_col, 
  col = ann_colors,
  show_legend = TRUE,
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 8)
)

ht <- Heatmap(
  plot_matrix, 
  name = "Expression",          
  top_annotation = col_ann,     
  show_row_names = FALSE,       
  show_column_names = FALSE,    
  cluster_columns = FALSE,      
  cluster_rows = TRUE,          
  col = colorRampPalette(c("blue", "white", "red"))(100),
  width = unit(14, "cm")
)

ht_opt(
  legend_title_gp = gpar(fontsize = 8, fontface = "bold"),
  legend_labels_gp = gpar(fontsize = 7),
  legend_grid_height = unit(3, "mm"),
  legend_grid_width = unit(3, "mm")
)

# Landscape A4: 11.69 x 8.27 in
png("/mnt/petasan_ccb/alessandro/SCANB/plots/comparisons/heatmap_rna_metadata.png", 
    width = 11.69, height = 8.27, units = "in", res = 300)
draw(ht, 
     merge_legends = TRUE, 
     heatmap_legend_side = "right", 
     annotation_legend_side = "right")
dev.off()

ht_opt(RESET = TRUE)


#mofa
mofa_clusters_2 <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_mofa/mofa_km_clusters_2.csv")
mofa_clusters_3 <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_mofa/mofa_km_clusters_3.csv")
mofa_clusters_4 <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_mofa/mofa_km_clusters_4.csv")

#icluster
icluster_2 <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_icluster/iclusters_clusters_2.csv")
icluster_3 <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_icluster/iclusters_clusters_3.csv")                                
icluster_4 <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_icluster/iclusters_clusters_4.csv")
icluster_5 <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_icluster/iclusters_clusters_5.csv")
                                 
#snf
SNF_clusters_2 <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_snf/snf_clusters_2.csv")
SNF_clusters_3 <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_snf/snf_clusters_3.csv")                                
SNF_clusters_4 <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_snf/snf_clusters_4.csv")
                                 
#xintnmf                                 
XintNMF_clusters_2_reg <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_xintnmf/xintNMF_clusters_k2_reg.csv")                    
XintNMF_clusters_3_reg <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_xintnmf/xintNMF_clusters_k3_reg.csv")
XintNMF_clusters_4_reg <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_xintnmf/xintNMF_clusters_k4_reg.csv")                                

# Vectors to map cluster to samples
methyl_map_probes <- setNames(methyl_probes$Cluster, methyl_probes$SampleID)
methyl_map_probes_4 <- setNames(methyl_probes_4$Cluster, methyl_probes_4$SampleID)
cnv_map <- setNames(cnv_data_clusters$Cluster, cnv_data_clusters$SampleID)
mofa_map_2 <- setNames(mofa_clusters_2$Cluster, mofa_clusters_2$SampleID)
mofa_map_3 <- setNames(mofa_clusters_3$Cluster, mofa_clusters_3$SampleID)
mofa_map_4 <- setNames(mofa_clusters_4$Cluster, mofa_clusters_4$SampleID)
SNF_map_2 <- setNames(SNF_clusters_2$Cluster, SNF_clusters_2$SampleID)
SNF_map_3 <- setNames(SNF_clusters_3$Cluster, SNF_clusters_3$SampleID)
SNF_map_4 <- setNames(SNF_clusters_4$Cluster, SNF_clusters_4$SampleID)
#xint_map_2 <- setNames(XintNMF_clusters_2$Cluster, XintNMF_clusters_2$SampleID)
#xint_map_3 <- setNames(XintNMF_clusters_3$Cluster, XintNMF_clusters_3$SampleID)
#xint_map_4 <- setNames(XintNMF_clusters_4$Cluster, XintNMF_clusters_4$SampleID)
xint_map_2_reg <- setNames(XintNMF_clusters_2_reg$Cluster, XintNMF_clusters_2_reg$SampleID)
xint_map_3_reg <- setNames(XintNMF_clusters_3_reg$Cluster, XintNMF_clusters_3_reg$SampleID)
xint_map_4_reg <- setNames(XintNMF_clusters_4_reg$Cluster, XintNMF_clusters_4_reg$SampleID)
icluster_map_2 <- setNames(icluster_2$Cluster, icluster_2$SampleID)
icluster_map_3 <- setNames(icluster_3$Cluster, icluster_3$SampleID)
icluster_map_4 <- setNames(icluster_4$Cluster, icluster_4$SampleID)
icluster_map_5 <- setNames(icluster_5$Cluster, icluster_5$SampleID)
                                 
# Allineiamo tutti i cluster alle colonne di plot_matrix
annotation_col_meth <- data.frame(
  RNAseq                 = factor(cluster_data$Cluster[sample_order]),
  #RNAseq_1000            = factor(rna_1000_map[colnames(plot_matrix)]),
  DNAm          = factor(methyl_map_probes[colnames(plot_matrix)]),
  DNAm_4          = factor(methyl_map_probes_4[colnames(plot_matrix)]),
  CNV                    = factor(cnv_map[colnames(plot_matrix)]),
  MOFA_2                 = factor(mofa_map_2[colnames(plot_matrix)]),
  MOFA_3                 = factor(mofa_map_3[colnames(plot_matrix)]),
  MOFA_4                 = factor(mofa_map_4[colnames(plot_matrix)]),
  SNF_2                  = factor(SNF_map_2[colnames(plot_matrix)]),
  SNF_3                  = factor(SNF_map_3[colnames(plot_matrix)]),
  SNF_4                  = factor(SNF_map_4[colnames(plot_matrix)]),
  icluster_2             = factor(icluster_map_2[colnames(plot_matrix)]),
  icluster_3             = factor(icluster_map_3[colnames(plot_matrix)]),
  icluster_4             = factor(icluster_map_4[colnames(plot_matrix)]),
#  icluster_5             = factor(icluster_map_5[colnames(plot_matrix)]),
  XintNMF_2                  = factor(xint_map_2_reg[colnames(plot_matrix)]),
  XintNMF_3                  = factor(xint_map_3_reg[colnames(plot_matrix)]),
  XintNMF_4                  = factor(xint_map_4_reg[colnames(plot_matrix)])
  )


# Impostiamo i nomi delle righe per farli coincidere con le colonne della heatmap
rownames(annotation_col_meth) <- colnames(plot_matrix)
                                 
ann_colors_meth = list(
  RNAseq          = c("1" = "tomato3", "2" = "#0984E3", "3" = "#00B894"),
  DNAm            = c("1" = "#ffb8b8", "2" = "#95afc0", "3" = "#badc58"), 
  DNAm_4          = c("1" = "#95afc0", "2" = "#ffb8b8", "3" = "#badc58", "4"= "#34495e"), 
  CNV             = c("1" = "#ff7675", "2" = "#74b9ff", "3" = "#fdcb6e"), 
  MOFA_2          = c("1" = "#ffeaa7", "2" = "#0984E3"),
  MOFA_3          = c("1" = "#0984E3", "2" = "#9b59b6", "3" = "#ffeaa7"),
  MOFA_4          = c("1" = "#ffeaa7", "2" = "#9b59b6", "3" = "#ff7f98", "4" = "#0984E3"),
  icluster_2      = c("1" = "#27ae33", "2" = "#e1b99c"),
  icluster_3      = c("1" = "#e1b99c", "2" = "#4c5ce1", "3" = "#27ae33"),
  icluster_4      = c("1" = "#ff6f11", "2" = "#e1b99c", "3" = "#27ae33", "4" = "#4c5ce1"),
  #icluster_5      = c("1" = "#ff6f11", "2" = "#e1b99c", "3" = "#27ae33", "4" = "#4c5ce1", "5" = "#5B8E8A"),
  SNF_2           = c("1" = "#ff7f50", "2" = "#008080"),
  SNF_3           = c("1" = "#34495e", "2" = "#ff7f50", "3" = "#008080"),
  SNF_4           = c("1" = "#34495e", "2" = "#008080", "3" = "#f1c40f", "4" = "#ff7f50"),
  XintNMF_2           = c("1" = "#27ae60", "2" = "#e1b12c"),
  XintNMF_3           = c("1" = "#27ae60", "2" = "#e1b12c", "3" = "#6c5ce7"),
  XintNMF_4           = c("1" = "#6c5ce7", "2" = "#ff6f99", "3" = "#27ae60", "4" = "#e1b12c")
)              

# HEATMAP 

col_ann_meth <- HeatmapAnnotation(
  df = annotation_col_meth, 
  col = ann_colors_meth,
  show_legend = TRUE,
  annotation_name_side = "left"
)

ht_meth <- Heatmap(
  plot_matrix, 
  name = "Expression",          
  top_annotation = col_ann_meth,     
  show_row_names = FALSE,       
  show_column_names = FALSE,    
  cluster_columns = FALSE,      # Mantiene l'ordine dei campioni NMF
  cluster_rows = TRUE,         
  col = colorRampPalette(c("blue", "white", "red"))(100),
  width = unit(12, "cm") 
)

png("/mnt/petasan_ccb/alessandro/SCANB/plots/comparisons/heatmap_omics.png", width = 10, height = 8, units = "in", res = 300)
draw(ht_meth, 
     heatmap_legend_side = "right", 
     annotation_legend_side = "right"
)
dev.off()
 
