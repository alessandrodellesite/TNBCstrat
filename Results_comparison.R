# RNAseq data loading and preprocessing


dt <- read.csv("matched_rna.tsv", sep="\t", row.names = 1)

tpm_matrix <- apply(dt, 2, function(x) x / sum(as.numeric(x))*10^6)
data_filtered <- tpm_matrix[rowSums(dt) > 0, ]

threshold_tpm <- 2
min_samples <- 0.05*dim(data_filtered)[2]
keep <- rowSums(data_filtered > threshold_tpm) >= min_samples
data_fil <- data_filtered[keep, ]

data_prepped <- log2(data_fil + 1)
```

## Feature selection (top 3000 genes by MAD)

```{r}
gene_mads <- apply(data_prepped, 1, mad)
ordered_mads <- order(gene_mads, decreasing = TRUE)

top_3000_indices <- ordered_mads[1:3000]
data_prepped <- data_prepped[top_3000_indices, ]

dt_matrix <- as.matrix(data_prepped)
if(any(is.na(dt_matrix))) {
  dt_matrix[is.na(dt_matrix)] <- 0 # Sostituisci i NA con 0
  }
```

# DEA of of genes' expression between NMF clusters from RNAseq 

```{r}
cluster_data <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/singleomic/nmf_rna/rna_nmf_clusters.csv")
cluster_data <- cluster_data[match(colnames(dt_matrix), cluster_data$SampleID), ]
# Create the factor and design matrix
groups <- factor(cluster_data$Cluster) 
design <- model.matrix(~0 + groups)
colnames(design) <- c("C1", "C2", "C3")
```

```{r}
library(limma)
fit <- lmFit(dt_matrix, design) 
cont <- makeContrasts(
  C1_vs_others = C1 - (C2 + C3)/2,  
  C2_vs_others = C2 - (C1 + C3)/2,
  C3_vs_others = C3 - (C1 + C2)/2,
  levels = design
)

fit2 <- contrasts.fit(fit, cont)
fit2 <- eBayes(fit2, trend = TRUE) 

results_1 <- topTable(fit2, coef="C1_vs_others", number=500, p.value=0.05, sort.by="logFC")
results_2 <- topTable(fit2, coef="C2_vs_others", number=500, p.value=0.05, sort.by="logFC")
results_3 <- topTable(fit2, coef="C3_vs_others", number=500, p.value=0.05, sort.by="logFC")

get_cluster_genes <- function(fit_obj, coef_name, n_genes = 500) {
  top_tab <- topTable(fit_obj, coef = coef_name, number = n_genes, p.value = 0.05, sort.by = "logFC") #sorts by absolute value of log2fc
  return(rownames(top_tab)) 
}

cluster1_genes <- get_cluster_genes(fit2, "C1_vs_others")
cluster2_genes <- get_cluster_genes(fit2, "C2_vs_others")
cluster3_genes <- get_cluster_genes(fit2, "C3_vs_others")
```

## Heatmap with ordered samples

```{r}
library(pheatmap)
c1_markers <- get_cluster_genes(fit2, "C1_vs_others", n_genes = 50)
c2_markers <- get_cluster_genes(fit2, "C2_vs_others", n_genes = 50)
c3_markers <- get_cluster_genes(fit2, "C3_vs_others", n_genes = 50)

# ORDER groups
ordered_genes <- c(c1_markers, c2_markers, c3_markers)

# Sort SAMPLES by their cluster group
sample_order <- order(cluster_data$Cluster)
plot_matrix <- dt_matrix[ordered_genes, sample_order]

# z-score to Center the data  
plot_matrix <- t(scale(t(plot_matrix)))

# Subtract the median of each row: plot_matrix <- plot_matrix - apply(plot_matrix, 1, median)

#Create the sorted annotation
annotation_col <- data.frame(Cluster = factor(cluster_data$Cluster[sample_order]))
rownames(annotation_col) <- colnames(plot_matrix)


pheatmap(plot_matrix, 
         annotation_col = annotation_col, 
         cluster_rows = TRUE,  # Keep genes grouped by cluster
         cluster_cols = FALSE,  # Keep samples grouped by cluster
         show_colnames = FALSE, 
         show_rownames = FALSE, 
         main = "Ordered NMF Cluster Markers",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         #breaks = seq(-2, 2, length.out = 101)
         ) 


```
# HO RIPETUTO IL CODICE DUE VOLTE: LA SECONDA è QUELLA ORIGINALE, LA PRIMA è QUELLA AGGIUSTATA X ASSICURARSI CHE PRENDA BENE LE COLONNE, MA I RISULTATI SONO UGUALI

```{r}

library(circlize)
library(ComplexHeatmap)
library(RColorBrewer)
library(readxl)

# Caricamento dei file di cluster alternativi (Methyl, CNV, ecc.)
methyl_clusters_df <- read.csv("methyl_nmf_clusters_genes.csv") 
methyl_en_df       <- read.csv("methyl_nmf_clusters_genes_en.csv") 
methyl_pr_df       <- read.csv("methyl_nmf_clusters_genes_pr.csv")
methyl_probes      <- read.csv("methyl_nmf_clusters_probes.csv") 
cnv_data_clusters  <- read.csv("cnv_nmf_clusters.csv")
cnv_rnaseq_clusters <- read.csv("cnv_rnaseq_nmf_clusters.csv")

# Caricamento del file Excel dei metadati clinici
dt_metadata <- read_excel("ids_cohorts_match.xlsx", sheet = "1a SCAN-B discovery") 
dt_metadata$ASCAT_PLOIDY
dt_metadata$ASCAT_TUM_FRAC

# METADATA ALIGNMENT & PREPROCESSING

# Specifica il nome esatto della colonna degli ID dei campioni dentro il file Excel
colonna_id_excel <- "PD_ID" 

metadata_matched <- dt_metadata[match(cluster_data$SampleID[sample_order], dt_metadata[[colonna_id_excel]]), ]

# Selezioniamo solo le colonne di interesse dal dataset correttamente allineato

#cols_to_keep <- c("TMB", "TumSize", "Grade", "LNbinary", "TILs", "PAM50_Basal_NCN", "PAM50_NCN", "TNBCtype4_n235_notPreCentered", "TNBCtype6_n235_notPreCentered", "CibersortX.Tcell", "CibersortX.endothelial", "CibersortX.Bcell", "CibersortX.stroma", "CibersortX.macrophage", "CibersortX.epithelial")

cols_to_keep <- c("TMB", "TILs", "PAM50_Basal_NCN", "PAM50_NCN", "TNBCtype4_n235_notPreCentered", "TNBCtype6_n235_notPreCentered", "CibersortX.Tcell", "CibersortX.endothelial", "CibersortX.Bcell", "CibersortX.stroma", "CibersortX.macrophage", "CibersortX.epithelial", "ASCAT_PLOIDY", "ASCAT_TUM_FRAC")

meta_sub <- metadata_matched[, cols_to_keep]

# Colonne numeriche
numeric_cols <- c("TMB", "TILs", "CibersortX.Tcell", "CibersortX.endothelial", "CibersortX.Bcell", "CibersortX.stroma", "CibersortX.macrophage", "CibersortX.epithelial","ASCAT_PLOIDY",  "ASCAT_TUM_FRAC")
meta_sub[numeric_cols] <- lapply(meta_sub[numeric_cols], function(x) as.numeric(as.character(x)))

# Colonne categoriali (Factor)
cat_cols <- setdiff(cols_to_keep, numeric_cols)
meta_sub[cat_cols] <- lapply(meta_sub[cat_cols], factor)

# Creazione del dataframe finale per le annotazioni dell'heatmap
annotation_col <- data.frame(
  Cluster = factor(cluster_data$Cluster[sample_order]),
  meta_sub  # meta_sub è già allineato e ordinato correttamente
)

# Impostiamo i row names di annotation_col affinché corrispondano al 100% alle colonne di plot_matrix
rownames(annotation_col) <- colnames(plot_matrix)


# Definizione dei colori per le annotazioni cliniche e i cluster
ann_colors = list(
  Cluster = c("1" = "tomato3", "2" = "#0984E3", "3" = "#00B894"),
  #Grade = c("2" = "#FAB1A0", "3" = "#E17055", "NA" = "#B2BEC3"),
  #LNbinary = c("positive" = "#636E72", "negative" = "#DFE6E9"),
  TNBCtype4_n235_notPreCentered = c("BL1" = "#A29BFE", "BL2" = "#74B9FF", 
                                    "M" = "#55E6C1", "LAR" = "#FDCB6E", "NA"= "#B2BEC3"),
  TNBCtype6_n235_notPreCentered = c("BL1" = "#A29BFE", "BL2" = "#74B9FF", "M" = "#55E6C1", 
                                    "LAR" = "#FDCB6E", "NA"= "#B2BEC3", "IM"= "#006266", 
                                    "MSL"= "#FFEAA7", "UNS"= "#FFADAD"),
  PAM50_NCN = c(Basal = "#D63031", Her2 = "#A29BFE", LumB = "#0984E3", 
                LumA = "#74B9FF", Normal = "#00B894", unclassified = "#B2BEC3"),
  PAM50_Basal_NCN = c(Basal = "#D63031", nonBasal = "#E0E0E0"),
  #TumSize = colorRamp2(c(0, max(annotation_col$TumSize, na.rm = TRUE)), c("#F5F6FA", "#273C75")),
  TMB     = colorRamp2(c(0, max_TMB     <- max(annotation_col$TMB, na.rm = TRUE)),     c("#F5F6FA", "#079992")),
  TILs    = colorRamp2(c(0, max_TILs    <- max(annotation_col$TILs, na.rm = TRUE)),    c("#F5F6FA", "#6C5CE7")),
  Cbx_epithelial    = colorRamp2(c(0, max_TILs    <- max(annotation_col$CibersortX.epithelial, na.rm = TRUE)),    c("#F5F6FA", "#6C5CE7")),
  Cbx_macrophage    = colorRamp2(c(0, max_TILs    <- max(annotation_col$CibersortX.macrophage, na.rm = TRUE)),    c("#F5F6FA", "#6C5CE7")),
  Cbx_stroma    = colorRamp2(c(0, max_TILs    <- max(annotation_col$CibersortX.stroma, na.rm = TRUE)),    c("#F5F6FA", "#6C5CE7")),
  Cbx_Bcell    = colorRamp2(c(0, max_TILs    <- max(annotation_col$CibersortX.Bcell, na.rm = TRUE)),    c("#F5F6FA", "#6C5CE7")),
  Cbx_endothelial    = colorRamp2(c(0, max_TILs    <- max(annotation_col$CibersortX.endothelial, na.rm = TRUE)),    c("#F5F6FA", "#6C5CE7")),
  Cbx_tcell    = colorRamp2(c(0, max_TILs    <- max(annotation_col$CibersortX.Tcell, na.rm = TRUE)),    c("#F5F6FA", "#6C5CE7")),
  ASCAT_PLOIDY    = colorRamp2(c(0, max_TILs    <- max(annotation_col$ASCAT_PLOIDY, na.rm = TRUE)),    c("#F5F6FA", "#6C5CE7")),
  ASCAT_TUM_FRAC    = colorRamp2(c(0, max_TILs    <- max(annotation_col$ASCAT_TUM_FRAC, na.rm = TRUE)),    c( "#6C5CE7", "#F5F6FA"))
)


# Costruzione dell'oggetto HeatmapAnnotation (la sidebar superiore)
col_ann <- HeatmapAnnotation(
  df = annotation_col, 
  col = ann_colors,
  show_legend = TRUE
)

# Generazione del Main Heatmap di espressione genica
ht <- Heatmap(
  plot_matrix, 
  name = "Expression",          
  top_annotation = col_ann,     
  show_row_names = FALSE,       
  show_column_names = FALSE,    
  cluster_columns = FALSE,      # I campioni rimangono raggruppati rigidamente per Cluster RNAseq
  cluster_rows = TRUE,          # I geni (i marker identificati prima) vengono clusterizzati tra loro
  col = colorRampPalette(c("blue", "white", "red"))(100),
  width = unit(12, "cm") 
)

# Rendering dell'heatmap a schermo
draw(ht, 
     heatmap_legend_side = "right", 
     annotation_legend_side = "right",
     legend_grouping = "adjusted" 
)
```



# Loading data for comparison

```{r}
library(circlize)
library(ComplexHeatmap)
library(RColorBrewer)
library(readxl)
```

```{r}
# metadata
dt_metadata <- read_excel("ids_cohorts_match.xlsx", sheet = "1a SCAN-B discovery") 

rna_1000 <- read.csv("rnaseq_nmf_clusters_1000.csv")

# NMF clusters from genes diff. expressed between samples with different methyl. status in probes inferred to be the enhancers or promoters of those genes (ELMER)
methyl_clusters_df <- read.csv("methyl_nmf_clusters_genes.csv") 

# NMF clusters from genes diff. expressed between samples with different methyl. status in probes inferred to be the enhancers of those genes (ELMER)
methyl_en_df <- read.csv("methyl_nmf_clusters_genes_en.csv") 

# NMF clusters from genes diff. expressed between samples with different methyl. status in probes inferred to be the promoters of those genes (ELMER)
methyl_pr_df <- read.csv("methyl_nmf_clusters_genes_pr.csv")

# NMF clusters from all probes linked to the genes diff. expressed between samples with different methyl. status (ELMER)
methyl_probes <- read.csv("methyl_nmf_clusters_probes.csv") 

# NMF clusters from CNV data on genes with amplification/deletion and overexpression/underexpression respectively
cnv_data_clusters <- read.csv("cnv_nmf_clusters.csv")

# same but on rnaseq data
cnv_rnaseq_clusters <- read.csv("cnv_rnaseq_nmf_clusters.csv")

#mofa
mofa_clusters_2 <- read.csv("mofa_km_clusters_2_no.csv")
mofa_clusters_3 <- read.csv("mofa_km_clusters_3_no.csv")
mofa_clusters_4 <- read.csv("mofa_km_clusters_4_no.csv")


#icluster
icluster_2 <- read.csv("iclusters_clusters_2.csv")
icluster_3 <- read.csv("iclusters_clusters_3.csv")
icluster_4 <- read.csv("iclusters_clusters_4.csv")

#iclusterbayes
ibayes_3 <- read.csv("iclusters_bayes_3.csv")
ibayes_4 <- read.csv("iclusters_bayes_4.csv")

#SNF
SNF_clusters_2 <- read.csv("SNF_clusters_2.csv")
SNF_clusters_4 <- read.csv("SNF_clusters_4.csv")
SNF_clusters_3 <- read.csv("SNF_clusters_3.csv")


#x-intNMF
XintNMF_clusters_2 <- read.csv("xintNMF_clusters_k2.csv")
XintNMF_clusters_3 <- read.csv("xintNMF_clusters_k3.csv")
XintNMF_clusters_4 <- read.csv("xintNMF_clusters_k4.csv")

XintNMF_clusters_2_reg <- read.csv("xintNMF_clusters_k2_reg.csv")
XintNMF_clusters_3_reg <- read.csv("xintNMF_clusters_k3_reg.csv")
XintNMF_clusters_4_reg <- read.csv("xintNMF_clusters_k4_reg.csv")

```


# METADATA

#columns of interest for metadata
# categorical
dt_metadata$Grade
dt_metadata$LNbinary 
dt_metadata$TNBCtype4_n235_notPreCentered
dt_metadata$TNBCtype6_n235_notPreCentered
dt_metadata$PAM50_Basal_NCN

# continuous
dt_metadata$TMB
dt_metadata$TumSize
dt_metadata$TILs

```{r}
cols_to_keep <- c("TMB", "TumSize", "Grade", "LNbinary", "TILs", "PAM50_Basal_NCN", "PAM50_NCN", "TNBCtype4_n235_notPreCentered", "TNBCtype6_n235_notPreCentered", "ASCAT_PLOIDY", "ASCAT_TUM_FRAC")
# Subset metadata
meta_sub <- dt_metadata[, cols_to_keep]

#Ensure numeric columns are numeric 
numeric_cols <- c("TMB", "TumSize", "TILs", "ASCAT_PLOIDY", "ASCAT_TUM_FRAC")
meta_sub[numeric_cols] <- lapply(meta_sub[numeric_cols], function(x) as.numeric(as.character(x)))

#Ensure categorical columns are factors
cat_cols <- setdiff(cols_to_keep, numeric_cols)
meta_sub[cat_cols] <- lapply(meta_sub[cat_cols], factor)
```
```{r}
# Create the data frame for the heatmap sidebar
# Use the sample_order from your previous code to keep them grouped by NMF cluster
annotation_col <- data.frame(
  Cluster = factor(cluster_data$Cluster[sample_order]),
  meta_sub[sample_order, ] # Align metadata to the sample order
)
# Set row names to match the plot_matrix columns
rownames(annotation_col) <- colnames(plot_matrix)


annotation_col$Cluster <- factor(cluster_data$Cluster[sample_order])

ann_colors = list(
  Cluster = c("1" = "tomato3", "2" = "#0984E3", "3" = "#00B894"),
  Grade = c("2" = "#FAB1A0", "3" = "#E17055", "NA" = "#B2BEC3"),
  LNbinary = c("positive" = "#636E72", "negative" = "#DFE6E9"),
  TNBCtype4_n235_notPreCentered = c("BL1" = "#A29BFE", "BL2" = "#74B9FF", 
                                    "M" = "#55E6C1", "LAR" = "#FDCB6E", "NA"= "#B2BEC3"),
  TNBCtype6_n235_notPreCentered = c("BL1" = "#A29BFE", "BL2" = "#74B9FF", "M" = "#55E6C1", 
                                    "LAR" = "#FDCB6E", "NA"= "#B2BEC3", "IM"= "#006266", 
                                    "MSL"= "#FFEAA7", "UNS"= "#FFADAD"),
  PAM50_NCN = c(Basal = "#D63031", Her2 = "#A29BFE", LumB = "#0984E3", 
                LumA = "#74B9FF", Normal = "#00B894", unclassified = "#B2BEC3"),
  PAM50_Basal_NCN = c(Basal = "#D63031", nonBasal = "#E0E0E0"),
  TumSize = colorRamp2(c(0, max(annotation_col$TumSize, na.rm = TRUE)), c("#F5F6FA", "#273C75")),
  TMB     = colorRamp2(c(0, max_TMB     <- max(annotation_col$TMB, na.rm = TRUE)),     c("#F5F6FA", "#079992")),
  TILs    = colorRamp2(c(0, max_TILs    <- max(annotation_col$TILs, na.rm = TRUE)),    c("#F5F6FA", "#6C5CE7")),
  ASCAT_PLOIDY = colorRamp2(c(0, max_TMB     <- max(annotation_col$ASCAT_PLOIDY, na.rm = TRUE)),     c("#F5F6FA", "#273C75")),
  ASCAT_TUM_FRAC = colorRamp2(c(0, max_TMB     <- max(annotation_col$ASCAT_TUM_FRAC, na.rm = TRUE)),     c("#F5F6FA", "#079992"))
)


col_ann <- HeatmapAnnotation(
  df = annotation_col, 
  col = ann_colors,
  show_legend = TRUE
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
  width = unit(12, "cm") 
)

draw(ht, 
     heatmap_legend_side = "right", 
     annotation_legend_side = "right",
     legend_grouping = "adjusted" 
)
```


# single omics and multiomics

```{r}

# Vectors to map cluster to samples
rna_1000_map  <- setNames(rna_1000$Cluster, rna_1000$SampleID)
methyl_map    <- setNames(methyl_clusters_df$Cluster, methyl_clusters_df$SampleID)
methyl_map_en <- setNames(methyl_en_df$Cluster,       methyl_en_df$SampleID)
methyl_map_pr <- setNames(methyl_pr_df$Cluster,       methyl_pr_df$SampleID)
methyl_map_probes <- setNames(methyl_probes$Cluster,       methyl_probes$SampleID)
cnv_map <- setNames(cnv_data_clusters$Cluster, cnv_data_clusters$SampleID)
cnv_rnaseq_map <- setNames(cnv_rnaseq_clusters$Cluster, cnv_rnaseq_clusters$SampleID)
mofa_map_2 <- setNames(mofa_clusters_2$Cluster, mofa_clusters_2$SampleID)
mofa_map_3 <- setNames(mofa_clusters_3$Cluster, mofa_clusters_3$SampleID)
mofa_map_4 <- setNames(mofa_clusters_4$Cluster, mofa_clusters_4$SampleID)
SNF_map_2 <- setNames(SNF_clusters_2$Cluster, SNF_clusters_2$SampleID)
SNF_map_3 <- setNames(SNF_clusters_3$Cluster, SNF_clusters_3$SampleID)
SNF_map_4 <- setNames(SNF_clusters_4$Cluster, SNF_clusters_4$SampleID)
xint_map_2 <- setNames(XintNMF_clusters_2$Cluster, XintNMF_clusters_2$SampleID)
xint_map_3 <- setNames(XintNMF_clusters_3$Cluster, XintNMF_clusters_3$SampleID)
xint_map_4 <- setNames(XintNMF_clusters_4$Cluster, XintNMF_clusters_4$SampleID)
xint_map_2_reg <- setNames(XintNMF_clusters_2_reg$Cluster, XintNMF_clusters_2_reg$SampleID)
xint_map_3_reg <- setNames(XintNMF_clusters_3_reg$Cluster, XintNMF_clusters_3_reg$SampleID)
xint_map_4_reg <- setNames(XintNMF_clusters_4_reg$Cluster, XintNMF_clusters_4_reg$SampleID)
icluster_map_2 <- setNames(icluster_2$Cluster, icluster_2$SampleID)
icluster_map_3 <- setNames(icluster_3$Cluster, icluster_3$SampleID)
icluster_map_4 <- setNames(icluster_4$Cluster, icluster_4$SampleID)
ibayes_map_3 <- setNames(ibayes_3$Cluster, ibayes_3$SampleID)
ibayes_map_4 <- setNames(ibayes_4$Cluster, ibayes_4$SampleID)

# Allineiamo tutti i cluster alle colonne di plot_matrix
annotation_col_meth <- data.frame(
  RNAseq                 = factor(cluster_data$Cluster[sample_order]),
  RNAseq_1000            = factor(rna_1000_map[colnames(plot_matrix)]),
  #Methyl_genes           = factor(methyl_map[colnames(plot_matrix)]),
  #Methyl_genes_en        = factor(methyl_map_en[colnames(plot_matrix)]),
  #Methyl_genes_pr        = factor(methyl_map_pr[colnames(plot_matrix)]),
  Methyl_probes          = factor(methyl_map_probes[colnames(plot_matrix)]),
  CNV                    = factor(cnv_map[colnames(plot_matrix)]),
  #CNV_rnaseq             = factor(cnv_rnaseq_map[colnames(plot_matrix)]),
  MOFA_2                 = factor(mofa_map_2[colnames(plot_matrix)]),
  MOFA_3                 = factor(mofa_map_3[colnames(plot_matrix)]),
  MOFA_4                 = factor(mofa_map_4[colnames(plot_matrix)]),
  SNF_2                  = factor(SNF_map_2[colnames(plot_matrix)]),
  SNF_3                  = factor(SNF_map_3[colnames(plot_matrix)]),
  SNF_4                  = factor(SNF_map_4[colnames(plot_matrix)]),
  icluster_2             = factor(icluster_map_2[colnames(plot_matrix)]),
  icluster_3             = factor(icluster_map_3[colnames(plot_matrix)]),
  icluster_4             = factor(icluster_map_4[colnames(plot_matrix)]),
  ibayes_3               = factor(ibayes_map_3[colnames(plot_matrix)]),
  ibayes_4               = factor(ibayes_map_4[colnames(plot_matrix)]),
  Xint2                  = factor(xint_map_2[colnames(plot_matrix)]),
  Xint3                  = factor(xint_map_3[colnames(plot_matrix)]),
  Xint4                  = factor(xint_map_4[colnames(plot_matrix)]),
  Xint2_reg                  = factor(xint_map_2_reg[colnames(plot_matrix)]),
  Xint3_reg                  = factor(xint_map_3_reg[colnames(plot_matrix)]),
  Xint4_reg                  = factor(xint_map_4_reg[colnames(plot_matrix)])
  )


# Impostiamo i nomi delle righe per farli coincidere con le colonne della heatmap
rownames(annotation_col_meth) <- colnames(plot_matrix)

ann_colors_meth = list(
  RNAseq          = c("1" = "tomato3", "2" = "#0984E3", "3" = "#00B894"),
  RNAseq_1000     = c("1" = "tomato3", "2" = "#0984E3", "3" = "#00B894"),
  #Methyl_genes    = c("1" = "#fab1a0", "2" = "#a29bfe", "3" = "#55efc4"),
  
  # Sostituiti i duplicati con nuove palette pastello ad alto contrasto:
  #Methyl_genes_en = c("1" = "#ff9ff3", "2" = "#feca57", "3" = "#54a0ff"), 
  #Methyl_genes_pr = c("1" = "#ffeaa7", "2" = "#d63031", "3" = "#2d3436"),
  Methyl_probes   = c("1" = "#ffb8b8", "2" = "#95afc0", "3" = "#badc58"), 
  CNV             = c("1" = "#ff7675", "2" = "#74b9ff", "3" = "#fdcb6e"), 
  #CNV_rnaseq      = c("1" = "#e67e22", "2" = "#2ecc71", "3" = "#9b59b6"),
  MOFA_2          = c("1" = "#ffeaa7", "2" = "#0984E3"),
  MOFA_3          = c("1" = "#0984E3", "2" = "#9b59b6", "3" = "#ffeaa7"),
  MOFA_4          = c("1" = "#9b59b6", "2" = "#ffeaa7", "3" = "#0984E3", "4" = "#ff7f98"),
  icluster_2      = c("1" = "#27ae33", "2" = "#e1b99c"),
  icluster_3      = c("1" = "#27ae33", "2" = "#e1b99c", "3" = "#4c5ce1"),
  icluster_4      = c("1" = "#4c5ce1", "2" = "#e1b99c", "3" = "#ff6f11", "4" = "#27ae33"),
  ibayes_3        = c("1" = "#27ae33", "2" = "#e1b99c", "3" = "#4c5ce1"),
  ibayes_4        = c("1" = "#4c5ce1", "2" = "#e1b99c", "3" = "#ff6f11", "4" = "#27ae33"),
  SNF_2           = c("1" = "#ff7f50", "2" = "#008080"),
  SNF_3           = c("1" = "#ff7f50", "2" = "#34495e", "3" = "#008080"),
  SNF_4           = c("1" = "#ff7f50", "2" = "#008080", "3" = "#34495e", "4" = "#f1c40f"),
  Xint2           = c("1" = "#27ae60", "2" = "#e1b12c"),
  Xint3           = c("1" = "#e1b12c", "2" = "#27ae60", "3" = "#6c5ce7"),
  Xint4           = c("1" = "#6c5ce7", "2" = "#27ae60", "3" = "#ff6f99", "4" = "#e1b12c"),
  Xint2_reg           = c("1" = "#27ae60", "2" = "#e1b12c"),
  Xint3_reg           = c("1" = "#27ae60", "2" = "#e1b12c", "3" = "#6c5ce7"),
  Xint4_reg           = c("1" = "#27ae60", "2" = "#e1b12c", "3" = "#6c5ce7", "4" = "#ff6f99")
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

draw(ht_meth, 
     heatmap_legend_side = "right", 
     annotation_legend_side = "right"
)
```

```{r}
# Assicurati che plot_matrix sia già ordinata secondo sample_order come nei tuoi blocchi precedenti

# 1. Estraiamo i metadati e i cluster di metilazione per i campioni presenti in plot_matrix
# Usiamo i nomi delle colonne di plot_matrix per garantire l'allineamento perfetto

cols_to_keep <- c("TNBCtype4_n235_notPreCentered", "TNBCtype6_n235_notPreCentered", "ASCAT_PLOIDY", "ASCAT_TUM_FRAC")
# Subset metadata
meta_sub <- dt_metadata[, cols_to_keep]

#Ensure numeric columns are numeric 
numeric_cols <- c("ASCAT_PLOIDY", "ASCAT_TUM_FRAC")
meta_sub[numeric_cols] <- lapply(meta_sub[numeric_cols], function(x) as.numeric(as.character(x)))

#Ensure categorical columns are factors
cat_cols <- c("TNBCtype4_n235_notPreCentered", "TNBCtype6_n235_notPreCentered")
meta_sub[cat_cols] <- lapply(meta_sub[cat_cols], factor)

annotation_selected <- data.frame(
  RNAseq           = factor(cluster_data$Cluster[sample_order]),
  RNAseq_1000      = factor(rna_1000_map[colnames(plot_matrix)]),
  Methylation      = factor(methyl_map_probes[colnames(plot_matrix)]),
  CNV              = factor(cnv_map[colnames(plot_matrix)]),
  #CNV_rnaseq       = factor(cnv_rnaseq_map[colnames(plot_matrix)]),
  MOFA_2                 = factor(mofa_map_2[colnames(plot_matrix)]),
  MOFA_3                 = factor(mofa_map_3[colnames(plot_matrix)]),
  #SNF_2                  = factor(SNF_map_2[colnames(plot_matrix)]),
  SNF_4                  = factor(SNF_map_4[colnames(plot_matrix)]),
  TNBC_type4       = meta_sub[sample_order, ]$TNBCtype4_n235_notPreCentered,
  TNBC_type6       = meta_sub[sample_order, ]$TNBCtype6_n235_notPreCentered,
  ASCAT_PLOIDY       = meta_sub[sample_order, ]$ASCAT_PLOIDY,
  ASCAT_TUM_FRAC       = meta_sub[sample_order, ]$ASCAT_TUM_FRAC
)

# Impostiamo i nomi delle righe per ComplexHeatmap
rownames(annotation_selected) <- colnames(plot_matrix)
```

```{r}
ann_colors_selected = list(
  RNAseq       = c("1" = "tomato3", "2" = "#0984E3", "3" = "#00B894"),
  RNAseq_1000     = c("1" = "#e65055", "2" = "#06cce7", "3" = "#00B666"),
  Methylation  = c("1" = "#e17055", "2" = "#6c5ce7", "3" = "#00cec9"),
  CNV             =  c("1" = "#e67e22", "2" = "#2ecc71", "3" = "#9b59b6"), 
  #Methylation  = c("1" = "#00B894", "2" = "tomato3", "3" = "#0984E3"),
  #CNV_rnaseq      = c("1" = "#ff7675", "2" = "#74b9ff", "3" = "#2ecc71"),
  MOFA_2          = c("1" = "#ffeaa7", "2" = "#0984E3"),
  MOFA_3          = c("1" = "#ffeaa7", "2" = "#0984E3", "3" = "#9b59b6"),
  SNF_2           = c("1" = "#ff7f50", "2" = "#008080"),
  SNF_3           = c("1" = "#ff7f50", "2" = "#008080", "3" = "#34495e"),
  SNF_4           = c("1" = "#ff7f50", "2" = "#008080", "3" = "#f1c40f", "4" = "#34495e"),
  
  
  TNBC_type4   = c("BL1" = "#A29BFE", "BL2" = "#74B9FF", 
                                    "M" = "#55E6C1", "LAR" = "#FDCB6E", "NA"= "#B2BEC3"),
  TNBC_type6   = c("BL1" = "#A29BFE", "BL2" = "#74B9FF", "M" = "#55E6C1", 
                                    "LAR" = "#FDCB6E", "NA"= "#B2BEC3", "IM"= "#006266", 
                                    "MSL"= "#FFEAA7", "UNS"= "#FFADAD"),
  ASCAT_PLOIDY = colorRamp2(c(0, max_TMB     <- max(annotation_col$ASCAT_PLOIDY, na.rm = TRUE)),     c("#F5F6FA", "#273C75")),
  ASCAT_TUM_FRAC = colorRamp2(c(0, max_TMB     <- max(annotation_col$ASCAT_TUM_FRAC, na.rm = TRUE)),     c("#F5F6FA", "#079992"))
)
```

```{r}
# Creazione dell'annotazione superiore
col_ann_final <- HeatmapAnnotation(
  df = annotation_selected,
  col = ann_colors_selected,
  show_legend = TRUE,
  annotation_name_side = "left"
)

# Generazione della Heatmap
ht_final <- Heatmap(
  plot_matrix, 
  name = "Expression",          
  top_annotation = col_ann_final,     
  show_row_names = FALSE,        
  show_column_names = FALSE,    
  cluster_columns = FALSE,      # Fondamentale: mantiene l'ordine basato sui cluster RNA
  cluster_rows = TRUE,          
  col = colorRampPalette(c("blue", "white", "red"))(100),
  #column_title = "Samples aligned by RNA Cluster, TNBC types, and Methyl Probes",
  width = unit(10, "cm") 
)

# Disegno del plot
draw(ht_final, 
     heatmap_legend_side = "right", 
     annotation_legend_side = "right"
)
```
