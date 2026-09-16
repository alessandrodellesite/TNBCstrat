# RNAseq data loading and preprocessing

library(limma)
library(pheatmap)

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

# Differential expression analysis

#create the factor and design matrix
groups <- factor(cluster_data$Cluster) 
design <- model.matrix(~0 + groups)
colnames(design) <- c("C1", "C2", "C3")
fit <- lmFit(dt_matrix, design) 

#to compute log2FC of each gene in each cluster compared to the average of the gene in the other 2 clusters 
cont <- makeContrasts(
  C1_vs_others = C1 - (C2 + C3)/2,  
  C2_vs_others = C2 - (C1 + C3)/2,
  C3_vs_others = C3 - (C1 + C2)/2,
  levels = design
)
fit2 <- contrasts.fit(fit, cont)
fit2 <- eBayes(fit2, trend = TRUE) 

#function to pull the top significant genes
get_cluster_genes <- function(fit_obj, coef_name, n_genes = 500) {
  top_tab <- topTable(fit_obj, coef = coef_name, number = n_genes, p.value = 0.05, sort.by = "logFC")
  return(rownames(top_tab)) 
}

# Heatmap

c1_markers <- get_cluster_genes(fit2, "C1_vs_others", n_genes = 50)
c2_markers <- get_cluster_genes(fit2, "C2_vs_others", n_genes = 50)
c3_markers <- get_cluster_genes(fit2, "C3_vs_others", n_genes = 50)
ordered_genes <- c(c1_markers, c2_markers, c3_markers)

# Relabel old cluster IDs so sorting gives visual order 3,1,2 
# old "3" -> new "1", old "1" -> new "2", old "2" -> new "3"
relabel_map <- c("3" = "1", "1" = "2", "2" = "3")
cluster_data$Cluster_relabelled <- factor(relabel_map[as.character(cluster_data$Cluster)],
                                           levels = c("1", "2", "3"))
sample_order <- order(cluster_data$Cluster_relabelled)

plot_matrix <- dt_matrix[ordered_genes, sample_order]
plot_matrix <- t(scale(t(plot_matrix)))

annotation_col <- data.frame(RNAseq = cluster_data$Cluster_relabelled[sample_order])
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


