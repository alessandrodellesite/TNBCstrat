# RNAseq data loading and preprocessing

library(limma)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

#library(circlize)
#library(ComplexHeatmap)
#library(RColorBrewer)
#library(readxl)

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

results_1 <- topTable(fit2, coef="C1_vs_others", number=500, p.value=0.05, sort.by="logFC")
results_2 <- topTable(fit2, coef="C2_vs_others", number=500, p.value=0.05, sort.by="logFC")
results_3 <- topTable(fit2, coef="C3_vs_others", number=500, p.value=0.05, sort.by="logFC")

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

ph <- pheatmap(plot_matrix, 
         annotation_col = annotation_col, 
         cluster_rows = TRUE,  # Keep genes grouped by cluster
         cluster_cols = FALSE,  # Keep samples grouped by cluster
         show_colnames = FALSE, 
         show_rownames = FALSE, 
         main = "Ordered NMF Cluster Markers",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         #breaks = seq(-2, 2, length.out = 101)
         ) 

png("/mnt/petasan_ccb/alessandro/SCANB/plots/comparisons/heatmap_rna.png", 
    width = 11.69, height = 8.27, units = "in", res = 300)
print(ph)
dev.off()




# Annotation

annotate_with_stats <- function(results_df) {
  clean_ids <- gsub("\\..*", "", rownames(results_df))
  
  tryCatch({
    #Add ENTREZID back to the toType list
    anno <- bitr(clean_ids, 
                 fromType = "ENSEMBL", 
                 toType = c("SYMBOL", "GENENAME", "ENTREZID"), 
                 OrgDb = org.Hs.eg.db,
                 drop = TRUE)
    
    results_df$ENSEMBL <- clean_ids
    merged <- merge(anno, results_df, by = "ENSEMBL")
    
    # Sort by absolute logFC
    merged$absLogFC <- abs(merged$logFC) #PRENDE ABS VALUES
    final_table <- merged[order(-merged$absLogFC), ]
    
    # Reset indices
    rownames(final_table) <- NULL
    
    #keep ENTREZID in the dataframe now so enrichGO can use it
    final_table <- final_table[, c("SYMBOL", "GENENAME", "logFC", "adj.P.Val", "ENSEMBL", "ENTREZID")]
    
    return(final_table)
  }, error = function(e) {
    message("Error during annotation: ", e)
    return(NULL)
  })
}

#Run the annotation for clusters
anno_c1 <- annotate_with_stats(results_1)
anno_c2 <- annotate_with_stats(results_2)
anno_c3 <- annotate_with_stats(results_3)

print(anno_c1[1:20,])
print(anno_c2[1:20,])
print(anno_c3[1:20,])



# Enrichment of each cluster with gene ontology

#POSITIVE
ego1_up <- enrichGO(gene          = anno_c1$ENTREZID[anno_c1$logFC > 0],
                     OrgDb         = org.Hs.eg.db,
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE) #converts IDs back to Symbols in the results

ego2_up <- enrichGO(gene          = anno_c2$ENTREZID[anno_c2$logFC > 0],
                 OrgDb         = org.Hs.eg.db,
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05,
                 readable      = TRUE)

ego3_up <- enrichGO(gene          = anno_c3$ENTREZID[anno_c3$logFC > 0],
                     OrgDb         = org.Hs.eg.db,
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)

#NEGATIVE

ego1_down <- enrichGO(gene          = anno_c1$ENTREZID[anno_c1$logFC < 0],
                      OrgDb         = org.Hs.eg.db,
                      keyType       = "ENTREZID", 
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05, 
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

ego2_down <- enrichGO(gene          = anno_c2$ENTREZID[anno_c2$logFC < 0],
                      OrgDb         = org.Hs.eg.db,
                      keyType       = "ENTREZID", 
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05, 
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

ego3_down <- enrichGO(gene          = anno_c3$ENTREZID[anno_c3$logFC < 0],
                      OrgDb         = org.Hs.eg.db,
                      keyType       = "ENTREZID", 
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05, 
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

merged_ego1 <- merge_result(list(Upregulated = ego1_up, Downregulated = ego1_down))
merged_ego2 <- merge_result(list(Upregulated = ego2_up, Downregulated = ego2_down))
merged_ego3 <- merge_result(list(Upregulated = ego3_up, Downregulated = ego3_down))

# plot

cl1 <- dotplot(merged_ego1, x = "Cluster", showCategory = 10) + 
  ggtitle("Cluster 1: Directional Pathways") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 7) 
  )+
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 50))

cl2 <- dotplot(merged_ego2, x = "Cluster", showCategory = 10) + 
  ggtitle("Cluster 2: Directional Pathways") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 7) 
  )+
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 50))

cl3 <- dotplot(merged_ego3, x = "Cluster", showCategory = 10) + 
  ggtitle("Cluster 3: Directional Pathways") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 7) 
  )+
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 50))

                   
png(file.path("/mnt/petasan_ccb/alessandro/SCANB/plots/comparisons",
               "rna_enrichment_cl1.png"),
    width = 1400, height = 1000, res = 150)
print(cl1)
dev.off()

png(file.path("/mnt/petasan_ccb/alessandro/SCANB/plots/comparisons",
               "rna_enrichment_cl2.png"),
    width = 1400, height = 1000, res = 150)
print(cl2)
dev.off()

png(file.path("/mnt/petasan_ccb/alessandro/SCANB/plots/comparisons",
               "rna_enrichment_cl3.png"),
    width = 1400, height = 1000, res = 150)
print(cl3)
dev.off()                   
