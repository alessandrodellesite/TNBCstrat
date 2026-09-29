# RNAseq data loading and preprocessing

library(limma)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

rna <- readRDS("/mnt/petasan_ccb/alessandro/SCANB/rna_logtransformed.rds")
gene_mads <- apply(rna, 1, mad)
ordered_mads <- order(gene_mads, decreasing = TRUE)
top_3000_indices <- ordered_mads[1:3000]
rna <- rna[top_3000_indices, ]
dt_matrix <- as.matrix(rna)

cluster_data <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/singleomic/nmf_rna/rna_nmf_clusters.csv")

# FIX: verify the sample match actually worked before proceeding.
# NMF's predict() writes SampleID as names(sample_groups_rna), which come from
# colnames(rna_matrix) at clustering time -- these should match colnames(dt_matrix)
# here, but we check rather than assume.
idx <- match(colnames(dt_matrix), cluster_data$SampleID)
if (any(is.na(idx))) {
  stop("Some samples in dt_matrix were not found in cluster_data$SampleID: ",
       paste(colnames(dt_matrix)[is.na(idx)], collapse = ", "))
}
cluster_data <- cluster_data[idx, ]
stopifnot(identical(colnames(dt_matrix), cluster_data$SampleID))

# FIX: relabel the RAW NMF cluster IDs to the desired display order ONCE,
# here, before any DE/design/contrasts/GO code runs -- rather than relabeling
# only the heatmap at the end. This guarantees "Cluster 1/2/3" means the same
# patients everywhere downstream: the design matrix, the contrasts, the
# results tables, the GO enrichment, and the heatmap all key off this single
# relabeled Cluster column, so numbering can never drift out of sync again.
# raw "3" -> display "1" (leftmost), raw "1" -> display "2" (center),
# raw "2" -> display "3" (rightmost)
relabel_map <- c("3" = "1", "1" = "2", "2" = "3")
cluster_data$Cluster <- factor(relabel_map[as.character(cluster_data$Cluster)],
                                levels = c("1", "2", "3"))

# Differential expression analysis

#create the factor and design matrix
groups <- cluster_data$Cluster

# FIX: don't blindly hardcode colnames(design) <- c("C1","C2","C3").
# model.matrix names columns after levels(groups) in whatever order factor()
# produces. NMF's predict() returns cluster indices 1..rank, so this should be
# "1","2","3" sorted lexically -- but we verify instead of assuming, and derive
# the column names from the actual levels rather than hardcoding them.
stopifnot(setequal(levels(groups), c("1", "2", "3")))
design <- model.matrix(~0 + groups)
colnames(design) <- paste0("C", levels(groups))  # -> C1, C2, C3, matched to actual levels
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

# FIX: topTable(sort.by="logFC") sorts by SIGNED logFC, not |logFC|. Combined
# with number=500, this truncates to the 500 MOST POSITIVE logFC genes and can
# silently drop every significant down-regulated gene before it ever reaches
# annotation or enrichGO(logFC < 0). Use number=Inf to get all significant
# genes, sorted by significance (p-value / B-stat), not by logFC direction.
results_1 <- topTable(fit2, coef = "C1_vs_others", number = Inf, p.value = 0.05, sort.by = "P")
results_2 <- topTable(fit2, coef = "C2_vs_others", number = Inf, p.value = 0.05, sort.by = "P")
results_3 <- topTable(fit2, coef = "C3_vs_others", number = Inf, p.value = 0.05, sort.by = "P")

#function to pull the top significant marker genes (for the heatmap only)
# FIX: get ALL significant genes first (number = Inf), then rank by |logFC|
# and take the top n -- this avoids the signed-logFC truncation bug above.
get_cluster_genes <- function(fit_obj, coef_name, n_genes = 500) {
  top_tab <- topTable(fit_obj, coef = coef_name, number = Inf, p.value = 0.05, sort.by = "P")
  top_tab <- top_tab[order(-abs(top_tab$logFC)), ]
  top_tab <- head(top_tab, n_genes)
  return(rownames(top_tab))
}

# Heatmap

c1_markers <- get_cluster_genes(fit2, "C1_vs_others", n_genes = 50)
c2_markers <- get_cluster_genes(fit2, "C2_vs_others", n_genes = 50)
c3_markers <- get_cluster_genes(fit2, "C3_vs_others", n_genes = 50)
ordered_genes <- c(c1_markers, c2_markers, c3_markers)

# cluster_data$Cluster was already relabeled to the desired display order right
# after loading (raw 3/1/2 -> display 1/2/3), so a plain ascending sort here
# gives left-to-right = 1, 2, 3 -- and these are the SAME cluster identities
# used in the C1/C2/C3 contrasts and GO plots above, so numbering matches.
sample_order <- order(cluster_data$Cluster)

plot_matrix <- dt_matrix[ordered_genes, sample_order]
plot_matrix <- t(scale(t(plot_matrix)))

annotation_col <- data.frame(RNAseq = cluster_data$Cluster[sample_order])
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
  # Note: rownames are already clean ENSEMBL IDs (cleaned upstream when
  # rna_logtransformed.rds was created), so this gsub is a no-op safety net.
  clean_ids <- gsub("\\..*", "", rownames(results_df))

  tryCatch({
    #Add ENTREZID back to the toType list
    anno <- bitr(clean_ids,
                 fromType = "ENSEMBL",
                 toType = c("SYMBOL", "GENENAME", "ENTREZID"),
                 OrgDb = org.Hs.eg.db,
                 drop = TRUE)

    # FIX: bitr's ENSEMBL->SYMBOL/ENTREZID mapping is not 1:1; some ENSEMBL
    # IDs map to multiple ENTREZIDs (or vice versa), which duplicates rows
    # after merge(). Deduplicate on the join key to keep one row per gene.
    anno <- anno[!duplicated(anno$ENSEMBL), ]

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

print(head(anno_c1, 20))
print(head(anno_c2, 20))
print(head(anno_c3, 20))


# Enrichment of each cluster with gene ontology

# FIX: build the background universe from the genes actually tested in the DE
# analysis (the 3000 MAD-filtered genes), not clusterProfiler's default of
# every annotated gene in org.Hs.eg.db. Testing DE on 3000 genes but enriching
# against the whole genome as background biases the enrichment p-values.
universe_ensembl <- gsub("\\..*", "", rownames(dt_matrix))
universe_map <- bitr(universe_ensembl,
                      fromType = "ENSEMBL", toType = "ENTREZID",
                      OrgDb = org.Hs.eg.db, drop = TRUE)
universe_ids <- unique(universe_map$ENTREZID)

#POSITIVE
ego1_up <- enrichGO(gene          = anno_c1$ENTREZID[anno_c1$logFC > 0],
                     universe      = universe_ids,
                     OrgDb         = org.Hs.eg.db,
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE) #converts IDs back to Symbols in the results

ego2_up <- enrichGO(gene          = anno_c2$ENTREZID[anno_c2$logFC > 0],
                 universe      = universe_ids,
                 OrgDb         = org.Hs.eg.db,
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05,
                 readable      = TRUE)

ego3_up <- enrichGO(gene          = anno_c3$ENTREZID[anno_c3$logFC > 0],
                     universe      = universe_ids,
                     OrgDb         = org.Hs.eg.db,
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)

#NEGATIVE

ego1_down <- enrichGO(gene          = anno_c1$ENTREZID[anno_c1$logFC < 0],
                      universe      = universe_ids,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = "ENTREZID",
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

ego2_down <- enrichGO(gene          = anno_c2$ENTREZID[anno_c2$logFC < 0],
                      universe      = universe_ids,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = "ENTREZID",
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

ego3_down <- enrichGO(gene          = anno_c3$ENTREZID[anno_c3$logFC < 0],
                      universe      = universe_ids,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = "ENTREZID",
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

# FIX: guard against a cluster having too few (or zero) up/down genes to
# produce enrichment results, which would otherwise crash merge_result()/dotplot().
safe_merge_and_plot <- function(ego_up, ego_down, title, outfile) {
  has_up   <- !is.null(ego_up)   && nrow(as.data.frame(ego_up))   > 0
  has_down <- !is.null(ego_down) && nrow(as.data.frame(ego_down)) > 0

  if (!has_up && !has_down) {
    message("No enrichment results (up or down) for '", title, "'; skipping plot.")
    return(invisible(NULL))
  }

  res_list <- list()
  if (has_up)   res_list[["Upregulated"]]   <- ego_up
  if (has_down) res_list[["Downregulated"]] <- ego_down

  merged_ego <- merge_result(res_list)

  p <- dotplot(merged_ego, x = "Cluster", showCategory = 10) +
  ggtitle(title) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 12),   # <- bigger pathway names (was 7)
    plot.title  = element_text(size = 14)
  ) +
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 40))

  png(outfile, width = 1000, height = 1100, res = 150)   # narrower, a bit taller
  print(p)
  dev.off()                   

  invisible(p)
}

cl1 <- safe_merge_and_plot(ego1_up, ego1_down, "Cluster 1: Directional Pathways",
                            "/mnt/petasan_ccb/alessandro/SCANB/plots/comparisons/rna_enrichment_cl1.png")
cl2 <- safe_merge_and_plot(ego2_up, ego2_down, "Cluster 2: Directional Pathways",
                            "/mnt/petasan_ccb/alessandro/SCANB/plots/comparisons/rna_enrichment_cl2.png")
cl3 <- safe_merge_and_plot(ego3_up, ego3_down, "Cluster 3: Directional Pathways",
                            "/mnt/petasan_ccb/alessandro/SCANB/plots/comparisons/rna_enrichment_cl3.png")
