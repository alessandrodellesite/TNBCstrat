# NMF on feature selected single-omic data

#to avoid NMF forcing cores=2
Sys.unsetenv("_R_CHECK_TIMINGS_")
Sys.unsetenv("_R_CHECK_CRAN_INCOMING_")
library(NMF)


# RNAseq

#RNAseq feature selection (top 3000 genes by MAD)
print("Rna preprocessing")
rna_data <- readRDS("/mnt/petasan_ccb/alessandro/SCANB/rna_logtransformed.rds")

gene_mads <- apply(rna_data, 1, mad)
ordered_mads <- order(gene_mads, decreasing = TRUE)

#Select top 3000 most variable genes
top_3000_indices <- ordered_mads[1:3000]
rna_data <- rna_data[top_3000_indices, ]
rna_matrix <- as.matrix(rna_data)

if(any(is.na(rna_matrix))) {
  rna_matrix[is.na(rna_matrix)] <- 0 
}

#NMF on RNAseq data
print("Starting NMF Rank Estimation on RNA data...")
estim.r_rna <- nmfEstimateRank(rna_matrix, 
                           range = 2:6, 
                           nrun = 50,          
                           seed = 123456, 
                           .options = "vp30") 

saveRDS(estim.r_rna, "/mnt/petasan_ccb/alessandro/SCANB/plots/NMF_rank_estimation_RNAseq.rds")
print("Rank Estimation finished and saved!")


print("Starting Final NMF Execution on RNA data...")
res_rna <- nmf(rna_matrix, 
           rank = 3, 
           nrun = 200,                         
           seed = 123456,
           .options = "vp30")         

saveRDS(res_rna, "/mnt/petasan_ccb/alessandro/SCANB/plots/NMF_final_results_RNAseq.rds")
print("RNA NMF Final execution finished!")

print("Saving clustering results...")
sample_groups_rna <- predict(res_rna)
table(sample_groups_rna)

export_groups_rna <- data.frame(
  SampleID = names(sample_groups_rna),
  Cluster = as.vector(sample_groups_rna)
)
write.csv(export_groups_rna, "/mnt/petasan_ccb/alessandro/SCANB/clustering_results/single_omics/rna_nmf_clusters.csv", row.names = FALSE)




# Methylation 

# methylation feature selection
print("Methylation feature selection")

#methylation data
adjusted_data <- readRDS("/mnt/petasan_ccb/alessandro/SCANB/methylation_data/adjusted_data.rds")
# Filtering
library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
#library(sesameData)
ann <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
keep_autosomes <- !(ann$chr %in% c("chrX", "chrY"))
keep_cpg <- grepl("^cg", ann$Name)
probes_to_keep <- ann$Name[keep_autosomes & keep_cpg]
probes_filtered <- adjusted_data[rownames(adjusted_data) %in% probes_to_keep, ]
met_filtered <- as.matrix(probes_filtered)

common_samples <- intersect(colnames(met_filtered), colnames(rna_data))
met_filtered   <- met_filtered[, common_samples]

enhancers_pairs <- readRDS("/mnt/petasan_ccb/alessandro/SCANB/methylation_data/result_pairs_enhancer.rds")
promoters_pairs <- readRDS("/mnt/petasan_ccb/alessandro/SCANB/methylation_data/result_pairs_promoters.rds")

# Enhancers filtering
# Filter for Pairs with Raw.p < 1e-8)
top_pairs_en <- en_pairs[en_pairs$Raw.p < 1e-8 & en_pairs$Pe < 0.001, ]
#select only truly distal pairs
top_pairs_en <- top_pairs_en[abs(top_pairs_en$Distance) > 2000, ]
cat("Top enhancers pairs:", dim(top_pairs_en), "\n")

# Extract unique probes 
top_cpg_en <- unique(top_pairs_en$Probe)
cat("Unique enhancer probes", length(top_cpg_en), "\n")

# Promoters filtering
# Filter for Pairs with Raw.p < 1e-8)
top_pairs_pr <- pr_pairs[pr_pairs$Raw.p < 1e-8 & pr_pairs$Pe < 0.001, ]
#select only truly proximal pairs
top_pairs_en <- top_pairs_en[abs(top_pairs_en$Distance) <= 2000, ]
cat("Top promoters pairs:", dim(top_pairs_pr), "\n")

# Extract unique probes 
top_cpg_pr <- unique(top_pairs_pr$Probe)
cat("Unique promoters probes", length(top_cpg_pr), "\n")

top_cpg_combined <- union(top_cpg_en, top_cpg_pr)
cat("Final CpG amount:", length(top_cpg_combined), "\n")
met_matrix_filtered <- met_filtered[rownames(met_filtered) %in% top_cpg_combined, ]
cat("Final methylation matrix dimesions:", dim(met_matrix_filtered), "\n")


#NMF on methylation data
print("Starting NMF Rank Estimation on meth data...")
estim.r_meth <- nmfEstimateRank(met_matrix_filtered, 
                           range = 2:6, 
                           nrun = 50,          
                           seed = 123456, 
                           .options = "vp30") 

saveRDS(estim.r_meth, "/mnt/petasan_ccb/alessandro/SCANB/plots/NMF_rank_estimation_meth.rds")
print("Rank Estimation finished and saved!")


print("Starting Final NMF Execution on meth data...")
res_meth <- nmf(met_matrix_filtered, 
           rank = 3, 
           nrun = 200,                         
           seed = 123456,
           .options = "vp30")         

saveRDS(res_meth, "/mnt/petasan_ccb/alessandro/SCANB/plots/NMF_final_results_meth.rds")
print("RNA NMF Final execution finished!")

print("Saving clustering results...")
sample_groups_meth <- predict(res_meth)
table(sample_groups_meth)

export_groups_meth <- data.frame(
  SampleID = names(sample_groups_meth),
  Cluster = as.vector(sample_groups_meth)
)
write.csv(export_groups_meth, "/mnt/petasan_ccb/alessandro/SCANB/clustering_results/single_omics/methyl_nmf_clusters.csv", row.names = FALSE)



# Copy number variants

cnv_matrix <- readRDS("/mnt/petasan_ccb/alessandro/SCANB/cnv_processed.rds")

#NMF on cnv data
print("Starting NMF Rank Estimation on cnv data...")
estim.r_cnv <- nmfEstimateRank(cnv_matrix, 
                           range = 2:6, 
                           nrun = 50,          
                           seed = 123456, 
                           .options = "vp30") 

saveRDS(estim.r_cnv, "/mnt/petasan_ccb/alessandro/SCANB/plots/NMF_rank_estimation_cnv.rds")
print("Rank Estimation finished and saved!")


print("Starting Final NMF Execution on cnv data...")
res_cnv <- nmf(cnv_matrix, 
           rank = 3, 
           nrun = 200,                         
           seed = 123456,
           .options = "vp30")         

saveRDS(res_cnv, "/mnt/petasan_ccb/alessandro/SCANB/plots/NMF_final_results_cnv.rds")
print("CNV NMF Final execution finished!")

print("Saving clustering results...")
sample_groups_cnv <- predict(res_cnv)
table(sample_groups_cnv)

export_groups_cnv <- data.frame(
  SampleID = names(sample_groups_cnv),
  Cluster = as.vector(sample_groups_cnv)
)
write.csv(export_groups_cnv, "/mnt/petasan_ccb/alessandro/SCANB/clustering_results/single_omics/cnv_nmf_clusters.csv", row.names = FALSE)

print("All single omics analyses finished!")
