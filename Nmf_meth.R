# NMF on feature selected single-omic data

Sys.unsetenv("_R_CHECK_TIMINGS_")
Sys.unsetenv("_R_CHECK_CRAN_INCOMING_")
library(NMF)


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

rna_data <- readRDS("/mnt/petasan_ccb/alessandro/SCANB/rna_logtransformed.rds")
common_samples <- intersect(colnames(met_filtered), colnames(rna_data))
met_filtered   <- met_filtered[, common_samples]

en_pairs <- readRDS("/mnt/petasan_ccb/alessandro/SCANB/methylation_data/result_pairs_enhancer.rds")
pr_pairs <- readRDS("/mnt/petasan_ccb/alessandro/SCANB/methylation_data/result_pairs_promoters.rds")

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
top_pairs_pr <- top_pairs_pr[abs(top_pairs_pr$Distance) <= 2000, ]
cat("Top promoters pairs:", dim(top_pairs_pr), "\n")

# Extract unique probes 
top_cpg_pr <- unique(top_pairs_pr$Probe)
cat("Unique promoters probes", length(top_cpg_pr), "\n")

top_cpg_combined <- union(top_cpg_en, top_cpg_pr)
cat("Final CpG amount:", length(top_cpg_combined), "\n")
met_matrix_filtered <- met_filtered[rownames(met_filtered) %in% top_cpg_combined, ]
cat("Final methylation matrix dimesions:", dim(met_matrix_filtered), "\n")


print("Starting Final NMF Execution on meth data with rank=4...")
res_meth <- nmf(met_matrix_filtered, 
           rank = 4, 
           nrun = 200,                         
           seed = 123456,
           .options = "vp30")         

saveRDS(res_meth, "/mnt/petasan_ccb/alessandro/SCANB/singleomic/nmf_meth/NMF_final_results_meth_4_clusters.rds")

print("Saving clustering results...")
sample_groups_meth <- predict(res_meth)
table(sample_groups_meth)

export_groups_meth <- data.frame(
  SampleID = names(sample_groups_meth),
  Cluster = as.vector(sample_groups_meth)
)
write.csv(export_groups_meth, "/mnt/petasan_ccb/alessandro/SCANB/singleomic/nmf_meth/methyl_nmf_clusters_4_clusters.csv", row.names = FALSE)

