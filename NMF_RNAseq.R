# NMF on feature selected single-omic data

#RNAseq feature selection (top 3000 genes by MAD)
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

library(NMF)

n_cores <- 30 #DEVE ESSERE UGUALE A SLURM
options(mc.cores = 30)

print("Starting NMF Rank Estimation...")
estim.r <- nmfEstimateRank(dt_matrix, 
                           range = 2:6, 
                           nrun = 50,          
                           seed = 123456, 
                           .options = "v") 

saveRDS(estim.r, "/mnt/petasan_ccb/alessandro/SCANB/NMF_rank_estimation_RNAseq.rds")
print("Rank Estimation finished and saved!")


print("Starting Final NMF Execution...")
res <- nmf(dt_matrix, 
           rank = 3, 
           nrun = 200,                         
           seed = 123456,
           .options = "v")         

saveRDS(res, "/mnt/petasan_ccb/alessandro/SCANB/NMF_final_results_RNAseq.rds")
print("NMF Final execution finished!")




# ------

# methylation feature selection

#methylation data

adjusted_data <- readRDS("adjusted_data.rds")
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

common_samples <- intersect(colnames(met_filtered), colnames(exp_matrix))
met_filtered   <- met_filtered[, common_samples]
exp_matrix     <- exp_matrix[, common_samples]

enhancers_pairs <- readRDS("/mnt/petasan_ccb/alessandro/SCANB/methylation_data/result_pairs_enhancer.rds")
promoters_pairs <- readRDS("/mnt/petasan_ccb/alessandro/SCANB/methylation_data/result_pairs_promoter.rds")

# Enhancers filtering
# Filter for Pairs with Raw.p < 1e-8)
top_pairs_en <- en_pairs[en_pairs$Raw.p < 1e-8 & en_pairs$Pe < 0.001, ]
dim(top_pairs_en)

# Extract unique probes 
top_cpg_en <- unique(top_pairs_en$Probe)
length(top_cpg_en)

# Promoters filtering
# Filter for Pairs with Raw.p < 1e-8)
top_pairs_pr <- pr_pairs[pr_pairs$Raw.p < 1e-8 & pr_pairs$Pe < 0.001, ]
dim(top_pairs_pr)

# Extract unique probes 
top_cpg_pr <- unique(top_pairs_pr$Probe)
length(top_cpg_pr)


