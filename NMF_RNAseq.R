# NMF on 1000 most variable genes by MAD

library(NMF)

n_cores <- 30 #DEVE ESSERE UGUALE A SLURM
options(mc.cores = 30)

rna_data <- read.table("/mnt/petasan_ccb/alessandro/SCANB/rna_logtransformed.tsv", header=TRUE, sep="\t", row.names=1)

gene_mads <- apply(rna_data, 1, mad)
ordered_mads <- order(gene_mads, decreasing = TRUE)

#Select top 1000 most variable genes
top_1000_indices <- ordered_mads[1:1000]
rna_data <- rna_data[top_1000_indices, ]
dt_matrix <- as.matrix(rna_data)

if(any(is.na(dt_matrix))) {
  dt_matrix[is.na(dt_matrix)] <- 0 
}


print("Starting NMF Rank Estimation...")
estim.r <- nmfEstimateRank(dt_matrix, 
                           range = 2:6, 
                           nrun = 50,          
                           seed = 123456, 
                           .options = "v",
                           .pckg = "mclapply") 

saveRDS(estim.r, "/mnt/petasan_ccb/alessandro/SCANB/NMF_rank_estimation_RNAseq.rds")
print("Rank Estimation finished and saved!")


print("Starting Final NMF Execution...")
res <- nmf(dt_matrix, 
           rank = 3, 
           nrun = 200,                         
           seed = 123456,
           .options = "v",
           .pckg = "mclapply")         

saveRDS(res, "/mnt/petasan_ccb/alessandro/SCANB/NMF_final_results_RNAseq.rds")
print("NMF Final execution finished!")