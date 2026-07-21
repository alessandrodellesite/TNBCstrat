library(metasnf)
library(ComplexHeatmap)
library(SNFtool)
library(future)

# data loading and preprocessing

## RNAseq

rna_data <- read.table("/mnt/petasan_ccb/alessandro/SCANB/rna_logtransformed.tsv", header=TRUE, sep="\t", row.names=1)
# Feature selection
gene_mads <- apply(rna_data, 1, mad)
ordered_mads <- order(gene_mads, decreasing = TRUE)
top_3000_indices <- ordered_mads[1:3000]
rna_data <- rna_data[top_3000_indices, ]
rna_matrix <- as.matrix(rna_data)

## Mehtylation

adjusted_data <- readRDS("/mnt/petasan_ccb/alessandro/SCANB/adjusted_data.rds")

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

# Enhancers results
enhancers_pairs <- readRDS("/mnt/petasan_ccb/alessandro/SCANB/result_pairs_enhancer.rds")
en_pairs <- enhancers_pairs[order(enhancers_pairs$Pe), ]
top_pairs_en <- en_pairs[en_pairs$Raw.p < 1e-9, ]
top_cpg_en <- unique(top_pairs_en$Probe)

# Promoters results
promoters_pairs <- readRDS("/mnt/petasan_ccb/alessandro/SCANB/result_pairs_promoter.rds")
pr_pairs <- promoters_pairs[order(promoters_pairs$Pe), ]
top_pairs_pr <- pr_pairs[pr_pairs$Raw.p < 1e-9, ]
top_cpg_pr <- unique(top_pairs_pr$Probe)

# Final data
top_cpg_combined <- union(top_cpg_en, top_cpg_pr)
met_matrix_filtered <- met_filtered[rownames(met_filtered) %in% top_cpg_combined, ]

met_mvals <- met_matrix_filtered
# to prevent ±Infinity errors
met_mvals[met_mvals == 0] <- 0.001
met_mvals[met_mvals == 1] <- 0.999
# logit transformation to get M-values
met_mvals <- log2(met_mvals / (1 - met_mvals))


## CNV

dna_matrix_filtered <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/dna_matrix_filtered_2900genes.csv", header = TRUE, row.names = 1)
cnv_log_ratio <- dna_matrix_filtered
# homozygous deletions (0 alleles) to avoid log2(0) = -Inf
cnv_log_ratio[cnv_log_ratio == 0] <- 0.5  
#log2 scale
cnv_log_ratio <- log2(cnv_log_ratio / 2)


#Align samples
common_samples <- intersect(colnames(rna_matrix), intersect(colnames(met_mvals), colnames(cnv_log_ratio)))
rna_matrix <- rna_matrix[, common_samples]
met_mvals <- met_mvals[, common_samples]
cnv_log_ratio <- cnv_log_ratio[, common_samples]

#Transpose all matrices 
data_rna <- t(rna_matrix)      
data_meth <- t(met_mvals)            
data_cnv <- t(cnv_log_ratio)

# Standard normalization (makes each feature (column) have a mean of 0 and SD of 1)

data_rna_norm  <- standardNormalization(data_rna)
data_meth_norm <- standardNormalization(data_meth)
data_cnv_norm  <- standardNormalization(data_cnv)

#Convert matrices into data frames
df_rna  <- as.data.frame(data_rna_norm)     # 235x3000
df_meth <- as.data.frame(data_meth_norm)    # 235x3527
df_cnv  <- as.data.frame(data_cnv_norm)     # 235x2983

# Add patient_id column data frames using the row names
df_rna$patient_id  <- rownames(df_rna)
df_meth$patient_id <- rownames(df_meth)
df_cnv$patient_id  <- rownames(df_cnv)


# Construct the formal metasnf data list
# This automatically aligns sample IDs and drops missing data
my_data_list <- data_list(
  list(df_rna,  "rna_expression",    "transcriptomics", "continuous"),
  list(df_meth, "dna_methylation",    "epigenomics",     "continuous"),
  list(df_cnv,  "copy_number_variants", "genomics",        "continuous"),
  uid = "patient_id"
)



# SLURM Parallelization
# Automatically grabs the CPUS assigned by SLURM (--cpus-per-task)
n_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 4))
plan(multisession, workers = n_cores)
cat(sprintf("Running metasnf using %d cores...\n", n_cores))

#Run batch_snf

solutions_df <- batch_snf(
  dl = my_data_list, 
  sc = my_config,
  return_sim_mats = TRUE  # keeps similarity matrices for quality metrics
)

# Save Output Safely
# -------------------------------------------------------------------
output_path <- "/mnt/petasan_ccb/alessandro/SCANB/sulutions_metasnf.rds"

# Ensure output directory exists before saving
dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)

saveRDS(solutions_df, file = output_path)
cat(sprintf("Success! Results saved to %s\n", output_path))

