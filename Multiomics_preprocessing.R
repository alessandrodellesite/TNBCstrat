# Preprocessing of each omic's dataset specific for each multiomics model
# and loading on ondemand directory


# RNAseq

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


# CNV



# MOFA specific preprocessing -> Mean-center the rows (genes) so that each gene's average across samples is 0 (scale = FALSE to not scale the variance)

# icluster specific

# SNF specific preprocessing



