# Preprocessing of each omic's dataset specific for each multiomics model
# and loading on ondemand directory


# RNAseq

rna_data <- read.table("rna_logtransformed.tsv", header=TRUE, sep="\t", row.names=1)
# Feature selection
gene_mads <- apply(rna_data, 1, mad)
ordered_mads <- order(gene_mads, decreasing = TRUE)
top_3000_indices <- ordered_mads[1:3000]
rna_data <- rna_data[top_3000_indices, ]
rna_matrix <- as.matrix(rna_data)

# MOFA specific preprocessing -> Mean-center the rows (genes) so that each gene's average across samples is 0 (scale = FALSE to not scale the variance)

# icluster specific

# SNF specific preprocessing



