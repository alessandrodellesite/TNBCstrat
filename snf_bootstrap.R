# Fits SNF (normalize -> distance -> affinity -> fuse -> spectral cluster) on a
# subsample (without replacement) of the full cohort. Fixed K (number of
# clusters) and fixed SNF hyperparameters (neighbors, sigma, iterations) from
# the original full-cohort run. Subsample composition is read from the SAME
# shared subsample list used by every other benchmarked method, so replicate
# #i uses the exact same set of patients across iCluster, MOFA, SNF, xintNMF.
#
# Usage (called automatically by the sbatch script below):
#   Rscript snf_bootstrap.R <iter_id> <outdir>

## Parse arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript snf_bootstrap.R <iter_id> <outdir>")
}

iter_id <- as.integer(args[1])
outdir  <- args[2]

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cat("Iteration:", iter_id, "| Outdir:", outdir, "\n")
cat("Start time:", as.character(Sys.time()), "\n")

library(SNFtool)

## Data loading of SNF-specific preprocessed matrices 
## standardNormalization() runs per-replicate below, since it must be recomputed on each subsample, not reused from the full cohort
data_dir <- "/mnt/petasan_ccb/alessandro/SCANB/snf_inputdata/"
data_rna  <- readRDS(file.path(data_dir, "rna_snf.rds"))
data_meth <- readRDS(file.path(data_dir, "met_snf.rds"))
data_cnv  <- readRDS(file.path(data_dir, "cnv_snf.rds"))
stopifnot(nrow(data_meth) == nrow(data_rna), nrow(data_cnv) == nrow(data_rna))

## Load shared subsample 
shared_dir <- "/mnt/petasan_ccb/alessandro/SCANB/bootstrap_shared/"
subsample_list <- readRDS(file.path(shared_dir, "boot_subsamples.rds"))

if (!as.character(iter_id) %in% names(subsample_list)) {
  stop(sprintf("iter_id %d not found in boot_subsamples.rds (max = %d)",
               iter_id, length(subsample_list)))
}
sample_ids_sub <- subsample_list[[as.character(iter_id)]]

missing_ids <- setdiff(sample_ids_sub, rownames(data_rna))
if (length(missing_ids) > 0) {
  stop(sprintf("Iteration %d: %d sample IDs from the shared subsample are missing from data_rna (e.g. %s)",
               iter_id, length(missing_ids), paste(head(missing_ids, 3), collapse = ", ")))
}

cat("Subsample size:", length(sample_ids_sub), "\n")

## Subset by sample ID (rows) before  normalization
rna_sub  <- data_rna[sample_ids_sub, , drop = FALSE]
meth_sub <- data_meth[sample_ids_sub, , drop = FALSE]
cnv_sub  <- data_cnv[sample_ids_sub, , drop = FALSE]

## Fixed hyperparameters from original full-cohort run
fixed_k_clusters <- 2      # number of clusters 
K_neighbors <- 20
sigma <- 0.5
T_iter <- 20

cat("Using K_clusters =", fixed_k_clusters,
    "| K_neighbors =", K_neighbors, "| sigma =", sigma, "| T_iter =", T_iter, "\n")

## Standard normalization 
rna_norm  <- standardNormalization(rna_sub)
meth_norm <- standardNormalization(meth_sub)
cnv_norm  <- standardNormalization(cnv_sub)

## Pairwise distances
dist_rna  <- dist2(as.matrix(rna_norm),  as.matrix(rna_norm))
dist_meth <- dist2(as.matrix(meth_norm), as.matrix(meth_norm))
dist_cnv  <- dist2(as.matrix(cnv_norm),  as.matrix(cnv_norm))

## Affinity matrices
W_rna  <- affinityMatrix(dist_rna,  K_neighbors, sigma)
W_meth <- affinityMatrix(dist_meth, K_neighbors, sigma)
W_cnv  <- affinityMatrix(dist_cnv,  K_neighbors, sigma)

## Fuse networks
W_fused <- SNF(list(W_rna, W_meth, W_cnv), K_neighbors, T_iter)

## Spectral clustering at fixed K
clusters <- spectralClustering(W_fused, fixed_k_clusters)
names(clusters) <- sample_ids_sub

## Save outputs
result <- list(
  iter_id     = iter_id,
  sample_ids  = sample_ids_sub,
  cluster     = clusters,
  K           = fixed_k_clusters,
  K_neighbors = K_neighbors,
  sigma       = sigma,
  T_iter      = T_iter
)

saveRDS(result, file = file.path(outdir, paste0("snf_boot_", iter_id, ".rds")))

cat("End time:", as.character(Sys.time()), "\n")
cat("Done.\n")
