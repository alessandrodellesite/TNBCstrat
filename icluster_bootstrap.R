# Fits a single iClusterPlus model (fixed K, fixed lambda) on a subsample
# (without replacement) of the full cohort. Intended to be launched once per
# SLURM array task. Subsample composition is read from a shared, pre-generated
# subsample list (boot_subsamples.rds) so that replicate #i uses the exact
# same set of patients across every benchmarked method, not just iCluster.
#
# Usage (called automatically by the sbatch script below):
#   Rscript icluster_bootstrap.R <iter_id> <outdir>

## Parse arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript icluster_bootstrap.R <iter_id> <outdir>")
}

iter_id <- as.integer(args[1])   # SLURM_ARRAY_TASK_ID -- indexes into the shared subsample list
outdir  <- args[2]

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cat("Iteration:", iter_id, "| Outdir:", outdir, "\n")
cat("Start time:", as.character(Sys.time()), "\n")

library(iClusterPlus)

## Data loading of iCluster-specific preprocessed matrices
data_dir <- "/mnt/petasan_ccb/alessandro/SCANB/icluster_inputdata/"
rna_scaled  <- readRDS(file.path(data_dir, "rna_icluster.rds"))
meth_scaled <- readRDS(file.path(data_dir, "met_icluster.rds"))
cnv_scaled  <- readRDS(file.path(data_dir, "cnv_icluster.rds"))
stopifnot(nrow(meth_scaled) == nrow(rna_scaled), nrow(cnv_scaled) == nrow(rna_scaled))

## Load shared subsample definition 
shared_dir <- "/mnt/petasan_ccb/alessandro/SCANB/bootstrap_shared/"
subsample_list <- readRDS(file.path(shared_dir, "boot_subsamples.rds"))

if (!as.character(iter_id) %in% names(subsample_list)) {
  stop(sprintf("iter_id %d not found in boot_subsamples.rds (max = %d)",
               iter_id, length(subsample_list)))
}
sample_ids_sub <- subsample_list[[as.character(iter_id)]]

# Confirm this method's data actually contains every sampled ID -- fail loudly
# rather than silently subsetting to fewer samples than intended
missing_ids <- setdiff(sample_ids_sub, rownames(rna_scaled))
if (length(missing_ids) > 0) {
  stop(sprintf("Iteration %d: %d sample IDs from the shared subsample are missing from rna_scaled (e.g. %s)",
               iter_id, length(missing_ids), paste(head(missing_ids, 3), collapse = ", ")))
}

cat("Subsample size:", length(sample_ids_sub), "\n")

## Fixed K / lambda from original tuning
fixed_lambda_2_cl <- c(0.11621622, 0.08918919, 0.60810811) # output[[1]]$lambda[minBICid[1],] -> for k=1 (2 clusters)
fixed_lambda_3_cl <- c(0.31081081, 0.14864865, 0.06216216) # output[[2]]$lambda[minBICid[2],] -> for k=2 (3 clusters)
fixed_lambda_4_cl <- c(0.15405405, 0.07297297, 0.02972973) # output[[3]]$lambda[minBICid[3],] -> for k=3 (4 clusters)

fixed_k      <- 1
fixed_lambda <- fixed_lambda_2_cl
cat("Using K =", fixed_k, "| lambda =", paste(round(fixed_lambda, 4), collapse = ", "), "\n")

## Subset by sample ID 
rna_sub  <- rna_scaled[sample_ids_sub, , drop = FALSE]
meth_sub <- meth_scaled[sample_ids_sub, , drop = FALSE]
cnv_sub  <- cnv_scaled[sample_ids_sub, , drop = FALSE]

fit <- iClusterPlus(
  dt1 = rna_sub,
  dt2 = meth_sub,
  dt3 = cnv_sub,
  type = c("gaussian", "gaussian", "gaussian"),
  K = fixed_k,
  lambda = fixed_lambda,
  maxiter = 20
)

## Save outputs
result <- list(
  iter_id    = iter_id,
  sample_ids = sample_ids_sub,
  cluster    = fit$clusters,
  meanZ      = fit$meanZ,
  K          = fixed_k,
  lambda     = fixed_lambda
)

saveRDS(result, file = file.path(outdir, paste0("icluster_boot_", iter_id, ".rds")))

cat("End time:", as.character(Sys.time()), "\n")
cat("Done.\n")
