# Fits a single iClusterPlus model (fixed K, fixed lambda) on a subsample
# (without replacement) of the full cohort. Intended to be launched once per
# SLURM array task, with SLURM_ARRAY_TASK_ID used both as the iteration ID
# and as the RNG seed for reproducible subsampling.
#
# Usage (called automatically by the sbatch script below):
#   Rscript run_icluster_bootstrap.R <iter_id> <subsample_fraction> <outdir>

## Parse arguments 
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript icluster_bootstrap.R <iter_id> <subsample_fraction> <outdir>")
}
 
iter_id <- as.integer(args[1])   # SLURM_ARRAY_TASK_ID, also used as RNG seed
frac    <- as.numeric(args[2])   # e.g. 0.8 for 80% subsampling
outdir  <- args[3]               # where to save this iteration's result
 
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
 
cat("Iteration:", iter_id, "| Subsample fraction:", frac, "| Outdir:", outdir, "\n")
cat("Start time:", as.character(Sys.time()), "\n")


# Data loading

data_dir <- "/mnt/petasan_ccb/alessandro/SCANB/icluster_inputdata/"  
rna_scaled  <- readRDS(file.path(data_dir, "rna_icluster.rds"))
meth_scaled <- readRDS(file.path(data_dir, "meth_icluster.rds"))
cnv_scaled  <- readRDS(file.path(data_dir, "cnv_icluster.rds"))

n_total <- nrow(rna_scaled)
stopifnot(nrow(meth_scaled) == n_total, nrow(cnv_scaled) == n_total)

## Fixed K / lambda from original tuning 
fixed_k      <- 3
fixed_lambda <- c(0.31081081, 0.14864865, 0.06216216)  
 
cat("Using K =", fixed_k, "| lambda =", paste(round(fixed_lambda, 4), collapse = ", "), "\n")
 
## Subsample without replacement 
set.seed(iter_id)  # reproducible: each array task gets its own, fixed seed
n_sub <- floor(frac * n_total)
sub_idx <- sample(seq_len(n_total), size = n_sub, replace = FALSE)
 
rna_sub  <- rna_scaled[sub_idx, , drop = FALSE]
meth_sub <- meth_scaled[sub_idx, , drop = FALSE]
cnv_sub  <- cnv_scaled[sub_idx, , drop = FALSE]
 
sample_ids_sub <- rownames(rna_scaled)[sub_idx]



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
  sub_idx    = sub_idx,
  sample_ids = sample_ids_sub,
  cluster    = fit$clusters,
  meanZ      = fit$meanZ,
  K          = fixed_k,
  lambda     = fixed_lambda
)
 
saveRDS(result, file = file.path(outdir, paste0("icluster_boot_", iter_id, ".rds")))
 
cat("End time:", as.character(Sys.time()), "\n")
cat("Done.\n")

