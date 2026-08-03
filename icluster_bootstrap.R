# Fits a single iClusterPlus model (fixed K, fixed lambda) on a subsample
# (without replacement) of the full cohort. Intended to be launched once per
# SLURM array task, with SLURM_ARRAY_TASK_ID used both as the iteration ID
# and as the RNG seed for reproducible subsampling.
#
# Usage (called automatically by the sbatch script below):
#   Rscript run_icluster_bootstrap.R <iter_id> <subsample_fraction> <outdir>

suppressMessages(library(iClusterPlus))

## ---- 1. Parse arguments -----------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript run_icluster_bootstrap.R <iter_id> <subsample_fraction> <outdir>")
}

iter_id     <- as.integer(args[1])   # e.g. SLURM_ARRAY_TASK_ID, used as seed too
frac        <- as.numeric(args[2])   # e.g. 0.8 for 80% subsampling
outdir      <- args[3]

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cat("Iteration:", iter_id, "| Subsample fraction:", frac, "| Outdir:", outdir, "\n")
cat("Start time:", as.character(Sys.time()), "\n")

## ---- 2. Load full data (already scaled) --------------------------------
# Replace these paths with wherever your scaled matrices live.
# Rows = samples, columns = features, matching your original tuning script.
rna_scaled  <- readRDS("rna_scaled.rds")
meth_scaled <- readRDS("meth_scaled.rds")
cnv_scaled  <- readRDS("cnv_scaled.rds")

n_total <- nrow(rna_scaled)
stopifnot(nrow(meth_scaled) == n_total, nrow(cnv_scaled) == n_total)

## ---- 3. Load the fixed K / lambda chosen from tuning --------------------
# From your tuning step: K = 3 (i.e. output[[2]], 3-cluster solution),
# lambda = the row of output[[2]]$lambda that minimized BIC (minBICid[2]).
# Save these once from your tuning session, e.g.:
#   saveRDS(list(K = 3, lambda = output[[2]]$lambda[minBICid[2], ]),
#           file = "chosen_model_params.rds")
model_params <- readRDS("chosen_model_params.rds")
K      <- model_params$K
lambda <- model_params$lambda

cat("Using K =", K, "| lambda =", paste(round(lambda, 4), collapse = ", "), "\n")

## ---- 4. Subsample without replacement -----------------------------------
set.seed(iter_id)  # reproducible: each array task gets its own, fixed seed
n_sub <- floor(frac * n_total)
sub_idx <- sample(seq_len(n_total), size = n_sub, replace = FALSE)

rna_sub  <- rna_scaled[sub_idx, , drop = FALSE]
meth_sub <- meth_scaled[sub_idx, , drop = FALSE]
cnv_sub  <- cnv_scaled[sub_idx, , drop = FALSE]

sample_ids_sub <- rownames(rna_scaled)[sub_idx]

## ---- 5. Fit iClusterPlus at fixed K / lambda ----------------------------
fit <- iClusterPlus(
  dt1 = rna_sub,
  dt2 = meth_sub,
  dt3 = cnv_sub,
  type = c("gaussian", "gaussian", "gaussian"),
  K = K,
  lambda = lambda,
  maxiter = 20
)

## ---- 6. Save outputs ------------------------------------------------------
result <- list(
  iter_id        = iter_id,
  sub_idx        = sub_idx,
  sample_ids     = sample_ids_sub,
  cluster        = fit$clusters,
  meanZ          = fit$meanZ,
  K              = K,
  lambda         = lambda
)

saveRDS(result, file = file.path(outdir, paste0("icluster_boot_", iter_id, ".rds")))

cat("End time:", as.character(Sys.time()), "\n")
cat("Done.\n")
