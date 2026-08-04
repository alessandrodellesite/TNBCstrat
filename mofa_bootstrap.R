# Fits MOFA (train -> extract factors -> k-means at fixed K) on a subsample
# (without replacement) of the full cohort. Fixed number of MOFA factors used
# for clustering (matching the original full-cohort choice) and fixed K
# (cluster count) from the original analysis. Subsample composition is read
# from the SAME shared subsample list used by every other benchmarked method.
#
# Usage:
#   Rscript mofa_bootstrap.R <iter_id> <outdir>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript mofa_bootstrap.R <iter_id> <outdir>")
}

iter_id <- as.integer(args[1])
outdir  <- args[2]

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cat("Iteration:", iter_id, "| Outdir:", outdir, "\n")
cat("Start time:", as.character(Sys.time()), "\n")

library(MOFA2)

## Data loading -- MOFA-specific preprocessed matrices (features x samples,
## matching create_mofa()'s expected orientation)
data_dir <- "/mnt/petasan_ccb/alessandro/SCANB/mofa_inputdata/"
rna_cent <- readRDS(file.path(data_dir, "rna_mofa.rds"))
met_cent <- readRDS(file.path(data_dir, "met_mofa.rds"))
cnv_cent <- readRDS(file.path(data_dir, "cnv_mofa.rds"))

## Load SHARED subsample definition
shared_dir <- "/mnt/petasan_ccb/alessandro/SCANB/bootstrap_shared/"
subsample_list <- readRDS(file.path(shared_dir, "boot_subsamples.rds"))

if (!as.character(iter_id) %in% names(subsample_list)) {
  stop(sprintf("iter_id %d not found in boot_subsamples.rds (max = %d)",
               iter_id, length(subsample_list)))
}
sample_ids_sub <- subsample_list[[as.character(iter_id)]]

# note: rna_cent etc. are features x samples here, so check COLUMN names, not row names
missing_ids <- setdiff(sample_ids_sub, colnames(rna_cent))
if (length(missing_ids) > 0) {
  stop(sprintf("Iteration %d: %d sample IDs from the shared subsample are missing from rna_cent (e.g. %s)",
               iter_id, length(missing_ids), paste(head(missing_ids, 3), collapse = ", ")))
}

cat("Subsample size:", length(sample_ids_sub), "\n")

## Subset by sample ID (columns, since orientation is features x samples)
rna_sub <- rna_cent[, sample_ids_sub, drop = FALSE]
met_sub <- met_cent[, sample_ids_sub, drop = FALSE]
cnv_sub <- cnv_cent[, sample_ids_sub, drop = FALSE]

multiomics_data_sub <- list(
  RNAseq      = rna_sub,
  Methylation = met_sub,
  CNV         = cnv_sub
)

## Fixed parameters from original full-cohort tuning
fixed_num_factors   <- 20    # upper limit for MOFA's factor search, matches original
fixed_factors_to_use <- 6    # how many top factors feed into k-means -- from your >2.5% variance threshold
fixed_k_clusters     <- 3    # cluster count being benchmarked

cat("num_factors =", fixed_num_factors,
    "| factors_to_use =", fixed_factors_to_use,
    "| K_clusters =", fixed_k_clusters, "\n")

## Build and train MOFA object on the subsample
MOFAobject <- create_mofa(multiomics_data_sub)

data_opts  <- get_default_data_options(MOFAobject)
data_opts$scale_views <- TRUE

model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- fixed_num_factors
model_opts$likelihoods <- c(RNAseq = "gaussian", Methylation = "gaussian", CNV = "gaussian")

train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "slow"   
train_opts$seed <- iter_id              # vary by replicate for genuine MOFA-internal randomness, NOT sample draw
train_opts$verbose <- FALSE

MOFAobject <- prepare_mofa(MOFAobject, data_options = data_opts,
                            model_options = model_opts, training_options = train_opts)

mofa_outfile <- tempfile(fileext = ".hdf5")   # per-replicate scratch file, not kept
MOFAobject <- run_mofa(MOFAobject, outfile = mofa_outfile, use_basilisk = TRUE)

## Extract factors, cluster with k-means at fixed K
factors_matrix <- do.call(rbind, get_factors(MOFAobject, factors = "all"))
n_available_factors <- ncol(factors_matrix)

if (n_available_factors < fixed_factors_to_use) {
  warning(sprintf("Iteration %d: only %d factors available (wanted %d) -- using all available",
                   iter_id, n_available_factors, fixed_factors_to_use))
  factors_to_use <- factors_matrix
} else {
  factors_to_use <- factors_matrix[, 1:fixed_factors_to_use, drop = FALSE]
}

set.seed(iter_id)
km <- kmeans(factors_to_use, centers = fixed_k_clusters, nstart = 50)

clusters <- km$cluster
names(clusters) <- rownames(factors_to_use)

## Save outputs
result <- list(
  iter_id       = iter_id,
  sample_ids    = sample_ids_sub,
  cluster       = clusters,
  K             = fixed_k_clusters,
  factors_used  = fixed_factors_to_use,
  n_factors_available = n_available_factors
)

saveRDS(result, file = file.path(outdir, paste0("mofa_boot_", iter_id, ".rds")))

file.remove(mofa_outfile)   # clean up scratch HDF5 file

cat("End time:", as.character(Sys.time()), "\n")
cat("Done.\n")
