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

library(MOFA2)
library(tidyverse)

## Data loading of iCluster-specific preprocessed matrices
data_dir <- "/mnt/petasan_ccb/alessandro/SCANB/mofa_inputdata/"
rna_scaled  <- readRDS(file.path(data_dir, "rna_mofa.rds"))
meth_scaled <- readRDS(file.path(data_dir, "met_mofa.rds"))
cnv_scaled  <- readRDS(file.path(data_dir, "cnv_mofa.rds"))
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

## Subset by sample ID 
rna_sub  <- rna_scaled[sample_ids_sub, , drop = FALSE]
meth_sub <- meth_scaled[sample_ids_sub, , drop = FALSE]
cnv_sub  <- cnv_scaled[sample_ids_sub, , drop = FALSE]



multiomics_data <- list(
  RNAseq      = rna_cent,        
  Methylation = met_cent,        
  CNV         = cnv_cent         
)

MOFAobject <- create_mofa(multiomics_data)
data_opts <- get_default_data_options(MOFAobject)
model_opts <- get_default_model_options(MOFAobject)
data_opts$scale_views <- TRUE
model_opts$num_factors <- 20
model_opts$likelihoods <- c(
  RNAseq      = "gaussian",
  Methylation = "gaussian",
  CNV         = "gaussian"
)
train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "slow"  
train_opts$seed <- 123 #set seed for reproducibility

MOFAobject <- prepare_mofa(
  MOFAobject,
  data_options     = data_opts,
  model_options    = model_opts,
  training_options = train_opts
)
