# produces boot_subsamples.rds, referenced identically by every method's bootstrap script so that 
# replicate #i always means the same 80% of patients across all methods.

out_dir <- "/mnt/petasan_ccb/alessandro/SCANB/bootstrap_shared/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ground-truth sample IDs 
rna_scaled <- readRDS("/mnt/petasan_ccb/alessandro/SCANB/icluster_inputdata/rna_icluster.rds")
common_samples <- rownames(rna_scaled)

cat("Master sample set size:", length(common_samples), "\n")
stopifnot(!any(duplicated(common_samples)))  

n_reps <- 1000
frac   <- 0.8
n_sub  <- floor(frac * length(common_samples))

subsample_list <- vector("list", n_reps)

for (i in seq_len(n_reps)) {
  set.seed(i)   # iter_id doubles as seed -- single canonical source of randomness
  idx <- sample(seq_along(common_samples), size = n_sub, replace = FALSE) #draws size values from 235 (seq_along(common_samples)) without replacement
  subsample_list[[i]] <- common_samples[idx]
}
names(subsample_list) <- as.character(seq_len(n_reps))

saveRDS(subsample_list, file.path(out_dir, "boot_subsamples.rds"))
saveRDS(common_samples, file.path(out_dir, "boot_master_samples.rds"))

manifest <- list(
  n_reps = n_reps,
  frac = frac,
  n_sub_per_rep = n_sub,
  n_total_samples = length(common_samples),
  seed_scheme = "set.seed(i) for i in 1:n_reps, sample() without replacement",
  master_samples_source = "rna_icluster.rds rownames",
  generated_on = as.character(Sys.time())
)
saveRDS(manifest, file.path(out_dir, "boot_manifest.rds"))

cat("Saved", n_reps, "subsample definitions to", out_dir, "\n")
