outdir <- "/mnt/petasan_ccb/alessandro/SCANB/multiomics/bootstrap/mofa/mofa_bootstrap_k2"  # adjust to your path
n_iter <- 1000

f  <- file.path(outdir, paste0("mofa_boot_", 1:n_iter, ".rds"))
ok <- file.exists(f) & !is.na(file.size(f)) & file.size(f) > 0

cat(sprintf("Done: %d/%d | Missing: %d\n", sum(ok), n_iter, sum(!ok)))
missing <- which(!ok)
cat(paste(missing, collapse = ","), "\n")

subs <- readRDS("/mnt/petasan_ccb/alessandro/SCANB/bootstrap_shared/boot_subsamples.rds")

bad <- integer()
for (i in which(ok)) {
  r <- tryCatch(readRDS(f[i]), error = function(e) NULL)
  if (is.null(r) ||
      r$iter_id != i ||
      !setequal(r$sample_ids, subs[[as.character(i)]]) ||
      length(r$cluster) != length(subs[[as.character(i)]]) ||
      anyNA(r$cluster) ||
      length(unique(r$cluster)) != r$K) {
    bad <- c(bad, i)
  }
}
cat("Corrupt/inconsistent:", length(bad), "\n")
print(bad)

to_rerun <- sort(union(missing, bad))
cat(sprintf("sbatch --array=%s%%25 xintnmf_boot.sh\n", paste(to_rerun, collapse = ",")))
