ks     <- 2:4
n_iter <- 500
base   <- "/mnt/petasan_ccb/alessandro/SCANB/multiomics/bootstrap/xintnmf"

missing_by_k <- list()

for (k in ks) {
  outdir <- file.path(base, sprintf("xintnmf_bootstrap_k%d", k))
  f  <- file.path(outdir, paste0("iter_", 1:n_iter), "sample_factor.csv")
  ok <- file.exists(f) & !is.na(file.size(f)) & file.size(f) > 0

  missing_by_k[[as.character(k)]] <- which(!ok)

  cat(sprintf("K=%d: %d/%d completate, %d mancanti\n",
              k, sum(ok), n_iter, sum(!ok)))

  if (any(!ok)) {
    cat(sprintf(
      "  sbatch --job-name=xintnmf_boot_k%d --export=ALL,K=%d --array=%s%%20 xintnmf_boot.sh\n",
      k, k, paste(which(!ok), collapse = ",")))
  }
}
