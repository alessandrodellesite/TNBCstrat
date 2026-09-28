K <- 2
outdir <- sprintf("/mnt/petasan_ccb/alessandro/SCANB/multiomics/bootstrap/xintnmf/xintnmf_bootstrap_k%d", K)
n_iter <- 1000

marker <- file.path(outdir, paste0("iter_", 1:n_iter), "sample_factor.csv")
ok <- file.exists(marker) & file.size(marker) > 0
ok[is.na(ok)] <- FALSE   # file.size returns NA for missing files

cat("Done:", sum(ok), "/", n_iter, "\n")
missing <- which(!ok)
cat("Missing:", length(missing), "\n")
cat(paste(missing, collapse = ","), "\n")
