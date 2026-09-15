# script to run kmeans on the H matrices resulted from xintNMF on the bootstrapped samples

#!/usr/bin/env Rscript
## ============================================================
## Convert raw X-intNMF bootstrap outputs (sample_factor.csv per
## iter_N folder) into per-iteration .rds files with the same
## structure as the icluster/mofa/snf bootstrap files:
##   iter_id, sample_ids, cluster, K
##
## This performs the k-means step on H that was previously done
## in the "Final export" chunk of the original rank-selection
## script -- but now applied per-bootstrap-iteration, at a FIXED
## k (no re-selection of k per replicate).
## ============================================================

library(purrr)

convert_xintnmf_bootstrap <- function(k,
                                       n_iter = 500,
                                       raw_base = "/mnt/petasan_ccb/alessandro/SCANB/multiomics/bootstrap/xintnmf",
                                       subsamples_path = "/mnt/petasan_ccb/alessandro/SCANB/bootstrap_shared/boot_subsamples.rds",
                                       out_dir = NULL,
                                       seed = 123,
                                       nstart = 50) {

  raw_dir <- file.path(raw_base, sprintf("xintnmf_bootstrap_k%d", k))
  if (is.null(out_dir)) {
    out_dir <- file.path(raw_base, sprintf("xintnmf_bootstrap_k%d_clusters", k))
  }
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  subs <- readRDS(subsamples_path)

  set.seed(seed)  # reproducible k-means across the whole conversion run

  n_ok <- 0
  n_fail <- 0

  for (i in seq_len(n_iter)) {
    iter_folder <- file.path(raw_dir, sprintf("iter_%d", i))
    h_path <- file.path(iter_folder, "sample_factor.csv")

    if (!file.exists(h_path)) {
      warning(sprintf("k=%d iter=%d: sample_factor.csv missing, skipping", k, i))
      n_fail <- n_fail + 1
      next
    }

    H <- as.matrix(read.csv(h_path, header = TRUE))

    sample_ids_sub <- subs[[as.character(i)]]
    if (is.null(sample_ids_sub)) {
      warning(sprintf("k=%d iter=%d: no entry in boot_subsamples.rds for this iter_id, skipping", k, i))
      n_fail <- n_fail + 1
      next
    }

    ## Integrity check before trusting row order 
    if (nrow(H) != length(sample_ids_sub)) {
      warning(sprintf(
        "k=%d iter=%d: row count mismatch (H has %d rows, subsample has %d ids) -- skipping",
        k, i, nrow(H), length(sample_ids_sub)
      ))
      n_fail <- n_fail + 1
      next
    }

    km <- kmeans(H, centers = k, nstart = nstart)
    cluster <- km$cluster
    names(cluster) <- sample_ids_sub

    out <- list(
      iter_id     = i,
      sample_ids  = sample_ids_sub,
      cluster     = as.integer(cluster),
      K           = k
    )

    saveRDS(out, file.path(out_dir, sprintf("xintnmf_boot_%d.rds", i)))
    n_ok <- n_ok + 1
  }

  message(sprintf("k=%d: converted %d/%d iterations (%d failed) -> %s",
                   k, n_ok, n_iter, n_fail, out_dir))
  invisible(list(n_ok = n_ok, n_fail = n_fail, out_dir = out_dir))
}

## Run for all three k's 
for (k in c(2, 3, 4)) {
  convert_xintnmf_bootstrap(k = k)
}


