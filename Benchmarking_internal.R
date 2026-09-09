library(cluster)
library(fpc)
library(clValid)
library(clusterSim)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)


## Internal validation



# function to validate a match() result instead of silently dropping samples
check_matches <- function(idx, sample_ids, source_name) {
  if (any(is.na(idx))) {
    missing <- sample_ids[is.na(idx)]
    stop(sprintf(
      "%s: %d SampleIDs in cluster_results have no match (e.g. %s)",
      source_name, length(missing), paste(head(missing, 3), collapse = ", ")
    ))
  }
  invisible(TRUE)
}

# Align each data source once

# MOFA: must use the SAME factor subset that was actually clustered on (factors_to_use <- factors_matrix[,1:6])
mofa_factors <- readRDS("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_mofa/mofa_factors.rds")
mofa_idx <- match(cluster_results$SampleID, rownames(mofa_factors))
check_matches(mofa_idx, cluster_results$SampleID, "MOFA")
mofa_factors_clean <- mofa_factors[mofa_idx, , drop = FALSE]

N_MOFA_FACTORS_USED <- 6  # must match factors_to_use in the MOFA k-means script
stopifnot(
  "mofa_factors.rds has fewer columns than N_MOFA_FACTORS_USED -- check the file" =
    ncol(mofa_factors_clean) >= N_MOFA_FACTORS_USED
)
if (ncol(mofa_factors_clean) != N_MOFA_FACTORS_USED) {
  warning(sprintf(
    "mofa_factors.rds has %d columns; subsetting to the first %d to match what k-means was run on.",
    ncol(mofa_factors_clean), N_MOFA_FACTORS_USED
  ))
}
mofa_factors_clean <- mofa_factors_clean[, 1:N_MOFA_FACTORS_USED, drop = FALSE]
mofa_dist_full <- dist(mofa_factors_clean, method = "euclidean")


# SNF

snf_w_fused <- readRDS("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_snf/snf_w_fused.rds")
snf_idx <- match(cluster_results$SampleID, rownames(snf_w_fused))
check_matches(snf_idx, cluster_results$SampleID, "SNF")
snf_w_clean <- snf_w_fused[snf_idx, snf_idx]
snf_dist_full <- as.dist(1 - snf_w_clean)

# SNF has no native feature-space embedding (only an affinity/distance matrix),
# but index.DB requires actual coordinates to compute centroid/medoid positions.
# Classical MDS gives a coordinate embedding that preserves the SNF distances.
snf_embedding <- cmdscale(snf_dist_full, k = min(10, nrow(snf_w_clean) - 2))

# iCluster 2/3/4
icluster_z2 <- readRDS("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_icluster/icluster_matrix_2.rds")
icluster_idx2 <- match(cluster_results$SampleID, rownames(icluster_z2))
check_matches(icluster_idx2, cluster_results$SampleID, "iCluster k=2")
icluster_z_clean2 <- icluster_z2[icluster_idx2, , drop = FALSE]
icluster_dist_full2 <- dist(icluster_z_clean2, method = "euclidean")

icluster_z3 <- readRDS("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_icluster/icluster_matrix_3.rds")
icluster_idx3 <- match(cluster_results$SampleID, rownames(icluster_z3))
check_matches(icluster_idx3, cluster_results$SampleID, "iCluster k=3")
icluster_z_clean3 <- icluster_z3[icluster_idx3, , drop = FALSE]
icluster_dist_full3 <- dist(icluster_z_clean3, method = "euclidean")

icluster_z4 <- readRDS("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_icluster/icluster_matrix_4.rds")
icluster_idx4 <- match(cluster_results$SampleID, rownames(icluster_z4))
check_matches(icluster_idx4, cluster_results$SampleID, "iCluster k=4")
icluster_z_clean4 <- icluster_z4[icluster_idx4, , drop = FALSE]
icluster_dist_full4 <- dist(icluster_z_clean4, method = "euclidean")


# xintNMF

# Recover patient IDs in the order X-intNMF was actually run on
input_dir <- "/mnt/petasan_ccb/alessandro/SCANB/multiomics/input_data/xintnmf_inputdata/"
rna_header  <- colnames(read.table(paste0(input_dir, "rna.tsv"), sep = "\t", header = TRUE, nrows = 1))
patient_ids <- rna_header[-1]  # drop the gene/feature ID column
cat("Number of patient IDs:", length(patient_ids), "\n")
cat("Any duplicates:", any(duplicated(patient_ids)), "\n")


xint_factors <- list()
base_dir <- "/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_xintnmf/"
for (k in c(2, 3, 4)) {
  file_path <- paste0(base_dir, "rankselect_k", k, "_graphreg/sample_factor.csv")
  mat <- read.csv(file_path, header = TRUE)
  stopifnot(
    "sample_factor.csv row count doesn't match patient_ids length" =
      nrow(mat) == length(patient_ids),
    "patient_ids has duplicates - can't be used as unique rownames" =
      !any(duplicated(patient_ids)),
    "sample_factor.csv contains non-numeric columns" =
      all(sapply(mat, is.numeric))
  )
  rownames(mat) <- patient_ids
  match_idx <- match(cluster_results$SampleID, rownames(mat))
  check_matches(match_idx, cluster_results$SampleID, sprintf("xintNMF k=%d", k))
  xint_factors[[as.character(k)]] <- mat[match_idx, , drop = FALSE]
}

# IMPORTANT: use the SAME distance metric that actually produced the xintNMF clusters
# (corr_hclust_cluster used 1 - abs(cor(t(H))), not Euclidean -- keep consistent here)
# QUA ME SA CHE DEVO CAMBIA

-----
xint_dist_full <- purrr::map(xint_factors, function(H) {
  as.dist(1 - abs(cor(t(as.matrix(H)))))
})

xint_embeddings <- purrr::map(xint_dist_full, function(d) {
  cmdscale(d, k = min(10, attr(d, "Size") - 2))
})

----



# One function for all metrics computing all the internal validation metrics

compute_clustering_metrics <- function(labels, dist_full,
                                        db_x, db_d = NULL, centrotypes = "centroids") {
  # samples valid = non-NA label AND non-NA row in the aligned distance matrix
  dist_mat <- as.matrix(dist_full)
  valid <- !is.na(labels) & complete.cases(dist_mat)
  labels_v <- as.integer(labels)[valid]
  dist_v   <- as.dist(dist_mat[valid, valid])

  # index.DB always needs an actual coordinate matrix (raw factors or an MDS embedding),
  # never a raw distance matrix -- see clusterSim::index.DB source, which computes
  # column-wise means of x directly, even for centrotypes = "medoids"
  db_x_v <- db_x[valid, , drop = FALSE]
  db_d_v <- if (!is.null(db_d)) as.dist(as.matrix(db_d)[valid, valid]) else NULL

  stats <- fpc::cluster.stats(d = dist_v, clustering = labels_v)

  # connectivity: pass the SAME aligned distance used everywhere else via the
  # `distance` argument (not `Data`) -- `Data` triggers an internal recomputation
  # of distances from raw features, which is wrong when we only have a distance/affinity
  # matrix (SNF), and inconsistent with dist_full for every other method too.
  conn <- clValid::connectivity(distance = dist_v, clusters = labels_v)

  db <- if (!is.null(db_d_v)) {
    clusterSim::index.DB(db_x_v, labels_v, d = db_d_v, centrotypes = centrotypes)$DB
  } else {
    clusterSim::index.DB(db_x_v, labels_v, centrotypes = centrotypes)$DB
  }

  c(avg_silwidth = stats$avg.silwidth,
    ch_index      = stats$ch,
    db_index      = db,
    dunn_index    = stats$dunn,
    connectivity  = conn)
}


# Config: which inputs each clustering solution needs

approach_configs <- list(
  MOFA_KM_2 = list(dist_full = mofa_dist_full, db_x = mofa_factors_clean, db_d = NULL, centrotypes = "centroids"),
  MOFA_KM_3 = list(dist_full = mofa_dist_full, db_x = mofa_factors_clean, db_d = NULL, centrotypes = "centroids"),
  MOFA_KM_4 = list(dist_full = mofa_dist_full, db_x = mofa_factors_clean, db_d = NULL, centrotypes = "centroids"),

  # icluster_2 = list(dist_full = icluster_dist_full2, db_x = icluster_z_clean2, db_d = NULL, centrotypes = "centroids"),
  icluster_3 = list(dist_full = icluster_dist_full3, db_x = icluster_z_clean3, db_d = NULL, centrotypes = "centroids"),
  icluster_4 = list(dist_full = icluster_dist_full4, db_x = icluster_z_clean4, db_d = NULL, centrotypes = "centroids"),

  ibayes_3 = list(dist_full = ibayes_dist_full3, db_x = ibayes_z_clean3, db_d = NULL, centrotypes = "centroids"),
  ibayes_4 = list(dist_full = ibayes_dist_full4, db_x = ibayes_z_clean4, db_d = NULL, centrotypes = "centroids"),

  # SNF: no raw feature space, so index.DB uses the MDS embedding derived from
  # the fused distance, with medoids identified via the actual SNF distance
  SNF_2 = list(dist_full = snf_dist_full, db_x = snf_embedding, db_d = snf_dist_full, centrotypes = "medoids"),
  SNF_3 = list(dist_full = snf_dist_full, db_x = snf_embedding, db_d = snf_dist_full, centrotypes = "medoids"),
  SNF_4 = list(dist_full = snf_dist_full, db_x = snf_embedding, db_d = snf_dist_full, centrotypes = "medoids"),

  xintNMF_2 = list(dist_full = xint_dist_full[["2"]], db_x = xint_embeddings[["2"]],
                  db_d = xint_dist_full[["2"]], centrotypes = "medoids"),
  xintNMF_3 = list(dist_full = xint_dist_full[["3"]], db_x = xint_embeddings[["3"]],
                  db_d = xint_dist_full[["3"]], centrotypes = "medoids"),
  xintNMF_4 = list(dist_full = xint_dist_full[["4"]], db_x = xint_embeddings[["4"]],
                  db_d = xint_dist_full[["4"]], centrotypes = "medoids")
)


### Run all metrics in one loop

all_metrics <- imap(approach_configs, function(cfg, name) {
  compute_clustering_metrics(
    labels      = cluster_results[[name]],
    dist_full   = cfg$dist_full,
    db_x        = cfg$db_x,
    db_d        = cfg$db_d,
    centrotypes = cfg$centrotypes
  )
})

results_table <- as.data.frame(all_metrics)
rownames(results_table) <- c("Avg_Silhouette", "Calinski_Harabasz", "Davies_Bouldin", "Dunn_Index", "Connectivity")
print(round(results_table, 3))


----

### Plot

```{r}
library(stringr)

plot_data <- results_table %>%
  tibble::rownames_to_column("Metric") %>%
  pivot_longer(-Metric, names_to = "Approach", values_to = "Value") %>%
  mutate(Method_Family = str_extract(Approach, "^[A-Za-z]+"))

ggplot(plot_data, aes(x = Approach, y = Value, fill = Method_Family)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  facet_wrap(~Metric, scales = "free_y") +
  theme_minimal() +
  labs(title = "Multi-Omics Internal Validation Comparison",
       x = NULL, y = "Metric Value", fill = "Method") +
  scale_fill_brewer(palette = "Set2") +
  theme(strip.text = element_text(face = "bold", size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1))
```
