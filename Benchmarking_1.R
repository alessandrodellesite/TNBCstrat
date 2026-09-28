library(dplyr)
library(purrr)
library(mclust)   # Adjusted Rand Index
library(aricode)  # NMI
library(mcclust)  # VI
library(pheatmap)
library(stringr)
library(cluster)
library(fpc)
library(clValid)
library(clusterSim)
library(tidyr)
library(ggplot2)

# load clustering results
#mofa
mofa2_results <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_mofa/mofa_km_clusters_2.csv")
mofa3_results <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_mofa/mofa_km_clusters_3.csv")
mofa4_results <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_mofa/mofa_km_clusters_4.csv")
#icluster
icluster2_results <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_icluster/iclusters_clusters_2.csv")
icluster3_results <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_icluster/iclusters_clusters_3.csv")
icluster4_results <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_icluster/iclusters_clusters_4.csv")

#snf
snf2_results <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_snf/snf_clusters_2.csv")
snf3_results <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_snf/snf_clusters_3.csv")
snf4_results <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_snf/snf_clusters_4.csv")

#xintnmf
XintNMF2_results_reg <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_xintnmf/xintNMF_clusters_k2_reg.csv")
XintNMF3_results_reg <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_xintnmf/xintNMF_clusters_k3_reg.csv")
XintNMF4_results_reg <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_xintnmf/xintNMF_clusters_k4_reg.csv")

# merge results into one dataframe

lista_dataframe <- list(
  "MOFA_2"      = mofa2_results,
  "MOFA_3"      = mofa3_results,
  "MOFA_4"      = mofa4_results,
  "iCluster+_2"  = icluster2_results,
  "iCluster+_3"  = icluster3_results,
  "iCluster+_4"  = icluster4_results,
  "SNF_2"       = snf2_results,
  "SNF_3"       = snf3_results,
  "SNF_4"       = snf4_results,
  "X-intNMF_2"   = XintNMF2_results_reg,
  "X-intNMF_3"   = XintNMF3_results_reg,
  "X-intNMF_4"   = XintNMF4_results_reg
)

cluster_results <- imap(lista_dataframe, function(df, nuovo_nome) {
  df <- dplyr::select(df, SampleID, Cluster)
  colnames(df)[colnames(df) == "Cluster"] <- nuovo_nome
  df
}) %>%
  purrr::reduce(dplyr::full_join, by = "SampleID")

head(cluster_results)


# Benchmarking separated by k

out_dir <- "/mnt/petasan_ccb/alessandro/SCANB/plots/benchmarking/concordance"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

all_cols <- colnames(cluster_results)[colnames(cluster_results) != "SampleID"]
col_k <- str_extract(all_cols, "(?<=_)[0-9]+$")
k_values <- sort(unique(as.numeric(col_k)))

ari_matrices <- list()
nmi_matrices <- list()
vi_matrices  <- list()

for (k in k_values) {

  cols_k <- all_cols[col_k == as.character(k)]
  method_names <- str_remove(cols_k, paste0("_", k, "$"))

  n_methods <- length(cols_k)

  ari_matrix <- matrix(1, nrow = n_methods, ncol = n_methods,
                        dimnames = list(method_names, method_names))
  nmi_matrix <- matrix(1, nrow = n_methods, ncol = n_methods,
                        dimnames = list(method_names, method_names))
  vi_matrix  <- matrix(0, nrow = n_methods, ncol = n_methods,
                        dimnames = list(method_names, method_names))

  if (n_methods > 1) {
    for (i in 1:(n_methods - 1)) {
      for (j in (i + 1):n_methods) {
        m1 <- cols_k[i]
        m2 <- cols_k[j]
        name1 <- method_names[i]
        name2 <- method_names[j]

        labels_1 <- cluster_results[[m1]]
        labels_2 <- cluster_results[[m2]]

        valid_idx <- !is.na(labels_1) & !is.na(labels_2)
        clean_l1  <- labels_1[valid_idx]
        clean_l2  <- labels_2[valid_idx]

        if (length(clean_l1) > 0) {
          ari_val <- mclust::adjustedRandIndex(clean_l1, clean_l2)
          ari_matrix[name1, name2] <- ari_val
          ari_matrix[name2, name1] <- ari_val

          nmi_val <- aricode::NMI(clean_l1, clean_l2)
          nmi_matrix[name1, name2] <- nmi_val
          nmi_matrix[name2, name1] <- nmi_val

          raw_vi <- mcclust::vi.dist(clean_l1, clean_l2, base = exp(1))
          max_possible_entropy <- log(length(clean_l1))
          normalized_vi <- raw_vi / max_possible_entropy

          vi_matrix[name1, name2] <- normalized_vi
          vi_matrix[name2, name1] <- normalized_vi
        } else {
          ari_matrix[name1, name2] <- NA; ari_matrix[name2, name1] <- NA
          nmi_matrix[name1, name2] <- NA; nmi_matrix[name2, name1] <- NA
          vi_matrix[name1, name2]  <- NA; vi_matrix[name2, name1]  <- NA
        }
      }
    }
  }

  ari_matrices[[as.character(k)]] <- ari_matrix
  nmi_matrices[[as.character(k)]] <- nmi_matrix
  vi_matrices[[as.character(k)]]  <- vi_matrix

  cat("k =", k, "\n")
  cat("ADJUSTED RAND INDEX (ARI) \n")
  print(round(ari_matrix, 3))
  cat("\nNORMALIZED MUTUAL INFORMATION (NMI) \n")
  print(round(nmi_matrix, 3))
  cat("\nNORMALIZED VARIATION OF INFORMATION (VI) \n")
  print(round(vi_matrix, 3))

  #  Save ARI heatmap 
  png(file.path(out_dir, paste0("ARI_k", k, ".png")),
      width = 1200, height = 1000, res = 150)
  pheatmap(ari_matrix,
           main = paste0("Pairwise Cluster Concordance (ARI), k = ", k),
           display_numbers = TRUE,
           color = colorRampPalette(c("white", "#E8F0FE", "#1A73E8"))(50),
           number_color = "black")
  dev.off()

  # Save NMI heatmap 
  png(file.path(out_dir, paste0("NMI_k", k, ".png")),
      width = 1200, height = 1000, res = 150)
  pheatmap(nmi_matrix,
           main = paste0("Normalized Mutual Information (NMI), k = ", k),
           display_numbers = TRUE,
           color = colorRampPalette(c("white", "#E8F0FE", "#2E6651"))(50),
           number_color = "black")
  dev.off()

  # Save VI heatmap 
  png(file.path(out_dir, paste0("VI_k", k, ".png")),
      width = 1200, height = 1000, res = 150)
  pheatmap(vi_matrix,
           main = paste0("Normalized Variation of Information (VI), k = ", k),
           display_numbers = TRUE,
           color = colorRampPalette(c("#1A7666", "#E8F0FE", "white"))(50),
           number_color = "black")
  dev.off()
}









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


# Replicate the exact embedding SNFtool::spectralClustering() computes internally 
# : row-normalized eigenvectors of the symmetric
# normalized graph Laplacian. This (not 1-W) is the space the clusters
# were actually formed in.
snf_spectral_embedding <- function(W, K) {
  d <- rowSums(W)
  d[d == 0] <- .Machine$double.eps
  D <- diag(d)
  L <- D - W
  Di <- diag(1 / sqrt(d))
  NL <- Di %*% L %*% Di
  NL <- (NL + t(NL)) / 2   # guard against floating-point asymmetry before eigen()

  eig <- eigen(NL)
  ord <- order(abs(eig$values))                      # smallest |eigenvalue| first, same as spectralClustering()
  U <- eig$vectors[, ord[1:K], drop = FALSE]
  U <- t(apply(U, 1, function(x) x / sqrt(sum(x^2)))) # row-normalize -- the "type = 3" step
  rownames(U) <- rownames(W)
  U
}

snf_embeddings <- purrr::map(c("2" = 2, "3" = 3, "4" = 4), function(k) {
  snf_spectral_embedding(snf_w_clean, K = k)
})

snf_dist_full <- purrr::map(snf_embeddings, function(U) dist(U, method = "euclidean"))




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

# use the SAME distance metric that  produced the xintNMF clusters

xint_dist_full <- purrr::map(xint_factors, function(H) {
  dist(as.matrix(H), method = "euclidean")
})



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
  MOFA_2 = list(dist_full = mofa_dist_full, db_x = mofa_factors_clean, db_d = NULL, centrotypes = "centroids"),
  MOFA_3 = list(dist_full = mofa_dist_full, db_x = mofa_factors_clean, db_d = NULL, centrotypes = "centroids"),
  MOFA_4 = list(dist_full = mofa_dist_full, db_x = mofa_factors_clean, db_d = NULL, centrotypes = "centroids"),

  "iCluster+_2" = list(dist_full = icluster_dist_full2, db_x = icluster_z_clean2, db_d = NULL, centrotypes = "centroids"),
  "iCluster+_3" = list(dist_full = icluster_dist_full3, db_x = icluster_z_clean3, db_d = NULL, centrotypes = "centroids"),
  "iCluster+_4" = list(dist_full = icluster_dist_full4, db_x = icluster_z_clean4, db_d = NULL, centrotypes = "centroids"),

  SNF_2 = list(dist_full = snf_dist_full[["2"]], db_x = snf_embeddings[["2"]], db_d = NULL, centrotypes = "centroids"),
  SNF_3 = list(dist_full = snf_dist_full[["3"]], db_x = snf_embeddings[["3"]], db_d = NULL, centrotypes = "centroids"),
  SNF_4 = list(dist_full = snf_dist_full[["4"]], db_x = snf_embeddings[["4"]], db_d = NULL, centrotypes = "centroids"),

  "X-intNMF_2" = list(dist_full = xint_dist_full[["2"]], db_x = as.matrix(xint_factors[["2"]]), db_d = NULL, centrotypes = "centroids"),
  "X-intNMF_3" = list(dist_full = xint_dist_full[["3"]], db_x = as.matrix(xint_factors[["3"]]), db_d = NULL, centrotypes = "centroids"),
  "X-intNMF_4" = list(dist_full = xint_dist_full[["4"]], db_x = as.matrix(xint_factors[["4"]]), db_d = NULL, centrotypes = "centroids")
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
val_dir <- "/mnt/petasan_ccb/alessandro/SCANB/plots/benchmarking/internal_validation"
write.csv(round(results_table, 3), file.path(val_dir, "internal_validation.csv"), row.names = FALSE)

# Plot
                            
metric_info <- tibble::tribble(
  ~Metric,               ~Direction, ~Label,
  "Avg_Silhouette",      "higher",   "Average Silhouette\n(higher is better \u2191)",
  "Calinski_Harabasz",   "higher",   "Calinski-Harabasz\n(higher is better \u2191)",
  "Dunn_Index",          "higher",   "Dunn Index\n(higher is better \u2191)",
  "Davies_Bouldin",      "lower",    "Davies-Bouldin\n(lower is better \u2193)",
  "Connectivity",        "lower",    "Connectivity\n(lower is better \u2193)"
)

plot_data <- results_table %>%
  tibble::rownames_to_column("Metric") %>%
  pivot_longer(-Metric, names_to = "Approach", values_to = "Value") %>%
  mutate(Method_Family = str_remove(Approach, "_\\d+$")) %>%
  left_join(metric_info, by = "Metric") %>%
  mutate(
    Label    = factor(Label, levels = metric_info$Label),   # fixes panel order
    Approach = factor(Approach, levels = names(approach_configs))
  )

# One background rectangle per panel tinted by direction
bg_data <- metric_info %>%
  mutate(Label = factor(Label, levels = metric_info$Label))

p <- ggplot(plot_data, aes(x = Approach, y = Value, fill = Method_Family)) +
  geom_rect(data = bg_data, inherit.aes = FALSE,
            aes(fill = NULL, alpha = NULL),
            xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
            fill = ifelse(bg_data$Direction == "higher", "#E8F4EA", "#FBE9E6")) +
  geom_col(position = "dodge", width = 0.6) +
  facet_wrap(~Label, ncol = 3, scales = "free_y") +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  labs(title = "Multi-omics internal validation comparison",
       subtitle = "Top row: higher is better (\u2191)   |   Bottom row: lower is better (\u2193)",
       x = NULL, y = "Metric Value", fill = "Method") +
  theme(strip.text = element_text(face = "bold", size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.spacing = unit(1, "lines"),
        legend.position = "bottom")

png(file.path("/mnt/petasan_ccb/alessandro/SCANB/plots/benchmarking/internal_validation",
              "internal_validation_comparison.png"),
    width = 1400, height = 1000, res = 150)
print(p)
dev.off()


                            

                  
# plot by number of clusters

plot_data <- results_table %>%
  tibble::rownames_to_column("Metric") %>%
  pivot_longer(-Metric, names_to = "Approach", values_to = "Value") %>%
  mutate(
    Method_Family   = str_extract(Approach, "^[A-Za-z]+"),
    Cluster_Version = str_extract(Approach, "\\d+")
  )

                            

# one plot per cluster solution (2, 3, 4)
plots_by_k <- plot_data %>%
  group_split(Cluster_Version) %>%
  set_names(purrr::map_chr(., ~ unique(.x$Cluster_Version))) %>%
  purrr::map(function(df) {
    k <- unique(df$Cluster_Version)
    ggplot(df, aes(x = Method_Family, y = Value, fill = Method_Family)) +
      geom_bar(stat = "identity", position = "dodge", width = 0.6) +
      facet_wrap(~ Metric, scales = "free_y") +
      theme_minimal() +
      labs(title = paste0("Internal validation (", k, " clusters)"),
           x = NULL, y = "Metric Value", fill = "Method") +
      scale_fill_brewer(palette = "Set2") +
      theme(strip.text = element_text(face = "bold", size = 11),
            axis.text.x = element_text(angle = 45, hjust = 1))
  })
                            
                            
out_dir_internal <- "/mnt/petasan_ccb/alessandro/SCANB/plots/benchmarking/internal_validation"
if (!dir.exists(out_dir_internal)) dir.create(out_dir_internal, recursive = TRUE)

purrr::iwalk(plots_by_k, function(p, k) {
  png(file.path(out_dir_internal, paste0("internal_validation_k", k, ".png")),
      width = 1400, height = 1000, res = 150)
  print(p)
  dev.off()
})
