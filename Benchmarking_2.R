library(dplyr)
library(purrr)
library(mclust)   # Adjusted Rand Index
library(aricode)  # NMI
library(mcclust)  # VI
library(pheatmap)
library(stringr)

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
  "icluster_2"  = icluster2_results,
  "icluster_3"  = icluster3_results,
  "icluster_4"  = icluster4_results,
  "SNF_2"       = snf2_results,
  "SNF_3"       = snf3_results,
  "SNF_4"       = snf4_results,
  "xintNMF_2"   = XintNMF2_results_reg,
  "xintNMF_3"   = XintNMF3_results_reg,
  "xintNMF_4"   = XintNMF4_results_reg
)

cluster_results <- imap(lista_dataframe, function(df, nuovo_nome) {
  df <- dplyr::select(df, SampleID, Cluster)
  colnames(df)[colnames(df) == "Cluster"] <- nuovo_nome
  df
}) %>%
  purrr::reduce(dplyr::full_join, by = "SampleID")

head(cluster_results)

# ----------------------------------------------------------------
# Benchmarking, separated by k
# ----------------------------------------------------------------
# columns other than SampleID, each named "<Method>_<k>"
all_cols <- colnames(cluster_results)[colnames(cluster_results) != "SampleID"]

# extract the k value from the column name suffix
col_k <- str_extract(all_cols, "(?<=_)[0-9]+$")
k_values <- sort(unique(as.numeric(col_k)))

# containers to keep results for all k, in case you need them later
ari_matrices <- list()
nmi_matrices <- list()
vi_matrices  <- list()

for (k in k_values) {

  cols_k <- all_cols[col_k == as.character(k)]
  # relabel columns with just the method name (strip "_k") for readability in plots
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
          # ARI
          ari_val <- mclust::adjustedRandIndex(clean_l1, clean_l2)
          ari_matrix[name1, name2] <- ari_val
          ari_matrix[name2, name1] <- ari_val

          # NMI
          nmi_val <- aricode::NMI(clean_l1, clean_l2)
          nmi_matrix[name1, name2] <- nmi_val
          nmi_matrix[name2, name1] <- nmi_val

          # Normalized VI
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

  cat("\n============================\n")
  cat("k =", k, "\n")
  cat("============================\n")
  cat("ADJUSTED RAND INDEX (ARI) \n")
  print(round(ari_matrix, 3))
  cat("\nNORMALIZED MUTUAL INFORMATION (NMI) \n")
  print(round(nmi_matrix, 3))
  cat("\nNORMALIZED VARIATION OF INFORMATION (VI) \n")
  print(round(vi_matrix, 3))

  # Plots for this k
  pheatmap(ari_matrix,
           main = paste0("Pairwise Cluster Concordance (ARI), k = ", k),
           display_numbers = TRUE,
           color = colorRampPalette(c("white", "#E8F0FE", "#1A73E8"))(50),
           number_color = "black")

  pheatmap(nmi_matrix,
           main = paste0("Normalized Mutual Information (NMI), k = ", k),
           display_numbers = TRUE,
           color = colorRampPalette(c("white", "#E8F0FE", "#2E6651"))(50),
           number_color = "black")

  pheatmap(vi_matrix,
           main = paste0("Normalized Variation of Information (VI), k = ", k),
           display_numbers = TRUE,
           color = colorRampPalette(c("#1A7666", "#E8F0FE", "white"))(50),
           number_color = "black")
}
