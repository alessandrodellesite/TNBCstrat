library(dplyr)
library(purrr)
library(mclust) #for Adjusted Rand Index
library(aricode) #for fast computation of Normalized Mutual Information, and it also computes Normalized Variation of Information
library(mcclust) #provides the original Variation of Information distance metric
library(pheatmap)

# load clustering results

#mofa
mofa2_results <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_mofa/mofa_km_clusters_2.csv")
mofa3_results <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_mofa/mofa_km_clusters_3.csv")
mofa4_results <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_mofa/mofa_km_clusters_4.csv")

#icluster
icluster2_results <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_icluster/iclusters_clusters_2.csv")
icluster3_results <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_icluster/iclusters_clusters_3.csv")                                
icluster4_results <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_icluster/iclusters_clusters_4.csv")
#icluster_5 <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_icluster/iclusters_clusters_5.csv")
                                 
#snf
snf2_results <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_snf/snf_clusters_2.csv")
snf3_results <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_snf/snf_clusters_3.csv")                                
snf4_results <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_snf/snf_clusters_4.csv")
                                 
#xintnmf                                 
XintNMF2_results_reg <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_xintnmf/xintNMF_clusters_k2_reg.csv")                    
XintNMF3_results_reg <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_xintnmf/xintNMF_clusters_k3_reg.csv")
XintNMF4_results_reg <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_xintnmf/xintNMF_clusters_k4_reg.csv")                                


# merge results into one dataframe

# change cluster column name
lista_dataframe <- list(
  "MOFA_2"  = mofa2_results,
  "MOFA_3"  = mofa3_results,
  "MOFA_4"  = mofa4_results,
  "icluster_2" = icluster2_results,
  "icluster_3" = icluster3_results,
  "icluster_4" = icluster4_results,
  "SNF_2"      = snf2_results,
  "SNF_3"      = snf3_results,
  "SNF_4"      = snf4_results,
  "xintNMF_2"  = XintNMF2_results_reg,
  "xintNMF_3"  = XintNMF3_results_reg,
  "xintNMF_4"  = XintNMF4_results_reg
)

#full_join using sampleID
cluster_results <- imap(lista_dataframe, function(df, nuovo_nome) {
  df <- dplyr::select(df, SampleID, Cluster)
  colnames(df)[colnames(df) == "Cluster"] <- nuovo_nome
  df
}) %>%
  purrr::reduce(dplyr::full_join, by = "SampleID")

head(cluster_results)

# benchmarking

## Results concordance between different methods 
#Metrics:  Adjusted Rand Index (ARI), Normalized Mutual information (NMI), Jaccard variation of information


# exclude the first column to get only the cluster label columns
clustering_methods <- colnames(cluster_results)[colnames(cluster_results) != "SampleID"]
n_methods <- length(clustering_methods)
n_patients <- nrow(cluster_results)

# initialize empty matrices for the metrics
ari_matrix <- matrix(1, nrow = n_methods, ncol = n_methods, dimnames = list(clustering_methods, clustering_methods))
nmi_matrix <- matrix(1, nrow = n_methods, ncol = n_methods, dimnames = list(clustering_methods, clustering_methods))
vi_matrix  <- matrix(0, nrow = n_methods, ncol = n_methods, dimnames = list(clustering_methods, clustering_methods))


# compute pairwise metrics
for (i in 1:(n_methods - 1)) {
  for (j in (i + 1):n_methods) {
    m1 <- clustering_methods[i]
    m2 <- clustering_methods[j]
    
    # Extract cluster labels for both methods
    labels_1 <- cluster_results[[m1]]
    labels_2 <- cluster_results[[m2]]
    
    # Handle potential missing data (NA values) if any sample didn't get clustered by a method
    valid_idx <- !is.na(labels_1) & !is.na(labels_2)
    clean_l1  <- labels_1[valid_idx]
    clean_l2  <- labels_2[valid_idx]
    
    if (length(clean_l1) > 0) {
      # Calculate Adjusted Rand Index (ARI)
      ari_val <- mclust::adjustedRandIndex(clean_l1, clean_l2)
      ari_matrix[m1, m2] <- ari_val
      ari_matrix[m2, m1] <- ari_val
      
      # Calculate Normalized Mutual Information (NMI)
      nmi_val <- aricode::NMI(clean_l1, clean_l2)
      nmi_matrix[m1, m2] <- nmi_val
      nmi_matrix[m2, m1] <- nmi_val
      
      # Calculate Normalized Variation of Information (VI)
      # Bounded between 0 (identical) and 1 (different)
      raw_vi <- mcclust::vi.dist(clean_l1, clean_l2, base = exp(1))
      max_possible_entropy <- log(length(clean_l1)) 
      normalized_vi <- raw_vi / max_possible_entropy
      
      vi_matrix[m1, m2] <- normalized_vi
      vi_matrix[m2, m1] <- normalized_vi
    } else {
      ari_matrix[m1, m2] <- NA; ari_matrix[m2, m1] <- NA
      nmi_matrix[m1, m2] <- NA; nmi_matrix[m2, m1] <- NA
      vi_matrix[m1, m2]  <- NA; vi_matrix[m2, m1]  <- NA
    }
  }
}



cat("ADJUSTED RAND INDEX (ARI) \n")
print(round(ari_matrix, 3))

cat("\nNORMALIZED MUTUAL INFORMATION (NMI) \n")
print(round(nmi_matrix, 3))

cat("\nNORMALIZED VARIATION OF INFORMATION (VI) \n")
print(round(vi_matrix, 3))

# Plotting the ARI similarity matrix as a heatmap
pheatmap(ari_matrix, 
         main = "Pairwise Cluster Concordance (ARI)", 
         display_numbers = TRUE, 
         color = colorRampPalette(c("white", "#E8F0FE", "#1A73E8"))(50),
         number_color = "black")



pheatmap(nmi_matrix, 
         main = "NORMALIZED MUTUAL INFORMATION (NMI)", 
         display_numbers = TRUE, 
         color = colorRampPalette(c("white", "#E8F0FE", "#2E6651"))(50),
         number_color = "black")


pheatmap(vi_matrix, 
         main = "NORMALIZED VARIATION OF INFORMATION (VI)", 
         display_numbers = TRUE, 
         color = colorRampPalette(c("#1A7666", "#E8F0FE", "white"))(50),
         number_color = "black")

