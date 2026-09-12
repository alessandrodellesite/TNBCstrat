library(cluster)
library(mclust)

# Load and extract matrices for each k (2:6) 

base_dir <- "/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_xintnmf/"
k_range  <- 2:6

h_matrices_reg <- list()
w_matrices_reg <- list()  # nested: w_matrices[["3"]]$rna, $meth, $cnv

for (k in k_range) {
  folder <- paste0(base_dir, "rankselect_k", k, "_graphreg/")

  h_path    <- paste0(folder, "sample_factor.csv")
  w_rna_p   <- paste0(folder, "rna_factor.csv")
  w_meth_p  <- paste0(folder, "methylation_factor.csv")
  w_cnv_p   <- paste0(folder, "cnv_factor.csv")

  if (all(file.exists(c(h_path, w_rna_p, w_meth_p, w_cnv_p)))) {
    h_matrices_reg[[as.character(k)]] <- as.matrix(read.csv(h_path, header = TRUE))
    w_matrices_reg[[as.character(k)]] <- list(
      rna  = as.matrix(read.csv(w_rna_p,  header = TRUE)),
      meth = as.matrix(read.csv(w_meth_p, header = TRUE)),
      cnv  = as.matrix(read.csv(w_cnv_p,  header = TRUE))
    )
    cat("Loaded k =", k, "| H dim:", paste(dim(h_matrices_reg[[as.character(k)]]), collapse = " x "), "\n")
  } else {
    warning(paste("Missing files for k =", k))
  }
}

# Kmeans on H matrices + silhouette and WSS evaluation

set.seed(123)  # reproducibility
cluster_assignments_km_reg <- list()
mean_silhouettes_km_reg <- numeric(length(k_range))
names(mean_silhouettes_km_reg) <- k_range
wss_per_k <- numeric(length(k_range))
names(wss_per_k) <- k_range


for (k in k_range) {
  H <- h_matrices_reg[[as.character(k)]]
  km <- kmeans(H, centers = k, nstart = 50)   
  cluster_assignments_km_reg[[as.character(k)]] <- km$cluster

  #silhouette
  dist_mat <- dist(H, method = "euclidean")
  sil <- silhouette(km$cluster, dist_mat)
  mean_silhouettes_km_reg[as.character(k)] <- mean(sil[, "sil_width"])

  #wss
  wss_per_k[as.character(k)] <- km$tot.withinss
}

print("Silhouette for kmeans at different k's:")
print(round(mean_silhouettes_km_reg, 3))
print("WSS for kmeans at different k's:")
print(round(wss_per_k, 3))

#RSS (requires the input matrices, Reconstruction is W %*% t(H), compared against the real data)

data_dir <- "/mnt/petasan_ccb/alessandro/SCANB/multiomics/input_data/xintnmf_inputdata"
rna_matrix <- read.csv(file.path(data_dir, "rna.tsv"), sep = "\t", row.names = 1, check.names = FALSE)
met_matrix <- read.csv(file.path(data_dir, "methylation.tsv"), sep = "\t", row.names = 1, check.names = FALSE)
dna_matrix <- read.csv(file.path(data_dir, "cnv.tsv"), sep = "\t", row.names = 1, check.names = FALSE)

rss_per_k_reg <- numeric(length(k_range))
names(rss_per_k_reg) <- k_range

for (k in k_range) {
  H <- h_matrices_reg[[as.character(k)]]
  W <- w_matrices_reg[[as.character(k)]]

  rss_rna  <- sum((as.matrix(rna_matrix) - (W$rna  %*% t(H)))^2)
  rss_meth <- sum((as.matrix(met_matrix) - (W$meth %*% t(H)))^2)
  rss_cnv  <- sum((as.matrix(dna_matrix) - (W$cnv  %*% t(H)))^2)

  rss_per_k_reg[as.character(k)] <- rss_rna + rss_meth + rss_cnv
}

print("RSS values per k:")
print(round(rss_per_k_reg, 1))

png("/mnt/petasan_ccb/alessandro/SCANB/plots/multiomics/xintnmf/kmeans_wss_silh.png",
    width = 10, height = 5, units = "in", res = 300)
par(mfrow = c(1, 2))
plot(k_range, wss_per_k, type = "b", pch = 19, col = "darkorange", lwd = 2,
     xlab = "Rank (k)", ylab = "Within-cluster sum of squares", main = "WSS profile")
plot(k_range, mean_silhouettes_km_reg, type = "b", pch = 19, col = "royalblue", lwd = 2,
     xlab = "Rank (k)", ylab = "Mean Silhouette Width", main = "Silhouette profile")
dev.off()

#Testing hclust (with euclidian distance) - measure silhouette

set.seed(123)

cluster_assignments_hc_reg <- list()
mean_silhouettes_hc_reg <- numeric(length(k_range))
names(mean_silhouettes_hc_reg) <- k_range

for (k in k_range) {
  H <- h_matrices_reg[[as.character(k)]]
  dist_mat <- dist(H, method = "euclidean")
  
  hc <- hclust(dist_mat, method = "ward.D2")  # ward.D2 is the closest analogue to kmeans' objective
  clusters_hc <- cutree(hc, k = k)
  cluster_assignments_hc_reg[[as.character(k)]] <- clusters_hc
  
  sil <- silhouette(clusters_hc, dist_mat)
  mean_silhouettes_hc_reg[as.character(k)] <- mean(sil[, "sil_width"])
}
print("Silhouette for hclust at different k's:")
print(round(mean_silhouettes_hc_reg, 3))


png("/mnt/petasan_ccb/alessandro/SCANB/plots/multiomics/xintnmf/kmeans_vs_hclust_silh.png",
    width = 6, height = 5, units = "in", res = 300)
plot(k_range, mean_silhouettes_km_reg, type = "b", pch = 19, col = "royalblue", lwd = 2,
     ylim = range(c(mean_silhouettes_km_reg, mean_silhouettes_hc_reg)),
     xlab = "Rank (k)", ylab = "Mean Silhouette Width", main = "k-means vs Hierarchical")
lines(k_range, mean_silhouettes_hc_reg, type = "b", pch = 17, col = "firebrick", lwd = 2)
legend("topright", legend = c("k-means", "hierarchical (ward.D2)"),
       col = c("royalblue", "firebrick"), pch = c(19, 17))
dev.off()


#Measure ARI: how concordant kmeans and hclust results are (for each k)

ari_km_hc_reg <- numeric(length(k_range))
names(ari_km_hc_reg) <- k_range

for (k in k_range) {
  km_labels <- cluster_assignments_km_reg[[as.character(k)]]
  hc_labels <- cluster_assignments_hc_reg[[as.character(k)]]
  ari_km_hc_reg[as.character(k)] <- adjustedRandIndex(km_labels, hc_labels)
}
print("ARI values comapring kmeans and hclust clusterings")
print(round(ari_km_hc_reg, 3))

png("/mnt/petasan_ccb/alessandro/SCANB/plots/multiomics/xintnmf/kmeans_vs_hclust_ari.png",
    width = 6, height = 5, units = "in", res = 300)
plot(k_range, ari_km_hc_reg, type = "b", pch = 19, col = "darkgreen", lwd = 2,
     ylim = c(-0.1, 1),
     xlab = "Rank (k)", ylab = "Adjusted Rand Index",
     main = "Agreement: k-means vs Hierarchical (Ward.D2)")
abline(h = 0, lty = 2, col = "gray")
dev.off()

#Final export: only kmeans (consistent with other multiomics approach)

out_dir <- "/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_xintnmf"

for (k in 2:4) {
  final_cluster_df_reg <- data.frame(
    SampleID = colnames(rna_matrix),
    Cluster  = as.factor(cluster_assignments_km_reg[[as.character(k)]])
  )
  write.csv(
    final_cluster_df_reg,
    file = file.path(out_dir, paste0("xintNMF_clusters_k", k, "_reg.csv")),
    row.names = FALSE
  )
}
