library(SNFtool)


# SNF

# Specific SNF's data loading 

data_dir <- "/mnt/petasan_ccb/alessandro/SCANB/multiomics/input_data/snf_inputdata/"
data_rna  <- readRDS(file.path(data_dir, "rna_snf.rds"))
data_meth <- readRDS(file.path(data_dir, "met_snf.rds"))
data_cnv  <- readRDS(file.path(data_dir, "cnv_snf.rds"))

# Standard normalization (makes each feature (column) have a mean of 0 and SD of 1)

data_rna_norm  <- standardNormalization(data_rna)
data_meth_norm <- standardNormalization(data_meth)
data_cnv_norm  <- standardNormalization(data_cnv)

# Calculate pairwise distance matrices for each data type 
#(retrieves values for each patient pair that represent the distance between the 2 -> how similar their profile is for that data type)

dist_rna  <- dist2(as.matrix(data_rna_norm),  as.matrix(data_rna_norm))
dist_meth <- dist2(as.matrix(data_meth_norm), as.matrix(data_meth_norm))
dist_cnv  <- dist2(as.matrix(data_cnv_norm),  as.matrix(data_cnv_norm))

# Convert to affinity matrices (distance values between 0-1)

# Set hyperparameters
K <- 20       # Number of nearest neighbors (typically between 10 and 30)
sigma <- 0.5  # Hyperparameter scaling the variance (typically between 0.3 and 0.8)

# Construct individual similarity graphs
W_rna  <- affinityMatrix(dist_rna, K, sigma)
W_meth <- affinityMatrix(dist_meth, K, sigma)
W_cnv  <- affinityMatrix(dist_cnv, K, sigma)

# Fuse networks into 1 (SNF algorithm)
# Set the number of iterations for network fusion
T_iter <- 20
# Fuse the 3 networks into 1 single unified network
W_fused <- SNF(list(W_rna, W_meth, W_cnv), K, T_iter)

# Save W_fused to an RDS file
saveRDS(W_fused, file = "/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_snf/snf_w_fused.rds")


# Convergence check: does the fused network stabilize before T_iter = 20?
max_T <- 20
diffs <- numeric(max_T - 1)

W_prev <- SNF(list(W_rna, W_meth, W_cnv), K, 1)
for (t in 2:max_T) {
  W_curr <- SNF(list(W_rna, W_meth, W_cnv), K, t)
  diffs[t - 1] <- norm(W_curr - W_prev, type = "F")
  W_prev <- W_curr
}

png("/mnt/petasan_ccb/alessandro/SCANB/plots/multiomics/snf/snf_convergence.png",
    width = 8, height = 5, units = "in", res = 300)

plot(2:max_T, diffs, type = "b",
     xlab = "Iteration (T)", ylab = "Frobenius norm of change",
     main = "SNF convergence check")

dev.off()

# Choose the number of clusters based on the graph 

estimated_k <- estimateNumberOfClustersGivenGraph(W_fused, NUMC=2:10)
print(estimated_k)


# Perform Spectral Clustering on the fused network
fused_groups_2 <- spectralClustering(W_fused, 2) 
fused_groups_4 <- spectralClustering(W_fused, 4)

# Create a final data frame mapping sample names to their assigned cluster
patient_subtypes_2 <- data.frame(
  SampleID = rownames(data_rna),
  Cluster = fused_groups_2
)
patient_subtypes_4 <- data.frame(
  SampleID = rownames(data_rna),
  Cluster = fused_groups_4
)

write.csv(patient_subtypes_2, "/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_snf/SNF_clusters_2.csv", row.names = FALSE)
write.csv(patient_subtypes_4, "/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_snf/SNF_clusters_4.csv", row.names = FALSE)

#3 clusters
fused_groups_3 <- spectralClustering(W_fused, 3)
patient_subtypes_3 <- data.frame(
  SampleID = rownames(data_rna),
  Cluster = fused_groups_3
)
write.csv(patient_subtypes_3, "/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_snf/SNF_clusters_3.csv", row.names = FALSE)



# Visualize the clusters on the fused similarity matrix
# Cluster 2 - simple display
png("/mnt/petasan_ccb/alessandro/SCANB/plots/multiomics/snf/fused_clusters_2.png",
    width = 8, height = 8, units = "in", res = 300)
displayClusters(W_fused, fused_groups_2)
dev.off()

# Cluster 2 - heatmap version
png("/mnt/petasan_ccb/alessandro/SCANB/plots/multiomics/snf/fused_clusters_2_heatmap.png",
    width = 8, height = 8, units = "in", res = 300)
displayClustersWithHeatmap(W_fused, fused_groups_2)
dev.off()

# Cluster 4 - simple display
png("/mnt/petasan_ccb/alessandro/SCANB/plots/multiomics/snf/fused_clusters_4.png",
    width = 8, height = 8, units = "in", res = 300)
displayClusters(W_fused, fused_groups_4)
dev.off()

# Cluster 4 - heatmap version
png("/mnt/petasan_ccb/alessandro/SCANB/plots/multiomics/snf/fused_clusters_4_heatmap.png",
    width = 8, height = 8, units = "in", res = 300)
displayClustersWithHeatmap(W_fused, fused_groups_4)
dev.off()


# Alluvial plot
#visualizes the shifting of samples between your K=2 and K=4 solutions

palette_K2 <- c("lightblue", "tomato")          
palette_K4 <- c("#1b9e77", "#7570b9", "#2b5c8f", "#e7298a") 

colors_K2 <- palette_K2[fused_groups_2]
colors_K4 <- palette_K4[fused_groups_4]

color_matrix <- cbind(colors_K2, colors_K4)
colnames(color_matrix) <- c("2 Clusters", "4 Clusters")
plotAlluvial(W_fused, clust.range = c(2, 4), col = color_matrix[, "2 Clusters"])

png("/mnt/petasan_ccb/alessandro/SCANB/plots/multiomics/snf/alluvional_plot.png",
    width = 8, height = 8, units = "in", res = 300)
plotAlluvial(W_fused, clust.range = c(2, 4), col = color_matrix[, "2 Clusters"])
dev.off()


# Check the concordance (contribution) of each individual data type to the fused network

concordance <- concordanceNetworkNMI(list(W_fused, W_rna, W_meth, W_cnv), 2)
dimnames(concordance) <- list(
  c("W_fused", "W_rna", "W_meth", "W_cnv"),
  c("W_fused", "W_rna", "W_meth", "W_cnv")
)
write.csv(concordance, "/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_snf/concordance_nmi.csv")
print(concordance)

concordance3 <- concordanceNetworkNMI(list(W_fused, W_rna, W_meth, W_cnv), 3)
dimnames(concordance3) <- list(
  c("W_fused", "W_rna", "W_meth", "W_cnv"),
  c("W_fused", "W_rna", "W_meth", "W_cnv")
)
write.csv(concordance3, "/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_snf/concordance3_nmi.csv")
print(concordance3)

concordance4 <- concordanceNetworkNMI(list(W_fused, W_rna, W_meth, W_cnv), 4)
dimnames(concordance4) <- list(
  c("W_fused", "W_rna", "W_meth", "W_cnv"),
  c("W_fused", "W_rna", "W_meth", "W_cnv")
)
write.csv(concordance4, "/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_snf/concordance4_nmi.csv")
print(concordance4)
