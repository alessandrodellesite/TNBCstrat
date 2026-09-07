library(SNFtool)


# SNF

#Data loading 

data_dir <- "/mnt/petasan_ccb/alessandro/SCANB/multiomics/input_data/snf_inputdata/"
rna_t  <- readRDS(file.path(data_dir, "rna_snf.rds"))
meth_t <- readRDS(file.path(data_dir, "met_snf.rds"))
cnv_t  <- readRDS(file.path(data_dir, "cnv_snf.rds"))


#Align samples

common_samples <- intersect(colnames(rna_matrix), intersect(colnames(met_mvals), colnames(cnv_log_ratio)))
rna_matrix <- rna_matrix[, common_samples]
met_mvals <- met_mvals[, common_samples]
cnv_log_ratio <- cnv_log_ratio[, common_samples]


#Transpose all matrices 

```{r}
data_rna <- t(rna_matrix)      
data_meth <- t(met_mvals)            
data_cnv <- t(cnv_log_ratio)
```

# Standard normalization (makes each feature (column) have a mean of 0 and SD of 1)

```{r}
data_rna_norm  <- standardNormalization(data_rna)
data_meth_norm <- standardNormalization(data_meth)
data_cnv_norm  <- standardNormalization(data_cnv)
```

# Calculate pairwise distance matrices for each data type 
(retrieves values for each patient pair that represent the distance between the 2 -> how similar their profile is for that data type)

```{r}
dist_rna  <- dist2(as.matrix(data_rna_norm),  as.matrix(data_rna_norm))
dist_meth <- dist2(as.matrix(data_meth_norm), as.matrix(data_meth_norm))
dist_cnv  <- dist2(as.matrix(data_cnv_norm),  as.matrix(data_cnv_norm))
```

# Convert to affinity matrices (distance values between 0-1)

```{r}
# Set hyperparameters
K <- 20       # Number of nearest neighbors (typically between 10 and 30)
sigma <- 0.5  # Hyperparameter scaling the variance (typically between 0.3 and 0.8)

# Construct individual similarity graphs
W_rna  <- affinityMatrix(dist_rna, K, sigma)
W_meth <- affinityMatrix(dist_meth, K, sigma)
W_cnv  <- affinityMatrix(dist_cnv, K, sigma)
```

# Fuse networks into 1 (SNF algorithm)

```{r}
# Set the number of iterations for network fusion
T_iter <- 20
# Fuse the 3 networks into 1 single unified network
W_fused <- SNF(list(W_rna, W_meth, W_cnv), K, T_iter)

# Save W_fused to an RDS file
saveRDS(W_fused, file = "snf_w_fused.rds")
```

# Choose the number of clusters based on the graph 

```{r}
estimated_k <- estimateNumberOfClustersGivenGraph(W_fused, NUMC=2:10)
print(estimated_k)
```

```{r}
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

write.csv(patient_subtypes_2, "SNF_clusters_2.csv", row.names = FALSE)
write.csv(patient_subtypes_4, "SNF_clusters_4.csv", row.names = FALSE)

#3 clusters
fused_groups_3 <- spectralClustering(W_fused, 3)
patient_subtypes_3 <- data.frame(
  SampleID = rownames(data_rna),
  Cluster = fused_groups_3
)
write.csv(patient_subtypes_3, "SNF_clusters_3.csv", row.names = FALSE)
```


symmetric similarity matrix - an adjacency matrix of a network (patients sorted by their cluster label)
hypothesis from 4-cluster division: one of the two main groups remains stable, while the other main group splits into three distinct sub-types when you look closer.
```{r}
# Visualize the clusters on the fused similarity matrix
displayClusters(W_fused, fused_groups_2)
displayClustersWithHeatmap(W_fused, fused_groups_2)

displayClusters(W_fused, fused_groups_4)
displayClustersWithHeatmap(W_fused, fused_groups_4)
```


```{r}
# Assignments for both solutions
group_2 <- spectralClustering(W_fused, 2)
group_4 <- spectralClustering(W_fused, 4)

# 2. Map them to your actual sample IDs
comparison_df <- data.frame(
  Sample_ID = rownames(data_rna),
  Clusters_2 = paste0("Subtype_", group_2),
  Clusters_4 = paste0("Subtype_", group_4)
)

# 3. Create a cross-tabulation table
cross_tab <- table(comparison_df$Clusters_2, comparison_df$Clusters_4)
print(cross_tab)
```
rows: groups from 2-clusters solution
columns: groups from 4-cluster solution

i valori sono quanti dei pazienti dei subtypes della 4-cluster solution provengono da quale dei 2 gruppi

next to verify if the K=4 sub-clusters neatly subdivide your K=2 primary groups

```{r}
palette_K2 <- c("lightblue", "tomato")          
palette_K4 <- c("#1b9e77", "#7570b9", "#2b5c8f", "#e7298a") 

colors_K2 <- palette_K2[group_2]
colors_K4 <- palette_K4[group_4]

color_matrix <- cbind(colors_K2, colors_K4)
colnames(color_matrix) <- c("2 Clusters", "4 Clusters")

displayClustersWithHeatmap(W_fused, group_4, ColSideColors = color_matrix)
```

# Alluvial plot
```{r}
#visualizes the shifting of samples between your K=2 and K=4 solutions
plotAlluvial(W_fused, clust.range = c(2, 4), col = color_matrix[, "2 Clusters"])
```
# Check the concordance (contribution) of each individual data type to the fused network

it runs spectralClustering separately on the multimics graph and the on the single omics affinity matrix.
It compares the resulting patient groupings pairwise using Normalized Mutual Information (NMI).
NMI is a metric from 0 to 1 that measures overlap (1: the two networks grouped the exact same patients together)
So the values are the agreement measures between the 2 clusters sol using single omics vs using all omics

```{r}
concordance <- concordanceNetworkNMI(list(W_fused, W_rna, W_meth, W_cnv), 2)
print(concordance)

concordance3 <- concordanceNetworkNMI(list(W_fused, W_rna, W_meth, W_cnv), 3)
print(concordance3)

concordance4 <- concordanceNetworkNMI(list(W_fused, W_rna, W_meth, W_cnv), 4)
print(concordance4)
```
