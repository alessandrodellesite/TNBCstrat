#Analysis of icluster model tuning

# visualize model tuning results: how many clusters to choose
#BiocManager::install("iClusterPlus")
library(iClusterPlus)

# Read output files
output <- list()
dir <- "/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_icluster"

for (k in 1:5) {
  path <- file.path(dir, paste0("cv.fit.k", k, ".Rdata"))
  env <- new.env()
  load(path, envir = env)
  output[[k]] <- env$cv.fit
}

nLambda = nrow(output[[1]]$lambda) #nlambda tested: 185
nK = length(output) #nK tested: 5

#selection parameters for the model (each matrix has nlambda rows and nK cols)
BIC = getBIC(output) #Bayesian Information Criterion matrix
devR = getDevR(output) #Deviance Ratio matrix 

#for each K find lambda that minimizes BIC 
#(BIC measures the model's error penalized by its complexity: lowest BIC is the perfect trade-off where the models captures the biological signals without including too much noise)
minBICid = apply(BIC, 2, which.min)

#for each k take the deviance ratio value corresponding to the perfect lambda (lowest BIC value)
#Deviance Ratio: % of explained variation from the clusters, the higher the better the clusters separate patients
devRatMinBIC = rep(NA, nK)
for(i in 1:nK){
  devRatMinBIC[i] = devR[minBICid[i], i]
}


#Elbow Plot: look for the point where adding another k doesn't yield an improvement in % explained var.
#x-axis: k clusters; y-axis: % explained variation by each cluster
#pick the K after which the curve plateaus (adding an additional cluster doesn't lead to an improvement in % explained variation) 

plot(1:(nK + 1), c(0, devRatMinBIC), type = "b", pch = 19, col = "blue",
     xlab = "Number of Clusters (K + 1)",
     ylab = "% Explained Variation",
     main = "Model Selection")

#%explained variance table
var_table <- data.frame(
  Clusters = 2:(nK + 1),
  Percent_EV = devRatMinBIC * 100 #%
)
print(var_table)

#between 3 and 4 clusters there is a 3.2% increase in explained variation (not big but not even negligible). a suggestion for noisy data is to plot heatmaps (in this case for 3 and 4 clusters) to see which one shows clearer patterns.

clusters = getClusters(output)

#3 clusters
best_cluster = clusters[, 2]
best_fit = output[[2]]$fit[[minBICid[2]]]
cat("fit for 3 clusters (output[[2]]$fit[[minBICid[2]]])", best_fit)

#4 clusters
best_cluster_4 = clusters[, 3] 
best_fit_4 = output[[3]]$fit[[minBICid[3]]]
cat("fit for 4 clusters (output[[3]]$fit[[minBICid[3]]])", best_fit_4)

#2 clusters
best_cluster_2 = clusters[, 1] 
best_fit_2 = output[[1]]$fit[[minBICid[1]]]
cat("fit for 2 clusters (output[[1]]$fit[[minBICid[1]]])", best_fit_2)


#extract the matrix space for internal benchmarking 
# for the 3-cluster solutions
#icluster_z <- best_fit$meanZ # Rows = Patients, Columns = Latent Factors (Z1, Z2, ...)
#assign patient IDs
#rownames(icluster_z) <- rownames(rna_scaled)
#saveRDS(icluster_z, file= "icluster_matrix_3.rds")

#for the 4-cluster solutions
icluster_z4 <- best_fit_4$meanZ # Rows = Patients, Columns = Latent Factors (Z1, Z2, ...)
#assign patient IDs
rownames(icluster_z4) <- rownames(rna_scaled)
saveRDS(icluster_z4, file= "icluster_matrix_4.rds")

icluster_z2 <- best_fit_2$meanZ # Rows = Patients, Columns = Latent Factors (Z1, Z2, ...)
#assign patient IDs
rownames(icluster_z2) <- rownames(rna_scaled)
saveRDS(icluster_z2, file= "icluster_matrix_2.rds")
#```

```{r}
lambda_prova = output[[2]]$lambda[minBICid[2],]

```

# Heatmaps

```{r}
library(gplots)

col.scheme = alist()
col.scheme[[1]] = bluered(256) # RNA-seq 
col.scheme[[2]] = bluered(256) # Methylation 
col.scheme[[3]] = colorpanel(256, low="blue", mid="white", high="red") # CNV
```

```{r}
#librery requested by iClusterPlus
library(lattice)

hm1 <- plotHeatmap(
  fit = best_fit, 
  datasets = list(rna_scaled, meth_scaled, cnv_scaled), 
  type = c("gaussian", "gaussian", "gaussian"),
  sample.order = NULL,               
  sparse = c(TRUE, TRUE, TRUE),      
  threshold = c(0.25, 0.25, 0.25), 
  col.scheme = col.scheme,
  plot.chr = c(FALSE, FALSE, FALSE), 
  cap = c(0, 0.99, 0.99)
)
```
rows: features
columns: samples

clear cluster separation

RNAseq: 2  blocks (one red/high expression in Cluster 1, white/low in Clusters 2 & 3; the other blue/downregulated in Cluster 1, red/upregulated in Clusters 2 & 3) point to a binary transcriptomic switch that  isolates Cluster 1 from the other two.



# 4 clusters (to check if it is truly capturing a distinct subtype or just splitting noise)

```{r}
best_cluster_4 = clusters[, 3] 
best_fit_4 = output[[3]]$fit[[minBICid[3]]]
hm2 <- plotHeatmap(
  fit = best_fit_4, 
  datasets = list(rna_scaled, meth_scaled, cnv_scaled), 
  type = c("gaussian", "gaussian", "gaussian"),
  sample.order = NULL,               
  sparse = c(TRUE, TRUE, TRUE),      
  threshold = c(0.25, 0.25, 0.25), 
  col.scheme = col.scheme,
  plot.chr = c(FALSE, FALSE, FALSE), 
  cap = c(0, 0.99, 0.99)
)

hm2 <- plotHeatmap(
  fit = best_fit_2, 
  datasets = list(rna_scaled, meth_scaled, cnv_scaled), 
  type = c("gaussian", "gaussian", "gaussian"),
  sample.order = NULL,               
  sparse = c(TRUE, TRUE, TRUE),      
  threshold = c(0.25, 0.25, 0.25), 
  col.scheme = col.scheme,
  plot.chr = c(FALSE, FALSE, FALSE), 
  cap = c(0, 0.99, 0.99)
)
```
The 4th cluster comes from separating the big one into 2; for every data type there are two clusters that are similar, so it might be oversplitting the data.


# Feature extraction

```{r}
# Create a named list of your feature names
features_list <- list(
  RNAseq      = colnames(rna_scaled),
  Methylation = colnames(meth_scaled),
  CNV         = colnames(cnv_scaled)
)

# Safely extract significant features for all 3 layers
sigfeatures <- lapply(1:3, function(i) {
  # Calculate row sums of absolute lasso coefficients
  rowsum <- apply(abs(best_fit$beta[[i]]), 1, sum)
  
  # Determine the 75th percentile cutoff threshold --> magari si può aggiustare perche potrebbe essere troppo o troppo poco, quindi vdere la distribuzione e prendere l'elbow
  upper <- quantile(rowsum, prob = 0.75)
  
  # Return the names of features exceeding the threshold
  return(features_list[[i]][which(rowsum > upper)])
})

# Assign the names directly to the resulting list
names(sigfeatures) <- c("RNAseq", "Methylation", "CNV")
```




```{r}
sapply(sigfeatures, length)
```


# SAVE FILES

#```{r}
iclusters_3 <- data.frame(
  SampleID = common_samples,
  Cluster  = best_cluster
)
write.csv(iclusters_3, file = "iclusters_clusters_3.csv", row.names = FALSE)
#```

