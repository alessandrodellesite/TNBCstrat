#Analysis of icluster model tuning
#BiocManager::install("iClusterPlus")
library(iClusterPlus)
library(gplots)
library(lattice)
library(ComplexHeatmap)

data_dir <- "/mnt/petasan_ccb/alessandro/SCANB/multiomics/input_data/icluster_inputdata/"
rna_scaled  <- readRDS(file.path(data_dir, "rna_icluster.rds"))
meth_scaled <- readRDS(file.path(data_dir, "met_icluster.rds"))
cnv_scaled  <- readRDS(file.path(data_dir, "cnv_icluster.rds"))

# visualize model tuning results: how many clusters to choose

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


png("/mnt/petasan_ccb/alessandro/SCANB/plots/multiomics/icluster/elbow_plot.png",
    width = 10, height = 5, units = "in", res = 300)

plot(1:(nK + 1), 
     c(0, devRatMinBIC), 
     type = "b", pch = 19, 
     col = "blue", 
     xlab = "Number of Clusters (K + 1)", 
     ylab = "% Explained Variation", 
     main = "Model Selection")

dev.off()

#%explained variance table
var_table <- data.frame(
  Clusters = 2:(nK + 1),
  Percent_EV = devRatMinBIC * 100 #%
)
print(var_table)
write.csv(var_table, "/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_icluster/expl_variance.csv")

#between 3 and 4 clusters there is a 3.2% increase in explained variation (not big but not even negligible). a suggestion for noisy data is to plot heatmaps (in this case for 3 and 4 clusters) to see which one shows clearer patterns.

clusters = getClusters(output)

#2 clusters
best_cluster_2 = clusters[, 1] 
best_fit_2 = output[[1]]$fit[[minBICid[1]]]
#cat("fit for 2 clusters (output[[1]]$fit[[minBICid[1]]]):\n")
#print(best_fit_2)
lambda_2 <- output[[1]]$lambda[minBICid[1],]
print("best lambda values for 2 clusters (lowest BIC):\n")
print(lambda_2)

#3 clusters
best_cluster = clusters[, 2]
best_fit = output[[2]]$fit[[minBICid[2]]]
#cat("fit for 3 clusters (output[[2]]$fit[[minBICid[2]]]):\n")
#print(best_fit)
lambda_3 <- output[[2]]$lambda[minBICid[2],]
print("best lambda values for 3 clusters (lowest BIC)")
print(lambda_3)

#4 clusters
best_cluster_4 = clusters[, 3] 
best_fit_4 = output[[3]]$fit[[minBICid[3]]]
#cat("fit for 4 clusters (output[[3]]$fit[[minBICid[3]]]):\n")
#print(best_fit_4)
lambda_4 <- output[[3]]$lambda[minBICid[3],]
print("best lambda values for 4 clusters (lowest BIC)")
print(lambda_4)

#5 clusters
best_cluster_5 = clusters[, 4] 
best_fit_5 = output[[4]]$fit[[minBICid[4]]]
#cat("fit for 5 clusters (output[[4]]$fit[[minBICid[4]]]):\n")
#print(best_fit_5)
lambda_5 <- output[[4]]$lambda[minBICid[4],]
print("best lambda values for 5 clusters (lowest BIC)")
print(lambda_5)


# Heatmaps

col.scheme = alist()
col.scheme[[1]] = bluered(256) # RNA-seq 
col.scheme[[2]] = bluered(256) # Methylation 
col.scheme[[3]] = colorpanel(256, low="blue", mid="white", high="red") # CNV

# 2 clusters
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

png("/mnt/petasan_ccb/alessandro/SCANB/plots/multiomics/icluster/heatmap_omics_2.png", width = 10, height = 8, units = "in", res = 300)
draw(hm2, 
     heatmap_legend_side = "right", 
     annotation_legend_side = "right"
)
dev.off()

# 3 clusters
hm3 <- plotHeatmap(
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

png("/mnt/petasan_ccb/alessandro/SCANB/plots/multiomics/icluster/heatmap_omics_3.png", width = 10, height = 8, units = "in", res = 300)
draw(hm3, 
     heatmap_legend_side = "right", 
     annotation_legend_side = "right"
)
dev.off()

# 4 clusters
hm4 <- plotHeatmap(
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

png("/mnt/petasan_ccb/alessandro/SCANB/plots/multiomics/icluster/heatmap_omics_4.png", width = 10, height = 8, units = "in", res = 300)
draw(hm4, 
     heatmap_legend_side = "right", 
     annotation_legend_side = "right"
)
dev.off()

# 5 clusters
hm5 <- plotHeatmap(
  fit = best_fit_5, 
  datasets = list(rna_scaled, meth_scaled, cnv_scaled), 
  type = c("gaussian", "gaussian", "gaussian"),
  sample.order = NULL,               
  sparse = c(TRUE, TRUE, TRUE),      
  threshold = c(0.25, 0.25, 0.25), 
  col.scheme = col.scheme,
  plot.chr = c(FALSE, FALSE, FALSE), 
  cap = c(0, 0.99, 0.99)
)

png("/mnt/petasan_ccb/alessandro/SCANB/plots/multiomics/icluster/heatmap_omics_5.png", width = 10, height = 8, units = "in", res = 300)
draw(hm5, 
     heatmap_legend_side = "right", 
     annotation_legend_side = "right"
)
dev.off()


# Feature extraction for each result

# Create a named list of your feature names
features_list <- list(
  RNAseq      = colnames(rna_scaled),
  Methylation = colnames(meth_scaled),
  CNV         = colnames(cnv_scaled)
)

# 2-clusters: extract significant features for all 2 layers
sigfeatures_2 <- lapply(1:2, function(i) {
  # Calculate row sums of absolute lasso coefficients
  rowsum <- apply(abs(best_fit_2$beta[[i]]), 1, sum)
  
  # Determine the 75th percentile cutoff threshold --> magari si può aggiustare perche potrebbe essere troppo o troppo poco, quindi vdere la distribuzione e prendere l'elbow
  upper <- quantile(rowsum, prob = 0.75)
  
  # Return the names of features exceeding the threshold
  return(features_list[[i]][which(rowsum > upper)])
})

# Assign the names directly to the resulting list
names(sigfeatures_2) <- c("RNAseq", "Methylation", "CNV")
print("Feature extraction for 2-clusters solution")
print(sapply(sigfeatures_2, length))

# 3-clusters: extract significant features for all 3 layers
sigfeatures_3 <- lapply(1:3, function(i) {
  # Calculate row sums of absolute lasso coefficients
  rowsum <- apply(abs(best_fit$beta[[i]]), 1, sum)
  
  # Determine the 75th percentile cutoff threshold --> magari si può aggiustare perche potrebbe essere troppo o troppo poco, quindi vdere la distribuzione e prendere l'elbow
  upper <- quantile(rowsum, prob = 0.75)
  
  # Return the names of features exceeding the threshold
  return(features_list[[i]][which(rowsum > upper)])
})

# Assign the names directly to the resulting list
names(sigfeatures_3) <- c("RNAseq", "Methylation", "CNV")
print("Feature extraction for 3-clusters solution")
print(sapply(sigfeatures_3, length))

# 4-clusters: extract significant features for all 4 layers
sigfeatures_4 <- lapply(1:4, function(i) {
  # Calculate row sums of absolute lasso coefficients
  rowsum <- apply(abs(best_fit_4$beta[[i]]), 1, sum)
  
  # Determine the 75th percentile cutoff threshold --> magari si può aggiustare perche potrebbe essere troppo o troppo poco, quindi vdere la distribuzione e prendere l'elbow
  upper <- quantile(rowsum, prob = 0.75)
  
  # Return the names of features exceeding the threshold
  return(features_list[[i]][which(rowsum > upper)])
})

# Assign the names directly to the resulting list
names(sigfeatures_4) <- c("RNAseq", "Methylation", "CNV")
print("Feature extraction for 4-clusters solution")
print(sapply(sigfeatures_4, length))

# 5-clusters: extract significant features for all 5 layers
sigfeatures_5 <- lapply(1:5, function(i) {
  # Calculate row sums of absolute lasso coefficients
  rowsum <- apply(abs(best_fit_5$beta[[i]]), 1, sum)
  
  # Determine the 75th percentile cutoff threshold --> magari si può aggiustare perche potrebbe essere troppo o troppo poco, quindi vdere la distribuzione e prendere l'elbow
  upper <- quantile(rowsum, prob = 0.75)
  
  # Return the names of features exceeding the threshold
  return(features_list[[i]][which(rowsum > upper)])
})

# Assign the names directly to the resulting list
names(sigfeatures_5) <- c("RNAseq", "Methylation", "CNV")
print("Feature extraction for 5-clusters solution")
print(sapply(sigfeatures_5, length))



# EXTRACT MATRIX SPACE for internal benchmarking
# for the 3-cluster solutions
icluster_z <- best_fit$meanZ # Rows = Patients, Columns = Latent Factors (Z1, Z2, ...)
#assign patient IDs
rownames(icluster_z) <- rownames(rna_scaled)
saveRDS(icluster_z, file= "/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_icluster/icluster_matrix_3.rds")

#for the 2-cluster solutions
icluster_z2 <- best_fit_2$meanZ # Rows = Patients, Columns = Latent Factors (Z1, Z2, ...)
#assign patient IDs
rownames(icluster_z2) <- rownames(rna_scaled)
saveRDS(icluster_z2, file= "/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_icluster/icluster_matrix_2.rds")

#for the 4-cluster solutions
icluster_z4 <- best_fit_4$meanZ # Rows = Patients, Columns = Latent Factors (Z1, Z2, ...)
#assign patient IDs
rownames(icluster_z4) <- rownames(rna_scaled)
saveRDS(icluster_z4, file= "/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_icluster/icluster_matrix_4.rds")

#for the 5-cluster solutions
icluster_z5 <- best_fit_5$meanZ # Rows = Patients, Columns = Latent Factors (Z1, Z2, ...)
#assign patient IDs
rownames(icluster_z5) <- rownames(rna_scaled)
saveRDS(icluster_z5, file= "/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_icluster/icluster_matrix_5.rds")


# Save clustering results

data_dir <- "/mnt/petasan_ccb/alessandro/SCANB/multiomics/input_data/icluster_inputdata/"
rna_scaled  <- readRDS(file.path(data_dir, "rna_icluster.rds"))
sample_ids <- rownames(rna_scaled)

export_2 <- data.frame(
  SampleID = sample_ids,
  Cluster  = best_cluster_2
)
write.csv(export_2, "/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_icluster/iclusters_clusters_2.csv", row.names = FALSE)

export_3 <- data.frame(
  SampleID = sample_ids,
  Cluster  = best_cluster
)
write.csv(export_3, "/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_icluster/iclusters_clusters_3.csv", row.names = FALSE)

export_4 <- data.frame(
  SampleID = sample_ids,
  Cluster  = best_cluster_4
)
write.csv(export_4, "/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_icluster/iclusters_clusters_4.csv", row.names = FALSE)

export_5 <- data.frame(
  SampleID = sample_ids,
  Cluster  = best_cluster_5
)
write.csv(export_5, "/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_icluster/iclusters_clusters_5.csv", row.names = FALSE)

