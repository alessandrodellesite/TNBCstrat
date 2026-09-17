library(readxl)
library(dplyr)
library(cluster)
library(ggplot2)
#reticulate::py_install("mofapy2", pip = TRUE)
#reticulate::py_module_available("mofapy2")
library(MOFA2)

# load data
data_dir <- "/mnt/petasan_ccb/alessandro/SCANB/multiomics/input_data/mofa_inputdata/"
rna_cent <- readRDS(file.path(data_dir, "rna_mofa.rds"))
met_cent <- readRDS(file.path(data_dir, "met_mofa.rds"))
cnv_cent <- readRDS(file.path(data_dir, "cnv_mofa.rds"))

multiomics_data <- list(
  RNAseq      = rna_cent,        
  Methylation = met_cent,        
  CNV         = cnv_cent         
)
lapply(multiomics_data, dim)

# order columns in the same way
common_samples <- intersect(
  intersect(colnames(multiomics_data$RNAseq), colnames(multiomics_data$Methylation)),
  colnames(multiomics_data$CNV)
)

multiomics_data$RNAseq      <- multiomics_data$RNAseq[, common_samples]
multiomics_data$Methylation <- multiomics_data$Methylation[, common_samples]
multiomics_data$CNV         <- multiomics_data$CNV[, common_samples]


# Create the untrained object
MOFAobject <- create_mofa(multiomics_data)

print(MOFAobject)

# Plot data overview to make sure overlapping samples are aligned properly
plot_data_overview(MOFAobject)

ggsave("/mnt/petasan_ccb/alessandro/SCANB/plots/multiomics/mofa/overlapping_samples.png",
       plot = plot_data_overview(MOFAobject),
       width = 10, height = 5, dpi = 300)

#All three data are continuous
#Default Data Options
data_opts <- get_default_data_options(MOFAobject)
model_opts <- get_default_model_options(MOFAobject)

# SET SCALE TO TRUE --> can have very different total variances. Without scaling each view to unit variance, MOFA's inference is easier to dominate by whichever view happens to have larger overall variance
data_opts$scale_views <- TRUE

# 20 defined as the upper limit of factors MOFA will test
model_opts$num_factors <- 20

#likelihoods set to gaussian for all three layers
model_opts$likelihoods <- c(
  RNAseq      = "gaussian",
  Methylation = "gaussian",
  CNV         = "gaussian"
)

#default training options
train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "slow"  
train_opts$seed <- 123 #set seed for reproducibility

#prepare the object
MOFAobject <- prepare_mofa(
  MOFAobject,
  data_options     = data_opts,
  model_options    = model_opts,
  training_options = train_opts
)


#Train the model

MOFAobject <- run_mofa(MOFAobject, outfile = "/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_mofa/mofafirst_mofa_model.hdf5")

# Save your R object for easier reloading later
#saveRDS(MOFAobject, "my_mofa_model.rds")

# Variance decomposition to see how much each factor explains per view
var_dec <- plot_variance_explained(MOFAobject, max_r2 = 100)

ggsave("/mnt/petasan_ccb/alessandro/SCANB/plots/multiomics/mofa/variance_decomposition.png",
       plot = var_dec,
       width = 10, height = 5, dpi = 300)

#Factors should be largely uncorrelated (orthogonal): correlation suggests poor model fit

png("/mnt/petasan_ccb/alessandro/SCANB/plots/multiomics/mofa/factor_correlation.png",
    width = 10, height = 5, units = "in", res = 300)
plot_factor_cor(MOFAobject)
dev.off()

# Total variance explained per view (using all factors)
tot_var <- plot_variance_explained(MOFAobject, plot_total = T)[[2]]

ggsave("/mnt/petasan_ccb/alessandro/SCANB/plots/multiomics/mofa/tot_variance_per_view.png",
       plot = tot_var,
       width = 10, height = 5, dpi = 300)

#How much variance each factor captures across the different data types, to choose which factors to consider for the k-means. 

# Use the official MOFA2 function to calculate R2
r2_list <- get_variance_explained(MOFAobject)

# Access the factor-specific matrix safely
# MOFA2 already converts these to actual percentages (0 to 100) 
factors_r2 <- r2_list$r2_per_factor[[1]]

# View the result
fact_file <- round(factors_r2, 2)
write.csv(fact_file, "/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_mofa/factors_values.csv")

#extract the factor values for each sample 
factors_matrix <- do.call(rbind, get_factors(MOFAobject, factors = "all")) #to compress into a single matrix
factors_to_use <- factors_matrix[, 1:6] # choose the factors to include based previous chunck (at least one variance over 2.5)



out_dir <- "/mnt/petasan_ccb/alessandro/SCANB/plots/multiomics/mofa"

# Elbow (WSS) method
set.seed(123)
k_range <- 1:10
wss <- sapply(k_range, function(k) {
  kmeans(factors_to_use, centers = k, nstart = 50, iter.max = 100)$tot.withinss
})

png(file.path(out_dir, "elbow_plot.png"), width = 10, height = 5, units = "in", res = 300)
plot(k_range, wss, type = "b", pch = 19, xaxt = "n",
     xlab = "Number of clusters k", ylab = "Total Within Sum of Square",
     main = "Optimal number of clusters")
axis(1, at = k_range)
dev.off()

# Silhouette method
set.seed(123)
k_range_sil <- 2:10
sil_width <- sapply(k_range_sil, function(k) {
  km <- kmeans(factors_to_use, centers = k, nstart = 50, iter.max = 100)
  mean(silhouette(km$cluster, dist(factors_to_use))[, "sil_width"])
})

png(file.path(out_dir, "silhouette.png"), width = 10, height = 5, units = "in", res = 300)
plot(k_range_sil, sil_width, type = "b", pch = 19, xaxt = "n",
     xlab = "Number of clusters k", ylab = "Average silhouette width",
     main = "Optimal number of clusters")
axis(1, at = k_range_sil)
dev.off()

# Gap Statistic method
set.seed(123)
gap_stat <- clusGap(factors_to_use, FUN = kmeans, nstart = 50, iter.max = 100, K.max = 10, B = 500)
gap_df <- as.data.frame(gap_stat$Tab)
gap_df$k <- seq_len(nrow(gap_df))

png(file.path(out_dir, "gap_stat.png"), width = 10, height = 5, units = "in", res = 300)
plot(gap_df$k, gap_df$gap, type = "b", pch = 19, xaxt = "n",
     xlab = "Number of clusters k", ylab = "Gap statistic",
     main = "Optimal number of clusters",
     ylim = range(c(gap_df$gap - gap_df$SE.sim, gap_df$gap + gap_df$SE.sim)))
axis(1, at = gap_df$k)
arrows(gap_df$k, gap_df$gap - gap_df$SE.sim,
       gap_df$k, gap_df$gap + gap_df$SE.sim,
       angle = 90, code = 3, length = 0.05)
dev.off()

optimal_k <- maxSE(gap_df$gap, gap_df$SE.sim, method = "Tibs2001SE")
print("Optimal k according to gap stat: ")
print(optimal_k)

saveRDS(factors_to_use, file = "/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_mofa/mofa_factors.rds")


# Run k-means for all 3 options options
set.seed(123)
km2 <- kmeans(factors_to_use, centers = 2, nstart = 50, iter.max = 100)
km3 <- kmeans(factors_to_use, centers = 3, nstart = 50, iter.max = 100)
km4 <- kmeans(factors_to_use, centers = 4, nstart = 50, iter.max = 100)

# Plot Factor 1 vs Factor 2 colored by the 2/3-cluster solution
png("/mnt/petasan_ccb/alessandro/SCANB/plots/multiomics/mofa/factor1_2_km2.png",
    width = 10, height = 8, units = "in", res = 300)
plot(factors_to_use[,1], factors_to_use[,2], col = km2$cluster, 
     pch = 19, xlab = "Factor 1", ylab = "Factor 2", 
     main = "K-means with K=2")
dev.off()

png("/mnt/petasan_ccb/alessandro/SCANB/plots/multiomics/mofa/factor1_2_km3.png",
    width = 10, height = 8, units = "in", res = 300)
plot(factors_to_use[,1], factors_to_use[,2], col = km3$cluster, 
     pch = 19, xlab = "Factor 1", ylab = "Factor 2", 
     main = "K-means with K=3")
dev.off()

png("/mnt/petasan_ccb/alessandro/SCANB/plots/multiomics/mofa/factor1_2_km4.png",
    width = 10, height = 8, units = "in", res = 300)
plot(factors_to_use[,1], factors_to_use[,2], col = km4$cluster, 
     pch = 19, xlab = "Factor 1", ylab = "Factor 2", 
     main = "K-means with K=4")
dev.off()



#Convert the k-means cluster vector into data frame
export_mofa_clusters_2 <- data.frame(
  SampleID = names(km2$cluster),
  Cluster  = as.vector(km2$cluster)
)
write.csv(export_mofa_clusters_2, "/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_mofa/mofa_km_clusters_2.csv", row.names = FALSE)

export_mofa_clusters_3 <- data.frame(
  SampleID = names(km3$cluster),
  Cluster  = as.vector(km3$cluster)
)
write.csv(export_mofa_clusters_3, "/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_mofa/mofa_km_clusters_3.csv", row.names = FALSE)

export_mofa_clusters_4 <- data.frame(
  SampleID = names(km4$cluster),
  Cluster  = as.vector(km4$cluster)
)
write.csv(export_mofa_clusters_4, "/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_mofa/mofa_km_clusters_4.csv", row.names = FALSE)

