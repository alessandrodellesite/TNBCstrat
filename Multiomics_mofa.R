library(factoextra)
reticulate::py_install("mofapy2", pip = TRUE)
reticulate::py_module_available("mofapy2")
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
plot_variance_explained(MOFAobject, max_r2 = 100)

ggsave("/mnt/petasan_ccb/alessandro/SCANB/plots/multiomics/mofa/variance_decomposition.png",
       plot = plot_variance_explained(MOFAobject, max_r2 = 100),
       width = 10, height = 5, dpi = 300)

#Factors should be largely uncorrelated (orthogonal): correlation suggests poor model fit
plot_factor_cor(MOFAobject)

ggsave("/mnt/petasan_ccb/alessandro/SCANB/plots/multiomics/mofa/factor_correlation.png",
       plot = plot_factor_cor(MOFAobject),
       width = 10, height = 5, dpi = 300)

# Total variance explained per view (using all factors)
plot_variance_explained(MOFAobject, plot_total = T)[[2]]

ggsave("/mnt/petasan_ccb/alessandro/SCANB/plots/multiomics/mofa/tot_variance_per_view.png",
       plot = plot_variance_explained(MOFAobject, plot_total = T)[[2]],
       width = 10, height = 5, dpi = 300)

#How much variance each factor captures across the different data types, to choose which factors to consider for the k-means. 

# Use the official MOFA2 function to calculate R2
r2_list <- get_variance_explained(MOFAobject)

# Access the factor-specific matrix safely
# MOFA2 already converts these to actual percentages (0 to 100) 
factors_r2 <- r2_list$r2_per_factor[[1]]

# View the result
round(factors_r2, 2)

#extract the factor values for each sample 
factors_matrix <- do.call(rbind, get_factors(MOFAobject, factors = "all")) #to compress into a single matrix
factors_to_use <- factors_matrix[, 1:6] # choose the factors to include based previous chunck (at least one variance over 2.5)


set.seed(123) 
elbow_plot <- factoextra::fviz_nbclust(factors_to_use, kmeans, method = "wss")
#print(elbow_plot)
ggsave("/mnt/petasan_ccb/alessandro/SCANB/plots/multiomics/mofa/elbow_plot.png",
       plot = elbow_plot,
       width = 10, height = 5, dpi = 300)


#Silhouette method on MOFA factors
set.seed(123)
silhouette_plot <- factoextra::fviz_nbclust(factors_to_use, kmeans, method = "silhouette", nstart = 50, iter.max = 100)

ggsave("/mnt/petasan_ccb/alessandro/SCANB/plots/multiomics/mofa/silhouette.png",
       plot = silhouette_plot,
       width = 10, height = 5, dpi = 300)

#Gap Statistic method
set.seed(123)
gap_stat <- factoextra::fviz_nbclust(factors_to_use, kmeans, method = "gap_stat", nstart = 50, nboot = 500)

ggsave("/mnt/petasan_ccb/alessandro/SCANB/plots/multiomics/mofa/gap_stat.png",
       plot = gap_stat,
       width = 10, height = 5, dpi = 300)

saveRDS(factors_to_use, file = "/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_mofa/mofa_factors.rds")


# Run k-means for all 3 options options
set.seed(123)
km2 <- kmeans(factors_to_use, centers = 2, nstart = 50)
km3 <- kmeans(factors_to_use, centers = 3, nstart = 50)
km4 <- kmeans(factors_to_use, centers = 4, nstart = 50)

# Plot Factor 1 vs Factor 2 colored by the 2/3-cluster solution
plot(factors_to_use[,1], factors_to_use[,2], col = km2$cluster, 
     pch = 19, xlab = "Factor 1", ylab = "Factor 2", 
     main = "K-means with K=2")

plot(factors_to_use[,1], factors_to_use[,2], col = km3$cluster, 
     pch = 19, xlab = "Factor 1", ylab = "Factor 2", 
     main = "K-means with K=3")

plot(factors_to_use[,1], factors_to_use[,2], col = km4$cluster, 
     pch = 19, xlab = "Factor 1", ylab = "Factor 2", 
     main = "K-means with K=4")




```{r}
#k-means clustering on the MOFA factors (3 subtype specified)
cluster_results_2 <- km2

#sdd cluster assignments to metadata 
MOFAobject@samples_metadata$Subtype <- as.factor(cluster_results_2$cluster)

#plot factors colored by clusters
plot_factor(MOFAobject, 
            factors = c(1:15), 
            color_by = "Subtype")
```

#```{r}
#Convert the k-means cluster vector into data frame
export_mofa_clusters <- data.frame(
  SampleID = names(km4$cluster),
  Cluster  = as.vector(km4$cluster)
)
write.csv(export_mofa_clusters, "/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_mofa/mofa_km_clusters_4.csv", row.names = FALSE)
#```

```{r}
cluster_results_3 <- km4

#sdd cluster assignments to metadata 
MOFAobject@samples_metadata$Subtype <- as.factor(cluster_results_3$cluster)

#plot factors colored by clusters
plot_factor(MOFAobject, 
            factors = c(1:15), 
            color_by = "Subtype")
```





# Add sample metadata (da ricontrollare se corretto)

```{r}
library(readxl)
library(dplyr)

dt_metadata <- read_excel("ids_cohorts_match.xlsx", sheet = "1a SCAN-B discovery")

# MOFA strictly requires sample column named 'sample' 
dt_metadata <- dt_metadata %>% 
  rename(sample = PD_ID)

# check if all samples in mofa object == metadata file
mofa_samples <- unlist(samples_names(MOFAobject))
missing_metadata <- setdiff(mofa_samples, dt_metadata$sample)

# filter the metadata so it only includes the samples present in mofa model
dt_metadata_cleaned <- dt_metadata %>% 
  filter(sample %in% mofa_samples)

# inject the metadata into the model
samples_metadata(MOFAobject) <- as.data.frame(dt_metadata_cleaned)
```

# Association analysis
test association between MOFA factors and some metadata

```{r}
correlate_factors_with_covariates(MOFAobject, 
  covariates = c("TMB","TILs","Age","ASCAT_PLOIDY","CibersortX.Tcell", 
  "CibersortX.Bcell",
  "CibersortX.macrophage",
  "CibersortX.stroma", 
  "CibersortX.endothelial", 
  "CibersortX.epithelial"), 
  plot="log_pval"
)
```

strong association → the factor explains a certain amount of variance that is dependent on the metadata
