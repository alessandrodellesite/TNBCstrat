library(iClusterPlus)

data_dir <- "/mnt/petasan_ccb/alessandro/SCANB/multiomics/input_data/icluster_inputdata/"
rna_scaled  <- readRDS(file.path(data_dir, "rna_icluster.rds"))
meth_scaled <- readRDS(file.path(data_dir, "met_icluster.rds"))
cnv_scaled  <- readRDS(file.path(data_dir, "cnv_icluster.rds"))

# iCluster model tuning

# Run the function for a range of K latent variables (Clusters = K + 1)
args <- commandArgs(trailingOnly = TRUE)
k <- as.integer(args[1])

cat("Starting model tuning for K =", k, "at", as.character(Sys.time()), "\n")

set.seed(123) #seed for reproducibility
cv.fit <- tune.iClusterPlus(
  cpus = 16, 
  dt1 = rna_scaled, dt2 = meth_scaled, dt3 = cnv_scaled,
  type = c("gaussian", "gaussian", "gaussian"), 
  K = k, 
  n.lambda = 185,
  scale.lambda = c(1, 1, 1),
  maxiter = 20
)

save(cv.fit, file = paste0("cv.fit.k", k, ".Rdata"))
