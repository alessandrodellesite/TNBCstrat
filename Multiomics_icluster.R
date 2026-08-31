

data_dir <- /mnt/petasan_ccb/alessandro/SCANB/multiomics/input_data/icluster_inputdata/"
rna_scaled  <- readRDS(file.path(data_dir, "rna_icluster.rds"))
meth_scaled <- readRDS(file.path(data_dir, "met_icluster.rds"))
cnv_scaled  <- readRDS(file.path(data_dir, "cnv_icluster.rds"))

