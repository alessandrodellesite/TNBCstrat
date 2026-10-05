library(NMF)
library(ggplot2)

# Helper function to save base/grid heatmaps cleanly
save_nmf_map <- function(filename, plot_fun, width = 10, height = 5, res = 300) {
  # Ensure destination directory exists
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  png(filename, width = width, height = height, units = "in", res = res)
  plot_fun()
  dev.off()
}

# RNA
estim.r_rna <- readRDS("/mnt/petasan_ccb/alessandro/SCANB/singleomic/nmf_rna/NMF_rank_estimation_RNAseq.rds")
res_rna     <- readRDS("/mnt/petasan_ccb/alessandro/SCANB/singleomic/nmf_rna/NMF_final_results_RNAseq.rds")

# estim plot returns a ggplot object, so ggsave works fine
ggsave("/mnt/petasan_ccb/alessandro/SCANB/plots/singleomic_nmf/rna_nmf/estim.png",
       plot = plot(estim.r_rna),
       width = 10, height = 5, dpi = 300)

# Heatmaps must be rendered to a png device
save_nmf_map("/mnt/petasan_ccb/alessandro/SCANB/plots/singleomic_nmf/rna_nmf/consensus.png",
             function() consensusmap(res_rna))

save_nmf_map("/mnt/petasan_ccb/alessandro/SCANB/plots/singleomic_nmf/rna_nmf/coef.png",
             function() coefmap(res_rna))


# METH
estim.r_meth <- readRDS("/mnt/petasan_ccb/alessandro/SCANB/singleomic/nmf_meth/NMF_rank_estimation_meth.rds")
res_meth     <- readRDS("/mnt/petasan_ccb/alessandro/SCANB/singleomic/nmf_meth/NMF_final_results_meth.rds")

ggsave("/mnt/petasan_ccb/alessandro/SCANB/plots/singleomic_nmf/meth_nmf/estim.png",
       plot = plot(estim.r_meth),
       width = 10, height = 5, dpi = 300)

save_nmf_map("/mnt/petasan_ccb/alessandro/SCANB/plots/singleomic_nmf/meth_nmf/consensus.png",
             function() consensusmap(res_meth))

save_nmf_map("/mnt/petasan_ccb/alessandro/SCANB/plots/singleomic_nmf/meth_nmf/coef.png",
             function() coefmap(res_meth))



# CNV
estim.r_cnv <- readRDS("/mnt/petasan_ccb/alessandro/SCANB/singleomic/nmf_cnv/NMF_rank_estimation_cnv.rds")
res_cnv     <- readRDS("/mnt/petasan_ccb/alessandro/SCANB/singleomic/nmf_cnv/NMF_final_results_cnv.rds")

ggsave("/mnt/petasan_ccb/alessandro/SCANB/plots/singleomic_nmf/cnv_nmf/estim.png",
       plot = plot(estim.r_cnv),
       width = 10, height = 5, dpi = 300)

save_nmf_map("/mnt/petasan_ccb/alessandro/SCANB/plots/singleomic_nmf/cnv_nmf/consensus.png",
             function() consensusmap(res_cnv))

save_nmf_map("/mnt/petasan_ccb/alessandro/SCANB/plots/singleomic_nmf/cnv_nmf/coef.png",
             function() coefmap(res_cnv))
