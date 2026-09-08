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


# METH
res_meth     <- readRDS("/mnt/petasan_ccb/alessandro/SCANB/singleomic/nmf_meth/	NMF_final_results_meth_4.rds")

save_nmf_map("/mnt/petasan_ccb/alessandro/SCANB/plots/singleomic_nmf/meth_nmf/consensus_4.png",
             function() consensusmap(res_meth))

save_nmf_map("/mnt/petasan_ccb/alessandro/SCANB/plots/singleomic_nmf/meth_nmf/coef_4.png",
             function() coefmap(res_meth))


# CNV
res_cnv     <- readRDS("/mnt/petasan_ccb/alessandro/SCANB/singleomic/nmf_cnv/NMF_final_results_cnv_6.rds")

save_nmf_map("/mnt/petasan_ccb/alessandro/SCANB/plots/singleomic_nmf/cnv_nmf/consensus_6.png",
             function() consensusmap(res_cnv))

save_nmf_map("/mnt/petasan_ccb/alessandro/SCANB/plots/singleomic_nmf/cnv_nmf/coef_6.png",
             function() coefmap(res_cnv))
