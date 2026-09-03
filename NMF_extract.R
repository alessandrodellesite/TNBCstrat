library(NMF)
library(ggplot2)

#rna
estim.r_rna <- readRDS("/mnt/petasan_ccb/alessandro/SCANB/singleomic/nmf_rna/NMF_rank_estimation_RNAseq.rds")
res_rna     <- readRDS("/mnt/petasan_ccb/alessandro/SCANB/singleomic/nmf_rna/NMF_final_results_RNAseq.rds")

ggsave("/mnt/petasan_ccb/alessandro/SCANB/plots/singleomic_nmf/rna_nmf/estim.png",
      plot= plot(estim.r_rna),
      width=10, height=5, dpi= 300 )

ggsave("/mnt/petasan_ccb/alessandro/SCANB/plots/singleomic_nmf/rna_nmf/consensus.png",
      plot= consensusmap(res_rna),
      width=10, height=5, dpi= 300 )
ggsave("/mnt/petasan_ccb/alessandro/SCANB/plots/singleomic_nmf/rna_nmf/coef.png",
      plot= coefmap(res_rna),
      width=10, height=5, dpi= 300 )


#----

#meth
estim.r_meth <- readRDS("/mnt/petasan_ccb/alessandro/SCANB/singleomic/nmf_meth/NMF_rank_estimation_meth.rds")
res_meth     <- readRDS("/mnt/petasan_ccb/alessandro/SCANB/singleomic/nmf_meth/NMF_final_results_meth.rds")

ggsave("/mnt/petasan_ccb/alessandro/SCANB/plots/singleomic_meth/meth_nmf/estim.png",
      plot= plot(estim.r_meth),
      width=10, height=5, dpi= 300 )

ggsave("/mnt/petasan_ccb/alessandro/SCANB/plots/singleomic_meth/meth_nmf/consensus.png",
      plot= consensusmap(res_meth),
      width=10, height=5, dpi= 300 )
ggsave("/mnt/petasan_ccb/alessandro/SCANB/plots/singleomic_meth/meth_nmf/coef.png",
      plot= coefmap(res_meth),
      width=10, height=5, dpi= 300 )

#---

#cnv
estim.r_cnv <- readRDS("/mnt/petasan_ccb/alessandro/SCANB/singleomic/nmf_cnv/NMF_rank_estimation_cnv.rds")
res_cnv     <- readRDS("/mnt/petasan_ccb/alessandro/SCANB/singleomic/nmf_cnv/NMF_final_results_cnv.rds")

ggsave("/mnt/petasan_ccb/alessandro/SCANB/plots/singleomic_cnv/cnv_nmf/estim.png",
      plot= plot(estim.r_cnv),
      width=10, height=5, dpi= 300 )

ggsave("/mnt/petasan_ccb/alessandro/SCANB/plots/singleomic_cnv/cnv_nmf/consensus.png",
      plot= consensusmap(res_cnv),
      width=10, height=5, dpi= 300 )
ggsave("/mnt/petasan_ccb/alessandro/SCANB/plots/singleomic_cnv/cnv_nmf/coef.png",
      plot= coefmap(res_cnv),
      width=10, height=5, dpi= 300 )

