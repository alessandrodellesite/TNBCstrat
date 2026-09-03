library(NMF)
library(ggplot2)

dir_rna <- "/mnt/petasan_ccb/alessandro/SCANB/singleomic/nmf_rna"
dir_meth <- "/mnt/petasan_ccb/alessandro/SCANB/singleomic/nmf_meth"
dir_cnv <- "/mnt/petasan_ccb/alessandro/SCANB/singleomic/nmf_cnv"

#rna
estim.r_rna <- readRDS("NMF_rank_estimation_RNAseq.rds")
res_rna     <- readRDS("NMF_final_results_RNAseq.rds")

ggsave("/mnt/petasan_ccb/alessandro/SCANB/plots/singleomic_nmf/rna_nmf/estim.png",
      plot= plot(estim.r_rna),
      width=10, height=5, dpi= 300 )

ggsave("/mnt/petasan_ccb/alessandro/SCANB/plots/singleomic_nmf/rna_nmf/consensus.png",
      plot= consensusmap(res_rna),
      width=10, height=5, dpi= 300 )
ggsave("/mnt/petasan_ccb/alessandro/SCANB/plots/singleomic_nmf/rna_nmf/coef.png",
      plot= coefmap(res_rna),
      width=10, height=5, dpi= 300 )

sample_groups_rna <- predict(res_rna)
table(sample_groups_rna)



