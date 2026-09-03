library(NMF)
dir <- ""
estim.r <- readRDS("NMF_rank_estimation_RNAseq.rds")
res     <- readRDS("NMF_final_results_RNAseq.rds")

library(NMF)
plot(estim.r)

consensusmap(res)
coefmap(res)

sample_groups <- predict(res)
table(sample_groups)


# Convert to a data frame for saving
export_groups <- data.frame(
  SampleID = names(sample_groups),
  Cluster = as.vector(sample_groups)
)
write.csv(export_groups, "rnaseq_nmf_clusters.csv", row.names = FALSE)

