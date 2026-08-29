# Methylation data preprocessing

library(flexmix)
library(parallel)
library(lattice)
library(readxl)
library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(ELMER)
library(SummarizedExperiment)
library(MultiAssayExperiment)
library(GenomeInfoDb)
library(GenomicRanges)
source("/home/alessandrodelle@vhio.org/ondemand/TNBCstrat/function_correctBetas.R")

# Point to the root mount where the data is stored
probes <- read.table("/mnt/petasan_ccb/alessandro/SCANB/methylation_data/matched_methyl_normalized.tsv", header=TRUE, sep="\t", row.names=1)
dt_metadata <- read_excel("/mnt/petasan_ccb/alessandro/SCANB/ids_cohorts_match.xlsx", sheet = "1a SCAN-B discovery")

## Purity adjustment using ASCAT purity estimates
purity_vector <- as.numeric(dt_metadata$ASCAT_TUM_FRAC)
names(purity_vector) <- dt_metadata$PD_ID
samples <- intersect(colnames(probes), names(purity_vector))
probes <- probes[, samples]
purity_vector <- purity_vector[samples]

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
if (is.na(n_cores) || n_cores < 1) n_cores <- 1
cat("Using", n_cores, "cores (from SLURM_CPUS_PER_TASK)\n")

probes_mat <- as.matrix(probes)  # single copy, no per-row duplication

results_list <- mclapply(
  seq_len(nrow(probes_mat)),
  FUN = function(i) {
    adjustBeta(methylation = probes_mat[i, ],
               purity = purity_vector,
               snames = samples,
               seed = FALSE)
  },
  mc.cores = n_cores,
  mc.preschedule = TRUE   # chunks work into n_cores static blocks — far fewer forks
)
names(results_list) <- rownames(probes_mat)

print("Iteration for purity adjustment finished!")
adjusted_data <- do.call(rbind, lapply(results_list, function(x) x$y.tum))
colnames(adjusted_data) <- samples
print("Purity adjustment finished!")
print("Saving preprocessed data to Petasan...")
saveRDS(adjusted_data, "/mnt/petasan_ccb/alessandro/SCANB/methylation_data/adjusted_data.rds")    
