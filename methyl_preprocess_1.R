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

n_cores <- detectCores()  # or set explicitly, e.g. 64, 128...
cat("Using", n_cores, "cores\n")

# convert to a list of rows once (apply() does this internally each time, mclapply needs a list)
row_list <- split(probes, seq(nrow(probes)))
row_list <- lapply(row_list, function(x) setNames(as.numeric(x), colnames(probes)))
names(row_list) <- rownames(probes)

results_list <- mclapply(
  row_list,
  FUN = function(cpg_row) {
    adjustBeta(methylation = cpg_row,
               purity = purity_vector,
               snames = samples,
               seed = FALSE)
  },
  mc.cores = n_cores,
  mc.preschedule = FALSE  # dispatches rows one at a time
)

print("Iteration for purity adjustment finished!")

# extract the correct tumor values frm the list (y.tum)
adjusted_data <- do.call(rbind, lapply(results_list, function(x) x$y.tum))
colnames(adjusted_data) <- samples

print("Purity adjustement finished!")

print("Saving preprocessed data to Petasan...")
saveRDS(adjusted_data, "/mnt/petasan_ccb/alessandro/SCANB/methylation_data/adjusted_data.rds")    
