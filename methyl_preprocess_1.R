# Methylation data preprocessing

library(flexmix)
library(parallel)
library(lattice)
library(readxl)
library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(sesameData)
library(ELMER)
library(SummarizedExperiment)
library(MultiAssayExperiment)
library(GenomeInfoDb)
library(GenomicRanges)
source("/home/alessandrodelle@vhio.org/ondemand/TNBCstrat/function_correctBetas.R")

# Point to the root mount where the data is stored
probes <- read.table("/mnt/petasan_ccb/alessandro/SCANB/matched_methyl.tsv", header=TRUE, sep="\t", row.names=1)
exp.data <- read.table("/mnt/petasan_ccb/alessandro/SCANB/rna_logtransformed.tsv", header=TRUE, sep="\t", row.names=1)
dt_metadata <- read_excel("/mnt/petasan_ccb/alessandro/SCANB/ids_cohorts_match.xlsx", sheet = "1a SCAN-B discovery")

## Purity adjustment using ASCAT purity estimates

purity_vector <- as.numeric(dt_metadata$ASCAT_TUM_FRAC)
names(purity_vector) <- dt_metadata$PD_ID

samples <- intersect(colnames(probes), names(purity_vector))
probes <- probes[, samples]
purity_vector <- purity_vector[samples]

print("Starting purity adjustment for every cpg")
# Apply adjustBeta() on every row (cpg)
results_list <- apply(probes, MARGIN = 1, FUN = function(cpg_row) {
  adjustBeta(methylation = cpg_row, 
             purity = purity_vector,
             snames = samples, 
             seed = FALSE)
})

print("Iteration for purity adjustment finished!")

# extract the correct tumor values frm the list (y.tum)
adjusted_data <- do.call(rbind, lapply(results_list, function(x) x$y.tum))
colnames(adjusted_data) <- samples

print("Purity adjustement finished!")

print("Saving preprocessed data to Petasan...")
saveRDS(adjusted_data, "/mnt/petasan_ccb/alessandro/SCANB/adjusted_data.rds")    
